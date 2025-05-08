#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"

using namespace mfem;

namespace pulmtln {

Driver Driver::loadFromFile(const std::string& fn)
{
	Parser p{ fn };
	return Driver{
		p.readModel(),
		p.readDriverOptions()
	};
}

Driver::Driver(Model&& model, const DriverOptions& opts) :
	model_{ std::move(model) },
	opts_{ opts }
{}

AttrToValueMap buildAttrToValueMap(
	const NameToAttrMap& matToAtt, double value)
{
	AttrToValueMap bcs;
	for (const auto& [mat, att] : matToAtt) {
		bcs[att] = value;
	}
	return bcs;
}

std::vector<int> getAttributesInMap(const NameToAttrMap& m)
{
	std::vector<int> res;
	for (const auto& [k, v] : m) {
		res.push_back(v);
	}
	return res;
}

void exportFieldSolutions(
	const DriverOptions& opts,
	ElectrostaticSolver& s,
	const std::string name,
	bool ignoreDielectrics)
{
	if (opts.exportParaViewSolution) {
		std::string outputName{ opts.exportFolder + "/" + "ParaView/" + name };
		if (ignoreDielectrics) {
			outputName += "_no_dielectrics";
		}
		ParaViewDataCollection pd{ outputName, s.getMesh() };
		s.writeParaViewFields(pd);
	}

	if (opts.exportVisItSolution) {
		std::string outputName{ opts.exportFolder + "/" + "VisIt/" + name };
		if (ignoreDielectrics) {
			outputName += "_no_dielectrics";
		}
		VisItDataCollection dC{ outputName, s.getMesh() };
		s.writeVisItFields(dC);
	}
}


SolverParameters buildSolverParameters(
	const Model& model,
	bool ignoreDielectrics)
{
	const Materials& mats = model.getMaterials();

	SolverParameters parameters;

	auto domainToEpsr{
		buildAttrToValueMap(mats.buildNameToAttrMapFor<Dielectric>(), 1.0)
	};
	if (!ignoreDielectrics) {
		for (const auto& d : mats.dielectrics) {
			domainToEpsr.at(d.attribute) = d.relativePermittivity;
		}
	}
	parameters.domainPermittivities = domainToEpsr;
	parameters.openBoundaries = getAttributesInMap(mats.buildNameToAttrMapFor<OpenBoundary>());
	parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMapFor<PEC>(), 0.0);

	return parameters;
}

mfem::Vector getCapacitancesWithOpenBoundary(const Model& m, const DriverOptions& opts, bool ignoreDielectrics)
{
	auto parameters{ buildSolverParameters(m, ignoreDielectrics) };
	const double VPrescribed{ 1.0 };
	for (auto& [k, v] : parameters.dirichletBoundaries) {
		v = VPrescribed;
	}

	auto mesh{ *m.getMesh() };
	ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
	s.Solve();

	auto openBdrs = m.getMaterials().buildNameToAttrMapFor<OpenBoundary>();
	if (openBdrs.size() != 1) {
		throw std::runtime_error("Only one open boundary is allowed.");
	}
	auto Vb = s.averagePotentialInBoundary(openBdrs.begin()->second);
	auto Vd = VPrescribed - Vb;

	auto conductors = m.getMaterials().buildNameToAttrMapFor<PEC>();

	mfem::Vector res(conductors.size());
	for (const auto& [name, attr] : conductors) {
		auto condId = Materials::getNumberContainedInName(name);
		auto Q = s.chargeInBoundary(attr);
		res[condId] = Q / Vd;
	}

	return res;
}

mfem::DenseMatrix getCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false)
{
	// PUL capacitance matrix as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// Contains N-1 x N-1 entries for a closed problem of N conductors.
	// Contains N x N entries for a open problem of N conductors.

	// Preconditions. 
	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	const auto openness{ model.determineOpenness() };
	if (conductors.size() == 1 
		&& openness == Model::OpennessType::open) {
		throw std::runtime_error(
			"The number of conductors must be at least 2."
		);
	}

	mfem::Vector Cib(0);
	if (openness == Model::OpennessType::open) {
		Cib.SetSize(conductors.size());
		Cib = getCapacitancesWithOpenBoundary(model, opts, ignoreDielectrics);
	}

	int CSize{ 0 };
	if (openness == Model::OpennessType::closed 
		|| openness == Model::OpennessType::semiopen) {
		CSize = (int)conductors.size() - 1;
	}
	else {
		CSize = (int) conductors.size();
	}
	mfem::DenseMatrix C(CSize);

	// Solves a electrostatic problem for each conductor besides the
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		if (Materials::getNumberContainedInName(nameI) == model.getGroundConductorId() &&
			openness != Model::OpennessType::open) {
			continue;
		}

		auto parameters{ buildSolverParameters(model, ignoreDielectrics) };
		parameters.dirichletBoundaries[bdrAttI] = 1.0;

		Mesh mesh{ *model.getMesh() };

		ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
		s.Solve();

		double Vb = 0.0;
		if (openness == Model::OpennessType::open) {
			int openBoundaryTag = *parameters.openBoundaries.begin();
			Vb = s.averagePotentialInBoundary(openBoundaryTag);
		}

		// Fills row
		int col{ 0 };
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			if (Materials::getNumberContainedInName(nameJ) == model.getGroundConductorId()) {
				continue;
			}
			
			double Qi = s.chargeInBoundary(conductors.at(nameJ));
			switch (openness) {
			case Model::OpennessType::closed:
			case Model::OpennessType::semiopen:
				// C_ij = Q_i / V_j. V_j is always 1.0
				C(row, col) = Qi;
				break;
			case Model::OpennessType::open:
				// C_ij = Q_i / V_j - Cib*(V_j -V_b).
				if (row == col) {
					C(row, col) = 0.0;
				}
				else {
					C(row, col) = Qi - Cib(row) * (1.0 - Vb);
				}
				break;
			default:
				throw std::runtime_error(
					"C matrix calculation not implemented for this kind of openness."
				);
			}
			col++;
		}

		exportFieldSolutions(opts, s, nameI, ignoreDielectrics);
		row++;
	}

	return C;
}

DenseMatrix getLMatrix(const Model& model, const DriverOptions& opts)
{
	// PUL inductance matrix as defined in:
	//   Clayton Paul's book: Analysis of Multiconductor Transmission Lines
	// Contains N-1 x N-1 entries for a problem of N conductors.
	// Inductance matrix can be computed from the 
	// capacitance obtained ignoring dielectrics as
	//          L = mu0 * eps0 * C^{-1}
	auto res{ getCMatrix(model, opts, true) };
	res.Invert();
	res *= MU0_NATURAL * EPSILON0_NATURAL;
	return res;
}

PULParameters buildPULParametersForModel(const Model& m, const DriverOptions& opts)
{
	PULParameters res;

	res.C = getCMatrix(m, opts);
	res.C *= EPSILON0_SI;

	res.L = getLMatrix(m, opts);
	res.L *= MU0_SI;

	if (opts.makeMatricesSymmetric) {
		res.C.Symmetrize();
		res.L.Symmetrize();
	}

	return res;
}

PULParameters Driver::getMTLPUL() const
{
	auto res{ buildPULParametersForModel(model_, opts_) };

	if (opts_.exportMatrices) {
		res.saveToJSONFile(opts_.exportFolder + "matrices.pulmtln.out.json");
	}

	return res;
}

PULParametersByDomain Driver::getMTLPULByDomains() const
{
	PULParametersByDomain res;

auto idToDomain{ Domain::buildDomains(model_) };

for (const auto& [id, domain] : idToDomain) {
	auto globalMesh{ *model_.getMesh() };
	res.domainToPUL[id] = buildPULParametersForModel(
		Domain::buildModelForDomain(globalMesh, model_.getMaterials(), domain),
		opts_
	);
}

res.domainTree = DomainTree{ idToDomain };

return res;
}

mfem::DenseMatrix floatingPotentialsForClosedCases(const mfem::DenseMatrix& C)
{
	// Forms system of equations to determine floating potentials. 
	// Q1 = C11*V1 + C21*V2 + ....
	// For prescribed V_1 = 1.0 we have C V = Q
	//    [ C ] [1.0, V_2, ...]^T = [Q_1, 0.0, ...]
	// 
	// which can be converted to A x = b with unknowns x = [Q_1, V_2, ...]^T
	auto N = C.NumRows();
	mfem::DenseMatrix res(N, N);
	for (int i{ 0 }; i < N; ++i) {
		mfem::DenseMatrix A{ C };
		mfem::Vector negativeQ(N);
		negativeQ = 0.0;
		negativeQ(i) = -1.0;
		
		A.SetCol(i, negativeQ);
		
		mfem::Vector b(N);
		b = C.GetColumn(i);
		b *= -1.0;

		mfem::Vector x(N);
		mfem::DenseMatrixInverse Ainv(A);
		Ainv.Mult(b, x);

		for (int j{ 0 }; j < res.NumCols(); ++j) {
			if (i == j) {
				res(i, i) = 1.0;
			}
			else {
				res(i, j) = x(j);
			}
		}
	}

	return res;
}

mfem::DenseMatrix floatingPotentialsForOpenCases(const mfem::DenseMatrix& C, const mfem::Vector& Cib)
{
	auto N = C.NumRows()+1; // N is the total number of conductors.
	if (Cib.Size() != N) {
		throw std::runtime_error("Cib should have the size of the number of conductors.");
	}

	// Build expanded capacitance matrix including conductor 0 and capacitance with boundary. 
	mfem::DenseMatrix Cexp(N+1, N+1);
	// Fill upper triangular part.
	{
		mfem::DenseMatrix U(N + 1, N + 1);
		for (int r{ 0 }; r < N; ++r) {
			for (int c{ r }; c < N; ++c) {
				if (r == c) {
					U(r, c) = 0.0;
				}
				else {
					U(r, c) = C(r, c-1);
				}
			}
			U(r, N) = Cib(r);

		}
		Cexp = U;
		U.Transpose();
		Cexp += U;
	}

	//  We must solve a system with N+1 unknowns once for each of the N conductors.
	mfem::DenseMatrix res(N, N);
	for (int i{ 0 }; i < N; ++i) {
		mfem::DenseMatrix D(N + 1, N + 1);
		{
			mfem::DenseMatrix E(N + 1, N + 1);
			E.Diag(-1.0, N + 1);
			for (int r{ 0 }; r < N + 1; ++r) {
				E(r, i) += 1.0;
			}

			mfem::Mult(Cexp, E, D);
		}

		mfem::Vector minusQ(N + 1); // Q is the vector of charges.
		minusQ = 0.0;
		minusQ(i) =  1.0;
		minusQ(N) =  1.0;
		minusQ *= -1.0; 


 		mfem::DenseMatrix A(N+1, N+1);
		A = D;
		A.SetCol(i, minusQ);

		mfem::Vector b(N + 1);
		b = D.GetColumn(i);
		b *= -1.0;

		// Solves system and fills result
		mfem::Vector x(N+1);
		mfem::DenseMatrixInverse Ainv(A);
		Ainv.Mult(b, x);

		for (int j{ 0 }; j < res.NumCols(); ++j) {
			if (i == j) {
				res(i, i) = 1.0;
			}
			else {
				res(i, j) = x(j);
			}
		}
	}

	return res;
}


mfem::DenseMatrix Driver::getFloatingPotentials(const bool ignoreDielectrics) const
{
	// For an open-problem with N conductors, returns a NxN matrix which has: 
	// - a main diagonal of 1s, representing a prescribed voltage of 1 in the n-th conductor.
	// - the off-diagonal terms are the voltages at the other conductors when they are assumed to be floating.

	// For closed and semi-open problems with N conductors, returns a N-1 x N-1 matrix and assumes that conductor 0 has 
	// alway a prescribed potential of zero.

	const auto conductors{ model_.getMaterials().buildNameToAttrMapFor<PEC>() };
	const int N = conductors.size();
	if (N == 1) {
		mfem::DenseMatrix res(1,1);
		res = 1.0;
		return res;
	}

	
	auto C{ getCMatrix(model_, opts_, ignoreDielectrics) };
	C.Symmetrize();
	
	switch (model_.determineOpenness()) {
	case Model::OpennessType::closed:
	{
		return floatingPotentialsForClosedCases(C);
	}
	case Model::OpennessType::open:
	{
		auto Cib{ getCapacitancesWithOpenBoundary(model_, opts_, ignoreDielectrics) };

		return floatingPotentialsForOpenCases(C, Cib);
	}
	default:
		throw std::runtime_error("Floating potentials not implemented for this kind of opennnes.");
	}

}

}