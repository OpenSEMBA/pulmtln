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

mfem::DenseMatrix getGeneralizedCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false)
{
	// PUL generalized capacitance matrix as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// For a problem with N conductors.
	// Only valid for open problems.
	// Contains N x N entries.

	// Preconditions. 
	if (model.determineOpenness() == Model::OpennessType::closed) {
		throw std::runtime_error("Generalized capacitance only implemented for open/semiopen problems");
	}

	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	int CSize = (int)conductors.size();
	mfem::DenseMatrix C(CSize);

	// Solves a electrostatic problem for each conductor.
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		auto parameters{ buildSolverParameters(model, ignoreDielectrics) };
		parameters.dirichletBoundaries[bdrAttI] = 1.0;

		Mesh mesh{ *model.getMesh() };
		ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
		s.Solve();

		//int openBoundaryTag = *parameters.openBoundaries.begin();
		//auto Vb = s.averagePotentialInBoundary(openBoundaryTag);
		//auto Vd = 1.0 - Vb;
		
		// Fills row
		int col{ 0 };
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			// C_ji = Q_j / V_i
			double Qj = s.chargeInBoundary(conductors.at(nameJ));
			C(col, row) = Qj;
			col++;
		}

		exportFieldSolutions(opts, s, nameI, ignoreDielectrics);
		row++;
	}

	return C;
}

mfem::DenseMatrix getCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false)
{
	// PUL capacitance matrix as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// Contains N-1 x N-1 entries for a problem of N conductors.

	// Preconditions. 
	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	const auto openness{ model.determineOpenness() };
	if (conductors.size() == 1) {
		throw std::runtime_error("The number of conductors must be at least 2.");
	}

	mfem::Vector Cib(0);
	if (openness == Model::OpennessType::open) {
		Cib.SetSize(conductors.size());
		Cib = getCapacitancesWithOpenBoundary(model, opts, ignoreDielectrics);
	}

	int CSize = (int)conductors.size() - 1;
	mfem::DenseMatrix C(CSize);

	// Solves a electrostatic problem for each conductor besides the
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		int condI = Materials::getNumberContainedInName(nameI);
		if (condI == model.getGroundConductorId()) {
			continue;
		}

		auto parameters{ buildSolverParameters(model, ignoreDielectrics) };
		parameters.dirichletBoundaries[bdrAttI] = 1.0;

		Mesh mesh{ *model.getMesh() };

		ElectrostaticSolver s(mesh, parameters, opts.solverOptions);
		s.Solve();

		// Fills row
		int col{ 0 };
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			int condJ = Materials::getNumberContainedInName(nameJ);
			if (condJ == model.getGroundConductorId()) {
				continue;
			}
			
			// C_ij = Q_j / V_i. V_i is always 1.0
			double Qj = s.chargeInBoundary(conductors.at(nameJ));
			C(row, col) = Qj;
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
			opts_);
		}

	res.domainTree = DomainTree{ idToDomain };

	return res;
}

mfem::DenseMatrix Driver::getFloatingPotentials(const bool ignoreDielectrics) const
{
	// For an open-problem with N conductors, returns a NxN matrix which has: 
	// - a main diagonal of 1s, representing a prescribed voltage of 1 in the n-th conductor.
	// - the off-diagonal terms are the voltages at the other conductors when they are assumed to be floating.

	// For closed and semi-open problems with N conductors, returns a N-1 x N-1 matrix and assumes that conductor 0 has 
	// alway a prescribed potential of zero.

	if (model_.getMaterials().buildNameToAttrMapFor<PEC>().size() == 1) {
		mfem::DenseMatrix res(1,1);
		res = 1.0;
		return res;
	}

	// Determine C matrix.
	mfem::DenseMatrix C;
	switch (model_.determineOpenness()) {
	case Model::OpennessType::closed:
		C = getCMatrix(model_, opts_, ignoreDielectrics);
		break;
	case Model::OpennessType::open:
		C = getGeneralizedCMatrix(model_, opts_, ignoreDielectrics);
		break;
	default:
		throw std::runtime_error(
			"Floating potentials not implemented for this kind of opennnes.");
	}
	C.Symmetrize();
	
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

}