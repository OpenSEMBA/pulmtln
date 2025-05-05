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

mfem::DenseMatrix getCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false,
	bool includeGround = false)
{
	// PUL capacitance matrix as defined in:
	//   Clayton Paul's book: Analysis of Multiconductor Transmission Lines
	// If ground is not included, contains N-1 x N-1 entries for a problem of N conductors.
	// if ground is included, then N x N.

	// Preconditions. 
	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	if (conductors.size() == 0) {
		throw std::runtime_error(
			"At least one conductor needed."
		);
	}
	if (conductors.size() == 1 && !includeGround) {
		throw std::runtime_error(
			"The number of conductors must be at least 2 if not including ground."
		);
	}

	int CSize;
	if (includeGround) {
		CSize = (int)conductors.size();
	}
	else {
		CSize = (int)conductors.size() - 1;
	}
	mfem::DenseMatrix C(CSize);

	// Solves a electrostatic problem for each conductor besides the
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		auto condId = Materials::getNumberContainedInName(nameI);
		if (condId == model.getGroundConductorId() 
			&& !includeGround) {
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
			if (Materials::getNumberContainedInName(nameJ) == model.getGroundConductorId() 
				&& !includeGround) {
				continue;
			}
			
			// C_ij = Q_i / V_j, V_j is always 1.0
			C(row, col) = s.chargeInBoundary(conductors.at(nameJ)); 
			
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

DenseMatrix symmetrizeMatrix(const DenseMatrix& m)
{
	assert(m.NumRows() == m.NumCols());
	
	DenseMatrix res(m.NumRows(), m.NumCols());
	for (int i{ 0 }; i < m.NumRows(); i++) {
		for (int j{ i }; j < m.NumCols(); j++) {
			if (i == j) {
				res(i, j) = m(i, j);
			}
			auto v{ (m(i,j) + m(j,i)) / 2.0 };
			res(i, j) = v;
			res(j, i) = v;
		}
	}

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
		res.C = symmetrizeMatrix(res.C);
		res.L = symmetrizeMatrix(res.L);
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


mfem::DenseMatrix Driver::getFloatingPotentials(const bool ignoreDielectrics) const
{
	// For an open-problem with N conductors, returns a NxN matrix which has: 
	// - a main diagonal of 1s, representing a prescribed voltage of 1 in the n-th conductor.
	// - the off-diagonal terms are the voltages at the other conductors when they are assumed to be floating.

	// For closed and semi-open problems with N conductors, returns a N-1 x N-1 matrix and assumes that conductor 0 has 
	// alway a prescribed potential of zero.

	const auto conductors{ model_.getMaterials().buildNameToAttrMapFor<PEC>() };
	if (conductors.size() == 1) {
		mfem::DenseMatrix res(1,1);
		res = 1.0;
		return res;
	}


	bool includeGround;
	auto openness{ model_.determineOpenness() };
	if (openness == Model::OpennessType::closed) {
		includeGround = false;
	}
	else {
		includeGround = true;
	}

	mfem::DenseMatrix C = symmetrizeMatrix(
		getCMatrix(model_, opts_, ignoreDielectrics, includeGround));
	const int N = C.NumRows();
	
	mfem::DenseMatrix res(N,N);
	
	for (int i{ 0 }; i < N; ++i) {
		// Forms system of equations to determine floating potentials. 
		// Q1 = C11*V1 + C21*V2 + ....
		// For prescribed V_1 = 1.0 we have C V = Q
		//    [ C ] [1.0, V_2, ...]^T = [Q_1, 0.0, ...]
		// 
		// which can be converted to A x = b with unknowns x = [Q_1, V_2, ...]^T

		mfem::DenseMatrix A{ C };
		for (int k{ 0 }; k < N; k++) {
			if (i == k) {
				A(k, i) = -1.0;
			}
			else {
				A(k, i) = 0.0;
			}
		}

		mfem::Vector b(N);
		for (int k{ 0 }; k < N; k++) {
			b(k) = -C(k,i);
		}

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