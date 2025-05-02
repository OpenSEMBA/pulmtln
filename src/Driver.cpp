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

	if (model.determineOpenness() == Model::OpennessType::open) {
		parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMapFor<PEC>(), -1.0);
	}
	else {
		parameters.dirichletBoundaries = buildAttrToValueMap(mats.buildNameToAttrMapFor<PEC>(), 0.0);
	}

	return parameters;
}

mfem::DenseMatrix getCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false)
{
	// PUL capacitance matrix as defined in:
	//   Clayton Paul's book: Analysis of Multiconductor Transmission Lines
	// Contains N-1 x N-1 entries for a problem of N conductors.

	// Preconditions. 
	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	if (conductors.size() < 2) {
		throw std::runtime_error(
			"The number of conductors must be greater than 2."
		);
	}

	// Solves a electrostatic problem for each conductor besides the
	// ground conductor
	mfem::DenseMatrix C((int)conductors.size() - 1);
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		if (Materials::getNumberContainedInName(nameI) == model.getGroundConductorId()) {
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
			if (Materials::getNumberContainedInName(nameJ) == model.getGroundConductorId()) {
				continue;
			}
			auto charge{ s.chargeInBoundary(conductors.at(nameJ)) };
			if (model.determineOpenness() == Model::OpennessType::open) {
				C(row, col) = charge / 2.0;
			}
			else {
				C(row, col) = charge;
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

	mfem::DenseMatrix C{ getCMatrix(model_, opts_, ignoreDielectrics) };
	const auto N = conductors.size();
	
	switch (model_.determineOpenness()) {
	case Model::OpennessType::closed:
		{
			mfem::DenseMatrix res(N - 1, N - 1);

			for (int i{ 0 }; i < res.NumRows(); ++i) {
				// Forms system of equations to determine floating potentials. 
				// 
				// For prescribed V_1 = 1.0 we have C V = Q
				//    [ C ] [1.0, V_2, V_3, ...]^T = [Q_1, 0.0, 0.0, ...]
				// 
				// which converts can be converted to A x = b as:
				//    [ -1.0, C_{1,2}, ... ] [] = [-C_{1,1}, 0.0, 0.0, ...] 

				mfem::DenseMatrix A{ C };
				A(i, i) = -1.0;
				mfem::Vector b(N - 1), x(N - 1);
				b = 0.0;
				b(i) = -C(i, i);

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
	default:
		{
			throw std::runtime_error("Not implemented.");
		}
	}
}

}