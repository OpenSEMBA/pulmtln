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

mfem::DenseMatrix solveCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics = false)
{
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

DenseMatrix solveLMatrix(const Model& model, const DriverOptions& opts)
{
	// Inductance matrix can be computed from the 
	// capacitance obtained ignoring dielectrics as
	//          L = mu0 * eps0 * C^{-1}
	auto res{ solveCMatrix(model, opts, true) };
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

	res.C = solveCMatrix(m, opts);
	res.C *= EPSILON0_SI;

	res.L = solveLMatrix(m, opts);
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

}