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
	bool ignoreDielectrics = false,
	bool generalized = false)
{
	// PUL capacitance matrix as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// Standard C contains N-1 x N-1 entries for a problem of N conductors.

	// Preconditions. 
	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	const auto openness{ model.determineOpenness() };
	if (conductors.size() == 1 && openness == Model::Openness::closed) {
		throw std::runtime_error(
			"The number of conductors must be at least 2 for closed problems.");
	}

	int CSize;
	if (generalized) {
		CSize = (int) conductors.size();
	}
	else {
		CSize = (int) conductors.size() - 1;
	}
	mfem::DenseMatrix C(CSize);

	const auto baseParameters{ buildSolverParameters(model, ignoreDielectrics) };
	Mesh mesh{ *model.getMesh() };
	ElectrostaticSolver s(mesh, baseParameters, opts.solverOptions);
	
	// Solves a electrostatic problem for each conductor besides the
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		int condI = Materials::getNumberContainedInName(nameI);
		if (condI == model.getGroundConductorId() && !generalized) {
			continue;
		}
		
		auto dbcs = baseParameters.dirichletBoundaries;
		dbcs[bdrAttI] = 1.0;
		s.setDirichletConditions(dbcs);
		s.Solve();

		// Fills row
		int col{ 0 };
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			int condJ = Materials::getNumberContainedInName(nameJ);
			if (condJ == model.getGroundConductorId() && !generalized) {
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

mfem::DenseMatrix getFloatingPotentialsMatrix(
	const Model& model,
	const DriverOptions& opts,
	const bool ignoreDielectrics)
{
	// For an open-problem with N conductors, returns a NxN matrix which has: 
	// - a main diagonal of 1s, representing a prescribed voltage of 1 in the n-th conductor.
	// - the off-diagonal terms are the voltages at the other conductors when they are assumed to be floating.

	// For closed and semi-open problems with N conductors, returns a N-1 x N-1 matrix and assumes that conductor 0 has 
	// alway a prescribed potential of zero.

	if (model.getMaterials().buildNameToAttrMapFor<PEC>().size() == 1) {
		mfem::DenseMatrix res(1,1);
		res = 1.0;
		return res;
	}

	// Determine C matrix.
	bool useGeneralizedCMatrix;
	switch (model.determineOpenness()) {
	case Model::Openness::closed:
		useGeneralizedCMatrix = false;
		break;
	case Model::Openness::open:
		useGeneralizedCMatrix = true;
		break;
	default:
		throw std::runtime_error(
			"Floating potentials not implemented for this kind of openness.");
	}
	mfem::DenseMatrix C{ getCMatrix(model, opts, ignoreDielectrics, useGeneralizedCMatrix) };
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

FloatingPotentials Driver::getFloatingPotentials() const
{
	FloatingPotentials res;

	res.electric = getFloatingPotentialsMatrix(model_, opts_, false);
	res.magnetic = getFloatingPotentialsMatrix(model_, opts_, true);

	return res;
}

std::list<std::string> listMaterialsInInnerRegion(const Model& m)
{
	std::list<std::string> res;
	res.push_back("Vacuum_0");
	for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<Dielectric>()) {
		res.push_back(name);
	}
	for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<PEC>()) {
		res.push_back(name);
	}
	return res;
}

double getInnerRegionArea(const Model& m)
{
	double res = 0.0;
	for (auto name : listMaterialsInInnerRegion(m)) {
		res += m.getAreaOfMaterial(name);
	}
	return res;
}

double getInnerRegionAveragePotential(const Model& m, const ElectrostaticSolver& s)
{
	double res = 0.0;

	// TODO

	return res;
}

std::map<std::string, InCellParameters::FieldParameters> getFieldParameters(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics)
{
	std::map<std::string, InCellParameters::FieldParameters> res;

	auto fp = getFloatingPotentialsMatrix(model, opts, ignoreDielectrics);
	const auto baseParameters{ buildSolverParameters(model, ignoreDielectrics) };
	Mesh mesh{ *model.getMesh() };
	ElectrostaticSolver s(mesh, baseParameters, opts.solverOptions);

	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	for (const auto& [nameI, bdrAttI] : conductors) {
		auto condI = Materials::getNumberContainedInName(nameI);

		auto dbcs = baseParameters.dirichletBoundaries;
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			auto condJ = Materials::getNumberContainedInName(nameJ);
			dbcs[bdrAttJ] = fp(condI, condJ);
		}
		s.setDirichletConditions(dbcs);
		s.Solve();

		res[nameI].innerRegionAveragePotential = getInnerRegionAveragePotential(model, s);
		res[nameI].expansionCenter = s.getCenterOfCharge();
		res[nameI].ab = s.getMultipolarCoefficients(opts.multipolarExpansionOrder);
	}

	return res;
}

InCellParameters Driver::getInCellParameters() const
{
	InCellParameters res;

	if (model_.determineOpenness() != Model::Openness::open) {
		throw std::runtime_error("In cell paramters can only be computed for open problems.");
	}

	const double innerRegionArea{ getInnerRegionArea(model_) };
	res.innerRegionRadius = std::sqrt(innerRegionArea / M_PI);

	res.electric = getFieldParameters(model_, opts_, false);
	res.magnetic = getFieldParameters(model_, opts_, true);

	return res;
}

}