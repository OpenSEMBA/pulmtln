#include "Driver.h"

#include "ElectrostaticSolver.h"
#include "Parser.h"
#include "multipolarExpansion.h"

using namespace mfem;

namespace pulmtln {

const std::string INNER_VACUUM_DEFAULT_NAME = "Vacuum_0";

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

SolverInputs Driver::buildSolverInputsFromModel(
	const Model& model,
	bool ignoreDielectrics)
{
	const Materials& mats = model.getMaterials();

	SolverInputs parameters;

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

DenseMatrix Driver::getCMatrix(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics,
	bool generalized)
{
	// PUL capacitance matrix for a N conductors system as defined in:
	// "Clayton Paul's book: Analysis of Multiconductor Transmission Lines"
	// - Standard C contains N-1 x N-1 entries
	// - Generalized C contains N x N entries.

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

	const auto baseParameters{ buildSolverInputsFromModel(model, ignoreDielectrics) };
	Mesh mesh{ *model.getMesh() };
	ElectrostaticSolver s(mesh, baseParameters, opts.solverOptions);
	
	// Solves a electrostatic problem for each conductor besides the
	int row{ 0 };
	for (const auto& [nameI, bdrAttI] : conductors) {
		int condI = Materials::getMaterialIdFromName(nameI);
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
			int condJ = Materials::getMaterialIdFromName(nameJ);
			if (condJ == model.getGroundConductorId() && !generalized) {
				continue;
			}
			
			// C_ij = Q_j / V_i. V_i is always 1.0
			double Qj = s.getChargeInBoundary(conductors.at(nameJ));
			C(row, col) = Qj;
			col++;
		}

		exportFieldSolutions(opts, s, nameI, ignoreDielectrics);
		row++;
	}

	return C;
}

DenseMatrix Driver::getLMatrix(const Model& model, const DriverOptions& opts)
{
	// PUL inductance matrix as defined in:
	//   Clayton Paul's book: Analysis of Multiconductor Transmission Lines
	// Contains N-1 x N-1 entries for a problem of N conductors.
	// Inductance matrix can be computed from the 
	// capacitance obtained ignoring dielectrics as
	//          L = mu0 * eps0 * C^{-1}
	auto res{ Driver::getCMatrix(model, opts, true) };
	res.Invert();
	res *= MU0_NATURAL * EPSILON0_NATURAL;
	return res;
}

PULParameters buildPULParametersForModel(const Model& m, const DriverOptions& opts)
{
	PULParameters res;

	res.C = Driver::getCMatrix(m, opts);
	res.C *= EPSILON0_SI;

	res.L = Driver::getLMatrix(m, opts);
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

DenseMatrix Driver::getFloatingPotentialsMatrix(
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

std::list<std::string> listMaterialsInInnerRegion(
	const Model& m,
	bool includeConductors = true)
{
	std::list<std::string> res;
	res.push_back(INNER_VACUUM_DEFAULT_NAME);
	for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<Dielectric>()) {
		if (name.find("Vacuum_") != std::string::npos) {
			continue;
		}
		res.push_back(name);
	}
	if (includeConductors) {
		for (auto [name, tag] : m.getMaterials().buildNameToAttrMapFor<PEC>()) {
			res.push_back(name);
		}
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

double getInnerRegionAveragePotential(
	const Model& m, 
	const ElectrostaticSolver& s,
	bool includeConductors)
{

	double totalPotential = 0.0;
	double totalArea = 0.0;
	
	auto innerRegionMaterials = listMaterialsInInnerRegion(m, includeConductors);
	auto materials = m.getMaterials().buildNameToAttrMap();
	
	for (const auto& name: innerRegionMaterials ) {
		auto tag = materials.at(name);
		double area = m.getAreaOfMaterial(name);
		if (m.getMaterials().isDomainMaterial(name)) {
			totalPotential += s.getAveragePotentialInDomain(tag) * area;
		} else {
			totalPotential += s.getAveragePotentialInBoundary(tag) * area;
		}
		totalArea += area;
	}

	return totalPotential / totalArea;
}

std::map<MaterialId, FieldReconstruction> getFieldParameters(
	const Model& model,
	const DriverOptions& opts,
	bool ignoreDielectrics)
{
	std::map<MaterialId, FieldReconstruction> res;

	auto fp = Driver::getFloatingPotentialsMatrix(model, opts, ignoreDielectrics);
	const auto baseParameters{ Driver::buildSolverInputsFromModel(model, ignoreDielectrics) };
	Mesh mesh{ *model.getMesh() };
	ElectrostaticSolver s(mesh, baseParameters, opts.solverOptions);

	const auto conductors{ model.getMaterials().buildNameToAttrMapFor<PEC>() };
	for (const auto& [nameI, bdrAttI] : conductors) {
		auto condI = Materials::getMaterialIdFromName(nameI);

		auto dbcs = baseParameters.dirichletBoundaries;
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			auto condJ = Materials::getMaterialIdFromName(nameJ);
			dbcs[bdrAttJ] = fp(condI, condJ);
		}
		s.setDirichletConditions(dbcs);
		s.Solve();

		exportFieldSolutions(opts, s, 
			nameI + "_prescribed_and_others_floating", ignoreDielectrics);

		res[condI].innerRegionAveragePotential = 
			getInnerRegionAveragePotential(model, s, true);
		res[condI].expansionCenter = s.getCenterOfCharge();
		res[condI].ab = s.getMultipolarCoefficients(opts.multipolarExpansionOrder);
		for (const auto& [nameJ, bdrAttJ] : conductors) {
			auto condJ = Materials::getMaterialIdFromName(nameJ);
			res[condI].conductorPotentials[condJ] = fp(condI, condJ);
		}
	}

	return res;
}

InCellPotentials Driver::getInCellPotentials() const
{
	InCellPotentials res;

	if (model_.determineOpenness() != Model::Openness::open) {
		throw std::runtime_error("In cell parameters can only be computed for open problems.");
	}

	const double innerRegionArea{ getInnerRegionArea(model_) };
	res.innerRegionRadius = std::sqrt(innerRegionArea / M_PI);

	res.electric = getFieldParameters(model_, opts_, false);
	res.magnetic = getFieldParameters(model_, opts_, true);


	return res;
}


double Driver::getInCellCapacitanceUsingInnerRegion(
	const InCellPotentials& potential, int i, int j)
{
	double Qj = potential.electric.at(j).ab[0].first;
	double avVj = potential.electric.at(j).innerRegionAveragePotential;
	double ViWhenPrescribedVj = potential.electric.at(j).conductorPotentials.at(i);
	avVj = -avVj + ViWhenPrescribedVj;
	return Qj / avVj * EPSILON0_SI;
}

double Driver::getInCellInductanceUsingInnerRegion(
	const InCellPotentials& potential, int i, int j)
{
	double Ij = potential.magnetic.at(j).ab[0].first;
	double avAj = potential.magnetic.at(j).innerRegionAveragePotential;
	double AiWhenPrescribedAj = potential.magnetic.at(j).conductorPotentials.at(i);
	avAj = -avAj + AiWhenPrescribedAj;
	return avAj / Ij * MU0_SI;
}


}