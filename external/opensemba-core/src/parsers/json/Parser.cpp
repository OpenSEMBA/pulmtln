#include "Parser.h"

#include "parsers/stl/Parser.h"

#include "core/math/function/Gaussian.h"
#include "core/math/function/BandLimited.h"
#include "core/model/Model.h"
#include "core/geometry/element/Line2.h"
#include "core/geometry/element/Triangle3.h"
#include "core/geometry/element/Quadrilateral4.h"
#include "core/geometry/element/Tetrahedron4.h"
#include "core/geometry/element/Hexahedron8.h"
#include "core/physicalModel/physicalModels.h"
#include "core/source/port/WaveguideRectangular.h"
#include "core/source/port/TEMCoaxial.h"
#include "core/source/PlaneWave.h"
#include "core/source/Generator.h"
#include "core/source/OnLine.h"
#include "core/outputRequest/BulkCurrent.h"
#include "core/outputRequest/FarField.h"
#include "core/outputRequest/OnPoint.h"
#include "core/outputRequest/OnLine.h"
#include "core/outputRequest/OnSurface.h"
#include "core/outputRequest/OnLayer.h"

using namespace semba;
using namespace geometry;
using namespace math;

using json = nlohmann::json;


namespace semba::parsers::JSON {

CVecI3 strToCVecI3(std::string str);
CVecR3 strToCVecR3(std::string str);
Constants::CartesianAxis strToCartesianAxis(std::string);
std::pair<CVecR3, CVecR3> strToBox(const std::string& str);
LocalAxis strToLocalAxes(const std::string& str);
const ElemR* boxToElemGroup(mesh::Unstructured& mesh, const std::string& line);

LayerGroup readLayers(const json&);
CoordR3Group readCoordinates(const json&);
ElemRGroup readElements(const PMGroup&, LayerGroup&, CoordR3Group&, const json&, const std::string& folder);
ElemRGroup readElementsFromFile(const PMGroup&, LayerGroup&, CoordR3Group&, const json&, const std::string& folder);
ElemRGroup readElementsFromSTLFile(const PMGroup&, LayerGroup&, CoordR3Group&, const json&, const std::string& folder);

std::unique_ptr<physicalModel::surface::Multilayer> readMultilayerSurface(const json& layers);
std::unique_ptr<PM> readPhysicalModel(const json& material);
physicalModel::PhysicalModel::Type strToMaterialType(std::string label);
physicalModel::multiport::Multiport::Type strToMultiportType(std::string label);
physicalModel::volume::Anisotropic::Model strToAnisotropicModel(std::string label);

std::unique_ptr<source::PlaneWave> readPlanewave(mesh::Unstructured& mesh, const json&);
std::unique_ptr<source::port::Waveguide> readPortWaveguide(mesh::Unstructured& mesh, const json&);
std::unique_ptr<source::port::TEM> readPortTEM(mesh::Unstructured& mesh, const json&);
std::unique_ptr<source::Generator> readGenerator(mesh::Unstructured& mesh, const json&);
std::unique_ptr<source::OnLine> readSourceOnLine(mesh::Unstructured& mesh, const json&);
std::unique_ptr<source::Magnitude::Magnitude> readMagnitude(const json&);
source::Generator::Type strToGeneratorType(std::string label);
source::Generator::Hardness strToGeneratorHardness(std::string str);
source::OnLine::Type strToNodalType(std::string label);
source::OnLine::Hardness strToNodalHardness(std::string label);
source::port::TEM::ExcitationMode strToTEMMode(std::string);
source::port::Waveguide::ExcitationMode strToWaveguideMode(std::string);

void checkVersionCompatibility(const std::string& version)
{
    if (version != std::string(VERSION)) {
        throw std::logic_error("File version " + version + " is not supported.");
    }
}

std::vector<ElemId> readElemIds(const json& j)
{
    std::vector<ElemId> r;
    for (const auto& eId : j) {
        r.push_back(ElemId{ std::size_t(eId.get<int>()) });
    }
    return r;
}


template<typename T>
element::Group<ElemR> readElemStrAs(
    const physicalModel::Group<>& mG,
    const LayerGroup& lG,
    const CoordR3Group& cG,
    const json& e) {
    element::Group<ElemR> res;

    for (auto it = e.begin(); it != e.end(); ++it) {

        ElemId elemId;
        MatId matId;
        LayerId layerId;
        std::vector<CoordId> vId;

        std::stringstream ss(it->get<std::string>());
        ss >> elemId >> matId >> layerId;
        vId.resize(T::sizeOfCoordinates);
        for (std::size_t j = 0; j < T::sizeOfCoordinates; j++) {
            ss >> vId[j];
        }

        const Layer* layerPtr;
        const physicalModel::PhysicalModel* matPtr;
        std::vector<const CoordR3*> vPtr;

        if (matId != MatId(0)) {
            matPtr = mG.getId(matId);
        }
        else {
            matPtr = nullptr;
        }
        if (layerId != LayerId(0)) {
            layerPtr = lG.getId(layerId);
        }
        else {
            layerPtr = nullptr;
        }
        vPtr.resize(vId.size(), nullptr);
        for (size_t i = 0; i < vId.size(); ++i) {
            vPtr[i] = cG.getId(vId[i]);
        }

        res.add(std::make_unique<T>(T(elemId, vPtr.data(), layerPtr, matPtr)));
    }

    return res;
}

double getProgressionStepByTotalNumber(const json& j, const std::string& jsonKey) {
	if (j.find("total" + jsonKey) != j.end()) {
		int total = j.at("total" + jsonKey).get<int>();

		// TODO: Yet to be clarify differences between UI / JSON
		auto step = j.at("final" + jsonKey).get<double>() - j.at("initial" + jsonKey).get<double>();

		if (total == 1 && step == 0.0) {
			// Need a silly value here for solver to run
			return 0.1;
		}

		if (total > 1) {
			step = step / (total - 1);
		}

		return step;

	}

	if (j.find("step" + jsonKey) != j.end()) {
		return j.at("step" + jsonKey).get<double>();
	}

	// TODO: Check better approach
	std::string toLower = jsonKey;
	for (auto& c : toLower) {
		c = tolower(c);
	}

	return j.at(toLower + "Step").get<double>();
}

std::vector<const CoordR3*> addAngGetCoordView(
    CoordR3Group& cG,
    const std::vector<CVecR3>& positions,
    const size_t& numberOfCoordinates
) {
    std::vector<const CoordR3*> coords;

    for (std::size_t i = 0; i < numberOfCoordinates; i++) {
        coords.push_back(
            cG.addAndAssignId(
                std::make_unique<CoordR3>(CoordId(), positions[i])
            )->get()
        );
    }

    return coords;
}

UnstructuredProblemDescription Parser::read() const 
{
	std::ifstream stream(this->filename);
	if (!stream.is_open()) {
		throw std::logic_error("Can not open file: " + this->filename);
	}

	json j;
	try {
		stream >> j;
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl;
	}

	checkVersionCompatibility(j.at("_version").get<std::string>());

    UnstructuredProblemDescription res;
    res.project = filename;
	res.analysis = readAnalysis(j);
	res.grids = readGrids(j);

    auto materialsGroup{ readMaterials(j.at("model")) };
	auto mesh = readUnstructuredMesh(materialsGroup, j.at("model"), this->filename.getFolder());

	res.sources = readSources(*mesh, j);
	res.outputRequests = readProbes(*mesh, j);

	readBoundary(*mesh, j, materialsGroup, res.grids);

    res.model = UnstructuredModel{*mesh, materialsGroup};

    postReadOperations(res);

	return res;
}

Parser::Parser(const std::string& fn) : semba::parsers::Parser(fn) 
{}

physicalModel::Bound::Type strToBoundType(const std::string& boundType) {
    if (boundType == "PEC") {
        return physicalModel::Bound::Type::pec;
    }

    if (boundType == "PMC") {
        return physicalModel::Bound::Type::pmc;
    }

    if (boundType == "PML") {
        return physicalModel::Bound::Type::pml;
    }

    if (boundType == "Periodic") {
        return physicalModel::Bound::Type::periodic;
    }

    if (boundType == "MUR1") {
        return physicalModel::Bound::Type::mur1;
    }

    if (boundType == "MUR2") {
        return physicalModel::Bound::Type::mur2;
    }

    throw std::logic_error("Unrecognized value in Bound ctor.");
}

void readBoundary(mesh::Unstructured& mesh, const json& j, PMGroup& physicalModelGroup, const Grid3& grid)  
{
    if (j.find("boundary") == j.end()) {
        return;
    }

    json lower = j.at("boundary").at("lower");
    json upper = j.at("boundary").at("upper");

    if (lower.size() != 3 || upper.size() != 3) {
        throw std::logic_error("Unexpected amount of materials for boundary layers specified. Three layers are expected");
    }

    const auto& box = grid.getFullDomainBoundingBox();

    for (const auto& bound : { Constants::CartesianBound::L, Constants::CartesianBound::U }) {
        const auto& boundaryResource = bound == Constants::CartesianBound::L ? lower : upper;

        for (const auto& axis : { Constants::CartesianAxis::x,Constants::CartesianAxis::y,Constants::CartesianAxis::z }) {
            physicalModel::Id id(0);

            for (const auto& boundI : physicalModelGroup.getOf<physicalModel::Bound>()) {
                if (boundI->getType() == strToBoundType(boundaryResource[axis].get<std::string>())) {
                    id = boundI->getId();
                    break;
                }
            }

            if (id == physicalModel::Id(0)) {
                auto it = physicalModelGroup.addAndAssignId(
                    std::make_unique<physicalModel::Bound>(physicalModel::Id(), strToBoundType(boundaryResource[axis].get<std::string>()))
                );

                id = (*it)->getId();
            }

            const auto& boxBoundary = box.getPosOfBound(axis, bound);
            std::array<const CoordR3*, 4> boundCoordinatesArray{};
            for (std::size_t boxBoundaryIndex = 0; boxBoundaryIndex < boxBoundary.size(); ++boxBoundaryIndex) {
                const auto& coordIt = mesh.coords().addAndAssignId(
                    std::make_unique<CoordR3>(boxBoundary[boxBoundaryIndex])
                );

                boundCoordinatesArray[boxBoundaryIndex] = coordIt->get();
            }

            mesh.elems().addAndAssignId(
                std::make_unique<QuaR4>(
                    ElemId(),
                    boundCoordinatesArray,
                    nullptr,
                    physicalModelGroup.getId(id)
                )
            );
        }
    }
}

json readAnalysis(const json& j)
{
    const std::string key{ "analysis" };
    return j.find(key) == j.end() ? json(): j.at(key).get<json>();
}

std::unique_ptr<mesh::Unstructured> readUnstructuredMesh(const PMGroup& physicalModels, const json& j, const std::string& folder)
{
    LayerGroup layers = readLayers(j);
	CoordR3Group coords = readCoordinates(j);
	return std::make_unique<mesh::Unstructured>(
		coords,
		readElements(physicalModels, layers, coords, j, folder),
		layers
	);
}

SourceGroup readSources(mesh::Unstructured& mesh, const json& j)
{
    auto sources = j.find("sources");
    
    SourceGroup res;
    if (sources == j.end()) {
        return res;
    }

    for (auto const& it: sources->get<json>()) {
        std::string sourceType = it.at("sourceType").get<std::string>();
        if (sourceType.compare("planewave") == 0) {
            res.addAndAssignId(readPlanewave(mesh, it));
            continue;
        } 
        
        if (sourceType.compare("generator") == 0) {
            res.addAndAssignId(readGenerator(mesh, it));
            continue;
        } 
        
        if (sourceType.compare("sourceOnLine") == 0) {
            res.addAndAssignId(readSourceOnLine(mesh, it));
            continue;
        } 
        
        if (sourceType.compare("waveguidePort") == 0) {
            res.addAndAssignId(readPortWaveguide(mesh, it));
            continue;
        } 
        
        if (sourceType.compare("temPort") == 0) {
            res.addAndAssignId(readPortTEM(mesh, it));
            continue;
        } 

        throw std::logic_error("Unrecognized source type: " + sourceType);
        
    }
    return res;
}

PMGroup readMaterials(const json& j)
{
    PMGroup res;
    for (auto const& mat: j.at("materials")) {
        res.addAndAssignId(readPhysicalModel( mat ));
    }
    return res;
}

std::unique_ptr<physicalModel::PhysicalModel> readPhysicalModel(const json& j)
{
    typedef physicalModel::PhysicalModel PM;
    
	MatId id(0);
	if (j.find("materialId") != j.end()) {
		id = MatId(j.at("materialId").get<int>());
	}

	std::string name;
	if (j.find("name") != j.end()) {
		name = j.at("name").get<std::string>();
	}

    auto type{ strToMaterialType( j.at("materialType").get<std::string>() ) };
    switch (type) {
    case PM::Type::PEC:
        return std::make_unique<physicalModel::PEC>(id, name);
    case PM::Type::PMC:
        return std::make_unique<physicalModel::PMC>(id, name);
    case PM::Type::SMA:
        return std::make_unique<physicalModel::SMA>(id, name);
    case PM::Type::vacuum:
        return std::make_unique<physicalModel::Vacuum>(id, name);
    case PM::Type::PML:
        return std::make_unique<physicalModel::volume::PML>(id, name, 
			strToLocalAxes(j.at("localAxes").get<std::string>()));

    case PM::Type::classic:
        return std::make_unique<physicalModel::volume::Classic>(
                id, 
                name,
                j.at("permittivity").get<double>(),
                j.at("permeability").get<double>(),
                j.at("electricConductivity").get<double>(),
                j.at("magneticConductivity").get<double>()
            );

    case PM::Type::elecDispersive:
        return std::make_unique <physicalModel::volume::Dispersive>(
                id, 
                name,
                j.at("filename").get<std::string>()
            );

    case PM::Type::wire:
    {
        std::string wireType = j.at("wireType").get<std::string>();
        if (wireType.compare("Dispersive") == 0) {
            return std::make_unique < physicalModel::wire::Wire>(id, name,
                    j.at("radius").get<double>(),
                    j.at("filename").get<std::string>());
        } else if(wireType.compare("SeriesParallel") == 0) {
            return std::make_unique < physicalModel::wire::Wire>(id, name,
                    j.at("radius").get<double>(),
                    j.at("resistance").get<double>(),
                    j.at("inductance").get<double>(),
                    j.at("capacitance").get<double>(),
                    j.at("parallelResistance").get<double>(),
                    j.at("parallelInductance").get<double>(),
                    j.at("parallelCapacitance").get<double>());
        } else if(wireType.compare("Standard") == 0) {
            return std::make_unique < physicalModel::wire::Wire>(id, name,
                    j.at("radius").get<double>(),
                    j.at("resistance").get<double>(),
                    j.at("inductance").get<double>());
        } else {
            throw std::logic_error("Unrecognized wire type" + wireType);
        }
    }

    case PM::Type::anisotropic:
    {
        std::string str = j.at("anisotropicModel").get<std::string>();
        if (str.compare("Crystal")==0) {
            return std::make_unique < physicalModel::volume::AnisotropicCrystal>(id, name,
                    strToLocalAxes(j.at("localAxes").get<std::string>()),
                    strToCVecR3(
                            j.at("relativePermittiviy").get<std::string>()),
                    j.at("crystalRelativePermeability").get<double>());
        } else if (str.compare("Ferrite")==0) {
            return std::make_unique < physicalModel::volume::AnisotropicFerrite>(id, name,
                    strToLocalAxes(j.at("localAxes").get<std::string>()),
                    j.at("kappa").get<double>(),
                    j.at("ferriteRelativePermeability").get<double>(),
                    j.at("ferriteRelativePermittivity").get<double>());
        } else {
            throw std::logic_error("Unrecognized Anisotropic Model: " + str);
        }
    }

    case PM::Type::isotropicsibc:
    {
        std::string sibcType = j.at("surfaceType").get<std::string>();
        if (sibcType.compare("File")==0) {
            return std::make_unique < physicalModel::surface::SIBCFile>(id, name,
                    j.at("filename").get<std::string>() );
        } else if (sibcType.compare("Layers")==0) {
            return readMultilayerSurface(j);
        } else {
            throw std::logic_error("Unrecognized SIBC type: " + sibcType);
        }
    }

    case PM::Type::gap:
        return std::make_unique < physicalModel::Gap>(id, name, j.at("width").get<double>());

    case PM::Type::multiport:
    {
        using namespace physicalModel::multiport;
        auto mpType = strToMultiportType(j.at("connectorType").get<std::string>());
        switch (mpType) {
        case Multiport::Type::shortCircuit:
            return std::make_unique<Predefined>(id, name, mpType);
        case Multiport::Type::openCircuit:
			return  std::make_unique<Predefined>(id, name, mpType);
        case Multiport::Type::dispersive:
            return  std::make_unique<Dispersive>(id, name, j.at("filename").get<std::string>());
        default:
            return  std::make_unique<RLC>(id, name, mpType,
                    j.at("resistance").get<double>(),
                    j.at("inductance").get<double>(),
                    j.at("capacitance").get<double>());
        }
    }

    default:
        throw std::logic_error("Material type not recognized for: " + name);
    }
}

outputRequest::OutputRequest::Type strToOutputType(std::string str) {
    str = trim(str);
    if (str.compare("electricField") == 0) {
        return outputRequest::OutputRequest::Type::electric;
    }
    else if (str.compare("magneticField") == 0) {
        return outputRequest::OutputRequest::Type::magnetic;
    }
    else if (str.compare("electricFieldNormals") == 0) {
        return outputRequest::OutputRequest::Type::electricFieldNormals;
    }
    else if (str.compare("magneticFieldNormals") == 0) {
        return outputRequest::OutputRequest::Type::magneticFieldNormals;
    }
    else if (str.compare("current") == 0) {
        return outputRequest::OutputRequest::Type::current;;
    }
    else if (str.compare("bulkCurrentElectric") == 0) {
        return outputRequest::OutputRequest::Type::bulkCurrentElectric;
    }
    else if (str.compare("bulkCurrentMagnetic") == 0) {
        return outputRequest::OutputRequest::Type::bulkCurrentMagnetic;
    }
    else if (str.compare("surfaceCurrentDensity") == 0) {
        return outputRequest::OutputRequest::Type::surfaceCurrentDensity;
    }
    else if (str.compare("farField") == 0) {
        return outputRequest::OutputRequest::Type::electric;
    }
    else {
        throw std::logic_error("Unrecognized output type: " + str);
    }
}

outputRequest::Domain readDomain(const json& j)
{
    bool timeDomain = false;
    Real initialTime = 0.0;
    Real finalTime = 0.0;
    Real samplingPeriod = 0.0;

    bool frequencyDomain = false;
    bool logFrequencySweep = false;
    bool usingTransferFunction = false;
    Real initialFrequency = 0.0;
    Real finalFrequency = 0.0;
    Real frequencyStep = 0.0;
    std::string transferFunctionFile;

    if (j.find("initialTime") != j.end()) {
        timeDomain = true;
        initialTime = j.at("initialTime").get<double>();
        finalTime = j.at("finalTime").get<double>();
        samplingPeriod = j.at("samplingPeriod").get<double>();
    }

    if (j.find("initialFrequency") != j.end()) {
        frequencyDomain = true;
        initialFrequency = j.at("initialFrequency").get<double>();
        finalFrequency = j.at("finalFrequency").get<double>();
        frequencyStep = getProgressionStepByTotalNumber(j, "Frequency");

        logFrequencySweep = j.at("logFrequencySweep").get<bool>();
        if (j.find("transferFunctionFile") != j.end()) {
            usingTransferFunction = true;
            transferFunctionFile = j.at(
                "transferFunctionFile"
            ).get<std::string>();
        }
    }

    return outputRequest::Domain(
        timeDomain, initialTime, finalTime, samplingPeriod,
        frequencyDomain, initialFrequency, finalFrequency,
        frequencyStep, logFrequencySweep,
        usingTransferFunction, transferFunctionFile);
}

std::unique_ptr<outputRequest::OutputRequest> readProbe(mesh::Unstructured& mesh, const json& j)
{
    std::string name = j.at("name").get<std::string>();
    outputRequest::OutputRequest::Type type = strToOutputType(j.at("type").get<std::string>());
    std::string gidOutputType = j.at("gidOutputType").get<std::string>();
    outputRequest::Domain domain = readDomain(j.at("domain").get<json>());

    if (type == outputRequest::OutputRequest::Type::bulkCurrentElectric ||
        type == outputRequest::OutputRequest::Type::bulkCurrentMagnetic) {
        if (gidOutputType.compare("OutRq_on_layer") == 0) {
            return std::make_unique<outputRequest::BulkCurrent>(
                outputRequest::BulkCurrent(
                    domain,
                    name,
                    { boxToElemGroup(mesh, j.at("box").get<std::string>())->getId() },
                    strToCartesianAxis(j.at("direction").get<std::string>()),
                    j.at("skip").get<int>()
                )
            );
        }

        if (gidOutputType.compare("OutRq_on_point") == 0) {
            return std::make_unique<outputRequest::BulkCurrent>(
                outputRequest::BulkCurrent(
                    domain,
                    name,
                    readElemIds(j.at("elemIds").get<json>()),
                    strToCartesianAxis(j.at("direction").get<std::string>()),
                    j.at("skip").get<int>()
                )
            );
        }

        return std::make_unique<outputRequest::BulkCurrent>(
            outputRequest::BulkCurrent(
                domain,
                name,
                readElemIds(j.at("elemIds").get<json>()),
                strToCartesianAxis(j.at("direction").get<std::string>()),
                j.at("skip").get<int>()
            )
        );
    }

    if (gidOutputType.compare("OutRq_on_point") == 0) {
        return std::make_unique<outputRequest::OnPoint>(
            type, domain, name, readElemIds(j.at("elemIds").get<json>())
        );
    }

    if (gidOutputType.compare("OutRq_on_line") == 0) {
        return std::make_unique<outputRequest::OnLine>(
            outputRequest::OnLine(
                type,
                domain,
                name,
                readElemIds(j.at("elemIds").get<json>())
            )
        );
    }

    if (gidOutputType.compare("OutRq_on_surface") == 0) {
        return std::make_unique<outputRequest::OnSurface>(
            outputRequest::OnSurface(
                type,
                domain,
                name,
                readElemIds(j.at("elemIds").get<json>())
            )
        );
    }

    if (gidOutputType.compare("OutRq_on_layer") == 0) {
        return std::make_unique<outputRequest::OnLayer>(
            outputRequest::OnLayer(
                type,
                domain,
                name,
                { boxToElemGroup(mesh, j.at("box").get<std::string>())->getId() }
            )
        );
    }

    if (gidOutputType.compare("Far_field") == 0) {
        static const Real degToRad = 2.0 * Constants::pi / 360.0;
        return std::make_unique<outputRequest::FarField>(
            outputRequest::FarField(
                domain,
                name,
                { boxToElemGroup(mesh, j.at(
                    j.find("layerBox") == j.end() ? "box" : "layerBox"
                ).get<std::string>())->getId() },
                j.at("initialTheta").get<double>() * degToRad,
                j.at("finalTheta").get<double>() * degToRad,
                getProgressionStepByTotalNumber(j, "Theta") * degToRad,
                j.at("initialPhi").get<double>() * degToRad,
                j.at("finalPhi").get<double>() * degToRad,
                getProgressionStepByTotalNumber(j, "Phi") * degToRad
            )
        );
    }

    throw std::logic_error("Unrecognized GiD Output request type: " + gidOutputType);
}

OutputRequestGroup readProbes(mesh::Unstructured& mesh, const json& j)
{
    OutputRequestGroup res;
    
    auto outs = j.find("probes");
    if (outs == j.end()) {
        return res;
    }

    for (auto const& out: outs->get<json>()) {
        res.addAndAssignId(readProbe(mesh, out));
    }

    return res;
}

const ElemR* boxToElemGroup(mesh::Unstructured& mesh, const std::string& line)
{
    BoxR3 box = strToBox(line);
    std::unique_ptr<ElemR> elem;

    auto posVec = box.getPos();

    if (box.isVolume()) {    
        elem = std::make_unique<HexR8>(ElemId(), addAngGetCoordView(mesh.coords(), posVec, 8).data());
    } else if (box.isSurface()) {
        elem = std::make_unique<QuaR4>(ElemId(), addAngGetCoordView(mesh.coords(), posVec, 4).data());
    } else if (box.isLine()) {
		elem = std::make_unique<LinR2>(ElemId(), addAngGetCoordView(mesh.coords(), posVec, 2).data());
	} else if (box.isPoint()) {
        elem = std::make_unique<NodR>(ElemId(), addAngGetCoordView(mesh.coords(), posVec, 1).data());
	} else {
        throw std::logic_error("Box to Elem Group only works for volumes and surfaces");
    }

    return mesh.elems().addAndAssignId(std::move(elem))->get();
}

Constants::CartesianAxis strToCartesianAxis(std::string str) {
    if (str.compare("x") == 0) {
        return Constants::x;
    } else if (str.compare("y") == 0) {
        return Constants::y;
    } else if (str.compare("z") == 0) {
        return Constants::z;
    } else {
        throw std::logic_error("Unrecognized cartesian axis label: " + str);
    }
}

LayerGroup readLayers(const json& j) 
{
    if (j.find("layers") == j.end()) {
        throw std::logic_error("layers object was not found.");
    }

    LayerGroup res;
    for (auto const& it: j.at("layers")) {
        res.add(
            std::make_unique<Layer>(
                LayerId(it.at("id").get<int>()),
                it.at("name").get<std::string>()
            )
        );
    }
    return res;
}

CoordR3Group readCoordinates(const json& j)
{

    if (j.find("coordinates") == j.end()) {
        throw std::logic_error("Coordinates label was not found.");
    }

    CoordR3Group res;
    const json& c = j.at("coordinates").get<json>();
    for (json::const_iterator it = c.begin(); it != c.end(); ++it) {
        CoordId id;
        CVecR3 pos;
        std::stringstream ss(it->get<std::string>());
        ss >> id >> pos(0) >> pos(1) >> pos(2);
        res.add(std::make_unique<CoordR3>(id, pos));
    }
    return res;
}

ElemRGroup readElements(
	const PMGroup& mG,
	LayerGroup& lG,
	CoordR3Group& cG,
	const json& j,
    const std::string& folder) 
{
	if (j.find("elements") == j.end()) {
		throw std::logic_error("Elements label was not found.");
	}

	ElemRGroup res;
	const json& elems = j.at("elements").get<json>();

	// Check how to iterate over `elems` keys, somehow
	if (elems.find("hexahedra") != elems.end()) {
		for (auto& element : readElemStrAs<HexR8>(mG, lG, cG, elems.at("hexahedra").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("tetrahedra") != elems.end()) {
		for (auto& element : readElemStrAs<Tet4>(mG, lG, cG, elems.at("tetrahedra").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("quadrilateral") != elems.end()) {
		for (auto& element : readElemStrAs<QuaR4>(mG, lG, cG, elems.at("quadrilateral").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("triangle") != elems.end()) {
		for (auto& element : readElemStrAs<Tri3>(mG, lG, cG, elems.at("triangle").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("line") != elems.end()) {
		for (auto& element : readElemStrAs<LinR2>(mG, lG, cG, elems.at("line").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("node") != elems.end()) {
		for (auto& element : readElemStrAs<NodR>(mG, lG, cG, elems.at("node").get<json>())) {
			res.add(std::move(element));
		}
	}

	if (elems.find("fromFile") != elems.end()) {
		res.addAndAssignIds(readElementsFromFile(mG, lG, cG, elems.at("fromFile").get<json>(), folder));
	}

    return res;
}

ElemRGroup readElementsFromSTLFile(
    const PMGroup& mG, LayerGroup& lG, CoordR3Group& cG, 
    const json& f, const std::string& folder)
{
    std::string fn = folder + f.at("file").get<std::string>();
    mesh::Unstructured m = parsers::STL::Parser(fn).readAsUnstructuredMesh();

    auto lay = lG.getId( LayerId(f.at("layerId").get<std::size_t>())  );
    auto mat = mG.getId( MatId(f.at("materialId").get<std::size_t>()) );
    
    std::map<const CoordR3*, const CoordR3*> readedToGlobal;
    for (const auto& c : m.coords()) {
        const auto it{ cG.addAndAssignId(std::make_unique<CoordR3>(*c)) };
        readedToGlobal.emplace(c.get(), it->get());
    }

	ElemRGroup res;
	for (const auto& elem : m.elems()) {
		std::vector<const CoordR3*> vs;
		for (std::size_t i = 0; i < elem->numberOfVertices(); i++) {
            auto it{ readedToGlobal.find(elem->getV(i)) };

            assert(it != readedToGlobal.end());

            vs.push_back(it->second);
		}

		res.addAndAssignId(std::make_unique<Tri3>(ElemId(0), vs.data(), lay, mat));
	}

    return res;
}

ElemRGroup readElementsFromFile(
    const PMGroup& mG,
    LayerGroup& lG,
    CoordR3Group& cG,
    const json& eFile,
    const std::string& folder) 
{
    ElemRGroup res;
    for (auto const& f : eFile) {
        if (f.find("format") == f.end()) {
            throw std::runtime_error("Format label not found when reading elements from file.");
        }

        std::string format = f.at("format").get<std::string>();
        if (format == "STL") {
            res.addAndAssignIds(readElementsFromSTLFile(mG, lG, cG, f, folder));

            continue;
		}
		
        if (format == "ignore") {
			// Ignoring file
            continue;
		}

		throw std::runtime_error("Unsupported file format when reading elements from file.");
	}
    return res;
}

std::unique_ptr<physicalModel::surface::Multilayer> 
readMultilayerSurface(const json& mat)
{
    MatId id = MatId(mat.at("materialId").get<int>());
    std::string name = mat.at("name").get<std::string>();

    std::vector<physicalModel::surface::Multilayer::Layer> layers;
    for (json::const_iterator it = mat.at("layers").begin(); it != mat.at("layers").end(); ++it) {
        layers.push_back(physicalModel::surface::Multilayer::Layer(
            it->at("thickness").get<double>(),
            it->at("permittivity").get<double>(),
            it->at("permeability").get<double>(),
            it->at("elecCond").get<double>() )
        );
    }

    if (mat.at("useSembaVectorFitting").get<bool>()) {
        physicalModel::surface::Multilayer::FittingOptions opts(
                std::make_pair(mat.at("freqMin").get<double>(),
                        mat.at("freqMax").get<double>()),
                        mat.at("numberOfPoles").get<int>());
        std::vector<physicalModel::surface::Multilayer::FittingOptions> optsVec = { opts };
        return std::make_unique<physicalModel::surface::Multilayer>(id, name, layers, optsVec);
    } else {
        return std::make_unique < physicalModel::surface::Multilayer>(id, name, layers);
    }
}

Grid3 readGrids(const json& j) 
{
    if (j.find("grids") == j.end()) {
        throw std::logic_error("Grids object not found.");
    }

    json g = j.at("grids").front();
    std::string gridType = g.at("gridType").get<std::string>();
    if (gridType.compare("gridCondition") == 0) {
        // Initializes basic grid.
        Grid3 res;
        if (g.at("type").get<std::string>().compare("Number_of_cells") == 0) {
            res = Grid3(
                    strToBox(g.at("layerBox").get<std::string>()),
                    strToCVecI3(g.at("numberOfCells").get<std::string>()));
        } else {
            BoxR3 box = strToBox(g.at("layerBox").get<std::string>());
            CVecR3 stepSize = strToCVecR3(g.at("stepSize").get<std::string>());
            if (g.at("fitSizeToBox").get<bool>()) {
                for (std::size_t i = 0; i < 3; i++) {
					std::size_t n = std::round(box.getLength()[i] / stepSize[i]);
					if (n == 0) {
						n = 1;
					}
                    stepSize[i] = box.getLength()[i] / n;
                }
            }
            res = Grid3(box, stepSize);
        }

        // Applies boundary padding operations.
        if (g.find("lowerPaddingMeshSize") != g.end()) {
            std::pair<CVecR3, CVecR3> boundaryMeshSize(
                    strToCVecR3(g.at("lowerPaddingMeshSize").get<std::string>()),
                    strToCVecR3(g.at("upperPaddingMeshSize").get<std::string>())
            );
            std::pair<CVecR3, CVecR3> boundaryPadding(
                    strToCVecR3(g.at("lowerPadding").get<std::string>()),
                    strToCVecR3(g.at("upperPadding").get<std::string>())
            );
            if (g.at("boundaryPaddingType").get<std::string>().compare(
                    "by_number_of_cells") == 0) {
                boundaryPadding.first  *= boundaryMeshSize.first;
                boundaryPadding.second *= boundaryMeshSize.second;
            }
            res.enlarge(boundaryPadding, boundaryMeshSize);
        }

		if (res.hasZeroSize()) {
			throw std::logic_error("Grid has zero size.");
		}
        return res;

    } 
    if (gridType.compare("nativeGiD") == 0) {
        CVecR3 corner = strToCVecR3(g.at("corner").get<std::string>());
        CVecR3 boxSize = strToCVecR3(g.at("boxSize").get<std::string>());
        CVecI3 nGridPoints = strToCVecI3(g.at("nGridPoints").get<std::string>());
        std::vector<Real> pos[3];
        pos[0] = g.at("xCoordinates").get<std::vector<double>>();
        pos[1] = g.at("yCoordinates").get<std::vector<double>>();
        pos[2] = g.at("zCoordinates").get<std::vector<double>>();
        if (!pos[0].empty()) {
            return Grid3(pos);
        }
        else {
            std::pair<CVecR3, CVecR3> box =
            { corner, corner + boxSize };
            return Grid3(box, nGridPoints);
        }
    } 
    if (gridType.compare("rectilinear") == 0) {
        std::map<std::string, std::vector<Real>> planes;
        for (const auto& label: { "xs", "ys", "zs"}) {
            auto dir{ g.find(label) };
            if (dir == g.end()) {
                throw std::logic_error("Missing position defintion in grid file in at least one direction");
            }
            auto plane{ dir->get<std::vector<double>>() };
            if (plane.empty()) {
                throw std::logic_error("Grid file had empty positions in at least one direction");
            }
            planes.emplace(label, plane);
        }
        std::vector<Real> pos[3]{ planes["xs"], planes["ys"],planes["zs"]};
        return Grid3(pos);
    } 
    
    throw std::logic_error("Unrecognized grid type: " + gridType);
    
}

std::unique_ptr<source::PlaneWave> readPlanewave(mesh::Unstructured& mesh, const json& j) {
    
    auto magnitude{ readMagnitude(j.at("magnitude").get<json>()) };
    auto elemId{ boxToElemGroup(mesh, j.at("layerBox").get<std::string>())->getId()};
    
    auto definitionMode = j.at("definitionMode").get<std::string>();
    if (definitionMode.compare("by_vectors")==0) {
		return std::make_unique<source::PlaneWave>(
			source::PlaneWave(
				magnitude,
                { elemId },
				strToCVecR3(j.at("directionVector").get<std::string>()),
				strToCVecR3(j.at("polarizationVector").get<std::string>())
			)
		);
    } else if (definitionMode.compare("by_angles")==0) {
        static const Real degToRad = 2.0 * Constants::pi / 360.0;
        std::pair<Real,Real> dirAngles, polAngles;
        dirAngles.first  = j.at("directionTheta").get<double>()    * degToRad;
        dirAngles.second = j.at("directionPhi").get<double>()      * degToRad;
        polAngles.first  = j.at("polarizationAlpha").get<double>() * degToRad;
        polAngles.second = j.at("polarizationBeta").get<double>()  * degToRad;
		return std::make_unique<source::PlaneWave>(
			source::PlaneWave(
				magnitude,
                { elemId },
				dirAngles,
				polAngles
			)
		);

    } else if (definitionMode.compare("randomized_multisource")==0) {
		return std::make_unique<source::PlaneWave>(
			source::PlaneWave(
				magnitude,
                { elemId },
				j.at("numberOfRandomPlanewaves").get<int>(),
				j.at("relativeVariationOfRandomDelay").get<double>()
			)
        );
    } else {
        throw std::logic_error("Unrecognized label: " + definitionMode);
    }
}

std::unique_ptr<source::port::Waveguide> readPortWaveguide(
    mesh::Unstructured& mesh, 
    const json& j
) {
	std::string shape = j.at("shape").get<std::string>();
	if (shape.compare("Rectangular") == 0) {
		return std::make_unique<source::port::WaveguideRectangular>(
			source::port::WaveguideRectangular(
				readMagnitude(j.at("magnitude").get<json>()),
				readElemIds(j.at("elemIds").get<json>()),
				strToWaveguideMode(j.at("excitationMode").get<std::string>()),
				{ j.at("firstMode").get<int>(), j.at("secondMode").get<int>() }
			)
		);
	}
	else {
		throw std::logic_error("Unrecognized waveguide port shape: " + shape);
	}
}

std::unique_ptr<source::port::TEM> readPortTEM(mesh::Unstructured& mesh, const json& j) 
{
	return std::make_unique<source::port::TEMCoaxial>(
		source::port::TEMCoaxial(
			readMagnitude(j.at("magnitude").get<json>()),
			readElemIds(j.at("elemIds").get<json>()),
			strToTEMMode(j.at("excitationMode").get<std::string>()),
			strToCVecR3(j.at("origin").get<std::string>()),
			j.at("innerRadius").get<double>(),
			j.at("outerRadius").get<double>()
		)
	);
}

std::unique_ptr<source::Generator> readGenerator(mesh::Unstructured& mesh, const json& j) 
{
	return std::make_unique<source::Generator>(
		source::Generator(
			readMagnitude(j.at("magnitude").get<json>()),
			readElemIds(j.at("elemIds").get<json>()),
			strToGeneratorType(j.at("type").get<std::string>()),
			source::Generator::soft
        )
	);
}

std::unique_ptr<source::OnLine> readSourceOnLine(mesh::Unstructured& mesh, const json& j) 
{
	return std::make_unique<source::OnLine>(
		source::OnLine(
			readMagnitude(j.at("magnitude").get<json>()),
			readElemIds(j.at("elemIds").get<json>()),
			strToNodalType(j.at("type").get<std::string>()),
			strToNodalHardness(j.at("hardness").get<std::string>())
        )
	);
}

source::Generator::Type strToGeneratorType(std::string str) 
{
    str = trim(str);
    if (str.compare("voltage")==0) {
        return source::Generator::voltage;
    }
    
    if (str.compare("current")==0) {
        return source::Generator::current;
    }

    throw std::logic_error("Unrecognized generator type: " + str);
}

source::Generator::Hardness strToGeneratorHardness(std::string str) 
{
    str = trim(str);
    if (str.compare("soft")==0) {
        return source::Generator::soft;
    } else if (str.compare("hard")==0) {
        return source::Generator::hard;
    } else {
        throw std::logic_error("Unrecognized generator hardness: " + str);
    }
}

physicalModel::PhysicalModel::Type strToMaterialType(std::string str) 
{
    using Type = semba::physicalModel::PhysicalModel::Type;

    str = trim(str);
    if (str.compare("PEC")==0) {
        return Type::PEC;
    } else if (str.compare("PMC")==0) {
        return Type::PMC;
    } else if (str.compare("PML")==0) {
        return Type::PML;
    } else if (str.compare("SMA")==0) {
        return Type::SMA;
    } else if (str.compare("Vacuum") == 0) {
        return Type::vacuum;
    } else if (str.compare("Classic")==0) {
        return Type::classic;
    } else if (str.compare("Dispersive")==0) {
        return Type::elecDispersive;
    } else if (str.compare("Anisotropic")==0) {
        return Type::anisotropic;
    } else if (str.compare("SIBC")==0) {
        return Type::isotropicsibc;
    } else if (str.compare("Wire")==0) {
        return Type::wire;
    } else if (str.compare("Connector")==0) {
        return Type::multiport;
    } else if (str.find("Thin_gap")==0) {
        return Type::gap;
    } else if (str.find("PriorityMaterial")==0) {
        return Type::priorityMaterial;
    } else {
        throw std::logic_error("Unrecognized material label: " + str);
    }
}

physicalModel::multiport::Multiport::Type strToMultiportType(std::string str) {
    using namespace physicalModel::multiport;

    str = trim(str);
    if (str.compare("Conn_short")==0) {
        return Multiport::Type::shortCircuit;
    } else if (str.compare("Conn_open")==0) {
        return Multiport::Type::openCircuit;
    } else if (str.compare("Conn_matched")==0) {
        return Multiport::Type::matched;
    } else if (str.compare("Conn_sRLC")==0) {
        return Multiport::Type::sRLC;
    } else if (str.compare("Conn_pRLC")==0) {
        return Multiport::Type::pRLC;
    } else if (str.compare("Conn_sLpRC")==0) {
        return Multiport::Type::sLpRC;
    } else if (str.compare("Conn_dispersive") == 0) {
        return Multiport::Type::dispersive;
    } else {
        throw std::logic_error("Unrecognized multiport label: " + str);
    }
}

std::pair<CVecR3, CVecR3> strToBox(
        const std::string& value) {
    std::size_t begin = value.find_first_of("{");
    std::size_t end = value.find_last_of("}");
    std::string aux = value.substr(begin+1,end-1);
    std::stringstream iss(aux);
    CVecR3 max, min;
    for (std::size_t i = 0; i < 3; i++) {
        iss >> max(i);
    }
    for (std::size_t i = 0; i < 3; i++) {
        iss >> min(i);
    }
    return {min, max};
}

CVecI3 strToCVecI3(std::string str) {
    str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
    std::stringstream ss(str);
    CVecI3 res;
    ss >> res(Constants::x)
       >> res(Constants::y)
       >> res(Constants::z);
    return res;
}

CVecR3 strToCVecR3(std::string str) {
    str.erase(std::remove(str.begin(), str.end(), '{'), str.end());
    str.erase(std::remove(str.begin(), str.end(), '}'), str.end());
    std::stringstream ss(str);
    CVecR3 res;
    ss >> res(Constants::x)
       >> res(Constants::y)
       >> res(Constants::z);
    return res;
}

source::OnLine::Type strToNodalType(std::string str) 
{
    str = trim(str);
    if (str.compare("electricField")==0) {
        return source::OnLine::Type::electric;
    } else if (str.compare("magneticField")==0) {
        return source::OnLine::Type::magnetic;
    } else {
        throw std::logic_error("Unrecognized nodal type: " + str);
    }
}

source::OnLine::Hardness strToNodalHardness(std::string str) 
{
    str = trim(str);
    if (str.compare("soft")==0) {
        return source::OnLine::Hardness::soft;
    } else if (str.compare("hard")==0) {
        return source::OnLine::Hardness::hard;
    } else {
        throw std::logic_error("Unrecognized nodal hardness: " + str);
    }
}

std::unique_ptr<source::Magnitude::Magnitude> readMagnitude(const json& j) 
{
    std::string type = j.at("type").get<std::string>();
    if (type.compare("File") == 0) {
		return std::make_unique<source::Magnitude::Numerical>(
			source::Magnitude::Numerical(
				j.at("filename").get<std::string>())
		);
    }

    if (type.compare("Gaussian") == 0) {
		return std::make_unique<source::Magnitude::Magnitude>(
			source::Magnitude::Magnitude(
				new function::Gaussian(
					function::Gaussian::buildFromMaximumFrequency(
						j.at("frequencyMaximum").get<double>(),
						1.0
					)
				)
			)
		);
    }

    if (type.compare("Band_limited") == 0) {
		return std::make_unique<source::Magnitude::Magnitude>(
			source::Magnitude::Magnitude(
				new function::BandLimited(
					j.at("frequencyMinimum").get<double>(),
					j.at("frequencyMaximum").get<double>()))
		);
    }

    throw std::logic_error("Unable to recognize magnitude type when reading excitation.");
}

LocalAxis strToLocalAxes(const std::string& str) {
    std::size_t begin = str.find_first_of("{");
    std::size_t end = str.find_first_of("}");
    CVecR3 eulerAngles = strToCVecR3(str.substr(begin+1,end-1));
    begin = str.find_last_of("{");
    end = str.find_last_of("}");
    CVecR3 origin = strToCVecR3(str.substr(begin+1,end-1));
    return LocalAxis(eulerAngles, origin);
}

source::port::TEM::ExcitationMode strToTEMMode(std::string str) {
    if (str.compare("Voltage") == 0) {
        return source::port::TEM::voltage;
    } else if (str.compare("Current") == 0) {
        return source::port::TEM::current;
    } else {
        throw std::logic_error("Unrecognized exc. mode label: " + str);
    }

}

source::port::Waveguide::ExcitationMode strToWaveguideMode(
        std::string str) {
    if (str.compare("TE") == 0) {
        return source::port::Waveguide::ExcitationMode::TE;
    } else if (str.compare("TM") == 0) {
        return source::port::Waveguide::ExcitationMode::TM;
    } else {
        throw std::logic_error("Unrecognized exc. mode label: " + str);
    }
}


}
