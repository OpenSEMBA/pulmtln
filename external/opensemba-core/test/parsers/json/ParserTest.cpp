#include <gtest/gtest.h>

#include "parsers/json/Parser.h"

#include "core/math/function/Gaussian.h"
#include "core/geometry/element/Line2.h"
#include "core/geometry/element/Triangle3.h"
#include "core/geometry/element/Tetrahedron4.h"
#include "core/outputRequest/FarField.h"
#include "core/outputRequest/OnPoint.h"
#include "core/physicalModel/wire/Wire.h"
#include "core/physicalModel/multiport/RLC.h"
#include "core/source/PlaneWave.h"
#include "core/source/Generator.h"

using namespace semba;
using namespace semba::parsers::JSON;
using namespace geometry::mesh;

class ParserJSONParserTest : public ::testing::Test {
protected:
    static std::string getFolder() {
        return "testData/";
    }

    static std::string getFilename(const std::string& name)
    {
        return getFolder() + name + "/" + name + ".smb.json";
    }
};

TEST_F(ParserJSONParserTest, dmcwf) 
{
    ASSERT_NO_THROW(Parser{ getFilename("dmcwf") }.read());
}

TEST_F(ParserJSONParserTest, sphere)
{
    ASSERT_NO_THROW(Parser{ getFilename("sphere") }.read());
}

TEST_F(ParserJSONParserTest, antennas)
{
    ASSERT_NO_THROW(Parser{ getFilename("antennas") }.read());
}

TEST_F(ParserJSONParserTest, wires)
{
    ASSERT_NO_THROW(Parser{ getFilename("wires") }.read());
}

TEST_F(ParserJSONParserTest, bowtie)
{
    ASSERT_NO_THROW(Parser{ getFilename("bowtie") }.read());
}

TEST_F(ParserJSONParserTest, b2)
{
    ASSERT_NO_THROW(Parser{ getFilename("b2") }.read());
}

TEST_F(ParserJSONParserTest, b2_detailed)
{
    auto data{ Parser{getFilename("b2") }.read() };

    EXPECT_EQ(353, data.model.mesh.coords().size());
    EXPECT_EQ(652, data.model.mesh.elems().sizeOf<geometry::Tri3>());
}

TEST_F(ParserJSONParserTest, bowtie_detailed)
{
    auto data{ Parser(getFilename("bowtie")).read() };

    EXPECT_EQ(422, data.model.mesh.coords().size());
    EXPECT_EQ(836, data.model.mesh.elems().size());
}

TEST_F(ParserJSONParserTest, sphere_detailed) 
{
    auto data{ Parser(getFilename("sphere")).read() };

    EXPECT_EQ(data.grids.getNumCells(), math::CVecR3(51, 23, 15));

    auto& sources = data.sources;
    EXPECT_EQ(sources.size(), 1);

    const source::PlaneWave* source = sources.get()[0]->castTo<source::PlaneWave>();
    EXPECT_EQ(
        source->getPolarization(),
        math::CVecR3(-0.4082482904638631, 0.8164965809277261, -0.4082482904638631)
    );
    EXPECT_EQ(
        source->getDirection(),
        math::CVecR3(1.0, 1.0, 1.0)
    );

    source::Magnitude::Magnitude magnitude = *source->getMagnitude();
    EXPECT_EQ(
        magnitude,
        source::Magnitude::Magnitude(
            new math::function::Gaussian(
                math::function::Gaussian::buildFromMaximumFrequency(
                    1000000000.0,
                    1.0
                )
            )
        )
    );

    auto& model = data.model;

    // 2 predefined + 6 bounds
    EXPECT_EQ(8, model.physicalModels.size());
    EXPECT_EQ("pec", model.physicalModels.get()[0]->getName());
    EXPECT_EQ("pmc", model.physicalModels.get()[1]->getName());

    const physicalModel::Bound::Type boundaryLowerUpperMaterials[6] = {
        physicalModel::Bound::Type::pml,
        physicalModel::Bound::Type::pec,
        physicalModel::Bound::Type::pmc,
        physicalModel::Bound::Type::mur1,
        physicalModel::Bound::Type::mur2,
        physicalModel::Bound::Type::periodic
    };

    auto boundaries = model.physicalModels.getOf<physicalModel::Bound>();
    for (auto& bound : boundaryLowerUpperMaterials) {
        EXPECT_TRUE(
            std::find_if(
                boundaries.cbegin(),
                boundaries.cend(),
                [&](const physicalModel::Bound* boundI) -> bool {return boundI->getType() == bound; }
            ) != boundaries.cend()
        );
    }

    // 384 coming from grid, 8 from source, 1 from point probe and 8 from far field
    // New elements added as part of Boundaries: 6 faces * 4 points/face
    ASSERT_EQ(384 + 8 + 1 + 8 + 24, model.mesh.coords().size());
    EXPECT_EQ(
        math::CVecR3(2.33333325, -5.71501865e-16, 1.66666663),
        model.mesh.coords().get()[0]->pos()
    );
    EXPECT_EQ(
        math::CVecR3(1.28204191, -1.31762123e+01, -1.70370862e-01),
        model.mesh.coords().get()[383]->pos()
    );

    auto& probes = data.outputRequests;

    ASSERT_EQ(2, probes.size());
    EXPECT_EQ(outputRequest::OutputRequest::Type::electric, probes.get()[0]->getType());
    EXPECT_EQ(outputRequest::OutputRequest::Type::electricFarField, probes.get()[1]->getType());

    ASSERT_EQ(1, probes.get()[0]->getTarget().size());

    auto recoveredNodeId = probes.get()[0]->getTarget().at(0);
    auto recoveredNode = model.mesh.elems().getId(recoveredNodeId)->castTo<geometry::NodR>();
    EXPECT_EQ(
        math::CVecR3(-0.8441360141053171, 12.017228978451016, 13.154724231963254),
        recoveredNode->getV(0)->pos()
    );
}

TEST_F(ParserJSONParserTest, sphere_rectilinear) 
{
    auto data{ Parser(getFolder() + "sphere/sphere-rectilinear.smb.json").read() };

    EXPECT_EQ(math::CVecI3(1,2,3), data.grids.getNumCells());
    EXPECT_EQ(std::vector<math::Real>({0, 1}), data.grids.getPos(0));
    EXPECT_EQ(std::vector<math::Real>({0, 1, 2}), data.grids.getPos(1));
    EXPECT_EQ(std::vector<math::Real>({0, 1, 2, 3}), data.grids.getPos(2));
}

TEST_F(ParserJSONParserTest, sphere_onePlaneFarField)
{
    auto data{ Parser{ getFolder() + "sphere/sphere-one-plane-farfield.smb.json" }.read() };

    EXPECT_EQ(data.outputRequests.sizeOf<outputRequest::FarField>(), 1);
    auto farFieldProbe{ data.outputRequests.getOf<outputRequest::FarField>().front() };
    EXPECT_EQ(farFieldProbe->initialPhi, 0.0);
    EXPECT_EQ(farFieldProbe->finalPhi, 0.0);
    EXPECT_EQ(farFieldProbe->stepPhi, 0.1 * 2 * math::Constants::pi / 360.0);
}

TEST_F(ParserJSONParserTest, antennas_detailed)
{
    auto data{ Parser{ getFilename("antennas") }.read() };

    EXPECT_EQ(data.outputRequests.sizeOf<outputRequest::OnPoint>(), 3);
    EXPECT_EQ(data.sources.sizeOf<source::Generator>(), 1);
    EXPECT_EQ(data.model.mesh.elems().sizeOf<geometry::NodR>(), 5);

    EXPECT_EQ(data.model.physicalModels.size(), 5); // Cable, 2 connector, 2 bounds (pec and pml)

    auto materialCableList = data.model.physicalModels.getOf<physicalModel::wire::Wire>();
    EXPECT_EQ(materialCableList.size(), 1);

    auto materialPortList = data.model.physicalModels.getOf<physicalModel::multiport::RLC>();
    EXPECT_EQ(materialPortList.size(), 1);

    auto& materialCable = materialCableList.front();
    auto& materialPort = materialPortList.front();

    geometry::ElemView elementsWithCableMaterial;
    for (auto& elem : data.model.mesh.elems()) {
        if (elem->getMatId() == materialCable->getId()) {
            elementsWithCableMaterial.push_back(elem.get());
        }
    }

    EXPECT_EQ(elementsWithCableMaterial.size(), 2);

    geometry::ElemView elementsWithPortMaterial;
    for (auto& elem : data.model.mesh.elems()) {
        if (elem->getMatId() == materialPort->getId()) {
            elementsWithPortMaterial.push_back(elem.get());
        }
    }

    EXPECT_EQ(elementsWithPortMaterial.size(), 2);
}

TEST_F(ParserJSONParserTest, readMaterials) 
{
    std::ifstream stream(getFolder() + "materials.smb.json");
    json j;
    stream >> j;
    
    auto materials{ parsers::JSON::readMaterials(j) };

    EXPECT_EQ(4, materials.size());
}