#include "MeshTest.h"
#include "core/geometry/mesh/Unstructured.h"

class MeshUnstructuredTest : public ::testing::Test,
    public MeshTest {
public:
    void SetUp() {
        MeshTest::SetUp();
        mesh_ = mesh::Unstructured(cG_, eG_, lG_);
    }

protected:
    mesh::Unstructured mesh_;
};
TEST_F(MeshUnstructuredTest, ctor) 
{
    EXPECT_EQ(cG_.size(), mesh_.coords().size());

    EXPECT_EQ(eG_.size(), mesh_.elems().size());
    EXPECT_EQ(eG_.sizeOf<Tet>(), mesh_.elems().getOf<Tet>().size());
    EXPECT_EQ(eG_.sizeOf<Tri>(), mesh_.elems().getOf<Tri>().size());

    EXPECT_EQ(lG_.size(), mesh_.layers().size());
}

TEST_F(MeshUnstructuredTest, copy_assignment)
{
    mesh::Unstructured mesh = mesh_;
    for (auto const& copiedCoordIt : mesh.coords()) {
        const CoordR3* copiedCoord = copiedCoordIt.get();

        auto id = copiedCoord->getId();
        auto coordIt = mesh_.coords().getId(id);
        const CoordR3* originalCoord = coordIt;

        EXPECT_EQ(*copiedCoord, *originalCoord);
    }
}

TEST_F(MeshUnstructuredTest, move_assignment)
{
    mesh::Unstructured mesh = std::move(mesh_);
    
    EXPECT_EQ(0, mesh_.coords().size());

}
