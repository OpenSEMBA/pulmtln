#include <gtest/gtest.h>

#include "core/util/GroupIdentifiableUnique.h"
#include "core/geometry/Layer.h"

using namespace semba;
using namespace geometry;
using namespace util;

class GroupIdentifiableUniqueTest : public ::testing::Test {
public:
    
protected:
    auto buildGroupOfThreeLayers() const 
    {
        GroupIdentifiableUnique<Layer> res;
        res.add(std::make_unique<Layer>(LayerId(1), "Patata"));
        res.add(std::make_unique<Layer>(LayerId(2), "Cebolla"));
        res.add(std::make_unique<Layer>(LayerId(3), "Huevos"));
        return res;
    }
};

TEST_F(GroupIdentifiableUniqueTest, CopyAndMoveCtor) 
{
    auto orig = buildGroupOfThreeLayers();
    std::size_t origSize = orig.size();

    auto copied( orig );
    EXPECT_EQ(origSize, copied.size());

    auto moved(std::move(orig));
    EXPECT_EQ(origSize, moved.size());
    EXPECT_EQ(origSize, copied.size());
    EXPECT_EQ(0, orig.size());
}

TEST_F(GroupIdentifiableUniqueTest, CopyAndMoveAssignment) 
{
    auto orig = buildGroupOfThreeLayers();
    std::size_t origSize = orig.size();
    
    auto copied = orig;
    EXPECT_EQ(origSize, copied.size());
    
    auto moved = std::move(orig);
    EXPECT_EQ(origSize, moved.size());
    EXPECT_EQ(origSize, copied.size());
    EXPECT_EQ(0, orig.size());
}

TEST_F(GroupIdentifiableUniqueTest, DeepCopyAndMoveAdd) 
{
    auto orig{ buildGroupOfThreeLayers() };
    
    orig.add(std::make_unique<Layer>(LayerId(5), "Melon"));

    // Keep building and moving as separate operations to avoid warnings by intel compiler.
    auto sandia{ std::make_unique<Layer>(LayerId(6), "Sandia") };
    orig.add(std::move(sandia));

    EXPECT_EQ(5, orig.size());
}
 
TEST_F(GroupIdentifiableUniqueTest, AddAndAssignId) 
{
    auto orig(buildGroupOfThreeLayers());
    auto melon = std::make_unique<Layer>("Melon");
    auto melonInGroup = orig.addAndAssignId(std::move(melon))->get();
    EXPECT_EQ(LayerId(4), melonInGroup->getId());
}

TEST_F(GroupIdentifiableUniqueTest, CopyAndAssignId)
{
    auto orig(buildGroupOfThreeLayers());
    auto melonInGroup = orig.copyAndAssignId(Layer{"Melon"})->get();
    EXPECT_EQ(LayerId(4), melonInGroup->getId());
}

TEST_F(GroupIdentifiableUniqueTest, AddAndAssignIds) 
{
	auto orig(buildGroupOfThreeLayers());

    // Keep building and moving as separate operations to avoid warnings by intel compiler.
    auto melon{ std::make_unique<Layer>("Melon") };
	orig.addAndAssignId(std::move(melon));

    auto nispora{ std::make_unique<Layer>("Nispora") };
	auto nisporaInGroup = orig.addAndAssignId(std::move(nispora))->get();

	EXPECT_EQ(LayerId(5), nisporaInGroup->getId());

	GroupIdentifiableUnique<Layer> another;
    auto naranja{ std::make_unique<Layer>("Naranja") };
	another.addAndAssignId(std::move(naranja));

	another.addAndAssignIds(orig);

	EXPECT_EQ(LayerId(5), nisporaInGroup->getId());

	EXPECT_EQ(6, another.size());
	EXPECT_EQ("Nispora", another.getId(LayerId(6))->getName());
}

TEST_F(GroupIdentifiableUniqueTest, DuplicatedIdShouldTriggerLogicException)
{
    auto orig(buildGroupOfThreeLayers());

    EXPECT_EQ(orig.getId(LayerId(1))->getName(), "Patata");

    orig.add(std::make_unique<Layer>(LayerId(1), "Another Patata"));

    EXPECT_EQ(orig.getId(LayerId(1))->getName(), "Another_Patata");

    for (const auto& item : orig) {
        EXPECT_TRUE(item->getName() != "Patata");
    }
}
