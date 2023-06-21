#include <gtest/gtest.h>

#include "core/physicalModel/Predefined.h"
#include "core/model/Model.h"

using ModelObject = semba::UnstructuredModel;

using namespace semba;
using namespace model;

TEST(ModelTest, CanCreate) {
	ModelObject model = ModelObject();

	EXPECT_NE(&model, nullptr);

	EXPECT_NE(&model.mesh, nullptr);
	EXPECT_TRUE(model.mesh.coords().empty());
	EXPECT_TRUE(model.mesh.elems().empty());
	EXPECT_TRUE(model.mesh.layers().empty());

	EXPECT_TRUE(model.physicalModels.empty());
}

TEST(ModelTest, CanInitializeGrid) {
	ModelObject model = ModelObject();

	CoordR3Group coordinatesGroup = CoordR3Group();
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(),
			math::CVecR3(0.0, 0.0, 0.0)
		)
	);
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(), 
			math::CVecR3(1.0, 2.0, 3.0)
		)
	);

	const CoordR3* coordinatesArgumentList[1] = {coordinatesGroup.getId(coordinate::Id(1))};
	ElemRGroup elementsGroup = ElemRGroup();
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			element::Id(), 
			coordinatesArgumentList
		)
	);

	mesh::Unstructured mesh(coordinatesGroup, elementsGroup);
	model.mesh = mesh;

	EXPECT_FALSE(model.mesh.coords().empty());
	EXPECT_EQ(2, model.mesh.coords().size());

	EXPECT_EQ(
		math::CVecR3(1.0, 2.0, 3.0), 
		(model.mesh.coords()).getId(coordinate::Id(2))->pos()
	);

	EXPECT_FALSE(model.mesh.elems().empty());
	EXPECT_EQ(1, model.mesh.elems().size());

	EXPECT_TRUE(model.mesh.layers().empty());
}

TEST(ModelTest, CanInitializePhysicalModels) {
	ModelObject model = ModelObject();

	PMGroup physicalModelsGroup = PMGroup();
	physicalModelsGroup.addAndAssignId(
		std::make_unique<physicalModel::PEC>(
			physicalModel::Id(),
			"Material PEC"
		)
	);

	model.physicalModels = physicalModelsGroup;

	EXPECT_EQ(1, model.physicalModels.size());
	EXPECT_EQ(
		"Material PEC",
		model.physicalModels.getId(physicalModel::Id(1))->getName()
	);
}

TEST(ModelTest, CanCopyConstructor) {
	CoordR3Group coordinatesGroup = CoordR3Group();
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(),
			math::CVecR3(0.0, 0.0, 0.0)
		)
	);
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(),
			math::CVecR3(1.0, 2.0, 3.0)
		)
	);

	const CoordR3* coordinatesArgumentList[1] = { coordinatesGroup.getId(coordinate::Id(1)) };
	ElemRGroup elementsGroup = ElemRGroup();
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			element::Id(),
			coordinatesArgumentList
		)
	);

	mesh::Unstructured mesh = mesh::Unstructured(coordinatesGroup, elementsGroup);

	PMGroup physicalModelsGroup = PMGroup();
	physicalModelsGroup.addAndAssignId(
		std::make_unique<physicalModel::PEC>(
			physicalModel::Id(),
			"Material PEC"
		)
	);

	ModelObject model(mesh, physicalModelsGroup);

	EXPECT_EQ(
		model.mesh.coords().getId(coordinate::Id(1)),
		model.mesh.elems().getId(ElemId(1))->getV(0)
	);

	ModelObject newModel(model);

	ASSERT_EQ(1, newModel.physicalModels.size());
	ASSERT_EQ(1, model.physicalModels.size());
	EXPECT_NE(
		model.physicalModels.get().front(),
		newModel.physicalModels.get().front()
	);

	EXPECT_EQ(
		model.physicalModels.get().front()->getName(),
		newModel.physicalModels.get().front()->getName()
	);

	auto newCoordinate1 = newModel.mesh.coords().getId(coordinate::Id(1));
	EXPECT_EQ(
		newCoordinate1,
		newModel.mesh.elems().getId(ElemId(1))->getV(0)
	);

	auto coordinate1 = model.mesh.coords().getId(coordinate::Id(1));
	EXPECT_NE(coordinate1, newCoordinate1);

	EXPECT_EQ(coordinate1->pos(), newCoordinate1->pos());
}

TEST(ModelTest, IsReassigningPhysicalGroupToMeshOnCopy) {
	PMGroup physicalModelsGroup = PMGroup();
	physicalModelsGroup.add(
		std::make_unique<physicalModel::PEC>(
			physicalModel::Id(18),
			"Material PEC"
		)
	);

	CoordR3Group coordinatesGroup = CoordR3Group();
	auto it = coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(),
			math::CVecR3(0.0, 0.0, 0.0)
		)
	);

	const CoordR3* coordinatesArgumentList[1] = { it->get()};
	ElemRGroup elementsGroup = ElemRGroup();
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			element::Id(),
			coordinatesArgumentList,
			nullptr,
			physicalModelsGroup.get().front()
		)
	);

	ModelObject newModel = ModelObject();
	{
		ModelObject model(
			mesh::Unstructured(coordinatesGroup, elementsGroup),
			physicalModelsGroup
		);

		newModel = model;
	}

	EXPECT_FALSE(newModel.physicalModels.empty());
	EXPECT_EQ(
		newModel.mesh.elems().getId(ElemId(1))->getModel()->getId(),
		physicalModel::Id(18)
	);

}

TEST(ModelTest, IsReassigningPhysicalGroupToMeshOnConstruct) {
	PMGroup physicalModelsGroup = PMGroup();
	physicalModelsGroup.add(
		std::make_unique<physicalModel::PEC>(
			physicalModel::Id(18),
			"Material PEC"
		)
	);

	CoordR3Group coordinatesGroup = CoordR3Group();
	auto it = coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			coordinate::Id(),
			math::CVecR3(0.0, 0.0, 0.0)
		)
	);

	const CoordR3* coordinatesArgumentList[1] = { it->get()};
	ElemRGroup elementsGroup = ElemRGroup();
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			element::Id(),
			coordinatesArgumentList,
			nullptr,
			physicalModelsGroup.get().front()
		)
	);

	ModelObject model(
		mesh::Unstructured(coordinatesGroup, elementsGroup), 
		physicalModelsGroup
	);

	physicalModelsGroup = PMGroup();

	EXPECT_FALSE(model.physicalModels.empty());
	EXPECT_TRUE(physicalModelsGroup.empty());
	EXPECT_EQ(
		model.mesh.elems().getId(ElemId(1))->getModel()->getId(),
		physicalModel::Id(18)
	);
}
