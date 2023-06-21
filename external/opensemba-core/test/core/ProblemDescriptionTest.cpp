#include <gtest/gtest.h>

#include "core/geometry/Grid.h"
#include "core/geometry/mesh/Unstructured.h"
#include "core/source/PlaneWave.h"
#include "core/math/function/Gaussian.h"
#include "core/physicalModel/Predefined.h"
#include "core/outputRequest/OnPoint.h"
#include "core/ProblemDescription.h"

using namespace semba;
using namespace geometry;

void fillProblemDescription(UnstructuredProblemDescription& pD)
{

	pD.project = util::ProjectFile("MyPath");
	
	{
		// Create mesh
		CoordR3Group coordinatesGroup;
		coordinatesGroup.copyAndAssignId(CoordR3{ CoordId{1}, { 0.0, 0.0, 0.0 } });
		coordinatesGroup.copyAndAssignId(CoordR3{ CoordId{2}, { 1.0, 2.0, 3.0 } });
	
		ElemRGroup elementsGroup;
		{
			const CoordR3* c[1] = { coordinatesGroup.getId(CoordId(1)) };
			elementsGroup.copyAndAssignId(NodR{ ElemId(1), c });
		}
		{
			const CoordR3* c[1] = { coordinatesGroup.getId(CoordId(2)) };
			elementsGroup.copyAndAssignId(NodR{ ElemId(2), c });
		}
		UnstructuredMesh m{coordinatesGroup, elementsGroup};

		// Create PMGroup.
		PMGroup physicalModelsGroup;
		physicalModelsGroup.copyAndAssignId(physicalModel::PEC{physicalModel::Id(), "Material PEC"});
		
		pD.model = UnstructuredModel{m, physicalModelsGroup};
	}
	
	// Create source
	pD.sources.copyAndAssignId(
		source::PlaneWave{
			std::make_unique<source::Magnitude::Magnitude>(
				new semba::math::function::Gaussian(0.5, 0.0, 1.0)
			),
			{ElemId{1}},
			math::CVecR3(1.0, 0.0, 0.0),
			math::CVecR3(0.0, 0.0, 1.0)
		}
	);

	// Create OutputRequest
	pD.outputRequests.copyAndAssignId(
		outputRequest::OnPoint{
			outputRequest::OutputRequest::Type::electric,
			outputRequest::Domain(),
			"My electric field point probe",
			{ElemId(2)}
		}
	);

}

auto returnProblemDescription()
{
	UnstructuredProblemDescription pD;
	fillProblemDescription(pD);
	return pD;
}

TEST(ProblemDescriptionTest, CanCreate) 
{
	UnstructuredProblemDescription problemDescription;

	EXPECT_NE(&problemDescription, nullptr);
}

TEST(ProblemDescriptionTest, CanInitializeProject) 
{
	const std::string path {"My/Project/Path/File.dat"};

	UnstructuredProblemDescription problemDescription;
	problemDescription.project = util::ProjectFile(path);

	EXPECT_EQ(problemDescription.project, path);
}

TEST(ProblemDescriptionTest, CanInitializeSources) {
	UnstructuredProblemDescription problemDescription;

	SourceGroup sources;

	source::PlaneWave planewave{
		std::make_unique<source::Magnitude::Magnitude>(
			new semba::math::function::Gaussian(0.5, 0.0, 1.0)
		),
		{},
		math::CVecR3(1.0, 0.0, 0.0),
		math::CVecR3(0.0, 0.0, 1.0)
	};

	sources.addAndAssignId(std::make_unique<source::PlaneWave>(planewave));

	problemDescription.sources = sources;

	EXPECT_EQ(problemDescription.grids, Grid3());

	EXPECT_EQ(problemDescription.sources.size(), 1);

	auto sourceInGroup = problemDescription.sources.getId(source::Id(1))->castTo<source::PlaneWave>();
	EXPECT_EQ(sourceInGroup->getDirection(), planewave.getDirection());
	EXPECT_EQ(sourceInGroup->getPolarization(), planewave.getPolarization());
	
	EXPECT_EQ(*sourceInGroup->getMagnitude(), *planewave.getMagnitude());
}

TEST(ProblemDescriptionTest, CanInitializeModel) {
	UnstructuredProblemDescription problemDescription = UnstructuredProblemDescription();

	PMGroup physicalModelsGroup = PMGroup();
	physicalModelsGroup.addAndAssignId(
		std::make_unique<physicalModel::PEC>(
			physicalModel::Id(),
			"Material PEC"
		)
	);

	const auto& physicalModelIt = physicalModelsGroup.addAndAssignId(
		std::make_unique<physicalModel::Bound>(physicalModel::Id(), physicalModel::Bound::Type::pml)
	);

	CoordR3Group coordinatesGroup = CoordR3Group();
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			CoordId(),
			math::CVecR3(0.0, 0.0, 0.0)
		)
	);
	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			CoordId(),
			math::CVecR3(1.0, 2.0, 3.0)
		)
	);

	coordinatesGroup.addAndAssignId(
		std::make_unique<CoordR3>(
			CoordId(),
			math::CVecR3(5.0, 5.0, 5.0)
			)
	);

	const CoordR3* coordinatesArgumentList[1] = { coordinatesGroup.getId(CoordId(1)) };
	ElemRGroup elementsGroup;
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			ElemId(), coordinatesArgumentList
		)
	);

	const CoordR3* coordinatesArgumentList2[1] = { coordinatesGroup.getId(CoordId(2)) };
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			ElemId(),
			coordinatesArgumentList2
		)
	);

	// Boundaries
	const CoordR3* coordinatesArgumentBoundaryList[1] = { coordinatesGroup.getId(CoordId(3)) };
	elementsGroup.addAndAssignId(
		std::make_unique<NodR>(
			ElemId(),
			coordinatesArgumentBoundaryList,
			nullptr,
			physicalModelIt->get()
		)
	);

	const UnstructuredModel model{
		mesh::Unstructured(coordinatesGroup, elementsGroup),
		physicalModelsGroup
	};

	problemDescription.model = model;

	ASSERT_FALSE(problemDescription.model.physicalModels.empty());
	EXPECT_EQ("Material PEC", problemDescription.model.physicalModels.get()[0]->getName());
	EXPECT_EQ("PML_Bound", problemDescription.model.physicalModels.get()[1]->getName());

	EXPECT_FALSE(problemDescription.model.mesh.coords().empty());
	EXPECT_EQ(
		math::CVecR3(1.0, 2.0, 3.0),
		(problemDescription.model.mesh.coords().get()[1])->pos()
	);

	const auto& element = problemDescription.model.mesh.elems().getId(ElemId(3));

	EXPECT_EQ(
		element->getMatId(),
		problemDescription.model.physicalModels.get()[1]->getId()
	);

	EXPECT_EQ(
		element->getCoordinates()[0]->pos(),
		math::CVecR3(5.0, 5.0, 5.0)
	);
}

TEST(ProblemDescriptionTest, CanCopyConstructor) 
{
	UnstructuredProblemDescription pD;
	fillProblemDescription(pD);
	
	UnstructuredProblemDescription copy{ pD };

	EXPECT_EQ(copy.project, pD.project);
}

TEST(ProblemDescriptionTest, CanAssign)
{
	UnstructuredProblemDescription pD;
	fillProblemDescription(pD);

	UnstructuredProblemDescription assignment;
	
	assignment = pD;

	EXPECT_EQ(assignment.project, pD.project);
}