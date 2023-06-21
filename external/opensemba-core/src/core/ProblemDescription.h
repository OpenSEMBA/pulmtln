#pragma once

#include <nlohmann/json.hpp>

#include "core/geometry/Grid.h"
#include "core/physicalModel/Bound.h"
#include "core/source/Group.h"
#include "core/outputRequest/Group.h"
#include "core/model/Model.h"
#include "core/util/ProjectFile.h"

namespace semba {

template<typename M = UnstructuredModel>
class ProblemDescriptionBase {
public:
	util::ProjectFile project;
	geometry::Grid3 grids;
	SourceGroup sources;
	nlohmann::json analysis;
	M model;
	OutputRequestGroup outputRequests;
};

typedef ProblemDescriptionBase<> UnstructuredProblemDescription;
typedef ProblemDescriptionBase<StructuredModel> StructuredProblemDescription;

}