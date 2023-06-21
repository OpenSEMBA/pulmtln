#pragma once

#include "PhysicalModel.h"
#include "core/util/GroupIdentifiableUnique.h"

namespace semba {
namespace physicalModel {

template<typename P = PhysicalModel>
class Group : public util::GroupIdentifiableUnique<P> {
};

} 

typedef physicalModel::Group<> PMGroup;

} 

