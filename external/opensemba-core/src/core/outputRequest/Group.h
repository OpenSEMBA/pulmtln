#pragma once

#include "OutputRequest.h"
#include "core/util/GroupIdentifiableUnique.h"

namespace semba {
namespace outputRequest {

template<typename O = OutputRequest>
class Group : public util::GroupIdentifiableUnique<O> {
};

} 

typedef outputRequest::Group<> OutputRequestGroup;

} 

