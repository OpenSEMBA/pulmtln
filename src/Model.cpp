#include "Model.h"

namespace pulmtln {

MatNameToAttribute Model::getMaterialsOfType(const MaterialType& type) const
{
	MatNameToAttribute res;
	for (const auto& m : materials_.pecs) {
		res[m.name] = m.tag;
	}
	return res;
}

}