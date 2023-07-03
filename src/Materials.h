#pragma once

#include "BdrConditionValues.h"
#include "Materials.h"

namespace pulmtln {

enum class MaterialType {
	PEC,
	Dielectric
};

struct MaterialPEC {
	std::string name;
	int tag;
};

struct Materials {
	std::vector<MaterialPEC> pecs;
};

}