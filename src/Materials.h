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


class Materials {
private:
	std::vector<MaterialPEC> PEC_;
};

}