#pragma once

#include "BdrConditionValues.h"
#include "Materials.h"

namespace pulmtln {

enum class MaterialType {
	PEC,
	Dielectric
};

struct Material {};

struct MaterialPEC : public Material {

};

class Materials {
private:
	std::vector<std::unique_ptr<Material>> v_;
};

}