#pragma once

#include "BdrConditionValues.h"
#include "Materials.h"

namespace pulmtln {

struct Material {
	std::string name;
	int tag;
};

struct PEC : public Material {
};

struct Dielectric : public Material {
	double relativePermittivity;
};

struct Materials {
	std::vector<PEC> pecs;
	std::vector<Dielectric> dielectrics;

	template <class T>
	MatNameToAttribute getMatNameToAttributeMap() const
	{
		MatNameToAttribute res;
		if constexpr (std::is_same<T, PEC>()) {
			for (const auto& m : pecs) {
				res[m.name] = m.tag;
			}
		}
		else if constexpr (std::is_same<T, Dielectric>()) {
			for (const auto& m : dielectrics) {
				res[m.name] = m.tag;
			}
		}
		return res;
	}
};

}