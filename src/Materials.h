#pragma once

#include "AttrToValueMap.h"
#include "Materials.h"

namespace pulmtln {

using MaterialId = int;
using Attribute = int;

struct Material {
	std::string name;
	Attribute attribute;
};

struct PEC : public Material {
};

struct OpenBoundary : public Material {
};

struct Vacuum : public Material {
};

struct Dielectric : public Material {
	double relativePermittivity;
};

struct Materials {
	std::vector<PEC> pecs;
	std::vector<OpenBoundary> openBoundaries;
	std::vector<Dielectric> dielectrics;
	std::vector<Vacuum> vacuums;

	template <class T>
	NameToAttrMap buildNameToAttrMap() const
	{
		NameToAttrMap res;
		if constexpr (std::is_same<T, PEC>()) {
			for (const auto& m : pecs) {
				res[m.name] = m.attribute;
			}
		}
		if constexpr (std::is_same<T, OpenBoundary>()) {
			for (const auto& m : openBoundaries) {
				res[m.name] = m.attribute;
			}
		}
		if constexpr (std::is_same<T, Vacuum>()) {
			for (const auto& m : vacuums) {
				res[m.name] = m.attribute;
			}
		}
		else if constexpr (std::is_same<T, Dielectric>()) {
			for (const auto& m : dielectrics) {
				res[m.name] = m.attribute;
			}
		}
		return res;
	}

	static int getNumberContainedInName(const std::string& name)
	{
		std::stringstream ss{ name.substr(name.find("_") + 1) };
		int res;
		ss >> res;
		return res;
	}
};

}