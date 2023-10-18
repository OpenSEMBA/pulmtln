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

struct Dielectric : public Material {
	double relativePermittivity;
};

struct Materials { 
	// TODO Convert to class to ensure no names are repeated.
	// TODO This does not scale well... should be converted to a polymorphic container.

	std::vector<PEC> pecs;
	std::vector<OpenBoundary> openBoundaries;
	std::vector<Dielectric> dielectrics;

	template <class T>
	NameToAttrMap buildNameToAttrMapFor() const
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
		else if constexpr (std::is_same<T, Dielectric>()) {
			for (const auto& m : dielectrics) {
				res[m.name] = m.attribute;
			}
		}
		return res;
	}
	NameToAttrMap buildNameToAttrMap() const
	{
		NameToAttrMap res{ buildNameToAttrMapFor<PEC>() };
		{
			auto aux{ buildNameToAttrMapFor<OpenBoundary>() };
			res.insert(aux.begin(), aux.end());
		}
		{
			auto aux{ buildNameToAttrMapFor<Dielectric>() };
			res.insert(aux.begin(), aux.end());
		}
		return res;
	}

	void removeMaterialsNotInList(const NameToAttrMap allowedMaterials)
	{

		auto condition = [&allowedMaterials](const Material& mat) {
			return allowedMaterials.count(mat.name) > 0;
			};

		{
			std::vector<PEC> vs;
			std::copy_if(pecs.begin(), pecs.end(), std::back_inserter(vs), condition);
			pecs = vs;
		}
		{
			std::vector<Dielectric> vs;
			std::copy_if(dielectrics.begin(), dielectrics.end(), std::back_inserter(vs), condition);
			dielectrics = vs;
		}
		{
			std::vector<OpenBoundary> vs;
			std::copy_if(openBoundaries.begin(), openBoundaries.end(), std::back_inserter(vs), condition);
			openBoundaries = vs;
		}
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