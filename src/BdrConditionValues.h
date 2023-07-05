#pragma once

#include <map>
#include <mfem.hpp>

namespace pulmtln {

class BdrConditionValues : public std::map<int,double> {
public:
	BdrConditionValues() = default;
	BdrConditionValues(const std::map<int,double>& attVals) :
		std::map<int,double>{attVals}
	{}
	
	mfem::Array<int> getAttributes() const
	{
		int dbcSize{ (int) size() };
		mfem::Array<int> bi(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bi[i] = it->first;
			++it;
		}
		return bi;
	}

	mfem::Vector getValues() const
	{
		int dbcSize{ (int) size() };
		mfem::Vector bv(dbcSize);
		auto it{ begin() };
		for (int i = 0; i < dbcSize; i++) {
			bv[i] = it->second;
			++it;
		}
		return bv;
	}
};

}