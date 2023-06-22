#pragma once

#include <map>
#include <mfem.hpp>

namespace pulmtln {

class BoundaryConditions {
public:
	BoundaryConditions() = default;
	BoundaryConditions(const std::map<int,double>& attVals) :
		attVals_{attVals}
	{}
	
	mfem::Array<int> getAttributes() const
	{
		int dbcSize{ (int)attVals_.size() };
		mfem::Array<int> bi(dbcSize);
		auto it{ attVals_.begin() };
		for (int i = 0; i < dbcSize; i++) {
			bi[i] = it->first;
			++it;
		}
		return bi;
	}

	mfem::Vector getValues() const
	{
		int dbcSize{ (int)attVals_.size() };
		mfem::Vector bv(dbcSize);
		auto it{ attVals_.begin() };
		for (int i = 0; i < dbcSize; i++) {
			bv[i] = it->second;
			++it;
		}
		return bv;
	}

private:
	std::map<int, double> attVals_;
};

}