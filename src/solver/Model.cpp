#include "Model.h"

using namespace mfem;

namespace pulmtln {

void convertToArrayAndVector(
	mfem::Array<int>& bi, 
	mfem::Vector& bv,
	const std::map<int,double>& bcs)
{
	int dbcSize{ (int) bcs.size() };
	bi = mfem::Array<int>(dbcSize);
	bv = mfem::Vector(dbcSize);
	auto it{ bcs.begin() };
	for (int i = 0; i < dbcSize; i++) {
		bi[i] = it->first;
		bv[i] = it->second;
		++it;
	}
}
Model::Model(
	Mesh& m,
	std::map<int, double> dirichletBCs,
	std::map<int, double> neumannBCs) :
	mesh(m)
{
	convertToArrayAndVector(dbcs, dbcv, dirichletBCs);
	convertToArrayAndVector(nbcs, nbcv, neumannBCs);
};
}