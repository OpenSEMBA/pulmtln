#pragma once

#include <mfem.hpp>
#include <map>

namespace pulmtln {

class Model {
public:
	using Mesh = mfem::Mesh;
	Model() = default;
	Model(
		Mesh& m,
		std::map<int, double> dirichletBCs
	);
	
	Mesh mesh;
	mfem::Array<int> dbcs;
	mfem::Vector dbcv;

};

}