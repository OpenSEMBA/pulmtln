#pragma once

#include "SolverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace pulmtln {

class Parser {
public:
	Parser(const std::string& filename);

	Model readModel() const;
	SolverOptions readSolverOptions() const;

private:
	nlohmann::json json_;
	std::string filename_;
};

}