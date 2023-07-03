#pragma once

#include "SolverOptions.h"
#include "Model.h"

#include <nlohmann/json.hpp>

namespace pulmtln {

class Parser {
public:
	Parser(const nlohmann::json&);

	Model readModel() const;
	SolverOptions readSolverOptions() const;

private:
	nlohmann::json json_;
};

}