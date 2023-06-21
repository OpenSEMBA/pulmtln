#pragma once

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include "core/geometry/mesh/Geometric.h"
#include "core/physicalModel/Predefined.h"

#include "parsers/Parser.h"

namespace semba::parsers::STL {

class Parser : public semba::parsers::Parser {
public:
    Parser(const std::string& fn);
    geometry::mesh::Unstructured readAsUnstructuredMesh() const;
};

}
