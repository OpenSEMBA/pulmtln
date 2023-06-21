#include "Parser.h"

namespace semba {
namespace parsers {

Parser::Parser(const std::string& fn) :
    filename(fn)
{
    std::ifstream ifs(fn);
    if (!ifs.is_open()) {
        throw std::runtime_error("Unable to open file: " + fn);
    }
};

bool strToBool(const std::string& value) {
    if (atoi(value.c_str()) == 1) {
        return true;
    } else {
        return false;
    }
}

std::string& trim(std::string& s) {
    return ltrim(rtrim(s));
}

std::string& ltrim(std::string& s) {
    s.erase(s.begin(),
        std::find_if(s.begin(), s.end(),
            [](int c) {return !std::isspace(c); })
    );
    return s;
}

std::string& rtrim(std::string& s) {
    s.erase(find_if(s.rbegin(), s.rend(),
        [](int c) {return !std::isspace(c); }).base(),
        s.end());
    return s;
}


bool toBool(const std::size_t param) {
    assert(param == 0 || param == 1);
    if (param == 1) {
        return true;
    }
    else {
        return false;
    }
}

void Parser::postReadOperations(UnstructuredProblemDescription& res)
{
    if (res.analysis.find("geometryScalingFactor") != res.analysis.end()) {
        math::Real scalingFactor{ 
            res.analysis.at("geometryScalingFactor").get<double>() };
        res.model.mesh.applyScalingFactor(scalingFactor);
        res.grids.applyScalingFactor(scalingFactor);
    }
}

}
} 
