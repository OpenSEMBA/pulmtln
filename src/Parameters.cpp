#include "Parameters.h"

namespace pulmtln {

using namespace mfem;
using json = nlohmann::json;

DenseMatrix Parameters::getCapacitiveCouplingCoefficients() const
{
    DenseMatrix r(C.NumRows(), C.NumCols());
    for (auto i{ 0 }; i < C.NumRows(); ++i) {
        auto selfC{ C(i,i) };
        for (auto j{ 0 }; j < C.NumCols(); ++j) {
            r(i, j) = C(i, j) / selfC;
        }
    }
    return r;
}

json Parameters::toJSON() const
{
    json res;

    std::vector<std::vector<int>> matrix = {
        {1, 2},
        {4, 5}
    };

    // Create an nlohmann::json object and store the matrix
    json jsonMatrix = matrix;

    // Print the JSON object
    std::cout << jsonMatrix.dump(4) << std::endl;

    return res;
}

void Parameters::saveToJSONFile(const std::string& filename) const
{

}

}