#include "Results.h"

#include "constants.h"

namespace pulmtln {

using namespace mfem;
using json = nlohmann::json;

std::vector<std::vector<double>> toVecVec(const DenseMatrix& m)
{
    std::vector<std::vector<double>> r(m.NumRows());
    for (auto i{ 0 }; i < m.NumRows(); ++i) {
        r[i].resize(m.NumCols());
        for (auto j{ 0 }; j < m.NumCols(); ++j) {
            r[i][j] = m(i, j);
        }
    }

    return r;
}

DenseMatrix toDenseMatrix(const std::vector<std::vector<double>>& v)
{
    
    if (v.size() == 0 || v[0].size() == 0) {
        return {};
    }

    for (const auto& row: v) {
        if (row.size() != v[0].size()) {
            throw std::runtime_error("vec-vec matrix is not square.");
        }
    }
    
    DenseMatrix r(int(v.size()), int(v[0].size()));
    for (auto i{ 0 }; i < r.NumRows(); ++i) {
        for (auto j{ 0 }; j < r.NumCols(); ++j) {
            r(i, j) = v[i][j];
        }
    }

    return r;
}

PULParameters::PULParameters(const json& j)
{
    if (!j.contains("C") || !j.contains("L")) {
        throw std::runtime_error("JSON does not contain C and/or L.");
    }
    C = toDenseMatrix(j["C"]);
    L = toDenseMatrix(j["L"]);
}

bool PULParameters::operator==(const PULParameters& rhs) const
{ 
    return
        toVecVec(C) == toVecVec(rhs.C) &&
        toVecVec(L) == toVecVec(rhs.L);
}

DenseMatrix PULParameters::getCapacitiveCouplingCoefficients() const
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

json PULParameters::toJSON() const
{
    json res;
    res["C"] = toVecVec(C);
    res["L"] = toVecVec(L);
    return res;
}

void PULParameters::saveToJSONFile(const std::string& filename) const
{
    std::ofstream ofs(filename);
    auto j{ toJSON() };
    ofs << j;
}

double InCellPotentials::getCapacitanceUsingInnerRegion(int i, int j) const
{
    double Qj = electric.at(j).ab[0].first;
    double avVj = electric.at(j).innerRegionAveragePotential;
    double ViWhenPrescribedVj = electric.at(j).conductorPotentials.at(i);
    avVj = -avVj + ViWhenPrescribedVj;
    return Qj / avVj * EPSILON0_SI;
}

double InCellPotentials::getInductanceUsingInnerRegion(int i, int j) const
{
    double Ij = magnetic.at(j).ab[0].first;
    double avAj = magnetic.at(j).innerRegionAveragePotential;
    double AiWhenPrescribedAj = magnetic.at(j).conductorPotentials.at(i);
    avAj = -avAj + AiWhenPrescribedAj;
    return avAj / Ij * MU0_SI;
}

double InCellPotentials::getCapacitanceOnBox(int i, int j, const Box& box) const
{
    return 0.0;
}

double InCellPotentials::getInductanceOnBox(int i, int j, const Box& box) const
{
    return 0.0;
}

}