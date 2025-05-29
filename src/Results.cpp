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

std::array<std::vector<double>, 2> buildIntegrationPlanesForBox(
    const Box& integrationBox, 
    const Box& innerRegionBox) 
{
    const int GRID_INTEGRATION_SAMPLING_POINTS = 100;

    for (int x = 0; x < 2; ++x) {
        if ((integrationBox.min(x) > innerRegionBox.min(x)) ||
            (integrationBox.max(x) < innerRegionBox.max(x))) {
            throw std::runtime_error("Integration box has to be larger than the  inner region.");
        }
    }

    std::array<std::set<double>, 2> planes;
    for (int x = 0; x < 2; ++x) {
        std::set<double> controlPoints;
        controlPoints.insert(integrationBox.min(x));
        controlPoints.insert(integrationBox.max(x));
        controlPoints.insert(innerRegionBox.min(x));
        controlPoints.insert(innerRegionBox.max(x));

        // Fills with equispaced points
        for (auto it = std::next(controlPoints.begin()); it != controlPoints.end(); ++it) {
            auto prev = std::prev(it);
            const double step = (*it - *prev) / GRID_INTEGRATION_SAMPLING_POINTS;
            for (int k = 0; k < GRID_INTEGRATION_SAMPLING_POINTS; ++k) {
                const double newPoint = *prev + k * step;
                planes[x].insert(newPoint);
            }
        }
    }

	std::array<std::vector<double>, 2> res;
	for (int x = 0; x < 2; ++x) {
		res[x].assign(planes[x].begin(), planes[x].end());
	}
	return res;
}

double InCellPotentials::getCapacitanceOnBox(int i, int j, const Box& box) const
{
	auto integrationPlanes{ buildIntegrationPlanesForBox(box, innerRegionBox) };
    double outerV = 0.0;
    for (int m = 1; m < integrationPlanes[0].size(); ++m) {
        for (int n = 1; n < integrationPlanes[1].size(); ++n) {
			mfem::Vector midPoint({
				(integrationPlanes[0][m] + integrationPlanes[0][m - 1]) / 2.0,
				(integrationPlanes[1][n] + integrationPlanes[1][n - 1]) / 2.0
			});
            if (innerRegionBox.isWithinBox(midPoint)) {
				// If the point is within the inner region, we do not integrate over it.
				continue;
            }
			double area = (integrationPlanes[0][m] - integrationPlanes[0][m - 1]) *
				(integrationPlanes[1][n] - integrationPlanes[1][n - 1]);
			outerV += area * multipolarExpansion(midPoint, electric.at(j).ab, electric.at(j).expansionCenter);
        }
    }

    double innerV = electric.at(j).innerRegionAveragePotential * innerRegionBox.area();

	double Qj = electric.at(j).ab[0].first;
	double avVj = (innerV + outerV) / box.area();
	double ViWhenPrescribedVj = electric.at(j).conductorPotentials.at(i);
	avVj = -avVj + ViWhenPrescribedVj;
	
	return Qj / avVj * EPSILON0_SI ;
}

double InCellPotentials::getInductanceOnBox(int i, int j, const Box& box) const
{
    throw std::runtime_error("Not implemented.");
    return 0.0;
}

}