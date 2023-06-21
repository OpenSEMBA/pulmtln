#pragma once

#include <vector>
#include <algorithm>
#include <array>

#include "core/geometry/Box.h"
#include "core/math/Real.h"

namespace semba {
namespace geometry {

template<class T, std::size_t D> class Box;

template<std::size_t D>
class Grid {
    typedef Box<math::Real,D> BoxRD;
    typedef Box<math::Int ,D> BoxID;
    typedef math::CartesianVector<math::Real,D> CVecRD;
    typedef math::CartesianVector<math::Int, D> CVecID;
public:
    static const math::Real tolerance;

    Grid() = default;
    Grid(const BoxRD&  boundingBox,
         const CVecRD& dxyz);
    Grid(const BoxRD&  boundingBox,
         const CVecID& dims);
    Grid(const std::vector<math::Real> positions[D]);
    Grid(const Grid& grid);
    virtual ~Grid () = default;


    bool operator==(const Grid& rhs) const;
    void setPos(const std::vector<math::Real> pos[D]);
    void setAdditionalSteps(const math::Constants::CartesianAxis d,
                            const math::Constants::CartesianBound b,
                            const std::vector<math::Real>& step);

    bool hasZeroSize() const;

    bool isInto(const CVecRD& pos) const;
    bool isInto(const std::size_t dir, const math::Real pos) const;
    bool isRegular() const;
    bool isRegular(const std::size_t d) const;
    bool isCartesian() const;
    bool isCell(const CVecRD& position,
                const math::Real tol = tolerance) const;
    bool isCell(const std::vector<CVecRD>& positions,
                const math::Real tol = tolerance) const;

    CVecID getNumCells() const;

    CVecID getOffset()   const; // DEPRECATED

    CVecRD getOrigin()   const;
    bool getNaturalCell(
            const math::Constants::CartesianAxis dir,
            const math::Real& x,
            long int& i,
            math::Real& relativeLen) const;

    std::vector<math::Real> getStep(const std::size_t dir) const;
    math::Real              getStep(const std::size_t dir,
                                    const math::Int& n) const;

    math::Real getMinimumSpaceStep() const;

    BoxRD getFullDomainBoundingBox() const;
    BoxID getFullDomainBoundingCellBox() const;
    BoxRD getBoundingBox(const BoxID& bound) const;
    BoxRD getBoxRContaining(const CVecRD& point) const;
    BoxID getBoxIContaining(const CVecRD& point) const;

    std::vector<CVecRD> getCenterOfCellsInside(const BoxRD& bound) const;
    std::vector<math::Real> getPosInRange(const std::size_t dir,
                                          const math::Real min,
                                          const math::Real max) const;

    std::vector<CVecRD>     getPos() const;
    std::vector<math::Real> getPos(const std::size_t dir) const;
    math::Real              getPos(const std::size_t dir,
                                   const math::Int i) const;
    CVecRD                  getPos(const CVecID& ijk) const;
    CVecRD                  getPos(const CVecRD& ijk) const { return ijk; }

    std::pair<math::Int, math::Real> getCellPair(
            const std::size_t       dir,
            const math::Real x,
            const bool approx = true,
            const math::Real tol = tolerance,
            bool* err = nullptr) const;
    std::pair<CVecID, CVecRD>        getCellPair(
            const CVecRD& pos,
            const bool approx = true,
            const math::Real tol = tolerance,
            bool* err = nullptr) const;

    math::Int getCell(const std::size_t dir,
                      const math::Real  x,
                      const bool  approx = true,
                      const math::Real  tol = tolerance,
                      bool* err = nullptr) const;
    CVecID    getCell(const CVecRD& pos,
                      const bool  approx = true,
                      const math::Real tol = tolerance,
                      bool* err = nullptr) const;
    CVecID    getCell(const CVecID& pos,
                      const bool approx = true,
                      const math::Real tol = tolerance) const { return pos; }

    void applyScalingFactor(const math::Real factor);

    void enlarge(const std::pair<CVecRD,CVecRD>& additionalCells,
                 const std::pair<CVecRD,CVecRD>& sizesOfNewCells);
    void enlargeBound(math::Constants::CartesianAxis d,
                      math::Constants::CartesianBound b,
                      math::Real pad, math::Real siz);

private:
    std::array<std::vector<math::Real>,D> pos_;
};


template<std::size_t D>
const math::Real Grid<D>::tolerance = 1e-2;

template<std::size_t D>
Grid<D>::Grid(const BoxRD& box,
    const CVecRD& dxyz) {
    CVecRD origin = box.getMin();
    for (std::size_t i = 0; i < D; i++) {
        math::Real boxLength = box.getMax()(i) - box.getMin()(i);
        std::size_t nCells;
        if (dxyz(i) == (math::Real)0.0) {
            nCells = 1;
        }
        else {
            nCells = math::ceil(boxLength / dxyz(i));
            if (math::greater(boxLength, nCells * dxyz(i),
                dxyz(i), tolerance)) {
                nCells++;
            }
        }
        pos_[i].resize(nCells + 1);
        for (std::size_t j = 0; j < nCells + 1; j++) {
            pos_[i][j] = origin(i) + j * dxyz(i);
        }
    }
}

template<std::size_t D>
Grid<D>::Grid(const BoxRD& boundingBox,
    const CVecID& dims) {
    CVecRD origin = boundingBox.getMin();
    for (std::size_t i = 0; i < D; i++) {
        math::Real step =
            (boundingBox.getMax()(i) - boundingBox.getMin()(i)) / dims(i);
        std::size_t nCells = dims(i);
        pos_[i].resize(nCells + 1);
        for (std::size_t j = 0; j < nCells + 1; j++) {
            pos_[i][j] = origin(i) + j * step;
        }
    }
}

template<std::size_t D>
Grid<D>::Grid(const std::vector<math::Real> pos[D]) {
    for (std::size_t d = 0; d < D; d++) {
        pos_[d] = pos[d];
    }
}


template<std::size_t D>
Grid<D>::Grid(const Grid<D>& grid) {
    for (std::size_t i = 0; i < D; i++) {
        pos_[i] = grid.pos_[i];
    }
}

template<std::size_t D>
bool Grid<D>::operator==(const Grid& rhs) const {
    return this->pos_ == rhs.pos_;
}

template<std::size_t D>
void Grid<D>::setPos(const std::vector<math::Real> pos[D]) {
    for (std::size_t d = 0; d < D; d++) {
        if (pos[d].size() == 0) {
            throw std::out_of_range(
                "Grid positions must contain at least one value");
        }
        pos_[d] = pos[d];
        if (pos_[d].size() == 1) {
            pos_[d].push_back(pos_[d][0]);
        }
    }
}

template<std::size_t D>
void Grid<D>::setAdditionalSteps(
    const math::Constants::CartesianAxis d,
    const math::Constants::CartesianBound b,
    const std::vector<math::Real>& step) {
    const std::size_t nCells = step.size();
    std::vector<math::Real> newPos(nCells);
    if (b == math::Constants::U) {
        newPos[0] = pos_[d].back() + step[0];
        for (std::size_t i = 1; i < nCells; i++) {
            newPos[i] = newPos[i - 1] + step[i];
        }
        pos_[d].insert(pos_[d].end(), newPos.begin(), newPos.end());
    }
    else {
        newPos[0] = pos_[d].front() - step[0];
        for (std::size_t i = 1; i < nCells; i++) {
            newPos[i] = newPos[i - 1] - step[i];
        }
        std::reverse(newPos.begin(), newPos.end());
        newPos.insert(newPos.end(), pos_[d].begin(), pos_[d].end());
        pos_[d] = newPos;
    }
}

template<std::size_t D>
bool Grid<D>::hasZeroSize() const {
    bool res = true;
    for (std::size_t i = 0; i < D; i++) {
        res &= (pos_[i].size() <= 1);
    }
    return res;
}

template<std::size_t D>
bool Grid<D>::isInto(const std::size_t dir, const math::Real pos) const {
    return pos >= getPos(dir)[0] && pos <= getPos(dir).back();
}

template<std::size_t D>
bool Grid<D>::getNaturalCell(
    const math::Constants::CartesianAxis dir,
    const math::Real& x,
    long int& i,
    math::Real& relativeLen) const {
    size_t n = 0;
    relativeLen = -1.0;
    if (x < getPos(dir, 0)) {
        i = 0;
        return false;
    }
    else if (getPos(dir, getNumCells()(dir)) <= x) {
        i = getNumCells()(dir);
        return false;
    }
    while (getPos(dir)[n] <= x) {
        n++;
    }  /*mod this: use sort*/
    i = n - 1;
    relativeLen = (x - getPos(dir)[i]) /
        getStep(dir, i);
    return true;
}

template<std::size_t D>
bool Grid<D>::isInto(const CVecRD& pos) const {
    for (std::size_t i = 0; i < D; i++) {
        if (!isInto(i, pos(i))) {
            return false;
        }
    }
    return true;
}

template<std::size_t D>
bool Grid<D>::isRegular() const {
    for (std::size_t i = 0; i < D; i++) {
        if (!isRegular(i)) {
            return false;
        }
    }
    return true;
}

template<std::size_t D>
bool Grid<D>::isRegular(const std::size_t d) const {
    std::vector<math::Real> step = getStep(d);
    for (std::size_t n = 1; n < step.size(); n++) {
        if (math::notEqual(step[n], step[0], step[0], tolerance)) {
            return false;
        }
    }
    return true;
}

template<std::size_t D>
bool Grid<D>::isCartesian() const {
    math::Real canon = getStep(math::Constants::x)[0];
    for (std::size_t i = 0; i < D; i++) {
        std::vector<math::Real> step = getStep(i);
        for (std::size_t n = 1; n < step.size(); n++) {
            if (math::notEqual(step[n], canon, canon, tolerance)) {
                return false;
            }
        }
    }
    return true;
}

template<std::size_t D>
bool Grid<D>::isCell(const CVecRD& position, const math::Real tol) const {
    std::pair<CVecID, CVecRD> natCell = getCellPair(position, true, tol);
    return natCell.second == CVecRD(0.0);
}

template<std::size_t D>
bool Grid<D>::isCell(const std::vector<CVecRD>& pos,
    const math::Real tol) const {
    for (std::size_t i = 0; i < pos.size(); i++) {
        if (!isCell(pos[i], tol)) {
            return false;
        }
    }
    return true;
}

template<std::size_t D>
math::CartesianVector<math::Int, D> Grid<D>::getNumCells() const {
    CVecID res;
    for (std::size_t d = 0; d < D; d++) {
        res(d) = getPos(d).size() - 1; // Minimum size of pos is 2.
    }
    return res;
}

template<std::size_t D>
math::CartesianVector<math::Int, D> Grid<D>::getOffset() const {
    return CVecID(0, 0, 0);
}

template<std::size_t D>
math::CartesianVector<math::Real, D> Grid<D>::getOrigin() const {
    CVecRD res;
    for (std::size_t d = 0; d < D; d++) {
        if (pos_[d].size() == 0) {
            throw std::out_of_range("Positions are not initialized.");
        }
        res(d) = pos_[d][0];
    }
    return res;
}

template<std::size_t D>
std::vector<math::Real> Grid<D>::getStep(const std::size_t dir) const {
    assert(dir >= 0 && dir < D);
    if (pos_[dir].size() == 0) {
        return std::vector<math::Real>();
    }
    std::vector<math::Real> res(pos_[dir].size() - 1);
    for (std::size_t i = 0; i < pos_[dir].size() - 1; i++) {
        res[i] = pos_[dir][i + 1] - pos_[dir][i];
    }
    return res;
}


template<std::size_t D>
math::Real Grid<D>::getStep(const std::size_t dir, const math::Int& n) const {
    assert(dir >= 0 && dir < D);
    assert(n >= 0 && n < (math::Int(pos_[dir].size()) - 1));

    if (pos_[dir].empty()) {
        return 0.0;
    }
    return pos_[dir][n + 1] - pos_[dir][n];
}



template<std::size_t D>
math::Real Grid<D>::getMinimumSpaceStep() const {
    math::Real minStep = std::numeric_limits<math::Real>::infinity();
    for (std::size_t i = 0; i < D; i++) {
        std::vector<math::Real> step = getStep(i);
        for (std::size_t n = 0; n < step.size(); n++) {
            if (step[n] < minStep) {
                minStep = step[n];
            }
        }
    }
    return minStep;
}

template<std::size_t D>
Box<math::Real, D> Grid<D>::getFullDomainBoundingBox() const {
    return getBoundingBox(
        std::pair<CVecID, CVecID>(CVecID(0, 0, 0), getNumCells()));
}

template<std::size_t D>
Box<math::Int, D> Grid<D>::getFullDomainBoundingCellBox() const {
    CVecID min, max, dims;
    for (std::size_t n = 0; n < D; n++) {
        dims[n] = pos_[n].size();
    }

    return BoxID(CVecID(0, 0, 0), dims);
}

template<std::size_t D>
Box<math::Real, D> Grid<D>::getBoundingBox(const BoxID& bound) const {
    BoxRD res(getPos(bound.getMin()), getPos(bound.getMax()));
    return res;
}

template<std::size_t D>
Box<math::Real, D> Grid<D>::getBoxRContaining(const CVecRD& point) const {
    BoxID boxI = getBoxIContaining(point);
    return getBoundingBox(boxI);
}


template<std::size_t D>
Box<math::Int, D> Grid<D>::getBoxIContaining(const CVecRD& point) const {
    CVecID min = getCell(point, false);
    CVecID max = min + 1;
    return BoxID(min, max);
}

template<std::size_t D>
std::vector< math::CartesianVector<math::Real, D> >
Grid<D>::getCenterOfCellsInside(const BoxRD& bound) const {
    // Determines centers of cells.
    std::vector<math::Real> center[D];
    for (std::size_t dir = 0; dir < D; dir++) {
        std::vector<math::Real> pos = getPosInRange(dir,
            bound.getMin()(dir),
            bound.getMax()(dir));
        if (pos.size() > 0) {
            center[dir].reserve(pos.size() - 1);
        }
        for (std::size_t i = 1; i < pos.size(); i++) {
            math::Real auxCenter = (pos[i - 1] + pos[i]) / 2.0;
            center[dir].push_back(auxCenter);
        }
    }
    // Combines centers in a std::vector of CVecRD positions.
    std::vector<CVecRD> res;
    res.reserve(center[math::Constants::x].size() *
        center[math::Constants::y].size() *
        center[math::Constants::z].size());
    for (std::size_t i = 0; i < center[math::Constants::x].size(); i++) {
        for (std::size_t j = 0; j < center[math::Constants::y].size(); j++) {
            for (std::size_t k = 0; k < center[math::Constants::z].size();
                k++) {
                res.push_back(CVecRD(center[math::Constants::x][i],
                    center[math::Constants::y][j],
                    center[math::Constants::z][k]));
            }
        }
    }
    return res;
}

template<std::size_t D>
std::vector<math::Real> Grid<D>::getPosInRange(const std::size_t dir,
    const math::Real min,
    const math::Real max) const {
    std::vector<math::Real> pos = getPos(dir);
    std::vector<math::Real> steps = getStep(dir);
    std::vector<math::Real> res;
    res.reserve(pos.size());
    for (std::size_t i = 0; i < pos.size(); i++) {
        math::Real step;
        if (i < steps.size()) {
            step = steps[i];
        }
        else {
            step = steps.back();
        }
        const bool inMin = math::equal(pos[i], min, step, tolerance);
        const bool inMax = math::equal(pos[i], max, step, tolerance);
        const bool inRange = (pos[i] >= min && pos[i] <= max);
        if (inMin || inMax || inRange) {
            res.push_back(pos[i]);
        }
    }
    return res;
}

template<std::size_t D>
std::vector< math::CartesianVector<math::Real, D> > Grid<D>::getPos() const {
    // Combines positions in a std::vector of CVecRD positions.
    std::vector<CVecRD> res;
    res.reserve(pos_[math::Constants::x].size() *
        pos_[math::Constants::y].size() *
        pos_[math::Constants::z].size());
    for (std::size_t i = 0; i < pos_[math::Constants::x].size(); i++) {
        for (std::size_t j = 0; j < pos_[math::Constants::y].size(); j++) {
            for (std::size_t k = 0; k < pos_[math::Constants::z].size(); k++) {
                res.push_back(CVecRD(pos_[math::Constants::x][i],
                    pos_[math::Constants::y][j],
                    pos_[math::Constants::z][k]));
            }
        }
    }
    return res;
}

template<std::size_t D>
std::vector<math::Real> Grid<D>::getPos(const std::size_t direction) const {
    assert(direction >= 0 && direction < D);
    return pos_[direction];
};

template<std::size_t D>
math::CartesianVector<math::Real, D> Grid<D>::getPos(
    const CVecID& ijk) const {
    CVecRD res;
    for (std::size_t i = 0; i < D; i++) {
        res(i) = pos_[i][ijk(i)];
    }
    return res;
};

template<std::size_t D>
math::Real Grid<D>::getPos(const std::size_t dir, const math::Int i) const {
    return  pos_[dir][i];
}

template<std::size_t D>
std::pair<math::Int, math::Real> Grid<D>::getCellPair(const std::size_t dir,
    const math::Real x,
    const bool approx,
    const math::Real tol,
    bool* err) const {
    if (err != nullptr) {
        *err = false;
    }

    math::Int  cell;
    math::Real dist;
    std::vector<math::Real> pos = getPos(dir);
    std::vector<math::Real> steps = getStep(dir);
    assert(pos_[dir].size() >= 1);
    // Checks if it is below the grid.
    if (math::lower(x, pos[0], steps[0], tol)) {
        cell = 0;
        dist = (x - pos[0]) / steps[0];
        if (err != nullptr) {
            *err = true;
        }
        return std::make_pair(cell, dist);
    }
    for (std::size_t i = 0; i < pos.size(); i++) {
        math::Real step;
        if (i == 0) {
            step = steps.front();
        }
        else if (i <= steps.size()) {
            step = steps[i - 1];
        }
        else {
            step = steps.back();
        }
        if (math::equal(x, pos[i], step, tol)) {
            cell = i;
            dist = 0.0;
            if (err != nullptr) {
                *err = false;
            }
            return std::make_pair(cell, dist);
        }
        else if (math::lower(x, pos[i], step, tol)) {
            cell = i - 1;
            dist = (x - pos[i - 1]) / step;
            if (math::equal(math::round(dist), 1.0) && approx) {
                cell++;
                dist -= 1.0;
            }
            if (err != nullptr) {
                *err = false;
            }
            return std::make_pair(cell, dist);
        }
    }
    cell = getNumCells()(dir);
    dist = (x - pos.back()) / steps.back();
    if (err != nullptr) {
        *err = true;
    }
    return std::make_pair(cell, dist);
}

template<std::size_t D>
std::pair<math::CartesianVector<math::Int, D>,
    math::CartesianVector<math::Real, D>>
    Grid<D>::getCellPair(const CVecRD& xyz,
        const bool approx,
        const math::Real tol,
        bool* err) const {
    if (err != nullptr) {
        *err = false;
    }
    bool stepErr = false;

    CVecID cell;
    CVecRD dist;
    for (std::size_t dir = 0; dir < D; dir++) {
        std::pair<math::Int, math::Real> res =
            getCellPair(dir, xyz(dir), approx, tol, &stepErr);
        cell(dir) = res.first;
        dist(dir) = res.second;
        if (err != nullptr) {
            *err = *err || stepErr;
        }
    }
    return std::make_pair(cell, dist);
}

template<std::size_t D>
math::Int Grid<D>::getCell(const std::size_t dir,
    const math::Real x,
    const bool approx,
    const math::Real tol,
    bool* err) const {
    return getCellPair(dir, x, approx, tol, err).first;
}

template<std::size_t D>
math::CartesianVector<math::Int, D> Grid<D>::getCell(const CVecRD& coords,
    const bool approx,
    const math::Real tol,
    bool* err) const {
    return getCellPair(coords, approx, tol, err).first;
}

template<std::size_t D>
void Grid<D>::applyScalingFactor(const math::Real factor) {
    for (std::size_t d = 0; d < D; d++) {
        for (std::size_t i = 0; i < pos_[d].size(); i++) {
            pos_[d][i] *= factor;
        }
    }
}

template<std::size_t D>
void Grid<D>::enlarge(const std::pair<CVecRD, CVecRD>& pad,
    const std::pair<CVecRD, CVecRD>& sizes) {
    for (std::size_t d = 0; d < D; d++) {
        for (std::size_t b = 0; b < 2; b++) {
            if (b == math::Constants::L) {
                enlargeBound(math::Constants::CartesianAxis(d),
                    math::Constants::CartesianBound(b),
                    pad.first(d), sizes.first(d));
            }
            else {
                enlargeBound(math::Constants::CartesianAxis(d),
                    math::Constants::CartesianBound(b),
                    pad.second(d), sizes.second(d));
            }
        }
    }
}

template<std::size_t D>
void Grid<D>::enlargeBound(math::Constants::CartesianAxis d,
    math::Constants::CartesianBound b,
    math::Real pad, math::Real siz) {
    assert(getNumCells()(d) > 0);
    if (std::abs(siz) > std::abs(pad)) {
        std::cerr << "WARNING @ Grid enlarge bound: "
            << "std::size_t was larger than padding. Ignoring padding in "
            << "axis " << d << " and bound " << b << "." << std::endl;
        return;
    }
    if (pad == 0.0) {
        return;
    }
    math::Int boundCell;
    if (b == math::Constants::L) {
        boundCell = 0;
    }
    else {
        boundCell = this->getNumCells()(d) - 1;
    }
    std::vector<math::Real> newSteps;
    if (math::greaterEqual(getStep(d, b), siz) || siz == 0.0) {
        siz = getStep(d, boundCell);
        // Computes enlargement for a padding with same size.
        math::Real nCellsFrac = std::abs(pad / siz);
        const math::Real tol = 0.01;
        std::size_t nCells = (std::size_t)math::ceil(nCellsFrac, tol);
        newSteps.resize(nCells, siz);
    }
    else {
        // Computes enlargement for padding with different size.
        // Taken from AutoCAD interface programmed in LISP (2001).
        math::Real d12 = getStep(d, boundCell);
        math::Real d14 = std::abs(pad) + d12 + std::abs(siz);
        math::Real d34 = std::abs(siz);
        math::Real d13 = d14 - d34;
        math::Real t0 = d12;
        math::Real r0 = (d14 - d12) / (d14 - d34);
        math::Real r = r0;
        math::Int n = math::ceil(log(d34 / d12) / log(r0),
            (math::Real)0.01) - 1;
        if (n > 1) {
            // Newton method to adjust the sum of available space.
            math::Real f = 1;
            const math::Real threshold =
                std::numeric_limits<math::Real>::epsilon() * 1.0e6;
            while (std::abs(f) >= threshold) {
                f = t0 * (1 - pow(r, n)) / (1 - r) - d13;
                math::Real df = t0 * (1 - pow(r, n)) / pow(1 - r, 2) -
                    t0 * n * pow(r, n - 1) / (1 - r);
                r = r - f / df;
            }
            newSteps.resize(n - 1);
            for (math::Int i = 0; i < n - 2; i++) {
                newSteps[i] = t0 * pow(r, (i + 1));
            }
            newSteps[n - 2] = d34;
        }
        else {
            newSteps.resize(1, d34);
        }
    }
    setAdditionalSteps(d, b, newSteps);
}

typedef Grid<3> Grid3;

} 
} 