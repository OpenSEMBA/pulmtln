#pragma once

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include "Function.h"

namespace semba {
namespace math {
namespace function {

template<class S, class T>
class LinearInterpolation : public Function<S,T> {
public:
    LinearInterpolation() = default;
    LinearInterpolation(const std::vector<std::pair<S,T>>& xy);
    LinearInterpolation(const std::string& file);
    virtual ~LinearInterpolation() = default;

    SEMBA_MATH_FUNCTION_DEFINE_CLONE(LinearInterpolation)

    T operator()(const S& arg) const;
    bool operator==(const Base& arg) const;
private:
    std::map<S,T> value_;

    void initFromFile(const std::string& file);
};

template<class S, class T>
inline LinearInterpolation<S, T>::LinearInterpolation(
        const std::vector<std::pair<S, T> >& xy) {
    for (std::size_t i = 0; i < xy.size(); i++) {
        value_.insert(xy[i]);
    }
}

template<class S, class T>
inline LinearInterpolation<S, T>::LinearInterpolation(
        const std::string& file) {
    initFromFile(file);
}

template<class S, class T>
inline T LinearInterpolation<S, T>::operator ()(const S& pos) const {
    if (value_.empty()) {
        throw std::out_of_range("Attempting to interpolate with no data.");
    }
    std::pair<typename std::map<S,T>::const_iterator,
              typename std::map<S,T>::const_iterator> range;
    range = value_.equal_range(pos);
    if (range.first == range.second) {
        if (range.first == value_.end()) {
            return value_.rbegin()->second;
        } else {
            if (pos <= range.first->first && range.first == value_.begin()) {
                return range.first->second;
            } else {
                typename std::map<S,T>::const_iterator prev = range.first;
                --prev;
                const Real x1 = prev->first;
                const Real y1 = prev->second;
                const Real x2 = range.second->first;
                const Real y2 = range.second->second;
                const Real m = (y2 - y1) / (x2 - x1);
                const Real n = (x2*y1 - x1*y2)/(x2 - x1);
                return (m * pos + n);
            }
        }
    } else {
        if (range.first->first == pos) {
            return range.first->second;
        } else {
            return range.first->second;
        }
    }
}

template<class S, class T>
inline bool LinearInterpolation<S, T>::operator ==(
        const Base& rhs) const {
    if (typeid(*this) != typeid(rhs)) {
        return false;
    }
    const LinearInterpolation<S,T>* rhsPtr =
        dynamic_cast<const LinearInterpolation<S,T>*>(&rhs);
    return value_ == rhsPtr->value_;
}

template<class S, class T>
void LinearInterpolation<S, T>::initFromFile(
        const std::string& file) {
    static_assert(std::is_same<S, Real>::value, "S must be Real");
    static_assert(std::is_same<T, Real>::value, "T must be Real");
    std::ifstream iStream;
    iStream.open(file, std::ifstream::in);
    if (iStream.fail()) {
        throw std::ios_base::failure(std::string("File: ") + file +
                                     std::string(" does not exist"));
    }

    while (!iStream.eof()) {
        std::pair<Real, Real> value;
        iStream >> value.first >> value.second;
        value_.insert(value);
    }

    if (value_.size() == 0) {
        throw std::ios_base::failure(std::string("No data in file: ") + file);
    }
}

} 
} 
} 


