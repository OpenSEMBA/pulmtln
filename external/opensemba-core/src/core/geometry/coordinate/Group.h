#pragma once

#include <type_traits>
#include <map>

#ifdef _MSC_VER
#pragma warning(disable:4250)
#endif

#include "core/util/GroupIdentifiableUnique.h"

#include "Coordinate.h"

namespace semba {
namespace geometry {
namespace coordinate {

struct CoordComparator {
    template<class T, std::size_t D>
    bool operator() (const Coordinate<T,D>* lhs,
                     const Coordinate<T,D>* rhs) const {
        return (lhs->pos() < rhs->pos());
    }
};

template<typename C = Coord>
class Group final : public util::GroupIdentifiableUnique<C> {
public:
    Group() = default;
    
    Group(const std::vector<math::CVecR3>&);
    Group(const std::vector<math::CVecI3>&);

    void applyScalingFactor(const math::Real factor);

    template<typename VEC>
    std::map<VEC, std::vector<const C*>> getIndex() const;

private:
    template<typename VEC>
    typename util::GroupIdentifiableUnique<C>::iterator  addPos(VEC);
};


template<typename C>
Group<C>::Group(const std::vector<math::CVecR3>& pos)
{
    addPos(pos);
}

template<typename C>
Group<C>::Group(const std::vector<math::CVecI3>& pos)
{
    addPos(pos);
}

template<typename C> template<typename VEC>
typename util::GroupIdentifiableUnique<C>::iterator
Group<C>::addPos(VEC newPosition)
{
    CoordId newId;
    if (this->size() == 0) {
        newId = CoordId(1);
    }
    else {
        auto backIt = --this->end();
        newId = ++(backIt->get()->getId());
    }
    return add(std::make_unique<C>(newId, newPosition));
}

template<typename C>
void Group<C>::applyScalingFactor(const math::Real factor)
{
    for (auto& c : *this) {
        if (c->template is<CoordR3>()) {
            *c *= factor;
        }
    }
}

template<typename C> template<typename VEC>
std::map<VEC, std::vector<const C*>> Group<C>::getIndex() const {
    std::map<VEC, std::vector<const C*>> index;

    for (auto& c : *this) {
        auto cIt = index.find(c->pos());
        if (cIt == index.end()) {
            index.emplace(c->pos(), std::vector<const C*>({ c.get() }));

            continue;
        }

        auto& vector = cIt->second;
        vector.push_back(c.get());
    }

    return index;
}

} 

typedef coordinate::Group<CoordR3> CoordR3Group;
typedef coordinate::Group<CoordI3> CoordI3Group;
typedef coordinate::Group<> CoordGroup;

} 
} 