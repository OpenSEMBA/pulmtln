#pragma once

#include "Class.h"

namespace semba {
namespace util {

template<class Id>
class Identifiable : public virtual Class {
public:
    Identifiable() = default;
    Identifiable(const Id& id);
    Identifiable(const Identifiable& rhs) {
        id_ = rhs.id_;
    }
    virtual ~Identifiable() = default;

    Id   getId() const;
    void setId(const Id& id);

    virtual bool operator==(const Identifiable& rhs) const;

private:
    Id id_{ 0 };
};

template<class Id>
Identifiable<Id>::Identifiable(const Id& id) {
    id_ = id;
}

template<class Id>
Id Identifiable<Id>::getId() const {
    return this->id_;
}

template<class Id>
void Identifiable<Id>::setId(const Id& id) {
    this->id_ = id;
}

template<class Id>
bool Identifiable<Id>::operator==(const Identifiable& rhs) const {
    return (this->id_ == rhs.id_);
}

}
} 