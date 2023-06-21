#include "Layer.h"

#include <cassert>
#include <iostream>
#include <sstream>

namespace semba {
namespace geometry {

Layer::Layer(const LayerId id, const std::string& name)
:   Identifiable<LayerId>(id) 
{
    name_ = spaceToUnderscore(name);
}

Layer::Layer(const std::string& name) 
{
    name_ = spaceToUnderscore(name);
}

Layer::Layer(const Layer& rhs)
:   Identifiable<LayerId>(rhs) 
{
    name_ = rhs.name_;
}

bool Layer::operator ==(const Layer& rhs) const {
    if (typeid(*this) != typeid(rhs)) {
        return false;
    }
    bool res = true;
    res &= (this->getId() == rhs.getId());
    res &= (this->getName() == rhs.getName());
    return res;
}

bool Layer::operator !=(const Layer& rhs) const {
    return !(*this == rhs);
}

std::string Layer::getName() const {
    return name_;
}

std::string Layer::getParentName() const {
   assert(false);
   return std::string();
}

std::string Layer::getChildName() const {
   assert(false);
   return std::string();
}

std::string Layer::toStr() const {
    std::stringstream ss;
    ss << "Layer. Id: " << getId() << " Name: " << getName();
    return ss.str();
}

std::string Layer::spaceToUnderscore(std::string rhs) {
    std::string str = rhs;
    for(std::string::iterator it = str.begin(); it != str.end(); ++it) {
        if(*it == ' ') {
            *it = '_';
        }
    }
    return str;
}

const Layer* LayerGroup::getName(const std::string name) const
{
    for (auto& it{ this->begin() }; it != this->end(); ++it) {
        if (it->get()->getName() == name) {
            return it->get();
        }
    }
    return nullptr;
}

} 
} 
