#pragma once

#include <string>
#include <memory>

#include "core/util/Identifiable.h"
#include "core/util/Identification.h"
#include "core/util/GroupIdentifiableUnique.h"

namespace semba {
namespace geometry {

class Layer;
using LayerId = util::Identification<Layer>;

class Layer final : public virtual util::Class,
              public virtual util::Identifiable<LayerId> {
public:
    Layer() = default;
    Layer(const LayerId id, const std::string& name);
    Layer(const std::string& name);
    Layer(const Layer& rhs);
    
    std::unique_ptr<Layer> clone() const {
        return std::make_unique<Layer>(*this);
    }
    
    virtual bool operator==(const Layer& rhs) const;
    virtual bool operator!=(const Layer& rhs) const;

    std::string getName() const;

    virtual std::string getParentName() const;
    virtual std::string getChildName() const;
    std::string toStr() const;

    friend std::ostream& operator<<(std::ostream& os, const Layer& lay) {
       return os << lay.toStr();
    }

private:
    std::string name_;

    static std::string spaceToUnderscore(std::string rhs);
};


class LayerGroup final : public util::GroupIdentifiableUnique<Layer>  {
public:
    const Layer* getName(const std::string name) const;
};

}
}
