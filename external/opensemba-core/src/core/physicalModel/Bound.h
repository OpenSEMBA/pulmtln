#pragma once

#include "PhysicalModel.h"

namespace semba {
namespace physicalModel {

class Bound : public virtual PhysicalModel {
public:
    enum class Type {
        mur1,
        mur2,
        pec,
        pmc,
        periodic,
        pml,
        sma
    };
    
    Bound(Id id, Type type);
    Bound(const Bound&);
    virtual ~Bound() = default;

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Bound>(*this);
    }

    Type getType() const;
    std::string getTypeName() const;

private:
    Type type;
};

} 
} 

