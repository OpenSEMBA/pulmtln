#pragma once

#include "PhysicalModel.h"

namespace semba {
namespace physicalModel {

class PredefinedGeneric : public virtual PhysicalModel {
};

template <class T>
class Predefined : public PredefinedGeneric {
public:   
    class PEC;
    class PMC;
    class SMA;
    class Vacuum;

    Predefined(Id id, const std::string& name) :
        Identifiable<Id>(id), 
        PhysicalModel(name) {
    }
    virtual ~Predefined() = default;
    
    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Predefined<T>>(*this);
    }
};

typedef Predefined<Predefined<void>::PEC> PEC;
typedef Predefined<Predefined<void>::PMC> PMC;
typedef Predefined<Predefined<void>::SMA> SMA;
typedef Predefined<Predefined<void>::Vacuum> Vacuum;

} 
} 

