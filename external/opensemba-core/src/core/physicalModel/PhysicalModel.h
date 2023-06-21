#pragma once

#include "core/math/Constants.h"
#include "core/util/Class.h"
#include "core/util/Identifiable.h"
#include "core/util/Identification.h"

namespace semba {
namespace physicalModel {

class PhysicalModel;
typedef util::Identification<PhysicalModel> Id;

class PhysicalModel : public virtual util::Identifiable<Id>,
                      public virtual util::Class {
public:
    enum class Type {
        PEC,
        PMC,
        SMA,
        vacuum,
        classic,
        elecDispersive,
        anisotropic,
        isotropicsibc,
        PML,
        wire,
        gap,
        multiport,
        priorityMaterial
    };

    virtual std::unique_ptr<PhysicalModel> clone() const = 0;
    
    PhysicalModel() = default;
    PhysicalModel(const PhysicalModel&);
    PhysicalModel(const std::string& name);
    
    
    virtual ~PhysicalModel() = default;

    const std::string& getName() const;
    void setName(const std::string& newName);

private:
    std::string name_ = "";
};


} 
} 

