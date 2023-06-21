#pragma once

#include "core/math/matrix/Static.h"
#include "core/util/ProjectFile.h"

#include "Surface.h"

namespace semba {
namespace physicalModel {
namespace surface {

class SIBCFile : public virtual Surface {
public:
    SIBCFile();
    SIBCFile(const Id id,
             const std::string& name,
             const util::ProjectFile& file);
    
    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<SIBCFile>(*this);
    }

    const util::ProjectFile getFile() const;

protected:
    util::ProjectFile file_;
};

} 
} 
} 

