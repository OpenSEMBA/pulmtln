

#include <physicalModel/surface/SIBCFile.h>

namespace semba {
namespace physicalModel {
namespace surface {

SIBCFile::SIBCFile(const Id id,
                         const std::string& name,
                         const util::ProjectFile& file)
:   Identifiable<Id>(id),
    PhysicalModel(name),
    file_(file) {
    std::string extension = file_.getExtension();
    if (extension.compare(".mibc") != 0) {
        throw std::logic_error("File extension must be .mibc in file: " + file_);
    }
}

const util::ProjectFile SIBCFile::getFile() const {
    return file_;
}

} 
} 
} 
