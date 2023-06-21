

#include <physicalModel/surface/Multilayer.h>
#include <vector>


namespace semba {
namespace physicalModel {
namespace surface {

Multilayer::Multilayer(
        const Id id,
        const std::string& name,
        const std::vector<Layer>& layers,
        const std::vector<FittingOptions>& options)
:   Identifiable<Id>(id),
    PhysicalModel(name),
    layers_(layers),
    options_(options) {

    if (options_.size() > 1) {
        throw std::runtime_error(
                "Multilayer Fitting options may contain up to one element");
    }
}

Multilayer::~Multilayer() {

}

std::string Multilayer::printLayer(const std::size_t i) const {
    assert(i < getNumberOfLayers());
    std::stringstream ss;
    ss << getElecCond(i) << " " << getPermittivity(i) <<
            " " << getPermeability(i) << " " << getThickness(i);
    return std::string(ss.str());
}

math::Real Multilayer::getThickness(const std::size_t i) const {
    return layers_[i].getThickness();
}

math::Real Multilayer::getPermittivity(const std::size_t i) const {
    return layers_[i].getRelPermittivity() * math::Constants::eps0;
}

math::Real Multilayer::getPermeability(const std::size_t i) const {
    return layers_[i].getRelPermeability() * math::Constants::mu0;
}

math::Real Multilayer::getElecCond(const std::size_t i) const {
    return layers_[i].getElecCond();
}

math::Real Multilayer::getMagnCond(const std::size_t i) const {
    return 0.0;
}

bool Multilayer::hasFittingOptions() const {
    if (options_.size() == 1) {
        return true;
    } else {
        return false;
    }
}

Multilayer::FittingOptions Multilayer::getFittingOptions() const {
    if (options_.size() != 1) {
        throw std::runtime_error(
                "Trying to access a corrupt fitting options value");
    }
    return options_.front();
}

} 
} 
} 


