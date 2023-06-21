

#pragma once

#include <cassert>
#include <exception>
#include <vector>

#include "Surface.h"

namespace semba {
namespace physicalModel {
namespace surface {

class Multilayer : public virtual Surface {
public:
    class Layer {
    public:
        Layer(math::Real thickness, math::Real relPermittivity,
                math::Real relPermeability, math::Real elecCond) :
                thickness_(thickness), relPermittivity_(relPermittivity),
                relPermeability_(relPermeability), elecCond_(elecCond) {
			if (relPermittivity_ == 0.0) {
				throw std::logic_error(
					"Layer relative permittivity must be greater than zero.");
			}
			if (relPermeability_ == 0.0) {
				throw std::logic_error(
					"Layer relative permeability must be greater than zero.");
			}
		}

        math::Real getThickness() const {return thickness_;}
        math::Real getRelPermittivity() const {return relPermittivity_;}
        math::Real getPermittivity() const {
            return relPermittivity_ * math::Constants::eps0;
        }
        math::Real getRelPermeability() const {return relPermeability_;}
        math::Real getPermeability() const {
            return relPermeability_ * math::Constants::mu0;
        }
        math::Real getElecCond() const {return elecCond_;}
    private:
        math::Real thickness_;
        math::Real relPermittivity_;
        math::Real relPermeability_;
        math::Real elecCond_;
    };

    class FittingOptions {
    public:
        FittingOptions(
                std::pair<math::Real, math::Real> minMaxFreq,
                size_t numberOfPoles ) :
                    minMaxFreq_(minMaxFreq), numberOfPoles_(numberOfPoles) { }

        std::pair<math::Real, math::Real> getMinMaxFreq() const {
            return minMaxFreq_;
        }

        std::size_t getNumberOfPoles() const {
            return numberOfPoles_;
        }

    private:
        std::pair<math::Real, math::Real> minMaxFreq_;
        std::size_t numberOfPoles_;
    };

    Multilayer(const Id id,
            const std::string& name,
            const std::vector<Layer>& layers,
            const std::vector<FittingOptions>& options = {});
    virtual ~Multilayer();

    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Multilayer>(*this);
    }

    Layer getLayer(size_t i) const {
        return layers_[i];
    }
    std::size_t getNumberOfLayers() const {
        return layers_.size();
    }

    math::Real getThickness   (const std::size_t i) const;
    math::Real getPermittivity(const std::size_t i) const;
    math::Real getPermeability(const std::size_t i) const;
    math::Real getElecCond    (const std::size_t i) const;
    math::Real getMagnCond    (const std::size_t i) const;

    bool hasFittingOptions() const;
    FittingOptions getFittingOptions() const;

    std::string printLayer(const std::size_t i) const;

private:
    std::vector<Layer> layers_;

    std::vector<FittingOptions> options_; // optional (0 or 1 size).
};

namespace Error {
namespace SurfaceMultilayer {

class IncompatibleSizes : public std::exception {
public:
    IncompatibleSizes() {}
    virtual ~IncompatibleSizes() throw() {}

    const char* what() const throw() {
        return "SurfaceMultilayer: "
                "Incompatible sizes of layers parameters.";
    }
};

}
} 
} 
} 
} 

