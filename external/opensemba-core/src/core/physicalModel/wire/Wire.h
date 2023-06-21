#pragma once

#include "core/physicalModel/PhysicalModel.h"

namespace semba {
namespace physicalModel {
namespace wire {

class Wire : public virtual PhysicalModel {
public:
    Wire(const Id id,
         const std::string name,
         const math::Real radius,
         const math::Real resistance,
         const math::Real inductance);
    Wire(const Id id,
         const std::string name,
         const math::Real radius,
         const math::Real resistance,
         const math::Real inductance,
         const math::Real capacitance,
         const math::Real pResistance,
         const math::Real pInductance,
         const math::Real pCapacitance);
    Wire(const Id id,
         const std::string name,
         const math::Real radius,
         const std::string filename);
    
    Wire(const Wire&);
    virtual ~Wire() = default;
    
    virtual std::unique_ptr<PhysicalModel> clone() const override {
        return std::make_unique<Wire>(*this);
    }

    math::Real getRadius() const;

    bool isSeriesParallel() const;
    bool isDispersive() const;

    math::Real getSeriesResistance() const;
    math::Real getSeriesInductance() const;
    math::Real getSeriesCapacitance() const;
    math::Real getParallelResistance() const;
    math::Real getParallelInductance() const;
    math::Real getParallelCapacitance() const;

    std::string getFilename() const;

private:
    math::Real radius_ = 0.0;
    bool isSeriesParallel_ = false;
    bool isDispersive_     = false;
    math::Real seriesResistance_   = 0.0;    // Resistance per meter.
    math::Real seriesInductance_   = 0.0;    // Inductance per meter.
    math::Real seriesCapacitance_  = 0.0;   // Capacitance per meter.
    math::Real parallelResistance_ = 0.0;  // Resistance per meter.
    math::Real parallelInductance_ = 0.0;  // Inductance per meter.
    math::Real parallelCapacitance_= 0.0; // Capacitance per meter.
    std::string filename_ = "";
};

}
} 
} 

