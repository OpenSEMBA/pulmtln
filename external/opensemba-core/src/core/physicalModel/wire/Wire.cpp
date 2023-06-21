#include "Wire.h"

namespace semba {
namespace physicalModel {
namespace wire {

Wire::Wire(const Id id,
           const std::string name,
           const math::Real radius,
           const math::Real resistance,
           const math::Real inductance) :   
    Identifiable<Id>(id),
    PhysicalModel(name)
{
    radius_ = radius;
    seriesResistance_ = resistance;
    seriesInductance_ = inductance;
}

Wire::Wire(const Id id,
           const std::string name,
           const math::Real radius,
           const math::Real resistance,
           const math::Real inductance,
           const math::Real capacitance,
           const math::Real pResistance,
           const math::Real pInductance,
           const math::Real pCapacitance) : 
    Identifiable<Id>(id),
    PhysicalModel(name)
{
    radius_ = radius;
    isSeriesParallel_ = true;
    isDispersive_ = false;
    seriesResistance_ = resistance;
    seriesInductance_ = inductance;
    seriesCapacitance_ = capacitance;
    parallelResistance_ = pResistance;
    parallelInductance_ = pInductance;
    parallelCapacitance_ = pCapacitance;
}

Wire::Wire(const Id id,
           const std::string name,
           const math::Real radius,
           const std::string filename) :
    Identifiable<Id>(id),
    PhysicalModel(name)
{
    radius_ = radius;
    filename_ = filename;
}

Wire::Wire(const Wire& rhs) : 
    Identifiable<Id>(rhs),
    PhysicalModel(rhs) 
{
    radius_ = rhs.radius_;
    isSeriesParallel_ = rhs.isSeriesParallel_;
    isDispersive_ = rhs.isDispersive_;
    seriesResistance_ = rhs.seriesResistance_;
    seriesInductance_ = rhs.seriesInductance_;
    seriesCapacitance_ = rhs.seriesCapacitance_;
    parallelResistance_ = rhs.parallelResistance_;
    parallelInductance_ = rhs.parallelInductance_;
    parallelCapacitance_ = rhs.parallelCapacitance_;
    filename_ = rhs.filename_;
}

bool Wire::isSeriesParallel() const 
{
    return isSeriesParallel_;
}

bool Wire::isDispersive() const 
{
    return isDispersive_;
}

math::Real Wire::getRadius() const 
{
    return radius_;
}

math::Real Wire::getSeriesResistance() const 
{
    return seriesResistance_;
}

math::Real Wire::getSeriesInductance() const 
{
    return seriesInductance_;
}

math::Real Wire::getSeriesCapacitance() const 
{
    return seriesCapacitance_;
}

math::Real Wire::getParallelResistance() const 
{
    return parallelResistance_;
}

math::Real Wire::getParallelInductance() const 
{
    return parallelInductance_;
}

math::Real Wire::getParallelCapacitance() const 
{
    return parallelCapacitance_;
}

std::string Wire::getFilename() const 
{
    return filename_;
}

}
} 
} 
