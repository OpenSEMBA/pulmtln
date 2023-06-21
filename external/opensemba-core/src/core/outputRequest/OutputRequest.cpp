#include "OutputRequest.h"

#include <iostream>

namespace semba {
namespace outputRequest {

OutputRequest::OutputRequest(
    const Type& type, const std::string& name, const Domain& domain, const Target& target) :
    name_(name),
    type_(type),
    domain_(domain),
    target_(target)
{
}

std::string OutputRequest::getTypeStr() const 
{
    switch (type_) {
    case Type::electric:
        return "Electric field";
    case Type::magnetic:
        return "Magnetic field";
	case Type::bulkCurrentElectric:
		return "Bulk current electric";
	case Type::bulkCurrentMagnetic:
		return "Bulk current magnetic";
	case Type::surfaceCurrentDensity:
		return "Surface current density";
    case Type::electricFieldNormals:
        return "Electric field normals";
    case Type::magneticFieldNormals:
        return "Magnetic field normals";
    case Type::current:
        return "Current";
    case Type::voltage:
        return "Voltage";
	case Type::electricFarField:
		return "Electric far field";
    default:
        throw std::logic_error("Unrecognized output type");
    }
}

} 
} 
