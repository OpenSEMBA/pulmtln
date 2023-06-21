#include "Port.h"

namespace semba {
namespace source {
namespace port {

Port::Port(const std::unique_ptr<Magnitude::Magnitude>& magnitude,
           const Target& elem) :   
    Source(magnitude, elem) 
{}

}
}
} 

