

#pragma once

#include <string>

#include "core/math/Types.h"

namespace semba {
namespace outputRequest {

class Domain {
public:
    typedef enum {
        NONE, TIME, FREQ, TRAN, TIFR, TITR, FRTR, ALL
    } Type;

    Domain();
    Domain(bool timeDomain,
           math::Real initialTime,
           math::Real finalTime,
           math::Real samplingPeriod,
           bool frequencyDomain,
           math::Real initialFrequency,
           math::Real finalFrequency,
           math::Real frequencyStep,
           bool logFrequencySweep,
           bool usingTransferFunction,
           std::string transferFunctionFile);
    Domain(const Domain& rhs);
    virtual ~Domain();

    Domain& operator=(const Domain& rhs);
    void setFinalTime(const math::Real finalTime);
    void setSamplingPeriod(const math::Real samplingPeriod);


    bool operator==(const Domain& rhs) const;

    bool       isTimeDomain() const;
    math::Real getInitialTime() const;
    math::Real getFinalTime() const;
    math::Real getSamplingPeriod() const;
    bool       isFrequencyDomain() const;
    math::Real getInitialFrequency() const;
    math::Real getFinalFrequency() const;
    math::Real getFrequencyStep() const;
    bool       isLogFrequencySweep() const;

    bool isUsingTransferFunction() const;
    const std::string& getTransferFunctionFile() const;

    Type getDomainType() const;

private:
    bool        timeDomain_;
    math::Real  initialTime_;
    math::Real  finalTime_;
    math::Real  samplingPeriod_;
    bool        frequencyDomain_;
    math::Real  initialFrequency_;
    math::Real  finalFrequency_;
    math::Real  frequencyStep_;
    bool        logFrequencySweep_;
    bool        usingTransferFunction_;
    std::string transferFunctionFile_;
};

} 
} 

