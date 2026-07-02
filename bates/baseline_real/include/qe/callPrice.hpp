#pragma once
#include <ql/quantlib.hpp>

namespace qe {
double batesCallPrice(
    double spot,
    double strike,
    double maturityYears,
    double r,
    double q,
    double v0,
    double kappa,
    double theta,
    double xi,
    double rho,
    double jumpIntensity,
    double jumpMean,
    double jumpVolatility,
    const QuantLib::Date& evalDate
);

double BatesCallPrice(
    double spot,
    double strike,
    double maturityYears,
    double r,
    double q,
    double v0,
    double kappa,
    double theta,
    double xi,
    double rho,
    const QuantLib::Date& evalDate
);
}
