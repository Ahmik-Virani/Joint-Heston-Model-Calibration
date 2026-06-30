#pragma once
#include <ql/quantlib.hpp>

namespace qe {
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