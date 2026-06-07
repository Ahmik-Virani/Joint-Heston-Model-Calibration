#pragma once
#include <ql/quantlib.hpp>
#include <vector>
#include <ostream>
#include "qe/params.hpp"

namespace qe{
    using QuantLib::Real;
    using QuantLib::Size;

    struct PPath{
        std::vector<Real> logS;
        std::vector<Real> v;
        std::vector<Real> returns;
    };

    std::ostream& operator<<(std::ostream& os, const PPath& p);
    PPath logReturns(const HestonPParams& P, Size steps, Size seed);
}

