#pragma once
#include <ql/quantlib.hpp>
#include <vector>
#include "qe/params.hpp"

namespace qe {
    using QuantLib::Date;
    using QuantLib::Real;
    using QuantLib::Rate;
    using QuantLib::Time;
    
    struct CallGrid{
        std::vector <Date> maturities;
        std::vector <int> strikes;
        std::vector<std::vector<Real>> C;
        std::vector<std::vector<Real>> ImpliedVol;
        Real S0;
        Real v0;
        Rate r;
        Rate q;
        Date evalDate;
    };

    void print_Callgrid(const CallGrid& grid);
    void print_Surfaces(const std::vector<CallGrid>& surfaces);
    void callgrid(const HestonQParams& Q,CallGrid& grid);
    Real impliedVolFromCallPrice(Real CallPrice,Real S0, Real K, Time T, Rate r, Rate q);
}
