#pragma once
#include <ql/quantlib.hpp>
#include <vector>
#include "qe/params.hpp"

using namespace std;

namespace qe {
    using QuantLib::Date;
    using QuantLib::Real;
    using QuantLib::Rate;
    using QuantLib::Time;
    
    struct CallGrid{
        // std::vector <Date> maturities;
        // std::vector <int> strikes;
        vector<vector<double>> C;
        // vector<vector<Real>> ImpliedVol;
        double S0;
        // Real v0;
        double r;
        double q;
        Date evalDate;
    };

    void print_Callgrid(const CallGrid& grid);
    void print_Surfaces(const std::vector<CallGrid>& surfaces);
    void callgrid(const HestonQParams& Q,CallGrid& grid);
    Real impliedVolFromCallPrice(Real CallPrice,Real S0, Real K, Time T, Rate r, Rate q);
}
