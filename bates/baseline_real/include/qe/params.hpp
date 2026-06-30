#pragma once
#include <ql/quantlib.hpp>
#include <ostream>
#include <vector>

// using namespace std;
// using namespace QuantLib;
namespace qe{
    using QuantLib::Date;
    using QuantLib::Real;
    using QuantLib::Rate;
    using QuantLib::DayCounter;
    

struct BatesPParams {
    QuantLib::Real S0;        // initial spot
    QuantLib::Real v0;        // initial variance
    QuantLib::Real mu;        // physical drift (P)

    QuantLib::Real kappaP;    // mean reversion speed (P)
    QuantLib::Real thetaP;    // long-run variance (P)
    QuantLib::Real xi;        // vol-of-vol (same in P/Q in your first model)
    QuantLib::Real rho;       // correlation (same in P/Q initially)

    QuantLib::Real JumpIntensityP;   // Jump Intensity (P)
    QuantLib::Real JumpMeanP;        // Jump Mean (P)
    QuantLib::Real JumpVolatility;  // Jump Volatility (P)
};

struct VRPParams {
    QuantLib::Real lambda;    // simplest: constant VRP for variance drift
    QuantLib::Real eta;       // [TODO] - define
};

struct BatesQParams {
    // State for pricing "today"
    QuantLib::Real S0;        // spot used for pricing (often S_t)
    QuantLib::Real v0;        // initial variance used for pricing (often v_t)

    // Market (Q drift uses r-q)
    QuantLib::Rate r;         // risk-free rate (or a curve later)
    QuantLib::Rate q;         // dividend yield (or carry)
    QuantLib::DayCounter dc;  // day count for time conversion

    // Q-model params
    QuantLib::Real kappaQ;
    QuantLib::Real thetaQ;
    QuantLib::Real xi;
    QuantLib::Real rho;

    QuantLib::Real JumpVolatility;
    QuantLib::Real JumpMeanQ;
    QuantLib::Real JumpIntensityQ;
};


std::ostream& operator<<(std::ostream& os, const BatesPParams& P);
std::ostream& operator<<(std::ostream& os, const BatesQParams& Q);
BatesQParams toQ(const BatesPParams& P, const VRPParams& V, Rate r, Rate q, const DayCounter& dc,
    const Real* S_today = nullptr, const Real* v_today = nullptr);
}