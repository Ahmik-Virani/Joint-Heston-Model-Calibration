#include "qe/params.hpp"

#include <iomanip>
#include <cmath>
#include <stdexcept>


using namespace QuantLib;
using namespace std;

using namespace qe;

namespace qe{
ostream& operator<<(ostream& os, const HestonPParams& P) {
    os << "=== Heston P Parameters ===\n";
    os << "S0      : " << P.S0 << "\n";
    os << "v0      : " << P.v0 << "\n";
    os << "mu      : " << P.mu << "\n";
    os << "kappaP  : " << P.kappaP << "\n";
    os << "thetaP  : " << P.thetaP << "\n";
    os << "xi      : " << P.xi << "\n";
    os << "rho     : " << P.rho << "\n";
    os << "===========================\n";
    return os;
}

ostream& operator<<(ostream& os, const HestonQParams& Q) {
    os << "=== Heston Q Parameters ===\n";
    os << "S0      : " << Q.S0 << "\n";
    os << "v0      : " << Q.v0 << "\n";
    os << "r       : " << Q.r << "\n";
    os << "q       : " << Q.q << "\n";
    os << "DayCount: " << Q.dc.name() << "\n";
    os << "kappaQ  : " << Q.kappaQ << "\n";
    os << "thetaQ  : " << Q.thetaQ << "\n";
    os << "xi      : " << Q.xi << "\n";
    os << "rho     : " << Q.rho << "\n";
    os << "===========================\n";
    return os;
}

HestonQParams toQ(const HestonPParams& P, const VRPParams& V, Rate r, Rate q, const DayCounter& dc,
    const Real* S_today, const Real* v_today){
        // Basic validation (keep it light but protective)
        if (P.kappaP <= 0.0) throw std::invalid_argument("kappaP must be > 0");
        if (P.thetaP <= 0.0) throw std::invalid_argument("thetaP must be > 0");
        if (P.xi < 0.0)      throw std::invalid_argument("xi must be >= 0");
        if (fabs(P.rho) > 1.0) throw std::invalid_argument("rho must be in [-1,1]");

        Real kappaQ = P.kappaP + V.lambda;
        if (kappaQ <= 0.0) throw invalid_argument("kappaQ must be > 0; Lambda is too negative?");
        Real thetaQ = (P.kappaP * P.thetaP) / kappaQ;

        HestonQParams Q;
        Q.r = r;
        Q.q = q;
        Q.dc = dc;

        Q.kappaQ = kappaQ;
        Q.thetaQ = thetaQ;
        Q.xi = P.xi;
        Q.rho = P.rho;

        Q.S0 = (S_today ? *S_today : P.S0);
        Q.v0 = (v_today ? *v_today : P.v0);

        if(Q.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
        if(Q.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");
        return Q;

    }
}//namespace qe
