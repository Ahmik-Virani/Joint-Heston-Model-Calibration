#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "qe/path.hpp"

using namespace std;
using namespace QuantLib;
using namespace qe;

namespace qe{
ostream& operator<<(ostream& os, const PPath& p) {
    os << "\n=========== PPath ===========\n";
    os << "Path length      : " << p.logS.size() << "\n";
    os << "Variance length  : " << p.v.size() << "\n";
    os << "Returns length   : " << p.returns.size() << "\n\n";

    auto print_sample = [&](const std::vector<Real>& vec,
                            const std::string& name) {
        os << name << " (first 5): ";
        Size n = std::min<Size>(5, vec.size());
        for (Size i = 0; i < n; ++i)
            os << std::fixed << std::setprecision(6) << vec[i] << " ";
        if (vec.size() > 5) os << "...";
        os << "\n";
    };

    print_sample(p.logS, "logS");
    print_sample(p.v, "v");
    print_sample(p.returns, "returns");

    os << "================================\n";
    return os;
}

PPath logReturns(const HestonPParams& P, Size steps, Size seed){
    PPath ppath;
    if (P.kappaP <= 0.0) throw invalid_argument("kappaP must be > 0");
    if (P.thetaP <= 0.0) throw invalid_argument("thetaP must be > 0");
    if (P.xi < 0.0)      throw invalid_argument("xi must be >= 0");
    if (fabs(P.rho) > 1.0) throw invalid_argument("rho must be in [-1,1]");
    if(P.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
    if(P.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");

    Calendar cal = TARGET();
    Date today(26,February,2026);
    Settings::instance().evaluationDate() = today;
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(today,P.mu,Actual365Fixed()));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(today, 0.0, Actual365Fixed()));
    Handle<Quote> s0(ext::make_shared<SimpleQuote>(P.S0));
    auto process = ext::make_shared<HestonProcess>(rTS,qTS,s0,P.v0,P.kappaP,P.thetaP,P.xi,P.rho,
                    HestonProcess::QuadraticExponential);

    Time T = 1.0;
    Time dt = T/steps;
    Size factors = 2;

    TimeGrid grid(T,steps);
    typedef PseudoRandom::rsg_type rsg_type;
    rsg_type rsg = PseudoRandom::make_sequence_generator(factors * grid.size()-factors, seed);
    MultiPathGenerator<rsg_type> gen(process,grid,rsg,false);
    const MultiPath path = gen.next().value;
    
    ppath.logS.resize(grid.size());
    ppath.v.resize(grid.size());
    ppath.returns.resize(grid.size() - 1);
    
    for(Size i = 0; i < grid.size();i++){
        ppath.logS[i] = log(path[0][i]);
        ppath.v[i] = path[1][i];
    }

    for(Size i = 1; i < grid.size();i++){
        ppath.returns[i-1] = ppath.logS[i] - ppath.logS[i-1]; 
    }
    return ppath;
}
} //namespace qe
