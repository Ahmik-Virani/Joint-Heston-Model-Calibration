#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "qe/surface.hpp"

using namespace QuantLib;
using namespace std;
using namespace qe;

namespace qe{
void print_Callgrid(const CallGrid& grid){

    cout <<"Starting S:"<<grid.S0<< endl;
    cout <<"Starting Vol:"<<grid.v0<< endl;
    cout <<"EvalDate:"<<grid.evalDate<< endl;
    
    
    cout << "printing inside the call grid"<<endl;
    cout<<"Call Price grid (rows = maturity, cols = strikes)"<<endl;
    cout << setw(20) << "T\\K";
    for (Real K : grid.strikes) cout << setw(12) << K;
    cout << "\n";

    for (Size iT = 0; iT < grid.maturities.size(); ++iT) {
        cout << setw(12) << grid.maturities[iT];
        for (Size iK = 0; iK < grid.strikes.size(); ++iK) {
            cout << setw(12) << fixed << setprecision(4) << grid.C[iT][iK];
        }
        cout << "\n";
    }

    cout<<"Implied Vol grid (rows = maturity, cols = strikes)"<<endl;
    cout << setw(20) << "T\\K";
    for (Real K : grid.strikes) cout << setw(12) << K;
    cout << "\n";

    for (Size iT = 0; iT < grid.maturities.size(); ++iT) {
        cout << setw(12) << grid.maturities[iT];
        for (Size iK = 0; iK < grid.strikes.size(); ++iK) {
            cout << setw(12) << fixed << setprecision(4) << grid.ImpliedVol[iT][iK];
        }
        cout << "\n";
    }

}

Real impliedVolFromCallPrice(Real CallPrice,Real S0, Real K, Time T, Rate r, Rate q){
    if (T < 0.0) throw invalid_argument("T > 0.0");
    Real discount = exp(-r * T);
    Real forward = S0 * exp((r-q)*T);
    Real guessStdDev = 0.2 * std::sqrt(T);
    // std::cout << "K=" << K
    //       << " T=" << T
    //       << " C=" << CallPrice
    //       << " F=" << forward
    //       << " D=" << discount
    //       << std::endl;

    Real stdDev = blackFormulaImpliedStdDev(
        Option::Call,
        K,
        forward,
        CallPrice,
        discount,
        0.0,
        guessStdDev,
        1e-8,     // accuracy
        200      // max iterations
    );
    return stdDev / sqrt(T);
}


void callgrid(const HestonQParams& Q,CallGrid& grid){
    if (Q.kappaQ <= 0.0) throw invalid_argument("kappaP must be > 0");
    if (Q.thetaQ <= 0.0) throw invalid_argument("thetaP must be > 0");
    if (Q.xi < 0.0)      throw invalid_argument("xi must be >= 0");
    if (fabs(Q.rho) > 1.0) throw invalid_argument("rho must be in [-1,1]");
    if(grid.S0 <= 0.0) throw invalid_argument("S0 must be >= 0");
    if(grid.v0 < 0.0) throw invalid_argument("V0 muste ve > 0");

    Calendar cal = TARGET();
    //Date today(25,February,2026);
    //Date today (grid.evalDate);
    DayCounter dc = Actual365Fixed();
    Settings::instance().evaluationDate() = grid.evalDate;
    
    Handle<Quote> spot(ext::make_shared<SimpleQuote>(grid.S0));
    Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(grid.evalDate, Q.r, Q.dc));
    Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(grid.evalDate, Q.q, Q.dc));

    auto process = ext::make_shared<HestonProcess>(rTS,qTS,spot,grid.v0,Q.kappaQ,Q.thetaQ,Q.xi,Q.rho,HestonProcess::QuadraticExponential);
    auto model = ext::make_shared<HestonModel>(process);
    auto engine =  ext::make_shared<AnalyticHestonEngine>(model);

    //vector <vector <Real>> C(grid.maturities.size(),vector<Real> (grid.strikes.size(),0));
    grid.C.assign(grid.maturities.size(),vector<Real> (grid.strikes.size(),0));
    for (Size iT = 0; iT < grid.maturities.size();iT++){
        auto exercise = ext::make_shared<EuropeanExercise>(grid.maturities[iT]);
        for (Size iK = 0; iK < grid.strikes.size();iK ++ ){
            auto payoff = ext::make_shared<PlainVanillaPayoff>(Option::Call,grid.strikes[iK]);
            VanillaOption opt(payoff,exercise);
            opt.setPricingEngine(engine);
            grid.C[iT][iK] = opt.NPV();
        }
    }
    grid.ImpliedVol.assign(grid.maturities.size(),vector<Real> (grid.strikes.size(),0));
    for(int iT = 0; iT < grid.maturities.size();iT++){
        for(int iK = 0; iK < grid.strikes.size();iK++){
            Time T = dc.yearFraction(grid.evalDate, grid.maturities[iT]);
            grid.ImpliedVol[iT][iK] = impliedVolFromCallPrice(grid.C[iT][iK],grid.S0,grid.strikes[iK],T,grid.r,grid.q);
        }
    }

    // cout<< "Printing inside the call surface function" <<endl;
    // cout<<"Call Price grid (rows = maturity, cols = strikes)"<<endl;
    // cout << setw(20) << "T\\K";
    // for (Real K : grid.strikes) cout << setw(12) << K;
    // cout << "\n";

    // for (Size iT = 0; iT < grid.maturities.size(); ++iT) {
    //     cout << setw(12) << grid.maturities[iT];
    //     for (Size iK = 0; iK < grid.strikes.size(); ++iK) {
    //         cout << setw(12) << fixed << setprecision(4) << C[iT][iK];
    //     }
    //     cout << "\n";
    // }

}   
} //namespace qe 
