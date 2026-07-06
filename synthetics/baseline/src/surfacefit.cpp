#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "qe/surface.hpp"
#include "qe/surfacefit.hpp"
#include <ql/math/solvers1d/brent.hpp>
#include <ql/math/functional.hpp>
// #include <ql/utilities/disposable.hpp>
#include <ql/math/optimization/costfunction.hpp>
#include <random>



using namespace QuantLib;
using namespace std;
using namespace qe;

namespace qe{

    ostream& operator<<(ostream& os, const HestonSurfaceFit& fit) {
        os << fixed << setprecision(6);
    
        os << "Heston Surface Fit\n";
        os << "-------------------\n";
        os << "v0       : " << fit.v0 << "\n";
        os << "kappaQ   : " << fit.kappaQ << "\n";
        os << "thetaQ   : " << fit.thetaQ << "\n";
        os << "xi       : " << fit.xi << "\n";
        os << "rho      : " << fit.rho << "\n";
        os << "\nFit Diagnostics\n";
        os << "----------------\n";
        os << "RMSE (IV)      : " << fit.rmseIv << "\n";
        os << "Max Abs IV Err : " << fit.maxAbsIvErr << "\n";
    
        return os;
    }


    ostream& operator<<(ostream& os, const HestonMultiSurfaceFit& fit)
    {
        os << "===== Heston Multi-Surface Fit =====\n";
        os << std::fixed << std::setprecision(6);

        os << "kappaQ   : " << fit.kappaQ   << "\n";
        os << "thetaQ   : " << fit.thetaQ   << "\n";
        os << "xi       : " << fit.xi       << "\n";
        os << "rho      : " << fit.rho      << "\n";
        os << "Total SSE: " << fit.TotalSSE << "\n";

        os << "v0 per surface:\n";
        for (Size i = 0; i < fit.v0_by_surface.size(); ++i) {
            os << "  Surface " << i
            << " : " << fit.v0_by_surface[i] << "\n";
        }

        os << "====================================";

        return os;
    }

    HestonSurfaceFit calibrateHestonQVolGrid(CallGrid& grid,HestonSurfaceFit& InitialGuess,const Date& today, const DayCounter& dc,const Calendar& cal){
        // Real v0_guess     = 0.06;   // true 0.04
        // Real kappaQ_guess = 1.0;    // true 2.0
        // Real thetaQ_guess = 0.05;   // true 0.03
        // Real xi_guess     = 0.45;   // true 0.30
        // Real rho_guess    = -0.30;  // true -0.70

        Real v0_guess     = InitialGuess.v0;
        Real kappaQ_guess = InitialGuess.kappaQ;
        Real thetaQ_guess = InitialGuess.thetaQ;
        Real xi_guess     = InitialGuess.xi;
        Real rho_guess    = InitialGuess.rho;

        
        Settings::instance().evaluationDate() = grid.evalDate;
        Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(grid.evalDate, grid.r,dc));
        Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(grid.evalDate, grid.q, dc));
        Handle<Quote> spot(ext::make_shared<SimpleQuote>(grid.S0));
        auto process = ext::make_shared<HestonProcess>(
            rTS, qTS, spot,
            v0_guess, kappaQ_guess, thetaQ_guess, xi_guess, rho_guess
        );
        auto model  = ext::make_shared<HestonModel>(process);
        auto engine = ext::make_shared<AnalyticHestonEngine>(model);
        
        // make the process -> model -> engine; helpers a pointer that will store the actual value and calibration error -> so
        // we need a vector of helpers for M*K (all maturity and grid values)
        
        vector<ext::shared_ptr<CalibrationHelper>> helpers;
        const Size M = grid.maturities.size();
        const Size K = grid.strikes.size();
        helpers.reserve(M * K);
    
        for (int iT = 0; iT < grid.maturities.size();iT++){
            Integer days = cal.businessDaysBetween(grid.evalDate, grid.maturities[iT]);
            QL_REQUIRE(days > 0, "Maturity must be after today");
            Period maturityPeriod(days,Days);
            for(int iK = 0; iK < grid.strikes.size();iK++){
                Real strike = grid.strikes[iK];
                Real iv_mkt = grid.ImpliedVol[iT][iK];
                Handle<Quote> volQuote(ext::make_shared<SimpleQuote>(iv_mkt));
                auto h = ext::make_shared<HestonModelHelper>(
                    maturityPeriod,         // maturity as Period
                    cal,
                    grid.S0,                   // spot
                    strike,
                    volQuote,               // market implied vol quote
                    rTS,
                    qTS
                );
                h->setPricingEngine(engine);
                helpers.push_back(h);
            }
        }
        //LevenbergMarquardt opt(1e-4, 1e-4, 1e-4);
        Simplex opt(0.05);
        EndCriteria ec(
            100, 20,   // max iterations, stationary iterations
            1e-4, 1e-4, 1e-4
        );
        // model is an object (HestonModel inherits from CalibrateModel)
        model->calibrate(helpers, opt, ec);
        Array params = model->params();
        Real sse = 0.0;
        Real maxAbs = 0.0;
        int n = 0;
        for(auto& h:helpers){
            Real err = h->calibrationError();
            sse = sse + err*err;
            maxAbs = max(maxAbs,fabs(err));
            n++;
        }
        HestonSurfaceFit out;
        out.thetaQ = params[0];
        out.kappaQ = params[1];
        out.xi     = params[2];
        out.rho    = params[3];
        out.v0     = params[4];
        out.rmseIv = sqrt(sse / max<Size>(1, n));
        out.maxAbsIvErr = maxAbs;
        return out;
    } 

    Real MultiSurfaceSSE(CallGrid& grid,HestonMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc,const Calendar& cal){
        
        //cout<<"MultiSurfaceSSE"<<endl;
        Real kappaQ_guess = phi.kappaQ;
        Real thetaQ_guess = phi.thetaQ;
        Real xi_guess     = phi.xi;
        Real rho_guess    = phi.rho;

        Settings::instance().evaluationDate() = grid.evalDate;
        Handle<YieldTermStructure> rTS(ext::make_shared<FlatForward>(grid.evalDate, grid.r,dc));
        Handle<YieldTermStructure> qTS(ext::make_shared<FlatForward>(grid.evalDate, grid.q, dc));
        Handle<Quote> spot(ext::make_shared<SimpleQuote>(grid.S0));

        // cout<<"phi.rho:"<<phi.rho<<endl;
        // cout<<"phi.kappaQ:"<<phi.kappaQ<<endl;
        // cout<<"phi.thetaQ:"<<phi.thetaQ<<endl;
        // cout<<"phi.xi:"<<phi.xi<<endl;
        // cout<<"v0_guess:"<<v0_guess<<endl;
        // cout<<"grid.S0:"<<grid.S0<<endl;
        
        
        // cout<<"Successfully spot."<<endl;
        auto process = ext::make_shared<HestonProcess>(
            rTS, qTS, spot,
            v0_guess, kappaQ_guess, thetaQ_guess, xi_guess, rho_guess
        );
        auto model  = ext::make_shared<HestonModel>(process);
        auto engine = ext::make_shared<AnalyticHestonEngine>(model);
        
        
        // make the process -> model -> engine; helpers a pointer that will store the actual value and calibration error -> so
        // we need a vector of helpers for M*K (all maturity and grid values)
        
        vector<ext::shared_ptr<CalibrationHelper>> helpers;
        const Size M = grid.maturities.size();
        const Size K = grid.strikes.size();
        helpers.reserve(M * K);

        Real sse = 0.0;
    
        for (int iT = 0; iT < grid.maturities.size();iT++){
            Integer days = cal.businessDaysBetween(grid.evalDate, grid.maturities[iT]);
            QL_REQUIRE(days > 0, "Maturity must be after today");
            Period maturityPeriod(days,Days);
            for(int iK = 0; iK < grid.strikes.size();iK++){
                Real strike = grid.strikes[iK];
                Real iv_mkt = grid.ImpliedVol[iT][iK];
                Handle<Quote> volQuote(ext::make_shared<SimpleQuote>(iv_mkt));
                auto h = ext::make_shared<HestonModelHelper>(
                    maturityPeriod,         // maturity as Period
                    cal,
                    grid.S0,                   // spot
                    strike,
                    volQuote,               // market implied vol quote
                    rTS,
                    qTS
                );
                h->setPricingEngine(engine);
                Real err = h->calibrationError();
                sse += (err*err);
            }
        }
        //cout<<"SSE:"<<sse<<endl;
        return sse;
    }

    Real bestv0ForSurface(CallGrid& grid, HestonMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc, const Calendar& cal){

        //cout<<"bestv0forSurface"<<endl;
        Real eps = 0.25;
        Real v0_min = max(0.001,v0_guess*(1-eps));
        Real v0_max = v0_guess*(1+eps);
        
        const int NScan = 30;
        Real bestv0 = v0_guess;
        Real bestVal = QL_MAX_REAL;
        for(int i = 0; i < NScan; i++){
            Real v0 = v0_min + (v0_max - v0_min)*Real(i)/Real(NScan);
            Real val = MultiSurfaceSSE(grid,phi,v0,dc,cal);
            if (val < bestVal){
                bestVal = val;
                bestv0 = v0;
            }
        }
        return bestv0;

    }

    HestonMultiSurfaceFit evaluatePhi(const vector<CallGrid>& surfaces,HestonMultiSurfaceFit& phi,
                                    const DayCounter& dc,const Calendar& cal){
        //cout<<"evaluatephi"<<endl;
        phi.TotalSSE = 0.0;
        random_device rd;                // seed source
        mt19937 gen(rd());               // Mersenne Twister RNG
        Real eps = 0.5;
        // possibly can be multi-threaded
        phi.v0_by_surface.resize(surfaces.size());
        for (int i = 0; i < surfaces.size(); i++){
            CallGrid grid = surfaces[i];
            uniform_real_distribution<> dist_v0(max(grid.v0 * (1-eps),0.001), grid.v0 * (1+eps));
            Real v0_init = dist_v0(gen);
            Real bestv0_surface = bestv0ForSurface(grid,phi,v0_init,dc,cal);
            phi.v0_by_surface[i] = bestv0_surface;
            Real t_sse = MultiSurfaceSSE(grid,phi,bestv0_surface,dc,cal);
            phi.TotalSSE += t_sse;
        }
        return phi;
    }

    HestonMultiSurfaceFit MultiSurfaceRandomSearch(const vector<CallGrid>& surfaces,const DayCounter& dc,const Calendar& cal,
                                                    const HestonMultiSurfaceFit& phi_init,int n_tries){
        //cout<<"multisurfacerandomsearch"<<endl;
        random_device rd;                // seed source
        mt19937 gen(rd());               // Mersenne Twister RNG
        Real eps = 0.5;
        uniform_real_distribution<> dist_kappaQ(phi_init.kappaQ*(1-eps),phi_init.kappaQ*(1+eps));
        uniform_real_distribution<> dist_thetaQ(phi_init.thetaQ*(1-eps),phi_init.thetaQ*(1+eps));
        uniform_real_distribution<> dist_xi(phi_init.xi*(1-eps),phi_init.xi*(1+eps));
        uniform_real_distribution<> dist_rho(max(-0.95, phi_init.rho - 0.4), min(-0.05, phi_init.rho + 0.4));
        
        HestonMultiSurfaceFit bestInitialFit;
        bestInitialFit.TotalSSE = QL_MAX_REAL;
        for (int i = 0; i < n_tries;i++){
            //cout<<"Try-"<< (i+1)<<endl;
            HestonMultiSurfaceFit candidate;
            candidate.kappaQ = dist_kappaQ(gen);
            candidate.thetaQ = dist_thetaQ(gen);
            candidate.xi = dist_xi(gen);
            candidate.rho = dist_rho(gen);
            auto newPhi = evaluatePhi(surfaces,candidate,dc,cal);
            if(newPhi.TotalSSE < bestInitialFit.TotalSSE){
                bestInitialFit = newPhi;
            }
        }
        return bestInitialFit;                                          
                                                    
    }

    struct MultiSurfaceObjective : public QuantLib::CostFunction{
    
    const vector<CallGrid>& surfaces;
    const DayCounter& dc;
    const Calendar& cal;

    // You may want these for the inner v0 scan bounds:
    Real v0_min = 1e-3;
    int  Nscan  = 30;   // inner scan resolution
    MultiSurfaceObjective(const std::vector<CallGrid>& s,const DayCounter& dc_,const Calendar& cal_):surfaces(s),dc(dc_),cal(cal_){}
    QuantLib::Real value(const QuantLib::Array& x) const override {
        Real a = x[0], b = x[1], c = x[2], d = x[3];
        HestonMultiSurfaceFit phi;
        phi.kappaQ = std::exp(a);
        phi.thetaQ = std::exp(b);
        phi.xi     = std::exp(c);
        phi.rho    = std::tanh(d); 
        
        auto cur = evaluatePhi(surfaces,phi,dc, cal);
        return cur.TotalSSE;
    }
    // Not needed for Simplex, but required by base class
    QuantLib::Array values(const QuantLib::Array& x) const override {
        QuantLib::Array y(1);
        y[0] = value(x);
        return y;
    }
};

HestonMultiSurfaceFit nedlerMeadMultiSurface(const vector<CallGrid>& surfaces,const DayCounter& dc,const Calendar& cal,const HestonMultiSurfaceFit& phi0){
    QuantLib::Array x0(4);
    x0[0] = std::log(std::max(Real(1e-12), phi0.kappaQ));
    x0[1] = std::log(std::max(Real(1e-12), phi0.thetaQ));
    x0[2] = std::log(std::max(Real(1e-12), phi0.xi));
    // atanh needs |rho|<1
    Real rho0 = std::min(Real(0.9999), std::max(Real(-0.9999), phi0.rho));
    x0[3] = 0.5 * std::log((1 + rho0) / (1 - rho0)); // atanh(rho0)
    
    MultiSurfaceObjective cost(surfaces, dc, cal);
    NoConstraint constraint;
    Problem problem(cost, constraint, x0);
    // Array step(4);
    // step[0] = 0.2;  // log-kappa
    // step[1] = 0.2;  // log-theta
    // step[2] = 0.2;  // log-xi
    // step[3] = 0.2;  // atanh-rho
    QuantLib::Real lambda = 0.1;
    Simplex solver(lambda);

    EndCriteria ec(
        300,    // maxIterations
        80,     // maxStationaryStateIterations
        1e-6,   // rootEpsilon (not super relevant here)
        1e-6,   // functionEpsilon
        1e-6    // gradientNormEpsilon (not used)
    );
    
    HestonMultiSurfaceFit out;

    EndCriteria::Type status = solver.minimize(problem, ec);
    QuantLib::Array xstar = problem.currentValue();
    out.kappaQ = exp(xstar[0]);
    out.thetaQ = exp(xstar[1]);
    out.xi     = exp(xstar[2]);
    out.rho    = tanh(xstar[3]);
    
    auto finalEval = evaluatePhi(surfaces,out,dc, cal);
    out.v0_by_surface = finalEval.v0_by_surface;
    out.TotalSSE      = finalEval.TotalSSE;
    return out;
}
}//namespace qe