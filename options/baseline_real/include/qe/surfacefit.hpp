#pragma once
#include <ql/quantlib.hpp>
#include <vector>
#include "qe/params.hpp"
#include "qe/surface.hpp"
#include <vector>

namespace qe{

    using QuantLib::Date;
    using QuantLib::DayCounter;
    using QuantLib::Calendar;
    using QuantLib::Real;
    using std::vector;

    struct HestonSurfaceFit{
        Real v0;
        Real kappaQ;
        Real thetaQ;
        Real xi;
        Real rho;
        Real rmseIv;
        Real maxAbsIvErr;
    
    };
    struct HestonMultiSurfaceFit{
        Real kappaQ;
        Real thetaQ;
        Real xi;
        Real rho;
        Real TotalSSE;
        vector<Real>v0_by_surface;
    };

    HestonSurfaceFit calibrateHestonQVolGrid(CallGrid& grid,HestonSurfaceFit& InitialGuess,const Date& today, const DayCounter& dc,const Calendar& cal);
    Real MultiSurfaceSSE(CallGrid& grid,HestonMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc,const Calendar& cal);
    Real bestv0ForSurface(CallGrid& grid, HestonMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc, const Calendar& cal);
    HestonMultiSurfaceFit evaluatePhi(const vector<CallGrid>& surfaces,HestonMultiSurfaceFit& phi,
        const DayCounter& dc,const Calendar& cal);
    HestonMultiSurfaceFit MultiSurfaceRandomSearch(const vector<CallGrid>& surfaces,const DayCounter& dc,const Calendar& cal,
            const HestonMultiSurfaceFit& phi_init,int n_tries);
            HestonMultiSurfaceFit nedlerMeadMultiSurface(const vector<CallGrid>& surfaces,const DayCounter& dc,
                const Calendar& cal,const HestonMultiSurfaceFit& phi0);
    std::ostream& operator<<(std::ostream& os, const HestonSurfaceFit& fit);
    std::ostream& operator<<(std::ostream& os, const HestonMultiSurfaceFit& fit);
    

}