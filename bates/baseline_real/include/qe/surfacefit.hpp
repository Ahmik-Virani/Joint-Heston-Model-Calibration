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

    struct BatesSurfaceFit{
        Real v0;
        Real kappaQ;
        Real thetaQ;
        Real xi;
        Real rho;
        Real rmseIv;
        Real maxAbsIvErr;

        Real JumpIntensityQ;
        Real JumpMeanQ;
        Real JumpVolatility;
    
    };
    struct BatesMultiSurfaceFit{
        Real kappaQ;
        Real thetaQ;
        Real xi;
        Real rho;
        Real TotalSSE;
        vector<Real>v0_by_surface;

        Real JumpIntensityQ;
        Real JumpMeanQ;
        Real JumpVolatility;
    };

    BatesSurfaceFit calibrateBatesQVolGrid(CallGrid& grid,BatesSurfaceFit& InitialGuess,const Date& today, const DayCounter& dc,const Calendar& cal);
    Real MultiSurfaceSSE(CallGrid& grid,BatesMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc,const Calendar& cal);
    Real bestv0ForSurface(CallGrid& grid, BatesMultiSurfaceFit& phi,Real v0_guess,const DayCounter& dc, const Calendar& cal);
    BatesMultiSurfaceFit evaluatePhi(const vector<CallGrid>& surfaces,BatesMultiSurfaceFit& phi,
        const DayCounter& dc,const Calendar& cal);
    BatesMultiSurfaceFit MultiSurfaceRandomSearch(const vector<CallGrid>& surfaces,const DayCounter& dc,const Calendar& cal,
            const BatesMultiSurfaceFit& phi_init,int n_tries);
            BatesMultiSurfaceFit nedlerMeadMultiSurface(const vector<CallGrid>& surfaces,const DayCounter& dc,
                const Calendar& cal,const BatesMultiSurfaceFit& phi0);
    std::ostream& operator<<(std::ostream& os, const BatesSurfaceFit& fit);
    std::ostream& operator<<(std::ostream& os, const BatesMultiSurfaceFit& fit);
    

}