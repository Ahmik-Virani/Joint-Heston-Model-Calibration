#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include "qe/pricing.hpp"
#include "qe/surface.hpp"
#include "qe/path.hpp"
#include "qe/params.hpp"

using namespace QuantLib;
using namespace std;

namespace qe{
void print_Surfaces(const vector<CallGrid>& surfaces) {

    std::cout << "\n========== Heston Call Surfaces ==========\n";
    std::cout << "Total surfaces: " << surfaces.size() << "\n\n";

    for (Size i = 0; i < 5; ++i) {
        std::cout << "------------------------------------------\n";
        std::cout << "Surface ID: " << i << "\n";
        std::cout << "------------------------------------------\n";
        print_Callgrid(surfaces[i]);
        std::cout << "\n";
    }
}

vector<CallGrid> build_surfaces_from_path(const HestonQParams& Q, const Date& today,
                                        const Calendar& cal, const PPath& ppath,int SurfaceFrequency,Rate r,Rate q){
    vector<CallGrid> surfaces;
    cout<<"logS size:"<<ppath.logS.size() << endl;
    for(int i = 0; i < ppath.logS.size();i = i+SurfaceFrequency){
        
        Date evalDate = cal.advance(today, i, Days);
        
        CallGrid grid;
        
        grid.evalDate = evalDate;
        grid.maturities = {
            cal.advance(evalDate,1,Months),
            cal.advance(evalDate,3,Months),
            cal.advance(evalDate,6,Months),
            cal.advance(evalDate,12,Months)
        };
    
        grid.strikes = {80,90,100,110,120};
        
        Real S_today = exp(ppath.logS[i]);
        Real v_today = ppath.v[i];
        
        grid.S0 = exp(ppath.logS[i]);
        grid.v0 = ppath.v[i];
        grid.r = r;
        grid.q = q;

        callgrid(Q,grid);
        surfaces.push_back(grid);
    }
    return surfaces;
}
} //namespace qe
