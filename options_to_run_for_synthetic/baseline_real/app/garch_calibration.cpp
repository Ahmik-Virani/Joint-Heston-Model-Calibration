#include <qe/params.hpp>
#include <qe/path.hpp>
// #include <qe/pricing.hpp>
#include <qe/surface.hpp>
#include <qe/surfacefit.hpp>
#include <qe/garch.hpp>
#include <qe/latent_path_mcmc.hpp>
#include <qe/particle_filters.hpp>

#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>

#include "../src/getRealData.cpp"

using namespace qe;
using namespace QuantLib;
using namespace std;
using namespace Eigen;

struct calParams{
    HestonSurfaceFit SingleSurfaceParams;
    HestonMultiSurfaceFit MultSurfaceParams;
    HestonPParams meanP_garch_mcmc,varP_garch_mcmc;
    HestonPParams meanP_pmcmc,varP_pmcmc;
};


void writePathToCSV(const PPath& ppath,const string& filename){
    ofstream file(filename);
    if(!file.is_open()){
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << fixed << setprecision(8);
    file << "log S,v,log Returns\n";
    int N = ppath.returns.size();
    for (int i = 1; i < N;i++){
        file<<ppath.returns[i]<<"\n";
    }
    file.close();
}


void writePathToCSV(const PPath& ppath,const vector<double>hPath,const string& filename){
    ofstream file(filename);
    if(!file.is_open()){
        cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    file << fixed << setprecision(8);
    file << "log S,v,log Returns,garch h\n";
    int N = ppath.logS.size();
    for (int i = 0; i < N;i++){
        file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<ppath.returns[i] <<"," << hPath[i]<<"\n";
    }
    file.close();
}


PPath buildExpandingPath(getRealData& data, int end_t) {
    PPath ppath;

    ppath.logS.resize(end_t + 1);
    ppath.returns.resize(end_t);

    for (int i = 0; i <= end_t; i++) {
        ppath.logS[i] = log(data.get_S(i));
    }

    for (int i = 1; i <= end_t; ++i) {
        ppath.returns[i - 1] = ppath.logS[i] - ppath.logS[i - 1];
    }

    return ppath;
}

int main() {
    
    getRealData Data;
    Calendar cal = TARGET();
    DayCounter dc = Actual365Fixed();
    Size steps = 252;

    random_device rd;                // seed source
    mt19937 gen(rd());               // Mersenne Twister RNG
    Real eps = 0.5;


    int no_of_timesteps = Data.get_time_steps();
    calParams PARAMS;
    ofstream garch_calibration_params("./Test_Output/garch_calibration_params.csv");
    garch_calibration_params << "date,mu_mean,mu_var,kappa_mean,kappa_var,theta_mean,theta_var,vol-of-vol_mean,vol-of-vol_var,rho_mean,rho_var" << '\n';
    cout<<"Garch Calibration Started."<<endl;

    Real best_rmseIv = 10000;    
    auto guess_P = Data.get_guess();
    const double v0_guess = guess_P[6] * guess_P[6];

    
    double dt = 1/double(steps);
    int n_iters = 5000;
    int num_particles = 1500;

    int starting_steps = 10;

    if (no_of_timesteps <= 2) {
        cout << "Not enough data for GARCH/MCMC. Need at least 3 price points, got "
             << no_of_timesteps << "." << endl;
        return 0;
    }
    if(starting_steps >= no_of_timesteps){
        const int end_t = no_of_timesteps - 1;
        PPath ppath = buildExpandingPath(Data, end_t);
        GarchParams gParams = garchPathFit(ppath);
        cout << "\nGARCH fit using full path through t=" << end_t << endl;
        cout << gParams;
        vector<double> daily_hPath = getGarchPath(gParams, ppath);
        vector<double> annual_hPath(daily_hPath.size(), 0.0);
        for (int i = 0; i < static_cast<int>(daily_hPath.size()); i++) {
            annual_hPath[i] = daily_hPath[i] * steps;
        }
        VectorXd x0(5);

        uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
        uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
        uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
        uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
        uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));

        HestonPParams P{
            Data.get_S(0),        // S0 for this expanding path
            annual_hPath.front(), // v0: variance proxy at path start
            guess_P[0],           // mu
            guess_P[2],           // kappaP
            guess_P[3],           // thetaP
            guess_P[1],           // xi
            guess_P[4]            // rho
        };

        double mu0 = dist_mu(gen);
        double kappa0 = dist_kappaP(gen);
        double theta0 = dist_thetaP(gen);
        double xi0 = dist_xi_mcmc(gen);
        double rho0 = dist_rho_mcmc(gen);

        x0[0] = log(mu0);
        x0[1] = log(kappa0);
        x0[2] = log(theta0);
        x0[3] = log(xi0);
        x0[4] = atanh(rho0);

        HestonPParams meanP,varP;
        vector<double>vProxy = annual_hPath;

        mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
        cout<<"Mean Statistics:"<<meanP<<endl;
        cout<<"Variance Statistics:"<<varP<<endl;
        PARAMS.meanP_garch_mcmc = meanP;
        PARAMS.varP_garch_mcmc = varP;

        // [TODO] - check if end_t or end_t+1
        garch_calibration_params << Data.get_date(end_t) << ',' << meanP.mu << ',' << varP.mu << ',' << meanP.kappaP << ',' << varP.kappaP << ',' << meanP.thetaP << ',' << varP.thetaP << ',' << meanP.xi << ',' << varP.xi << ',' << meanP.rho << ',' << varP.rho << ',' << '\n';
    }

    // [TODO] - can this be rolling window also?
    else{
        for(int end_t = starting_steps; end_t < no_of_timesteps; end_t ++){
            PPath ppath = buildExpandingPath(Data, end_t);
            GarchParams gParams = garchPathFit(ppath);
            cout << "\nGARCH fit using expanding path through t=" << end_t << endl;
            cout << gParams;
            vector<double> daily_hPath = getGarchPath(gParams, ppath);
            vector<double> annual_hPath(daily_hPath.size(), 0.0);
            for (int i = 0; i < static_cast<int>(daily_hPath.size()); i++) {
                annual_hPath[i] = daily_hPath[i] * steps;
            }
            VectorXd x0(5);

            uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
            uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
            uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
            uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
            uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));

            HestonPParams P{
                Data.get_S(0),        // S0 for this expanding path
                annual_hPath.front(), // v0: variance proxy at path start
                guess_P[0],           // mu
                guess_P[2],           // kappaP
                guess_P[3],           // thetaP
                guess_P[1],           // xi
                guess_P[4]            // rho
            };
            
            double mu0 = dist_mu(gen);
            double kappa0 = dist_kappaP(gen);
            double theta0 = dist_thetaP(gen);
            double xi0 = dist_xi_mcmc(gen);
            double rho0 = dist_rho_mcmc(gen);

            x0[0] = log(mu0);
            x0[1] = log(kappa0);
            x0[2] = log(theta0);
            x0[3] = log(xi0);
            x0[4] = atanh(rho0);
            HestonPParams meanP,varP;
            vector<double>vProxy = annual_hPath;
            mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
            cout<<"Mean Statistics:"<<meanP<<endl;
            cout<<"Variance Statistics:"<<varP<<endl;
            PARAMS.meanP_garch_mcmc = meanP;
            PARAMS.varP_garch_mcmc = varP;

            // [TODO] - same as before
            garch_calibration_params << Data.get_date(end_t) << ',' << meanP.mu << ',' << varP.mu << ',' << meanP.kappaP << ',' << varP.kappaP << ',' << meanP.thetaP << ',' << varP.thetaP << ',' << meanP.xi << ',' << varP.xi << ',' << meanP.rho << ',' << varP.rho << ',' << '\n';
        }
    }
    garch_calibration_params.close();
    return 0;

}
