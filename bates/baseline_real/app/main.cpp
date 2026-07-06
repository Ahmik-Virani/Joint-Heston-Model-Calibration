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
    BatesSurfaceFit SingleSurfaceParams;
    BatesMultiSurfaceFit MultSurfaceParams;
    BatesPParams meanP_garch_mcmc,varP_garch_mcmc;
    BatesPParams meanP_pmcmc,varP_pmcmc;
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
        // file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<ppath.returns[i] << "\n";
        file<<ppath.returns[i]<<"\n";
        // if (i == 0){
        //     file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<0.0 << "\n";
        // }
        // else{
        //     file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<ppath.returns[i] << "\n";
        // }
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
        // if (i == 0){
        //     file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<0.0 << "\n";
        // }
        // else{
        //     file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<ppath.returns[i] << "\n";
        // }
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
    
    // Date today(26, February, 2026);
    // Settings::instance().evaluationDate() = today;

    // [TODO] - check if this target is correct
    Calendar cal = TARGET();

    // BatesPParams P{
    //     100.0,   // S0
    //     0.04,    // v0
    //     0.05,    // mu
    //     1.5,     // kappaP
    //     0.04,    // thetaP
    //     0.3,     // xi
    //     -0.7     // rho
    // };

    // VRPParams V{0.5};

    // Rate r = 0.03;
    // Rate q = 0.00;
    DayCounter dc = Actual365Fixed();
    // BatesQParams Q = toQ(P, V, r, q, dc);

    // cout << P << endl;
    // cout << Q << endl;


    // [TODO] check if this is for 1 year, or steps of 
    Size steps = 252;
    // Size seed  = 42;

    // PPath ppath = logReturns(P, steps, seed);
    // string filename = "./logReturns.csv";
    // writePathToCSV(ppath,filename);
    // cout<<"Path Values written in file."<<endl;

    // cout << ppath << endl;
    // int SurfaceFrequency = 10;
    // // build_surfaces_from_path(const BatesQParams& Q, const Date& today,
    // //     const Calender& cal, const PPath& ppath,int SurfaceFrequency);
    
    // vector<CallGrid> surfaces = build_surfaces_from_path(Q,today,cal,ppath,SurfaceFrequency,r,q);
    // cout << "Built " << surfaces.size() << " surfaces\n";
    // print_Surfaces(surfaces);
    // cout<<"-----------------------------------------------------------"<<endl;

    calParams PARAMS;


    // Single Surface Calibration Starts
    // CallGrid CalibrationSurface = surfaces[10]; // or any surfaces[i]
    // cout << "The grid chosen to be calibrated is:"<<endl;
    // print_Callgrid(CalibrationSurface);

    int repeat_cal = 10;
    vector<BatesSurfaceFit>SurfaceFitsParams;
    vector<BatesSurfaceFit>SurfaceFitGuesses;
    
    random_device rd;                // seed source
    mt19937 gen(rd());               // Mersenne Twister RNG
    Real eps = 0.5;

    // uniform_real_distribution<> dist_v0(CalibrationSurface.v0 * (1-eps), CalibrationSurface.v0 * (1+eps));
    // uniform_real_distribution<> dist_kappaQ(Q.kappaQ * (1-eps), Q.kappaQ * (1+eps));
    // uniform_real_distribution<> dist_thetaQ(Q.thetaQ * (1-eps), Q.thetaQ * (1+eps));
    // uniform_real_distribution<> dist_xi(Q.xi * (1-eps), Q.xi * (1+eps));
    // uniform_real_distribution<> dist_rho(max(-0.95, Q.rho - 0.4),min(-0.05, Q.rho + 0.4));
    
    int best_idx = -1;
    Real best_rmseIv = 10000;    

    auto guess_P = Data.get_guess();

    int no_of_timesteps = Data.get_time_steps();

    ofstream single_state_calibration_errors("./Output/single_state_errors.csv");
    ofstream single_state_calibration_params("./Output/single_state_calibration_params.csv");
    single_state_calibration_errors << "date,strike,maturity,true_price,computed_price,abs_error" << '\n';
    single_state_calibration_params << "date,kappa,theta,vol-of-vol,rho,v0" << '\n';
    for(int t = 0 ; t < no_of_timesteps ; t++){
        CallGrid CalibrationSurface;
        CalibrationSurface.C = Data.get_grid(t);
        CalibrationSurface.S0 = Data.get_S(t);
        CalibrationSurface.r = Data.get_r(t);
        CalibrationSurface.q = Data.get_q(t);
        CalibrationSurface.v0 = Data.get_guess()[6];
        
        // [TODO] check if eval data is todays date
        CalibrationSurface.evalDate = Data.get_date(t);
        
        for(int i = 0;i<repeat_cal;i++){
            BatesSurfaceFit InitialGuess;
            InitialGuess.v0     = guess_P[6];
            InitialGuess.kappaQ = guess_P[2] + guess_P[5];
            InitialGuess.thetaQ = (double)(guess_P[2] * guess_P[3]) / (double)(guess_P[2] + guess_P[5]);
            InitialGuess.xi     = guess_P[1];
            InitialGuess.rho    = guess_P[4];
            // [TODO] - check if these are correct
            InitialGuess.JumpIntensityQ = guess_P[7] * exp(guess_P[10] * guess_P[8] + 0.5 * guess_P[10] * guess_P[10] * guess_P[9] * guess_P[9]);
            InitialGuess.JumpMeanQ = guess_P[8] + guess_P[10] * guess_P[9] * guess_P[9];
            InitialGuess.JumpVolatility = guess_P[9];
            SurfaceFitGuesses.push_back(InitialGuess);
            
            BatesSurfaceFit SingleSurfaceParams = calibrateBatesQVolGrid(CalibrationSurface,InitialGuess,Data.get_date(t),dc,cal);
            SurfaceFitsParams.push_back(SingleSurfaceParams);
            //cout<<"v0 guess in single surface:"<<InitialGuess.v0 <<endl;
            if (SingleSurfaceParams.rmseIv < best_rmseIv){
                best_idx = i;
                best_rmseIv = SingleSurfaceParams.rmseIv;
            }
        }

        Data.get_penalty(t, SurfaceFitsParams[best_idx].v0, SurfaceFitsParams[best_idx].kappaQ, SurfaceFitsParams[best_idx].thetaQ, SurfaceFitsParams[best_idx].xi, SurfaceFitsParams[best_idx].rho, SurfaceFitsParams[best_idx].JumpIntensityQ, SurfaceFitsParams[best_idx].JumpMeanQ, SurfaceFitsParams[best_idx].JumpVolatility, single_state_calibration_errors);
        single_state_calibration_params << Data.get_date(t) << ',' << SurfaceFitsParams[best_idx].kappaQ << ',' << SurfaceFitsParams[best_idx].thetaQ << ',' << SurfaceFitsParams[best_idx].xi << ',' << SurfaceFitsParams[best_idx].rho << ',' << SurfaceFitsParams[best_idx].v0 << '\n';
    }
    // cout << "Initial Guess" <<endl;
    // cout<<SurfaceFitGuesses[best_idx]<<endl;
    // cout << "Final Parameters" <<endl;
    // cout<<SurfaceFitsParams[best_idx]<<endl;
    // PARAMS.SingleSurfaceParams = SurfaceFitsParams[best_idx];
    // cout << Q << endl;
    single_state_calibration_errors.close();
    single_state_calibration_params.close();
    cout<<"Single Surface Calibration Done."<<endl;
    // Single Surface Calibration Ends

    // Multi Surface Calibration Starts
    ofstream multi_state_calibration_errors("./Output/multi_state_errors.csv");
    ofstream multi_state_calibration_params("./Output/multi_state_calibration_params.csv");
    multi_state_calibration_errors << "date,strike,maturity,true_price,computed_price,abs_error" << '\n';
    multi_state_calibration_params << "date,kappa,theta,vol-of-vol,rho,v0" << '\n';

    BatesMultiSurfaceFit phi_init;
    phi_init.kappaQ = guess_P[2] + guess_P[5];
    phi_init.thetaQ = (double)(guess_P[2] * guess_P[3]) / (double)(guess_P[2] + guess_P[5]);
    phi_init.xi = guess_P[1];
    phi_init.rho = guess_P[4];
    // [TODO] - check if these are correct
    phi_init.JumpIntensityQ = guess_P[7] * exp(guess_P[10] * guess_P[8] + 0.5 * guess_P[10] * guess_P[10] * guess_P[9] * guess_P[9]);
    phi_init.JumpMeanQ = guess_P[8] + guess_P[10] * guess_P[9] * guess_P[9];
    phi_init.JumpVolatility = guess_P[9];
    // [TODO] - can finetune
    int few = 3;
    for(int t = 0 ; t < Data.get_time_steps()-few+1 ; t+=few){
        vector<CallGrid> few_surfaces;
        for(int i = t; i < few+t; i++){
            CallGrid this_CalibrationSurface;
            this_CalibrationSurface.C = Data.get_grid(i);
            this_CalibrationSurface.S0 = Data.get_S(i);
            this_CalibrationSurface.r = Data.get_r(i);
            this_CalibrationSurface.q = Data.get_q(i);
            this_CalibrationSurface.v0 = Data.get_guess()[6];
            
            // [TODO] check if eval data is todays date
            this_CalibrationSurface.evalDate = Data.get_date(i);
            few_surfaces.push_back(this_CalibrationSurface);
        }
        BatesMultiSurfaceFit best_phi_init = MultiSurfaceRandomSearch(few_surfaces,dc,cal,phi_init,30);
        BatesMultiSurfaceFit best_phi = nedlerMeadMultiSurface(few_surfaces,dc,cal,best_phi_init);
        PARAMS.MultSurfaceParams = best_phi;

        for(int i = t ; i < few+t ; i++){
            Data.get_penalty(i, PARAMS.MultSurfaceParams.v0_by_surface[i-t], PARAMS.MultSurfaceParams.kappaQ, PARAMS.MultSurfaceParams.thetaQ, PARAMS.MultSurfaceParams.xi, PARAMS.MultSurfaceParams.rho, PARAMS.MultSurfaceParams.JumpIntensityQ, PARAMS.MultSurfaceParams.JumpMeanQ, PARAMS.MultSurfaceParams.JumpVolatility, multi_state_calibration_errors);
            multi_state_calibration_params << Data.get_date(i-t) << ',' << PARAMS.MultSurfaceParams.kappaQ << ',' << PARAMS.MultSurfaceParams.thetaQ << ',' << PARAMS.MultSurfaceParams.xi << ',' << PARAMS.MultSurfaceParams.rho << ',' << PARAMS.MultSurfaceParams.v0_by_surface[i-t] << '\n';
        }
    }

    // cout<<best_phi_init<<endl;
    // cout<<best_phi<<endl;
    // cout<<"Ground Truth v0:"<<endl;
    // for(int i = 0; i < few_surfaces.size();i++){
    //     cout<<few_surfaces[i].v0<<",  ";
    // }

    // cout<<"\n";
    // cout << Q << endl;
    // cout<<"\n";
    multi_state_calibration_errors.close();
    multi_state_calibration_params.close();
    cout<<"Multi Surface Calibration Done."<<endl;
    // Multi Surface Calibration Ends
    
    //GARCH Calibration Starts
    ofstream garch_errors("./Output/garch_errors.csv");
    garch_errors << "date,strike,maturity,true_price,computed_price,abs_error" << '\n';

    // [TODO] - print P space errors
    ofstream garch_calibration_params("./Output/garch_calibration_params.csv");
    garch_calibration_params << "date,mu_mean,mu_var,kappa_mean,kappa_var,theta_mean,theta_var,vol-of-vol_mean,vol-of-vol_var,rho_mean,rho_var,vproxy" << '\n';

    // number of time stamps required
    // [TODO] - tune
    int prev_path_steps = 5;

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

        VectorXd x0(8);

        uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
        uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
        uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
        uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
        uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
        uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
        uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
        uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

        BatesPParams P{
            // [TODO] change to index
            Data.get_S(0),   // S0
            guess_P[6],    // v0 // we do not have a value for this, using guess
            guess_P[0],    // mu
            guess_P[2],     // kappaP
            guess_P[3],    // thetaP
            guess_P[1],     // xi
            guess_P[4],     // rho
            // [TODO] - check if these indexes are correct
            guess_P[7],     // JumpIntensityP
            guess_P[8],     // JumpMeanP
            guess_P[9]      // JumpVolatility
        };

        double mu0 = dist_mu(gen);
        double kappa0 = dist_kappaP(gen);
        double theta0 = dist_thetaP(gen);
        double xi0 = dist_xi_mcmc(gen);
        double rho0 = dist_rho_mcmc(gen);
        double JumpIntensityP0 = dist_JumpIntensityP(gen);
        double JumpMeanP0 = dist_JumpMeanP(gen);
        double JumpVolatility0 = dist_JumpVolatility(gen);

        // [TODO] - check if removing log is fine for mu
        x0[0] = mu0;
        x0[1] = log(kappa0);
        x0[2] = log(theta0);
        x0[3] = log(xi0);
        x0[4] = atanh(rho0);
        // [TODO] - check if fine
        x0[5] = log(JumpIntensityP0);
        x0[6] = JumpMeanP0;
        x0[7] = log(JumpVolatility0);


        BatesPParams meanP,varP;
        vector<double>vProxy = annual_hPath;

        mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
        cout<<"Mean Statistics:"<<meanP<<endl;
        cout<<"Variance Statistics:"<<varP<<endl;
        PARAMS.meanP_garch_mcmc = meanP;
        PARAMS.varP_garch_mcmc = varP;

        // [TODO] - check if end_t or end_t+1
        garch_calibration_params << Data.get_date(end_t) << ',' << meanP.mu << ',' << varP.mu << ',' << meanP.kappaP << ',' << varP.kappaP << ',' << meanP.thetaP << ',' << varP.thetaP << ',' << meanP.xi << ',' << varP.xi << meanP.rho << ',' << varP.JumpMeanP << ',' << meanP.JumpIntensityP << ',' << varP.JumpIntensityP << ',' << meanP.JumpMeanP << ',' << varP.JumpMeanP << ',' << meanP.JumpVolatility << ',' << varP.JumpVolatility << '\n';
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
            VectorXd x0(8);

            uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
            uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
            uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
            uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
            uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
            uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
            uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
            uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

            BatesPParams P{
                // [TODO] change to index
                Data.get_S(0),   // S0
                guess_P[6],    // v0 // we do not have a value for this, using guess
                guess_P[0],    // mu
                guess_P[2],     // kappaP
                guess_P[3],    // thetaP
                guess_P[1],     // xi
                guess_P[4],     // rho
                // [TODO] - check if these indexes are correct
                guess_P[7],     // JumpIntensityP
                guess_P[8],     // JumpMeanP
                guess_P[9]      // JumpVolatility
            };
            
            double mu0 = dist_mu(gen);
            double kappa0 = dist_kappaP(gen);
            double theta0 = dist_thetaP(gen);
            double xi0 = dist_xi_mcmc(gen);
            double rho0 = dist_rho_mcmc(gen);
            double JumpIntensityP0 = dist_JumpIntensityP(gen);
            double JumpMeanP0 = dist_JumpMeanP(gen);
            double JumpVolatility0 = dist_JumpVolatility(gen);

            // [TODO] - check if removing log is fine for mu
            x0[0] = mu0;
            x0[1] = log(kappa0);
            x0[2] = log(theta0);
            x0[3] = log(xi0);
            x0[4] = atanh(rho0);
            // [TODO] - check if fine
            x0[5] = log(JumpIntensityP0);
            x0[6] = JumpMeanP0;
            x0[7] = log(JumpVolatility0);

            BatesPParams meanP,varP;
            vector<double>vProxy = annual_hPath;
            mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
            cout<<"Mean Statistics:"<<meanP<<endl;
            cout<<"Variance Statistics:"<<varP<<endl;
            PARAMS.meanP_garch_mcmc = meanP;
            PARAMS.varP_garch_mcmc = varP;

            // [TODO] - same as before
            garch_calibration_params << Data.get_date(end_t) << ',' << meanP.mu << ',' << varP.mu << ',' << meanP.kappaP << ',' << varP.kappaP << ',' << meanP.thetaP << ',' << varP.thetaP << ',' << meanP.xi << ',' << varP.xi << ','<< meanP.rho << ',' << varP.rho << ',' << vProxy.back() << '\n';
        }
    }

    
    garch_calibration_params.close();
    cout<<"GARCH Surface Calibration Done."<<endl;

    // for(int i = prev_path_steps ; i < no_of_timesteps ; i++){
    //     PPath ppath;
    //     ppath.returns.resize(prev_path_steps+1);
    //     for(int timestep = i-prev_path_steps ; timestep <= i ; timestep++){
    //         ppath.returns[timestep] = Data.get_log_return(timestep);
    //     }
    //     GarchParams gParams = garchPathFit(ppath);
    //     // cout<<gParams<<endl;
    //     vector<double>daily_hPath = getGarchPath(gParams,ppath);
    //     vector<double>annual_hPath(daily_hPath.size(),0.0);
    //     for (int i = 0; i < daily_hPath.size(); i++){
    //         annual_hPath[i] = daily_hPath[i] * steps;
    //     }

    //     // cout<<"v path length:"<<ppath.v.size()<<endl;
    //     // cout<<"h path length:"<<annual_hPath.size()<<endl;

    //     // string filename_garch = "./returns_with_garch.csv";
    //     // writePathToCSV(ppath,annual_hPath,filename_garch);
    //     // cout<<"GARCH Values written in file."<<endl;

    //     VectorXd x0(9);
        
    //         // 0 - mu
    //         // 1 - sigma (vol-of-vol)
    //         // 2 - kappa
    //         // 3 - theta
    //         // 4 - rho
    //         // 5 - lambda
    //         // 6 - Instantaneous volatility
    //     uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
    //     uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
    //     uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
    //     uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
    //     uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
    //     uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
    //     uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
    //     uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

    //     x0[0] = dist_mu(gen);
    //     x0[1] = dist_kappaP(gen);
    //     x0[2] = dist_thetaP(gen);
    //     x0[3] = dist_xi_mcmc(gen);
    //     x0[4] = dist_rho_mcmc(gen);
    //     x0[6] = dist_JumpIntensityP(gen);
    //     x0[7] = dist_JumpMeanP(gen);
    //     x0[8] = dist_JumpVolatility(gen);

    //     // mcmcOverLatent(BatesPParams& P, PPath ppath,vector<double>vProxy,VectorXd x0,double dt,BatesPParams& meanP,BatesPParams& varP)

    //     BatesPParams P{
    //         // [TODO] change to index
    //         Data.get_S(i-prev_path_steps),   // S0
    //         guess_P[6],    // v0 // we do not have a value for this, using guess
    //         guess_P[0],    // mu
    //         guess_P[2],     // kappaP
    //         guess_P[3],    // thetaP
    //         guess_P[1],     // xi
    //         guess_P[4],     // rho
    //         // [TODO] - check if these indexes are correct
    //         guess_P[7],     // JumpIntensityP
    //         guess_P[8],     // JumpMeanP
    //         guess_P[9]      // JumpVolatility
    //     };
    //     cout << P << endl;
    //     BatesPParams meanP,varP;
    //     vector<double>vProxy = annual_hPath;
    //     // [TODO] - check, passing guess, is it fine?
    //     mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
    //     cout<<"Mean Statistics:"<<meanP<<endl;
    //     cout<<"Variance Statistics:"<<varP<<endl;

    //     PARAMS.meanP_garch_mcmc = meanP;
    //     PARAMS.varP_garch_mcmc = varP;
        
    // }

    // garch_errors.close();
    //GARCH Calibration Ends

    //pmcmcOverLatent(BatesPParams& P, PPath& ppath,VectorXd x0,double dt,BatesPParams& meanP,BatesPParams& varP)
    //PMCMC Calibration Starts

    ofstream pmcmc_calibration_errors("./Output/pmcmc_errors.csv");
    pmcmc_calibration_errors << "date,strike,maturity,true_price,computed_price,abs_error" << '\n';

    // [TODO] - print P space errors
    ofstream pmcmc_calibration_params("./Output/pmcmc_calibration_params.csv");
    pmcmc_calibration_params << "date,mu_mean,mu_var,kappa_mean,kappa_var,theta_mean,theta_var,vol-of-vol_mean,vol-of-vol_var,rho_mean,rho_var,v_t_mean,v_t_var" << '\n';

    if(starting_steps >= no_of_timesteps){
        const int end_t = no_of_timesteps - 1;
        PPath ppath = buildExpandingPath(Data, end_t);
        VectorXd x0(8);

        uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
        uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
        uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
        uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
        uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
        uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
        uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
        uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

        BatesPParams P{
            // [TODO] change to index
            Data.get_S(0),   // S0
            guess_P[6],    // v0 // we do not have a value for this, using guess
            guess_P[0],    // mu
            guess_P[2],     // kappaP
            guess_P[3],    // thetaP
            guess_P[1],     // xi
            guess_P[4],     // rho
            // [TODO] - check if these indexes are correct
            guess_P[7],     // JumpIntensityP
            guess_P[8],     // JumpMeanP
            guess_P[9]      // JumpVolatility
        };
        
        double mu0 = dist_mu(gen);
        double kappa0 = dist_kappaP(gen);
        double theta0 = dist_thetaP(gen);
        double xi0 = dist_xi_mcmc(gen);
        double rho0 = dist_rho_mcmc(gen);
        double JumpIntensityP0 = dist_JumpIntensityP(gen);
        double JumpMeanP0 = dist_JumpMeanP(gen);
        double JumpVolatility0 = dist_JumpVolatility(gen);

        // [TODO] - check if removing log is fine for mu
        x0[0] = mu0;
        x0[1] = log(kappa0);
        x0[2] = log(theta0);
        x0[3] = log(xi0);
        x0[4] = atanh(rho0);
        // [TODO] - check if fine
        x0[5] = log(JumpIntensityP0);
        x0[6] = JumpMeanP0;
        x0[7] = log(JumpVolatility0);   

        BatesPParams meanP_pmcmc,varP_pmcmc;
        pmcmcOverLatent(P,ppath,x0,dt,meanP_pmcmc,varP_pmcmc,n_iters,num_particles);
        cout<<"Mean Statistics:"<<meanP_pmcmc<<endl;
        cout<<"Variance Statistics:"<<varP_pmcmc<<endl;
        filterValues finalFilter = ParticleFilter(meanP_pmcmc,P.v0,ppath,dt,num_particles);

        int N_pmcmc = finalFilter.particles[0].size();
        int T_pmcmc = finalFilter.particles.size();
        
        int numSampledPaths = 200;
        vector<vector<double>> sampledPaths(numSampledPaths,vector<double>(T_pmcmc,0));


        for(int i = 0; i < numSampledPaths;i++){
            vector<double>sample = ancestralSampling(finalFilter,gen);
            sampledPaths[i] = sample;
        }

        writeSampledPaths(sampledPaths);
        cout<<"Sampled Paths from PMCMC Written."<<endl; 

        PARAMS.meanP_pmcmc = meanP_pmcmc;
        PARAMS.varP_pmcmc = varP_pmcmc;

        // [TODO] - same end_t
        pmcmc_calibration_params << Data.get_date(end_t) << ',' << meanP_pmcmc.mu << ',' << varP_pmcmc.mu << ',' << meanP_pmcmc.kappaP << ',' << varP_pmcmc.kappaP << ',' << meanP_pmcmc.thetaP << ',' << varP_pmcmc.thetaP << ',' << meanP_pmcmc.xi << ',' << varP_pmcmc.xi << meanP_pmcmc.rho << ',' << varP_pmcmc.rho << ',' << '\n';

    }
    else{
        for(int end_t = starting_steps; end_t < no_of_timesteps; end_t ++){
            PPath ppath = buildExpandingPath(Data, end_t);
            VectorXd x0(8);

            uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
            uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
            uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
            uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
            uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
            uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
            uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
            uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

            BatesPParams P{
                // [TODO] change to index
                Data.get_S(0),   // S0
                guess_P[6],    // v0 // we do not have a value for this, using guess
                guess_P[0],    // mu
                guess_P[2],     // kappaP
                guess_P[3],    // thetaP
                guess_P[1],     // xi
                guess_P[4],     // rho
                // [TODO] - check if these indexes are correct
                guess_P[7],     // JumpIntensityP
                guess_P[8],     // JumpMeanP
                guess_P[9]      // JumpVolatility
            };
            
            double mu0 = dist_mu(gen);
            double kappa0 = dist_kappaP(gen);
            double theta0 = dist_thetaP(gen);
            double xi0 = dist_xi_mcmc(gen);
            double rho0 = dist_rho_mcmc(gen);
            double JumpIntensityP0 = dist_JumpIntensityP(gen);
            double JumpMeanP0 = dist_JumpMeanP(gen);
            double JumpVolatility0 = dist_JumpVolatility(gen);

            // [TODO] - check if removing log is fine for mu
            x0[0] = mu0;
            x0[1] = log(kappa0);
            x0[2] = log(theta0);
            x0[3] = log(xi0);
            x0[4] = atanh(rho0);
            // [TODO] - check if fine
            x0[5] = log(JumpIntensityP0);
            x0[6] = JumpMeanP0;
            x0[7] = log(JumpVolatility0);

            BatesPParams meanP_pmcmc,varP_pmcmc;
            pmcmcOverLatent(P,ppath,x0,dt,meanP_pmcmc,varP_pmcmc,n_iters,num_particles);
            cout<<"Mean Statistics:"<<meanP_pmcmc<<endl;
            cout<<"Variance Statistics:"<<varP_pmcmc<<endl;
            filterValues finalFilter = ParticleFilter(meanP_pmcmc,P.v0,ppath,dt,num_particles);

            int N_pmcmc = finalFilter.particles[0].size();
            int T_pmcmc = finalFilter.particles.size();

            int last_time = T_pmcmc - 1;
            double vt_mean = 0.0;
            for(double v:finalFilter.particles[last_time]){
                vt_mean += v;
            }
            vt_mean /= finalFilter.particles[last_time].size();
            double vt_var = 0.0;
            for(double v:finalFilter.particles[last_time]){
                vt_var += (v - vt_mean) * (v - vt_mean);
            }
            vt_var /= finalFilter.particles[last_time].size();
            
            int numSampledPaths = 200;
            vector<vector<double>> sampledPaths(numSampledPaths,vector<double>(T_pmcmc,0));


            for(int i = 0; i < numSampledPaths;i++){
                vector<double>sample = ancestralSampling(finalFilter,gen);
                sampledPaths[i] = sample;
            }

            writeSampledPaths(sampledPaths);
            cout<<"Sampled Paths from PMCMC Written."<<endl; 

            PARAMS.meanP_pmcmc = meanP_pmcmc;
            PARAMS.varP_pmcmc = varP_pmcmc;

            // [TODO] - same end_t
            pmcmc_calibration_params << Data.get_date(end_t) << ',' << meanP_pmcmc.mu << ',' << varP_pmcmc.mu << ',' << meanP_pmcmc.kappaP << ',' << varP_pmcmc.kappaP << ',' << meanP_pmcmc.thetaP << ',' << varP_pmcmc.thetaP << ',' << meanP_pmcmc.xi << ',' << varP_pmcmc.xi << ',' << meanP_pmcmc.rho << ',' << varP_pmcmc.rho << ',' <<vt_mean<<',' <<vt_var << '\n';
        }
    }

    pmcmc_calibration_params.close();
    cout<<"PMCMC Calibration Done."<<endl;
    
    // for(int i = prev_path_steps ; i < no_of_timesteps ; i++){
    //     PPath ppath;
    //     ppath.returns.resize(prev_path_steps+1);
    //     for(int timestep = i-prev_path_steps ; timestep <= i ; timestep++){
    //         ppath.returns[timestep] = Data.get_log_return(timestep);
    //     }

    //     VectorXd x0(9);
        
    //         // 0 - mu
    //         // 1 - sigma (vol-of-vol)
    //         // 2 - kappa
    //         // 3 - theta
    //         // 4 - rho
    //         // 5 - lambda
    //         // 6 - Instantaneous volatility
    //     uniform_real_distribution<> dist_mu(guess_P[0] * (1-eps), guess_P[0] * (1+eps));
    //     uniform_real_distribution<> dist_kappaP(guess_P[2] * (1-eps), guess_P[2] * (1+eps));
    //     uniform_real_distribution<> dist_thetaP(guess_P[3] * (1-eps), guess_P[3] * (1+eps));
    //     uniform_real_distribution<> dist_xi_mcmc(guess_P[1] * (1-eps), guess_P[1] * (1+eps));
    //     uniform_real_distribution<> dist_rho_mcmc(max(-0.95, guess_P[4] - 0.4),min(-0.05, guess_P[4] + 0.4));
    //     uniform_real_distribution<> dist_JumpIntensityP(guess_P[7] * (1-eps), guess_P[7] * (1+eps));
    //     uniform_real_distribution<> dist_JumpMeanP(guess_P[8] * (1-eps), guess_P[8] * (1+eps));
    //     uniform_real_distribution<> dist_JumpVolatility(guess_P[9] * (1-eps), guess_P[9] * (1+eps));

    //     x0[0] = dist_mu(gen);
    //     x0[1] = dist_kappaP(gen);
    //     x0[2] = dist_thetaP(gen);
    //     x0[3] = dist_xi_mcmc(gen);
    //     x0[4] = dist_rho_mcmc(gen);
    //     x0[6] = dist_JumpIntensityP(gen);
    //     x0[7] = dist_JumpMeanP(gen);
    //     x0[8] = dist_JumpVolatility(gen);

    //     BatesPParams P{
    //         // [TODO] change to index
    //         Data.get_S(i-prev_path_steps),   // S0
    //         guess_P[6],    // v0 // we do not have a value for this, using guess
    //         guess_P[0],    // mu
    //         guess_P[2],     // kappaP
    //         guess_P[3],    // thetaP
    //         guess_P[1],     // xi
    //         guess_P[4],     // rho
    //         // [TODO] - check if these indexes are correct
    //         guess_P[7],     // JumpIntensityP
    //         guess_P[8],     // JumpMeanP
    //         guess_P[9]      // JumpVolatility
    //     };

    //     BatesPParams meanP_pmcmc,varP_pmcmc;
    //     pmcmcOverLatent(P,ppath,x0,dt,meanP_pmcmc,varP_pmcmc,n_iters,num_particles);
    //     cout<<"Mean Statistics:"<<meanP_pmcmc<<endl;
    //     cout<<"Variance Statistics:"<<varP_pmcmc<<endl;
    //     filterValues finalFilter = ParticleFilter(meanP_pmcmc,P.v0,ppath,dt,num_particles);

    //     int N_pmcmc = finalFilter.particles[0].size();
    //     int T_pmcmc = finalFilter.particles.size();
        
    //     int numSampledPaths = 200;
    //     vector<vector<double>> sampledPaths(numSampledPaths,vector<double>(T_pmcmc,0));


    //     for(int i = 0; i < numSampledPaths;i++){
    //         vector<double>sample = ancestralSampling(finalFilter,gen);
    //         sampledPaths[i] = sample;
    //     }
        
    //     writeSampledPaths(sampledPaths);
    //     cout<<"Sampled Paths from PMCMC Written."<<endl; 

    //     PARAMS.meanP_pmcmc = meanP_pmcmc;
    //     PARAMS.varP_pmcmc = varP_pmcmc;
    // }

    pmcmc_calibration_errors.close();
    //PMCMC Calibration Ends
    // cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    // cout<<"Final Prints"<<endl;
    // cout<<P<<endl;
    // cout<<Q<<endl;
    // cout<<PARAMS.SingleSurfaceParams<<endl;
    // cout<<PARAMS.MultSurfaceParams<<endl;
    // cout<<PARAMS.meanP_garch_mcmc<<endl;
    // cout<<PARAMS.varP_garch_mcmc<<endl;
    // cout<<PARAMS.meanP_pmcmc<<endl;
    // cout<<PARAMS.varP_pmcmc<<endl;
    // cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

    return 0;
}
