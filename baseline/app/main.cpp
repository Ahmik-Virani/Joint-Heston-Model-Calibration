#include <qe/params.hpp>
#include <qe/path.hpp>
#include <qe/pricing.hpp>
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
    int N = ppath.logS.size();
    for (int i = 0; i < N;i++){
        file<<ppath.logS[i]<<","<<ppath.v[i]<<","<<ppath.returns[i] << "\n";
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

int main() {
    
    Date today(26, February, 2026);
    Settings::instance().evaluationDate() = today;
    Calendar cal = TARGET();

    HestonPParams P{
        100.0,   // S0
        0.04,    // v0
        0.05,    // mu
        1.5,     // kappaP
        0.04,    // thetaP
        0.3,     // xi
        -0.7     // rho
    };

    VRPParams V{0.5};

    Rate r = 0.03;
    Rate q = 0.00;
    DayCounter dc = Actual365Fixed();
    HestonQParams Q = toQ(P, V, r, q, dc);

    cout << P << endl;
    cout << Q << endl;


    Size steps = 252;
    Size seed  = 42;

    PPath ppath = logReturns(P, steps, seed);
    string filename = "./logReturns.csv";
    writePathToCSV(ppath,filename);
    cout<<"Path Values written in file."<<endl;

    cout << ppath << endl;
    int SurfaceFrequency = 10;
    // build_surfaces_from_path(const HestonQParams& Q, const Date& today,
    //     const Calender& cal, const PPath& ppath,int SurfaceFrequency);
    
    vector<CallGrid> surfaces = build_surfaces_from_path(Q,today,cal,ppath,SurfaceFrequency,r,q);
    cout << "Built " << surfaces.size() << " surfaces\n";
    print_Surfaces(surfaces);
    cout<<"-----------------------------------------------------------"<<endl;

    calParams PARAMS;


    // Single Surface Calibration Starts
    CallGrid CalibrationSurface = surfaces[10]; // or any surfaces[i]
    cout << "The grid chosen to be calibrated is:"<<endl;
    print_Callgrid(CalibrationSurface);

    int repeat_cal = 10;
    vector<HestonSurfaceFit>SurfaceFitsParams;
    vector<HestonSurfaceFit>SurfaceFitGuesses;
    
    random_device rd;                // seed source
    mt19937 gen(rd());               // Mersenne Twister RNG
    Real eps = 0.5;

    uniform_real_distribution<> dist_v0(CalibrationSurface.v0 * (1-eps), CalibrationSurface.v0 * (1+eps));
    uniform_real_distribution<> dist_kappaQ(Q.kappaQ * (1-eps), Q.kappaQ * (1+eps));
    uniform_real_distribution<> dist_thetaQ(Q.thetaQ * (1-eps), Q.thetaQ * (1+eps));
    uniform_real_distribution<> dist_xi(Q.xi * (1-eps), Q.xi * (1+eps));
    uniform_real_distribution<> dist_rho(max(-0.95, Q.rho - 0.4),min(-0.05, Q.rho + 0.4));
    
    int best_idx = -1;
    Real best_rmseIv = 10000;    
    for(int i = 0;i<repeat_cal;i++){
        HestonSurfaceFit InitialGuess;
        InitialGuess.v0     = dist_v0(gen);
        InitialGuess.kappaQ = dist_kappaQ(gen);
        InitialGuess.thetaQ = dist_thetaQ(gen);
        InitialGuess.xi     = dist_xi(gen);
        InitialGuess.rho    = dist_rho(gen);
        SurfaceFitGuesses.push_back(InitialGuess);
        HestonSurfaceFit SingleSurfaceParams = calibrateHestonQVolGrid(CalibrationSurface,InitialGuess,today,dc,cal);
        SurfaceFitsParams.push_back(SingleSurfaceParams);
        //cout<<"v0 guess in single surface:"<<InitialGuess.v0 <<endl;
        if (SingleSurfaceParams.rmseIv < best_rmseIv){
            best_idx = i;
            best_rmseIv = SingleSurfaceParams.rmseIv;
        }
    }
    cout << "Initial Guess" <<endl;
    cout<<SurfaceFitGuesses[best_idx]<<endl;
    cout << "Final Parameters" <<endl;
    cout<<SurfaceFitsParams[best_idx]<<endl;
    PARAMS.SingleSurfaceParams = SurfaceFitsParams[best_idx];
    cout << Q << endl;
    // Single Surface Calibration Ends

    // Multi Surface Calibration Starts
    HestonMultiSurfaceFit phi_init;
    phi_init.kappaQ = dist_kappaQ(gen);
    phi_init.thetaQ = dist_thetaQ(gen);
    phi_init.xi = dist_xi(gen);
    phi_init.rho = dist_rho(gen);
    vector<CallGrid> few_surfaces;
    int few = 5;
    for(int i = 0; i < few; i++){
        few_surfaces.push_back(surfaces[i]);
    }

    HestonMultiSurfaceFit best_phi_init = MultiSurfaceRandomSearch(few_surfaces,dc,cal,phi_init,30);
    cout<<best_phi_init<<endl;
    HestonMultiSurfaceFit best_phi = nedlerMeadMultiSurface(few_surfaces,dc,cal,best_phi_init);
    cout<<best_phi<<endl;
    PARAMS.MultSurfaceParams = best_phi;
    cout<<"Ground Truth v0:"<<endl;
    for(int i = 0; i < few_surfaces.size();i++){
        cout<<few_surfaces[i].v0<<",  ";
    }

    cout<<"\n";
    cout << Q << endl;
    cout<<"\n";
    // Multi Surface Calibration Ends
    
    //GARCH Calibration Starts
    GarchParams gParams = garchPathFit(ppath);
    cout<<gParams<<endl;
    vector<double>daily_hPath = getGarchPath(gParams,ppath);
    vector<double>annual_hPath(daily_hPath.size(),0.0);
    for (int i = 0; i < daily_hPath.size(); i++){
        annual_hPath[i] = daily_hPath[i] * steps;
    }

    // cout<<"v path length:"<<ppath.v.size()<<endl;
    // cout<<"h path length:"<<annual_hPath.size()<<endl;

    string filename_garch = "./returns_with_garch.csv";
    writePathToCSV(ppath,annual_hPath,filename_garch);
    cout<<"GARCH Values written in file."<<endl;

    VectorXd x0(5);
    uniform_real_distribution<> dist_mu(P.mu * (1-eps), P.mu * (1+eps));
    uniform_real_distribution<> dist_kappaP(P.kappaP * (1-eps), P.kappaP * (1+eps));
    uniform_real_distribution<> dist_thetaP(P.thetaP * (1-eps), P.thetaP * (1+eps));
    uniform_real_distribution<> dist_xi_mcmc(P.xi * (1-eps), P.xi * (1+eps));
    uniform_real_distribution<> dist_rho_mcmc(max(-0.95, P.rho - 0.4),min(-0.05, P.rho + 0.4));

    x0[0] = dist_mu(gen);
    x0[1] = dist_kappaP(gen);
    x0[2] = dist_thetaP(gen);
    x0[3] = dist_xi_mcmc(gen);
    x0[4] = dist_rho_mcmc(gen);
    double dt = 1/double(steps);
    int n_iters = 5000;
    int num_particles = 1500;

    // mcmcOverLatent(HestonPParams& P, PPath ppath,vector<double>vProxy,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP)
    HestonPParams meanP,varP;
    vector<double>vProxy = annual_hPath;
    mcmcOverLatent(P,ppath,vProxy,x0,dt,meanP,varP,n_iters);
    cout<<"Mean Statistics:"<<meanP<<endl;
    cout<<"Variance Statistics:"<<varP<<endl;

    PARAMS.meanP_garch_mcmc = meanP;
    PARAMS.varP_garch_mcmc = varP;


    //GARCH Calibration Ends

    //pmcmcOverLatent(HestonPParams& P, PPath& ppath,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP)
    //PMCMC Calibration Starts
    HestonPParams meanP_pmcmc,varP_pmcmc;
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
    //PMCMC Calibration Ends
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;
    cout<<"Final Prints"<<endl;
    cout<<P<<endl;
    cout<<Q<<endl;
    cout<<PARAMS.SingleSurfaceParams<<endl;
    cout<<PARAMS.MultSurfaceParams<<endl;
    cout<<PARAMS.meanP_garch_mcmc<<endl;
    cout<<PARAMS.varP_garch_mcmc<<endl;
    cout<<PARAMS.meanP_pmcmc<<endl;
    cout<<PARAMS.varP_pmcmc<<endl;
    cout<<"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"<<endl;

    return 0;
}

