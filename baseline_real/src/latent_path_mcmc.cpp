#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <Eigen/Core>
#include <iomanip>
#include <limits>
#include <Eigen/Cholesky>
#include "qe/params.hpp"
#include "qe/path.hpp"
#include "qe/latent_path_mcmc.hpp"
#include "qe/particle_filters.hpp"


using namespace std;
using namespace Eigen;


namespace qe {
double getPathLikelihood(const HestonPParams& P, const PPath& ppath,const vector<double>& vProxy,double dt){
    int n = ppath.returns.size();
    if (ppath.returns[n - 1] == 0 ) n = n - 1;
    //cout<<"Size of n:"<<n<<endl;
    vector<double> vDiff(vProxy.size() - 1,0);
    for(int i = 0; i < vProxy.size() - 1; i++){
        vDiff[i] = vProxy[i+1] - vProxy[i];
    }
    //cout<<"v diff size:"<<vDiff.size()<<endl;
    
    double mu = P.mu;
    double kappaP = P.kappaP;
    double thetaP = P.thetaP;
    double xi = P.xi;
    double rho = P.rho;
    const double pi = 3.141592;
    //double dt = 1.0/252;
    double mean_t_exp = 0.0;
    double logLikelihood  = 0.0;
    for (int i = 0; i < n; i++){
        Vector2d t_m;
        Matrix2d t_A_inverse;
        Vector2d del_t_y;

        t_m << (mu - 0.5 * vProxy[i]) * dt ,kappaP * (thetaP - vProxy[i]) * dt;
        double denom = vProxy[i] * (1 - rho*rho) * dt;
        t_A_inverse << 1/denom, -rho/(denom * xi),-rho/(denom * xi), 1/(xi * xi * denom);
        //t_A << vProxy[i],vProxy[i] * rho * xi,vProxy[i] * rho * xi,vProxy[i] * rho * rho;
        del_t_y << ppath.returns[i],vDiff[i];
        double t_exp = (del_t_y - t_m).dot(t_A_inverse * (del_t_y - t_m));
        double determ = (vProxy[i]*vProxy[i] * dt * dt) * (xi * xi) * (1 - rho*rho);
        logLikelihood += (-log(2*pi) - 0.5 * log(determ) - 0.5 * t_exp);
        mean_t_exp += t_exp;
    }
    mean_t_exp = mean_t_exp/n;
    //cout<<"Mean tExp:"<<mean_t_exp<<endl;
    return logLikelihood;
}

HestonPParams xToParamsMCMC(VectorXd& x){
    HestonPParams P;
    P.mu = exp(x[0]);
    P.kappaP = exp(x[1]);
    P.thetaP = exp(x[2]);
    P.xi = exp(x[3]);
    P.rho = tanh(x[4]);
    return P;
}

VectorXd ParamstoxMCMC(HestonPParams& P){
    VectorXd x;
    x[0] = log(P.mu);
    x[1] = log(P.kappaP);
    x[2] = log(P.thetaP);
    x[3] = log(P.xi);
    x[4] = atanh(P.rho);
    return x;
}

static inline double log_norm_prior(double x,double m, double s){
    double z = (m - x)/s;
    return -0.5 * log(2*M_PI) - log(s) - 0.5 * z*z; 

}

double logPrior(VectorXd x,double vbar){
    double lp = 0.0;
    lp += log_norm_prior(x[0], 0.0, 1.0);
    lp += log_norm_prior(x[1], std::log(1.0), 1.0);
    lp += log_norm_prior(x[2], std::log(vbar), 1.0);
    lp += log_norm_prior(x[3], std::log(0.3), 1.0);
    lp += log_norm_prior(x[4], 0.0, 1.5);
    return lp;
}

double logPosterior(VectorXd x,double vbar, const PPath& ppath,const vector<double>& vProxy,double dt){
    double lp = logPrior(x,vbar);
    HestonPParams P = xToParamsMCMC(x);
    double ll = getPathLikelihood(P,ppath,vProxy,dt);
    return ll + lp;  
}

double logPosterior(VectorXd x,double vbar, const PPath& ppath,double v0_actual,double dt,int N){
    double lp = logPrior(x,vbar);
    HestonPParams P = xToParamsMCMC(x);
    //double ll = getPathLikelihood(P,ppath);
    filterValues f = ParticleFilter(P,v0_actual,ppath,dt,N);
    double ll = f.likelihood;
    return ll + lp;  
}


VectorXd randn_vec(mt19937& rng, int d){
    static thread_local normal_distribution<double>nd(0.0,1.0);
    VectorXd z(d);
    for (int i = 0; i < d; i++) z[i] = nd(rng);
    return z;
}

// The metropolis needs to retur a chain

vector <VectorXd> AdaptiveMetropolis(const PPath& ppath, const vector<double>& vProxy,const VectorXd& x0, int n_iters, double dt,int adapt_start = 1000, int adapt_every = 50,double jitter = 1e-5){
    const int d = x0.size();
    mt19937 rng(123);
    double vbar = 0.0;
    for (int i = 0; i < ppath.v.size(); i++) vbar += vProxy[i];

    vbar = vbar/ppath.v.size();
    double s = 2.0/sqrt((double)d);
    MatrixXd C = MatrixXd::Identity(d,d);
    VectorXd x = x0;
    VectorXd mean = x;
    MatrixXd Sigma = MatrixXd::Identity(d,d);
    vector<VectorXd>chain;
    chain.reserve(n_iters);
    uniform_real_distribution<double> unif(0.0, 1.0);
    double logp = logPosterior(x,vbar,ppath,vProxy,dt);
    int accepted = 0;
    for(int it = 1; it < n_iters;it++){
        LLT<MatrixXd>llt(Sigma);
        if (llt.info() != Eigen::Success) {
            Sigma = Eigen::MatrixXd::Identity(d,d); // or reset to I
            llt.compute(Sigma);
        }
        MatrixXd L = llt.matrixL();
        VectorXd z = randn_vec(rng,d);
        VectorXd xnew = x + s*(L*z);
        double logpnew = logPosterior(xnew,vbar,ppath,vProxy,dt);
        double alpha_uniform = unif(rng);
        if (log(alpha_uniform) < (logpnew - logp)){
            x = xnew;
            logp = logpnew;
            chain.push_back(x);
            accepted ++;
        }
        if (it >= adapt_start){
            VectorXd x_it = x;
            VectorXd delta = x_it - mean;
            mean = mean + delta/(it - adapt_start + 1);
            VectorXd delta2 = x_it - mean; //we redeclare delta because of some precision issue
            C += delta*delta2.transpose();

            if((it % adapt_every) == 0){
                MatrixXd empCov = C/max(1,it - adapt_start);
                Sigma = empCov + jitter * MatrixXd::Identity(d,d);
            }

        }
    }
    double acRate = (double)accepted/(double)n_iters;
    cout<<"Acceptance Rate:"<<acRate<<endl;
    return chain;
}


vector <VectorXd> AdaptiveMetropolis(const PPath& ppath, const VectorXd& x0, int n_iters, double v0_actual,double dt,int num_particles,int adapt_start = 1000, int adapt_every = 50,double jitter = 1e-5){
    const int d = x0.size();
    mt19937 rng(123);
    double vbar = 0.0;
    for (int i = 0; i < ppath.v.size(); i++) vbar += ppath.v[i];
    vbar = vbar/ppath.v.size();
    double s = 2.0/sqrt((double)d);
    MatrixXd C = MatrixXd::Identity(d,d);
    VectorXd x = x0;
    VectorXd mean = x;
    MatrixXd Sigma = MatrixXd::Identity(d,d);
    vector<VectorXd>chain;
    chain.reserve(n_iters);
    uniform_real_distribution<double> unif(0.0, 1.0);
    double logp = logPosterior(x,vbar,ppath,v0_actual,dt,num_particles);
    int accepted = 0;
    for(int it = 1; it < n_iters;it++){
        if(it % 100 == 0) cout<<it<<" iters of MCMC is done."<<endl;
        LLT<MatrixXd>llt(Sigma);
        if (llt.info() != Eigen::Success) {
            Sigma = Eigen::MatrixXd::Identity(d,d); // or reset to I
            llt.compute(Sigma);
        }
        MatrixXd L = llt.matrixL();
        VectorXd z = randn_vec(rng,d);
        VectorXd xnew = x + s*(L*z);
        double logpnew = logPosterior(xnew,vbar,ppath,v0_actual,dt,num_particles);
        double alpha_uniform = unif(rng);
        if (log(alpha_uniform) < (logpnew - logp)){
            x = xnew;
            logp = logpnew;
            chain.push_back(x);
            accepted ++;
        }
        if (it >= adapt_start){
            VectorXd x_it = x;
            VectorXd delta = x_it - mean;
            mean = mean + delta/(it - adapt_start + 1);
            VectorXd delta2 = x_it - mean; //we redeclare delta because of some precision issue
            C += delta*delta2.transpose();

            if((it % adapt_every) == 0){
                MatrixXd empCov = C/max(1,it - adapt_start);
                Sigma = empCov + jitter * MatrixXd::Identity(d,d);
            }

        }
    }
    double acRate = (double)accepted/(double)n_iters;
    cout<<"Acceptance Rate:"<<acRate<<endl;
    return chain;
}



void chainStatistics(vector<VectorXd>& chain,int burn_in,HestonPParams& meanP,HestonPParams& varP){
    double mu_mean = 0.0;
    double kappaP_mean = 0.0;
    double thetaP_mean = 0.0;
    double xi_mean = 0.0;
    double rho_mean = 0.0;
    
    double mu_var = 0.0;
    double kappaP_var = 0.0;
    double thetaP_var = 0.0;
    double xi_var = 0.0;
    double rho_var = 0.0;
    
    for (int i = burn_in; i < chain.size();i++){
        HestonPParams tmp;
        tmp = xToParamsMCMC(chain[i]);
        mu_mean += tmp.mu;
        mu_var += (tmp.mu * tmp.mu);
        
        kappaP_mean += tmp.kappaP;
        kappaP_var += (tmp.kappaP * tmp.kappaP);

        thetaP_mean += tmp.thetaP;
        thetaP_var += (tmp.thetaP * tmp.thetaP);

        xi_mean += tmp.xi;
        xi_var += (tmp.xi * tmp.xi);

        rho_mean += tmp.rho;
        rho_var += (tmp.rho * tmp.rho);        
    }
    mu_mean = mu_mean/(chain.size() - burn_in);
    kappaP_mean = kappaP_mean/(chain.size() - burn_in);
    thetaP_mean = thetaP_mean/(chain.size() - burn_in);
    xi_mean = xi_mean/(chain.size() - burn_in);
    rho_mean = rho_mean/(chain.size() - burn_in);
    
    mu_var = mu_var/(chain.size() - burn_in) - mu_mean * mu_mean;
    kappaP_var = kappaP_var/(chain.size() - burn_in) - kappaP_mean * kappaP_mean;
    thetaP_var = thetaP_var/(chain.size() - burn_in) - thetaP_mean * thetaP_mean;
    xi_var = xi_var/(chain.size() - burn_in) - xi_mean * xi_mean;
    rho_var = rho_var/(chain.size() - burn_in) - rho_mean * rho_mean;

    meanP.mu = mu_mean;
    meanP.kappaP = kappaP_mean;
    meanP.thetaP = thetaP_mean;
    meanP.xi = xi_mean;
    meanP.rho = rho_mean;

    varP.mu = mu_var;
    varP.kappaP = kappaP_var;
    varP.thetaP = thetaP_var;
    varP.xi = xi_var;
    varP.rho = rho_var;
}


void mcmcOverLatent(HestonPParams& P, PPath& ppath, vector<double>& vProxy,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP,int n_iters = 5000){

    double logLikelihood = getPathLikelihood(P,ppath,vProxy,dt);
    cout<<"Likelihood for actual P:"<<logLikelihood<<endl;
    // mt19937 gen(123);
    // double eps = 0.5;
    // VectorXd x0(5);
    // uniform_real_distribution<> dist_mu(actualP.mu * (1-eps), actualP.mu * (1+eps));
    // uniform_real_distribution<> dist_kappaP(actualP.kappaP * (1-eps), actualP.kappaP * (1+eps));
    // uniform_real_distribution<> dist_thetaP(actualP.thetaP * (1-eps), actualP.thetaP * (1+eps));
    // uniform_real_distribution<> dist_xi(actualP.xi * (1-eps), actualP.xi * (1+eps));
    // uniform_real_distribution<> dist_rho(max(-0.95, actualP.rho - 0.4),min(-0.05, actualP.rho + 0.4));

    // x0[0] = dist_mu(gen);
    // x0[1] = dist_kappaP(gen);
    // x0[2] = dist_thetaP(gen);
    // x0[3] = dist_xi(gen);
    // x0[4] = dist_rho(gen);

    // int n_iters = 5000;
    vector<VectorXd>chain = AdaptiveMetropolis(ppath,vProxy,x0,n_iters,dt);
    cout<<"chain is built;"<<endl;

    chainStatistics(chain,1000,meanP,varP);
    // cout<<meanP<<endl;
    // cout<<varP<<endl;
    
}

void pmcmcOverLatent(HestonPParams& P, PPath& ppath,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP,int n_iters = 5000,int num_particles = 1500){

    //static thread_local std::mt19937 gen(69);
    // double logLikelihood = getPathLikelihood(actualP,ppath);
    // cout<<"Actual Likelihood with true P path:"<<logLikelihood<<endl;

    //mt19937 gen(123);
    double eps = 0.5;

    // VectorXd x0(5);
    // uniform_real_distribution<> dist_mu(actualP.mu * (1-eps), actualP.mu * (1+eps));
    // uniform_real_distribution<> dist_kappaP(actualP.kappaP * (1-eps), actualP.kappaP * (1+eps));
    // uniform_real_distribution<> dist_thetaP(actualP.thetaP * (1-eps), actualP.thetaP * (1+eps));
    // uniform_real_distribution<> dist_xi(actualP.xi * (1-eps), actualP.xi * (1+eps));
    // uniform_real_distribution<> dist_rho(max(-0.95, actualP.rho - 0.4),min(-0.05, actualP.rho + 0.4));

    // x0[0] = dist_mu(gen);
    // x0[1] = dist_kappaP(gen);
    // x0[2] = dist_thetaP(gen);
    // x0[3] = dist_xi(gen);
    // x0[4] = dist_rho(gen);

    double v0_actual = P.v0;
    // int n_iters = 5000;
    // int num_particles = 1500;
    vector<VectorXd>chain = AdaptiveMetropolis(ppath,x0,n_iters,v0_actual,dt,num_particles); // 5000 --> n_iters (no of MH steps)
    cout<<"chain is built;"<<endl;
    // HestonPParams meanP;
    // HestonPParams varP;
    
    chainStatistics(chain,1000,meanP,varP); // 1000 is burn in
    // cout<<meanP<<endl;
    // cout<<varP<<endl;
}
}//namespace qe;



