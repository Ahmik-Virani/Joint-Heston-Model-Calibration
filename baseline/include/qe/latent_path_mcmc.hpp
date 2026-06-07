#pragma once
#include <iomanip>
#include "qe/path.hpp"
#include <vector>
#include <Eigen/Core>
#include <Eigen/Cholesky>
#include "qe/particle_filters.hpp"

namespace qe{
    using Eigen::VectorXd;
    using std::vector;
    using std::mt19937;

    double getPathLikelihood(const HestonPParams& P, const PPath& ppath,const vector<double>& vProxy,double dt);
    HestonPParams xToParamsMCMC(VectorXd& x);
    VectorXd ParamstoxMCMC(HestonPParams& P);
    static inline double log_norm_prior(double x,double m, double s);
    double logPrior(VectorXd x,double vbar);
    double logPosterior(VectorXd x,double vbar, const PPath& ppath,const vector<double>& vProxy,double dt);
    double logPosterior(VectorXd x,double vbar, const PPath& ppath,double v0_actual,double dt,int N);
    VectorXd randn_vec(mt19937& rng, int d);
    vector <VectorXd> AdaptiveMetropolis(const PPath& ppath, const vector<double>& vProxy,const VectorXd& x0, int n_iters, double dt,int adapt_start, int adapt_every,double jitter);
    vector <VectorXd> AdaptiveMetropolis(const PPath& ppath, const VectorXd& x0, int n_iters, double v0_actual,double dt,int num_particles,int adapt_start, int adapt_every,double jitter);
    void chainStatistics(vector<VectorXd>& chain,int burn_in,HestonPParams& meanP,HestonPParams& varP);
    void mcmcOverLatent(HestonPParams& P, PPath& ppath,vector<double>& vProxy,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP,int n_iters);
    void pmcmcOverLatent(HestonPParams& P, PPath& ppath,VectorXd x0,double dt,HestonPParams& meanP,HestonPParams& varP,int n_iters,int num_particles);

}
