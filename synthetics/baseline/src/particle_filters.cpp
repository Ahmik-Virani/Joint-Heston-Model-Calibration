#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "qe/params.hpp"
#include "qe/particle_filters.hpp"

using namespace std;


namespace qe{

static inline double logNormal(double x, double mu, double sigma){
    double z = (mu - x)/sigma;
    return -0.5*log(2*M_PI) - log(sigma) - 0.5 * z * z;
}


double observationLikelihood(double r,double v,double mu,double dt){
    double vp = max(v,1e-6);
    double mean = (mu - 0.5 * vp) * dt;
    double sd = sqrt(vp * dt);
    double obsL = logNormal(r,mean,sd);
    return obsL;
}

double qLikelihood(const HestonPParams& P,double v,double r,double dt,double new_v){
    
    double kappaP = P.kappaP;
    double thetaP = P.thetaP;
    double xi = P.xi;
    double rho = P.rho;
    double mu = P.mu;
    double vp = max(v,1e-6);

    double mean = vp + kappaP * (thetaP - vp)*dt + xi*rho*(r - (mu - 0.5*v)*dt);
    double sd = xi * sqrt(vp*(1 - rho * rho)*dt);
    double qL = logNormal(new_v,mean,sd);
    return qL;

}

double priorLikelihood(const HestonPParams& P,double v,double dt,double new_v){
    double kappaP = P.kappaP;
    double thetaP = P.thetaP;
    double xi = P.xi;
    double rho = P.rho;
    double mu = P.mu;
    double vp = max(v,1e-6);

    double mean = vp + kappaP * (thetaP - vp)*dt;
    double sd = xi * sqrt(vp * dt);
    double pL = logNormal(new_v,mean,sd);
    return pL;

}



double transitionFunction(const HestonPParams& P, double v,double r, double dt){
    
    static thread_local std::mt19937 gen(std::random_device{}());
    static thread_local std::normal_distribution<double> normal(0.0, 1.0);

    double kappaP = P.kappaP;
    double thetaP = P.thetaP;
    double xi = P.xi;
    double rho = P.rho;
    double mu = P.mu;
    double vp = max(v,1e-6);
    double mean = vp + kappaP*(thetaP - vp)*dt + xi * rho * (r - (mu - 0.5*vp) * dt);
    double sd = xi * sqrt(vp * dt * (1 - rho*rho));
    double next_v = mean + sd * normal(gen);
    return next_v;  
}

vector<int> systematicResample(const vector<double>& W, std::mt19937& gen) {
    
    int N = W.size();
    vector<int> idx(N);

    vector<double> C(N);
    C[0] = W[0];
    for (int i = 1; i < N; i++)
        C[i] = C[i-1] + W[i];

    uniform_real_distribution<double> U(0.0, 1.0 / N);
    double u0 = U(gen);

    int i = 0;

    for (int j = 0; j < N; j++) {

        double u = u0 + (double)j / N;

        while (u > C[i])
            i++;

        idx[j] = i;
    }

    return idx;
}


filterValues ParticleFilter(const HestonPParams& P, double v0_actual,const PPath& ppath,double dt,int N){
    
    static thread_local std::mt19937 gen(std::random_device{}());
    double eps = 0.5;
    // double v0_actual = P.v0;
    int resample_count = 0;
    uniform_real_distribution<double> uniform(v0_actual*(1-eps), v0_actual*(1+eps));
    
    vector<double>logReturns = ppath.returns;
    int n = logReturns.size();
    //int N = 20; // no of particles
    
    if (logReturns[n-1] == 0) n = n-1;

    vector<double>logW (N,0.0);
    vector<double>v_chain (N,0.0);
    vector<double>W(N,0.0);
    vector<double>baseW(N,1.0);
    vector<double>new_v_chain(N,0.0);
    vector<vector<double>> particles(n + 1, vector<double>(N, 0.0));
    vector<vector<int>> ancestors(n + 1, vector<int>(N, -1));
    
    

    for (int j = 0; j < N; j++){
        v_chain[j] = uniform(gen); 
        ancestors[0][j] = -1;
        particles[0][j] = v_chain[j];
    }

    //double v = uniform(gen);
    double mu = P.mu;
    double logLik = 0.0;

    for (int i = 0; i < n; i++){
        double r = logReturns[i]; //logReturns[i] --> r_t
        for (int j = 0; j < N; j++){
            double v = v_chain[j]; // v_t for all chain

            double new_v = transitionFunction(P,v,r,dt); // get new transitions
            new_v_chain[j] = new_v; // store new transitions
            //logW[j] = (observationLikelihood(r,v,mu,dt) + priorLikelihood(P,v,dt,new_v) - qLikelihood(P,v,r,dt,new_v)); //calculate the weights
            logW[j] = observationLikelihood(r,v,mu,dt);
        }
        double maxLogW = *max_element(logW.begin(), logW.end());
        double sumW = 0.0;
        for(int j = 0; j < N; j++){
            logW[j] = logW[j] + log(baseW[j]);
            W[j] = exp(logW[j] - maxLogW);
            sumW += W[j];
        }
        logLik += (maxLogW + log(sumW) - log((double)N));

        double ess = 0.0;
        for(int j = 0; j < N; j++){
            W[j] = W[j]/sumW;
            ess += W[j] * W[j];
        }
        ess = 1/ess;
        if (ess < 0.3 * N){
            resample_count ++;
            //cout<<"Resampled at step:"<<i<<endl;
            vector<int> idx = systematicResample(W,gen);
            for(int j = 0; j < N;j++){
                v_chain[j] = new_v_chain[idx[j]];
                particles[i + 1][j] = v_chain[j];
                ancestors[i + 1][j] = idx[j];
                W[j] = 1.0/N;
            }    
        }
        else{
            for (int j = 0; j < N;j++){
                v_chain[j] = new_v_chain[j];
                particles[i + 1][j] = v_chain[j];
                ancestors[i + 1][j] = j;
            }
            
        }
        baseW = W;
    }
    //cout<<"Resample Count:"<<resample_count<<endl;
    filterValues filter;
    filter.resample_count = resample_count;
    filter.likelihood = logLik;
    filter.ancestors = ancestors;
    filter.particles = particles;
    filter.v_chain = v_chain;
    filter.W = W;
    return filter;
    
}


vector<double>  ancestralSampling(filterValues& filter, mt19937& gen){
    
    int N = filter.particles[0].size();
    int T = filter.particles.size();
    
    vector<double>sampledPath(T,0.0);

    discrete_distribution<>dist(filter.W.begin(),filter.W.end());
    int k = dist(gen);
    //cout<<"k:"<<k<<endl;
    for(int t = T-1; t >= 0; t--){
        sampledPath[t] = filter.particles[t][k];
        k = filter.ancestors[t][k];
        if ((k < 0) && (t > 0)) break;
    }
    return sampledPath;
}

void writeSampledPaths(vector<vector<double>>& sampledPaths,string filename){
    ofstream file(filename);
    int numSampledPaths = sampledPaths.size();
    int T = sampledPaths[0].size();
    for(int i = 0; i < numSampledPaths; i++){
        for(int t = 0; t < T; t++){
            file << sampledPaths[i][t];

            if(t < T-1)
                file << ",";
        }
        file << "\n";
    }

    file.close();
}
}//namespace qe
