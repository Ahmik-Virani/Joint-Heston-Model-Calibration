#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include "qe/params.hpp"
#include "qe/path.hpp"



namespace qe{

    using std::mt19937;
    using std::vector;
    using std::string;

    struct filterValues{
        vector<double>v_chain;
        vector<double>W;
        vector<vector<double>> particles;
        vector<vector<int>> ancestors;
        double likelihood;
        int resample_count;
    };

    static inline double logNormal(double x, double mu, double sigma);
    double observationLikelihood(double r,double v,double mu,double dt);
    double transitionFunction(const HestonPParams& P, double v,double r, double dt);
    vector<int> systematicResample(const vector<double>& W, std::mt19937& gen);
    filterValues ParticleFilter(const HestonPParams& P, double v0_actual,const PPath& ppath,double dt,int N);
    vector<double> ancestralSampling(filterValues& filter, mt19937& gen);
    void writeSampledPaths(vector<vector<double>>& sampledPaths,string filename = "sampled_paths.csv");
}