// This requires QuantLib, Eigen, Armadillo (optional), LBFGSPP
// g++ -std=c++17 app/main.cpp src/*.cpp -Iinclude -I$HOME/local/include -I$HOME/path/to//eigen-3.4.0 -I$HOME/path/to/LBFGSpp/include -L$HOME/local/lib -lQuantLib -o demo
 // Run this before running the ./demo (any executable). This is if you are installing QuantLib and Armadillo locally.
export LD_LIBRARY_PATH="$HOME/local/lib:$LD_LIBRARY_PATH"

This is mainly a Simulation and Calibration engine for stochastic volatility models (currently Heston) supporting bayesian calibration of model parameters.
It can also be used for Particle Filters and PMCMC techniques.