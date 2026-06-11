#pragma once

#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/constraint.hpp>
#include <ql/math/optimization/simplex.hpp>
#include <ql/math/optimization/endcriteria.hpp>

#include "callPrice.cpp"

using namespace std;

class surfaceCalibration {
private:
    // This is the vector we will be optimizing
    vector<double> result;
    vector<double> unconstrainedParams;

    vector<double> P_to_Q(const vector<double> &particle) const {
        vector<double> Q_space_params = particle;
        Q_space_params[2] = particle[2] + particle[5];
        Q_space_params[3] = (particle[2] * particle[3]) / (particle[2] + particle[5] + 1e-8);
        return Q_space_params;
    }

    vector<double> make_unconstrained(const vector<double> &params){
        vector<double> unconstrainedParams(params.size());
        
        unconstrainedParams[0] = params[0];
        unconstrainedParams[1] = log(params[1]);
        unconstrainedParams[2] = log(params[2]);
        unconstrainedParams[3] = log(params[3]);
        unconstrainedParams[4] = atanh(params[4]);
        unconstrainedParams[5] = params[5];
        unconstrainedParams[6] = log(params[6]);

        return unconstrainedParams;
    }

    vector<double> make_original(const vector<double> &unconstrained_params) const {
        vector<double> params(unconstrained_params.size());

        params[0] = unconstrained_params[0];
        params[1] = exp(unconstrained_params[1]);
        params[2] = exp(unconstrained_params[2]);
        params[3] = exp(unconstrained_params[3]);
        params[4] = tanh(unconstrained_params[4]);
        params[5] = unconstrained_params[5];
        params[6] = exp(unconstrained_params[6]);

        return params;
    }

    bool isValidQ(const vector<double> &qParams) const{
        // ensure that the new parameters in Q space are in limits
        if(qParams[2] <= 0) return false;
        if(qParams[3] <= 0) return false;

        return true;
    }

    // double loss_function(const vector<double> &unconstrained_params, const vector<double> &K, const vector<double> &T, const vector<vector<double>> &grid, const QuantLib::Date &today, const double &spot, const double &r, const double &q){
    double loss_function(const vector<double> &unconstrained_params) const {
        vector<double> qParams = P_to_Q(make_original(unconstrained_params));

        if(!isValidQ(qParams)){
            return 1e8;
        }

        // If the points are valid, let us compute sum of squared errors
        double SSE = 0.0;

        for(int i = 0 ; i < grid_.size() ; i++){
            double true_value = grid_[i][2];
            double predicted_value = qe::hestonCallPrice(
                spot_,
                grid_[i][0],
                grid_[i][1],
                r_,
                q_,
                qParams[6],
                qParams[2],
                qParams[3],
                qParams[1],
                qParams[4],
                today_
            );

            SSE += ((true_value - predicted_value) * (true_value - predicted_value));
        }
        return SSE;
    }

    vector<double> x0_;
    vector<vector<double>> grid_;
    QuantLib::Date today_;
    double spot_, r_, q_;

    void runSimplexCalibration() {
        struct Cost : public QuantLib::CostFunction {
            const surfaceCalibration& self;
            explicit Cost(const surfaceCalibration& s) : self(s) {}

            QuantLib::Real value(const QuantLib::Array& x) const override {
                vector<double> u(x.size());
                for (size_t i = 0; i < u.size(); ++i) u[i] = x[i];
                return self.loss_function(u);
            }

            QuantLib::Array values(const QuantLib::Array& x) const override {
                QuantLib::Array y(1);
                y[0] = value(x);
                return y;
            }
        };

        QuantLib::Array x0(x0_.size());
        for (size_t i = 0; i < x0_.size(); ++i) x0[i] = x0_[i];

        Cost cost(*this);
        QuantLib::NoConstraint constraint;
        QuantLib::Problem problem(cost, constraint, x0);

        QuantLib::Simplex solver(0.15);
        QuantLib::EndCriteria ec(100, 20, 1e-8, 1e-8, 1e-8);
        solver.minimize(problem, ec);

        QuantLib::Array xStar = problem.currentValue();
        vector<double> uStar(xStar.size());
        for (size_t i = 0; i < uStar.size(); ++i) uStar[i] = xStar[i];

        result = make_original(uStar); // return P-space params
    }

public:
    // we need to pass the C(K, T)
    // We also need to pass in the guesses

    // I am assuming initial guesses are good - not violating constraints
    surfaceCalibration(vector<vector<double>> grid, vector<double> guesses, QuantLib::Date today, double spot, double r, double q){
        // Note that we pass the guesses in P-space
        // This means that we have (index - name)
        // 0 - mu                   
        // 1 - sigma (vol-of-vol)       -> exp for > 0
        // 2 - kappa                    -> exp for > 0
        // 3 - theta                    -> exp for > 0
        // 4 - rho                      -> tanh for [-1,1]
        // 5 - lambda
        // 6 - Instantaneous volatility -> exp for > 0

        // First let us ensure that the guesses are valid
        guesses[1] = max(guesses[1], 1e-8);
        guesses[2] = max(guesses[2], 1e-8);
        guesses[3] = max(guesses[3], 1e-8);
        guesses[4] = clamp(guesses[4], -0.999, 0.999);
        guesses[6] = max(guesses[6], 1e-8);

        unconstrainedParams = make_unconstrained(guesses);

        // These lines are added for the simplex
        grid_ = grid; today_ = today; spot_ = spot; r_ = r; q_ = q;
        x0_ = make_unconstrained(guesses);
        // Let us do calibration here itself, and later return it 
        runSimplexCalibration();
    }

    vector<double> getCalibration(){
        return result;
    }
};