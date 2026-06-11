#pragma once

#include <ql/quantlib.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <vector>

#include <ql/math/optimization/costfunction.hpp>
#include <ql/math/optimization/endcriteria.hpp>
#include <ql/math/optimization/problem.hpp>
#include <ql/math/optimization/simplex.hpp>

#include "callPrice.cpp"

using namespace std;

class surfaceCalibrationMultiStartSimplex {
private:
    vector<double> result;
    vector<vector<double>> grid_;
    QuantLib::Date today_;
    double spot_ = 0.0, r_ = 0.0, q_ = 0.0;
    vector<double> x0_;
    mt19937 rng_{random_device{}()};

    vector<double> P_to_Q(const vector<double>& particle) const {
        vector<double> qParams = particle;
        qParams[2] = particle[2] + particle[5];
        qParams[3] = (particle[2] * particle[3]) / (particle[2] + particle[5] + 1e-8);
        return qParams;
    }

    vector<double> make_unconstrained(const vector<double>& params) const {
        vector<double> u(params.size());
        u[0] = params[0];
        u[1] = log(max(params[1], 1e-8));
        u[2] = log(max(params[2], 1e-8));
        u[3] = log(max(params[3], 1e-8));
        u[4] = atanh(clamp(params[4], -0.999, 0.999));
        u[5] = params[5];
        u[6] = log(max(params[6], 1e-8));
        return u;
    }

    vector<double> make_original(const vector<double>& u) const {
        vector<double> p(u.size());
        p[0] = u[0];
        p[1] = exp(u[1]);
        p[2] = exp(u[2]);
        p[3] = exp(u[3]);
        p[4] = tanh(u[4]);
        p[5] = u[5];
        p[6] = exp(u[6]);
        return p;
    }

    bool isValidQ(const vector<double>& qParams) const {
        return qParams[2] > 0.0 && qParams[3] > 0.0;
    }

    double loss_function(const vector<double>& unconstrained) const {
        const vector<double> qParams = P_to_Q(make_original(unconstrained));
        if (!isValidQ(qParams))
            return 1e8;

        double sse = 0.0;
        for (const auto& row : grid_) {
            const double true_value = row[2];
            const double predicted_value = qe::hestonCallPrice(
                spot_, row[0], row[1], r_, q_, qParams[6], qParams[2], qParams[3], qParams[1], qParams[4], today_);
            const double err = true_value - predicted_value;
            sse += err * err;
        }
        return sse;
    }

    vector<double> clampPParams(vector<double> p) const {
        p[1] = max(p[1], 1e-8);
        p[2] = max(p[2], 1e-8);
        p[3] = max(p[3], 1e-8);
        p[4] = clamp(p[4], -0.999, 0.999);
        p[6] = max(p[6], 1e-8);
        return p;
    }

    vector<double> jitterStart(const vector<double>& base, size_t restart) {
        normal_distribution<double> dist(0.0, 0.2 + 0.05 * static_cast<double>(restart));
        vector<double> start = base;
        for (double& v : start)
            v += dist(rng_);
        return make_unconstrained(clampPParams(make_original(start)));
    }

    vector<double> runSingleSimplex(const vector<double>& initial, double& bestLoss) {
        struct Cost : public QuantLib::CostFunction {
            const surfaceCalibrationMultiStartSimplex& self;
            explicit Cost(const surfaceCalibrationMultiStartSimplex& s) : self(s) {}

            QuantLib::Real value(const QuantLib::Array& x) const override {
                vector<double> u(x.size());
                for (size_t i = 0; i < u.size(); ++i)
                    u[i] = x[i];
                return self.loss_function(u);
            }

            QuantLib::Array values(const QuantLib::Array& x) const override {
                QuantLib::Array y(1);
                y[0] = value(x);
                return y;
            }
        };

        QuantLib::Array x(initial.size());
        for (size_t i = 0; i < initial.size(); ++i)
            x[i] = initial[i];

        Cost cost(*this);
        QuantLib::NoConstraint constraint;
        QuantLib::Problem problem(cost, constraint, x);

        QuantLib::Simplex solver(0.15);
        QuantLib::EndCriteria ec(10000, 20, 1e-8, 1e-8, 1e-8);
        solver.minimize(problem, ec);

        QuantLib::Array xStar = problem.currentValue();
        vector<double> uStar(xStar.size());
        for (size_t i = 0; i < uStar.size(); ++i)
            uStar[i] = xStar[i];

        bestLoss = loss_function(uStar);
        return make_original(uStar);
    }

    void runMultiStartSimplexCalibration() {
        const vector<double> base = make_unconstrained(x0_);
        const size_t restarts = 8;

        double bestLoss = numeric_limits<double>::infinity();
        vector<double> bestResult = x0_;

        for (size_t restart = 0; restart < restarts; ++restart) {
            const vector<double> start = (restart == 0) ? base : jitterStart(base, restart);
            double thisLoss = numeric_limits<double>::infinity();
            vector<double> candidate = runSingleSimplex(start, thisLoss);
            if (thisLoss < bestLoss) {
                bestLoss = thisLoss;
                bestResult = candidate;
            }
        }

        result = bestResult;
    }

public:
    surfaceCalibrationMultiStartSimplex(vector<vector<double>> grid, vector<double> guesses, QuantLib::Date today, double spot, double r, double q) {
        x0_ = clampPParams(std::move(guesses));
        grid_ = std::move(grid);
        today_ = today;
        spot_ = spot;
        r_ = r;
        q_ = q;
        runMultiStartSimplexCalibration();
    }

    vector<double> getCalibration() {
        return result;
    }
};
