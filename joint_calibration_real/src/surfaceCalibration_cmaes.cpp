#pragma once

#include <ql/quantlib.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "callPrice.cpp"

using namespace std;

namespace {

using Matrix = vector<vector<double>>;

Matrix identityMatrix(size_t n) {
    Matrix m(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i)
        m[i][i] = 1.0;
    return m;
}

Matrix outerProduct(const vector<double>& a, const vector<double>& b) {
    const size_t n = a.size();
    Matrix m(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j)
            m[i][j] = a[i] * b[j];
    }
    return m;
}

Matrix addScaled(const Matrix& a, const Matrix& b, double scale) {
    Matrix m = a;
    for (size_t i = 0; i < m.size(); ++i) {
        for (size_t j = 0; j < m[i].size(); ++j)
            m[i][j] += scale * b[i][j];
    }
    return m;
}

Matrix scalarMultiply(const Matrix& a, double scale) {
    Matrix m = a;
    for (auto& row : m) {
        for (double& x : row)
            x *= scale;
    }
    return m;
}

vector<double> matVec(const Matrix& a, const vector<double>& x) {
    vector<double> y(a.size(), 0.0);
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < x.size(); ++j)
            y[i] += a[i][j] * x[j];
    }
    return y;
}

Matrix choleskyDecompose(const Matrix& a) {
    const size_t n = a.size();
    Matrix l(n, vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double sum = a[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= l[i][k] * l[j][k];

            if (i == j) {
                l[i][j] = sqrt(max(sum, 1e-12));
            } else {
                l[i][j] = sum / max(l[j][j], 1e-12);
            }
        }
    }
    return l;
}

double squaredNorm(const vector<double>& x) {
    double s = 0.0;
    for (double v : x)
        s += v * v;
    return s;
}

vector<double> clampPParams(vector<double> p) {
    p[1] = max(p[1], 1e-8);
    p[2] = max(p[2], 1e-8);
    p[3] = max(p[3], 1e-8);
    p[4] = clamp(p[4], -0.999, 0.999);
    p[6] = max(p[6], 1e-8);
    return p;
}

} // namespace

class surfaceCalibrationCMAES {
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

    void runCMAESCalibration() {
        const vector<double> initial = make_unconstrained(x0_);
        const size_t n = initial.size();
        const size_t lambda = max<size_t>(8, 4 + static_cast<size_t>(3.0 * log(static_cast<double>(n) + 1.0)));
        const size_t mu = max<size_t>(2, lambda / 2);

        vector<double> weights(mu);
        for (size_t i = 0; i < mu; ++i)
            weights[i] = log(static_cast<double>(mu) + 0.5) - log(static_cast<double>(i) + 1.0);
        const double weightSum = accumulate(weights.begin(), weights.end(), 0.0);
        for (double& w : weights)
            w /= weightSum;

        double sumSquares = 0.0;
        for (double w : weights)
            sumSquares += w * w;
        const double mueff = 1.0 / sumSquares;

        vector<double> mean = initial;
        Matrix C = identityMatrix(n);
        double sigma = 0.5;

        vector<double> bestX = mean;
        double bestLoss = loss_function(bestX);

        const size_t maxIter = 120;
        const double chiN = sqrt(static_cast<double>(n)) * (1.0 - 1.0 / (4.0 * static_cast<double>(n)) + 1.0 / (21.0 * static_cast<double>(n) * static_cast<double>(n)));
        const double cmu = min(0.5, 0.4 + 1.0 / (static_cast<double>(n) + 1.0));

        for (size_t iter = 0; iter < maxIter; ++iter) {
            Matrix L = choleskyDecompose(C);

            struct Candidate {
                double loss;
                vector<double> x;
                vector<double> y;
            };

            vector<Candidate> population;
            population.reserve(lambda);

            normal_distribution<double> norm(0.0, 1.0);
            for (size_t k = 0; k < lambda; ++k) {
                vector<double> z(n, 0.0);
                for (double& zi : z)
                    zi = norm(rng_);

                vector<double> step = matVec(L, z);
                vector<double> x(n, 0.0);
                vector<double> y(n, 0.0);
                for (size_t j = 0; j < n; ++j) {
                    y[j] = step[j];
                    x[j] = mean[j] + sigma * step[j];
                }

                const double loss = loss_function(x);
                population.push_back({loss, x, y});

                if (loss < bestLoss) {
                    bestLoss = loss;
                    bestX = x;
                }
            }

            sort(population.begin(), population.end(), [](const Candidate& a, const Candidate& b) {
                return a.loss < b.loss;
            });

            vector<double> newMean(n, 0.0);
            vector<double> yw(n, 0.0);
            for (size_t i = 0; i < mu; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    newMean[j] += weights[i] * population[i].x[j];
                    yw[j] += weights[i] * population[i].y[j];
                }
            }

            Matrix newC(n, vector<double>(n, 0.0));
            for (size_t i = 0; i < mu; ++i)
                newC = addScaled(newC, outerProduct(population[i].y, population[i].y), weights[i]);

            C = addScaled(scalarMultiply(C, 1.0 - cmu), newC, cmu);
            for (size_t i = 0; i < n; ++i)
                C[i][i] = max(C[i][i], 1e-12);

            const double normYw = sqrt(squaredNorm(yw));
            sigma *= exp(0.2 * (normYw / max(chiN, 1e-12) - 1.0));
            sigma = clamp(sigma, 1e-4, 2.0);
            mean = newMean;

            if (sigma < 1e-4 || bestLoss < 1e-12)
                break;
        }

        result = make_original(bestX);
    }

public:
    surfaceCalibrationCMAES(vector<vector<double>> grid, vector<double> guesses, QuantLib::Date today, double spot, double r, double q) {
        x0_ = clampPParams(std::move(guesses));
        grid_ = std::move(grid);
        today_ = today;
        spot_ = spot;
        r_ = r;
        q_ = q;
        runCMAESCalibration();
    }

    vector<double> getCalibration() {
        return result;
    }
};
