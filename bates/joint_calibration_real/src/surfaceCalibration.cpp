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
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <Eigen/Dense>

#include "callPrice.cpp"

using namespace std;

class surfaceCalibration {
private:
    // This is the vector we will be optimizing
    vector<double> result;
    vector<double> unconstrainedParams;

    
    // [TODO] - check if these conversions are correct
    vector<double> P_to_Q(const vector<double> &particle){
        vector<double> Q_space_params = particle;
        Q_space_params[2] = particle[2] + particle[5];
        Q_space_params[3] = (particle[2] * particle[3]) / (particle[2] + particle[5] + 1e-8);
        Q_space_params[7] = particle[7] * exp(particle[10] * particle[8] + 0.5 * particle[10] * particle[10] * particle[9] * particle[9]);
        Q_space_params[8] = particle[8] + particle[10] * particle[9] * particle[9];
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
        unconstrainedParams[7] = log(params[7]);
        unconstrainedParams[8] = params[8];
        unconstrainedParams[9] = log(params[9]);
        unconstrainedParams[10] = params[10];

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
        params[7] = exp(unconstrained_params[7]);
        params[8] = unconstrained_params[8];
        params[9] = exp(unconstrained_params[9]);
        params[10] = unconstrained_params[10];

        return params;
    }

    // [TODO] - check if adding upper bounds is better
    bool isValidQ(const vector<double> &qParams) const{
        // ensure that the new parameters in Q space are in limits
        if (qParams[1] <= 1e-8) return false;
        if (qParams[2] <= 1e-8) return false;
        if (qParams[3] <= 1e-8) return false;
        if (qParams[4] <= -0.999 || qParams[4] >= 0.999) return false;
        if (qParams[6] <= 1e-10) return false;
        if (qParams[7] < 0.0) return false;
        if (qParams[9] <= 1e-8) return false;

        return true;
    }

    // double loss_function(const vector<double> &unconstrained_params, const vector<double> &K, const vector<double> &T, const vector<vector<double>> &grid, const QuantLib::Date &today, const double &spot, const double &r, const double &q){
    double loss_function(const vector<double> &unconstrained_params) const {
        vector<double> qParams = P_to_Q(make_original(unconstrained_params));

        if(!isValidQ(qParams)){
            cout << "---- I AM INVALID -----\n";
            return 1e8;
        }

        // If the points are valid, let us compute sum of squared errors
        double SSE = 0.0;

        for(int i = 0 ; i < grid_.size() ; i++){
            double true_value = grid_[i][2];
            double predicted_value = qe::batesCallPrice(
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
                qParams[7],
                qParams[8],
                qParams[9],
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
        // 7 - Jump Intensity           -> exp for > 0
        // 8 - Jump Mean            
        // 9 - Jump Volatility          -> exp for > 0
        // 10 - eta                     

        // First let us ensure that the guesses are valid
        guesses[1] = max(guesses[1], 1e-8);
        guesses[2] = max(guesses[2], 1e-8);
        guesses[3] = max(guesses[3], 1e-8);
        guesses[4] = clamp(guesses[4], -0.999, 0.999);
        guesses[6] = max(guesses[6], 1e-8);
        guesses[7] = max(guesses[7], 1e-8);
        guesses[9] = max(guesses[9], 1e-8);

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


class surfaceCalibrationLaplacian {
    private:
        // This is the vector we will be optimizing
        vector<double> unconstrainedParams;
        vector<double> prior_mean_;
        vector<double> prior_std_;
        double sigma_C_;

        // [TODO] - check if these conversions are correct
        vector<double> P_to_Q(const vector<double> &particle) const {
            vector<double> Q_space_params = particle;
            Q_space_params[2] = particle[2] + particle[5];
            Q_space_params[3] = (particle[2] * particle[3]) / (particle[2] + particle[5] + 1e-8);
            Q_space_params[7] = particle[7] * exp(particle[10] * particle[8] + 0.5 * particle[10] * particle[10] * particle[9] * particle[9]);
            Q_space_params[8] = particle[8] + particle[10] * particle[9] * particle[9];
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
            unconstrainedParams[7] = log(params[7]);
            unconstrainedParams[8] = params[8];
            unconstrainedParams[9] = log(params[9]);
            unconstrainedParams[10] = params[10];

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
            params[7] = exp(unconstrained_params[7]);
            params[8] = unconstrained_params[8];
            params[9] = exp(unconstrained_params[9]);
            params[10] = unconstrained_params[10];

            return params;
        }

        // [TODO] - Arka
        Eigen::MatrixXd make_original_covariance(Eigen::MatrixXd covariance_u,const vector<double> &unconstrained_params) const{
            Eigen::MatrixXd G = Eigen::MatrixXd::Zero(11, 11);
            G(0,0) = 1.0;
            G(1,1) = std::exp(unconstrained_params[1]);
            G(2,2) = std::exp(unconstrained_params[2]);
            G(3,3) = std::exp(unconstrained_params[3]);

            double rho = std::tanh(unconstrained_params[4]);
            G(4,4) = 1.0 - rho * rho;

            G(5,5) = 1.0;
            G(6,6) = std::exp(unconstrained_params[6]);
            G(7,7)   = std::exp(unconstrained_params[7]);  // jump intensity
            G(8,8)   = 1.0;                                // jump mean
            G(9,9)   = std::exp(unconstrained_params[9]);  // jump volatility
            G(10,10) = 1.0;                                // eta

            Eigen::MatrixXd covariance_p =
                G * covariance_u * G.transpose();
            return covariance_p;
        }
    
        // [TODO] - check if adding upper bounds is better
        bool isValidQ(const vector<double> &qParams) const{
            // ensure that the new parameters in Q space are in limits
            if (qParams[1] <= 1e-8) return false;
            if (qParams[2] <= 1e-8) return false;
            if (qParams[3] <= 1e-8) return false;
            if (qParams[4] <= -0.999 || qParams[4] >= 0.999) return false;
            if (qParams[6] <= 1e-10) return false;
            if (qParams[7] < 0.0) return false;
            if (qParams[9] <= 1e-8) return false;

            return true;
        }

        double prior_penalty(const vector<double>& pParams)const{
            double penalty = 0.0;
            for(int i = 0; i < pParams.size(); i++){
                double z = (pParams[i] - prior_mean_[i])/prior_std_[i];
                penalty += z * z;
            }
            return penalty * 0.5;
        }
    
        // double loss_function(const vector<double> &unconstrained_params, const vector<double> &K, const vector<double> &T, const vector<vector<double>> &grid, const QuantLib::Date &today, const double &spot, const double &r, const double &q){
        double loss_function(const vector<double> &unconstrained_params) const {
            vector<double>pParams = make_original(unconstrained_params);
            vector<double> qParams = P_to_Q(pParams);

    
            if(!isValidQ(qParams)){
                cout << "I am INVIALID qParams line 316\n";
                return 1e8;
            }
    
            // If the points are valid, let us compute sum of squared errors
            double surface_nll = 0.0;
    
            for(int i = 0 ; i < grid_.size() ; i++){
                double true_value = grid_[i][2];
                double predicted_value = qe::batesCallPrice(
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
                    qParams[7],
                    qParams[8],
                    qParams[9],
                    today_
                );
    
                surface_nll += ((true_value - predicted_value) * (true_value - predicted_value));
            }
            double prior_nll = prior_penalty(pParams);
            double SSE = prior_nll + surface_nll;
            return SSE;
        }

        vector<double> loss_function_lms(const vector<double> &unconstrained_params)const {
            vector<double> pParams = make_original(unconstrained_params);
            vector<double> qParams = P_to_Q(pParams);
        
            vector<double> residuals;
        
            if(!isValidQ(qParams)){
                residuals.assign(grid_.size() + pParams.size(), 1e4);
                return residuals;
            }
        
            for(int i = 0; i < grid_.size(); i++){
                double true_value = grid_[i][2];
        
                double predicted_value = qe::batesCallPrice(
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
                    qParams[7],
                    qParams[8],
                    qParams[9],
                    today_
                );
                residuals.push_back((true_value - predicted_value) / sigma_C_);
            }
        
            for(int i = 0; i < pParams.size(); i++){
                double z = (pParams[i] - prior_mean_[i]) / (prior_std_[i] + 1e-12);
                residuals.push_back(z);
            }
        
            return residuals;
        }
        
        vector<double> x0_;
        vector<vector<double>> grid_;
        QuantLib::Date today_;
        double spot_, r_, q_;
    
        void runSimplexCalibration() {
            struct Cost : public QuantLib::CostFunction {
                const surfaceCalibrationLaplacian& self;
                explicit Cost(const surfaceCalibrationLaplacian& s) : self(s) {}
    
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
            // [TODO] increase for bates
            QuantLib::EndCriteria ec(2000, 200, 1e-6, 1e-6, 1e-6);
            solver.minimize(problem, ec);
    
            QuantLib::Array xStar = problem.currentValue();
            vector<double> uStar(xStar.size());
            for (size_t i = 0; i < uStar.size(); ++i) uStar[i] = xStar[i];
    
            map_result.result = make_original(uStar); // return P-space params
        }

        void runLMSCalibration() {
            struct LMSFunctor {
                typedef double Scalar;
                enum {
                    InputsAtCompileTime = Eigen::Dynamic,
                    ValuesAtCompileTime = Eigen::Dynamic
                };
                typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
                typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
                typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

                const surfaceCalibrationLaplacian& self;
        
                LMSFunctor(const surfaceCalibrationLaplacian& s)
                    : self(s) {}
        
                int inputs() const {
                    return 11;   // number of parameters
                }
        
                int values() const {
                    return static_cast<int>(self.grid_.size()) + 11;
                }
        
                int operator()(const Eigen::VectorXd& x,
                               Eigen::VectorXd& fvec) const {
                    vector<double> u(x.size());
        
                    for(int i = 0; i < x.size(); i++) {
                        u[i] = x[i];
                    }
        
                    vector<double> r = self.loss_function_lms(u);
        
                    for(int i = 0; i < r.size(); i++) {
                        fvec[i] = r[i];
                    }
        
                    return 0;
                }
            };
        
            Eigen::VectorXd x(x0_.size());
        
            for(int i = 0; i < x0_.size(); i++) {
                x[i] = x0_[i];
            }
        
            LMSFunctor functor(*this);
            Eigen::NumericalDiff<LMSFunctor> numDiff(functor);
        
            Eigen::LevenbergMarquardt<
                Eigen::NumericalDiff<LMSFunctor>,
                double
            > lm(numDiff);
        
            lm.parameters.maxfev = 200;
            lm.parameters.xtol = 1e-8;
            lm.parameters.ftol = 1e-8;
            lm.parameters.gtol = 1e-8;
        
            int status = lm.minimize(x);
        
            vector<double> uStar(x.size());
        
            for(int i = 0; i < x.size(); i++) {
                uStar[i] = x[i];
            }
        
            map_result.result = make_original(uStar);
            Eigen::MatrixXd J(functor.values(),functor.inputs());
            numDiff.df(x,J);
            Eigen::MatrixXd H = J.transpose() * J;
            Eigen::MatrixXd H_reg = H + 1e-6 * Eigen::MatrixXd::Identity(H.rows(), H.cols());
            Eigen::MatrixXd covariance_u = H_reg.inverse();
            map_result.covariance_p = make_original_covariance(covariance_u,uStar);

        }

    public:
        struct MAP{
            vector<double> result;
            Eigen::MatrixXd covariance_p;

        };
        MAP map_result;
        // we need to pass the C(K, T)
        // We also need to pass in the guesses
    
        // I am assuming initial guesses are good - not violating constraints
        surfaceCalibrationLaplacian(vector<vector<double>> grid, vector<double> guesses, QuantLib::Date today, double spot, double r, double q){
            // [TODO] - note it is always good to bound them because they usually never explode
            // Note that we pass the guesses in P-space
            // This means that we have (index - name)
            // 0 - mu                   
            // 1 - sigma (vol-of-vol)       -> exp for > 0
            // 2 - kappa                    -> exp for > 0
            // 3 - theta                    -> exp for > 0
            // 4 - rho                      -> tanh for [-1,1]
            // 5 - lambda
            // 6 - Instantaneous volatility -> exp for > 0
            // 7 - Jump Intensity           -> exp for > 0
            // 8 - Jump Mean            
            // 9 - Jump Volatility          -> exp for > 0
            // 10 - eta                     
    
            // First let us ensure that the guesses are valid
            guesses[1] = max(guesses[1], 1e-8);
            guesses[2] = max(guesses[2], 1e-8);
            guesses[3] = max(guesses[3], 1e-8);
            guesses[4] = clamp(guesses[4], -0.999, 0.999);
            guesses[6] = max(guesses[6], 1e-8);
            guesses[7] = max(guesses[7], 1e-8);
            guesses[9] = max(guesses[9], 1e-8);
            prior_mean_ = guesses;
            // [TODO] - check if these are correct
            prior_std_ = {
                0.20,   // mu
                0.30,   // xi
                1.00,   // kappa
                0.05,   // theta
                0.30,   // rho
                0.50,   // lambda
                0.05,    // v0

                0.30,   // jump intensity
                0.05,   // jump mean
                0.10,   // jump volatility
                1.00    // eta
            };
            sigma_C_ = 1.0;

            unconstrainedParams = make_unconstrained(guesses);
            
            // These lines are added for the simplex
            grid_ = grid; today_ = today; spot_ = spot; r_ = r; q_ = q;
            x0_ = make_unconstrained(guesses);
            // Let us do calibration here itself, and later return it 
            //runSimplexCalibration();
            runLMSCalibration();
        }
    
        MAP getCalibration(){
            return map_result;
        }
    };
    
    
