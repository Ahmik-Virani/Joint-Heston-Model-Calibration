#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <random>
#include <string>

#include "callPrice.cpp"

using namespace std;

class generateData{
private:
    // number of days we will be going ahead
    int t;
    double dt;

    // first let us define the heston parameters and other necessary constants
    double S0, v0, mu, kappa_P, theta_P, xi, rho, r, q, lambda;

    vector<double> S, v;

    // strike prices
    vector<double> K;

    // time to maturity in years
    vector<double> T;

    // The date of maturity
    string maturity;
    string today;

    // Grid for each time steps
    // time, K, T
    vector<vector<vector<double>>> grid;

    // An array for the values of C we are finding (fixed maturity and strike)
    vector<double> C;

    double strike;

public:
    // Here we will initialize all the values
    generateData(){
        t = 20;
        dt = 1.0/252;
        S0 = 100.0, v0 = 0.04, mu = 0.05, kappa_P = 1.5, theta_P = 0.04, xi = 0.3, rho = -0.7, r = 0.03, q = 0.00, lambda = 0.50;

        K = {80.0, 90.0, 100.0, 110.0, 120.0};
        T = {1.0/12, 3.0/12, 6.0/12, 1.0};

        strike = 100.0;

        // push the current values of S and v
        S.resize(t+1);
        v.resize(t+1);

        S[0] = S0;
        v[0] = v0;

        // Generate the S path and v path using heston equation

        random_device rd{};
        mt19937 gen{rd()};
        normal_distribution<double> d{0.0, 1.0};
        for(int i = 1 ; i <= t ; i++){
            double dZ2 = d(gen) * sqrt(dt);
            v[i] = v[i-1] + kappa_P * (theta_P - v[i-1]) * dt + xi * sqrt(v[i-1]) * dZ2;
            v[i] = max(1e-8, v[i]);

            double dW = d(gen) * sqrt(dt);
            double dZ1 = rho * dZ2 + sqrt(1.0 - rho * rho) * dW;
            S[i] = S[i-1] + mu * S[i-1] * dt + sqrt(v[i]) * S[i-1] * dZ1;
        }

        maturity = "2026-12-31";
        today = "2026-03-09";

        QuantLib::Date todayDate = QuantLib::DateParser::parseISO(today);
        QuantLib::Date maturityDate = QuantLib::DateParser::parseISO(maturity);
        QuantLib::Calendar cal = QuantLib::TARGET();

        double kappa_Q = kappa_P + lambda;
        double theta_Q = kappa_P * theta_P / (kappa_P + lambda);

        grid.resize(t+1, vector<vector<double>> (K.size(), vector<double> (T.size())));

        // Make the grid C(K,T) for each day
        for(int i = 0 ; i <= t ; i++){
            for(int k = 0 ; k < K.size() ; k++){
                for(int tau = 0 ; tau < T.size() ; tau++){
                    QuantLib::Date evalDate = cal.advance(todayDate, i, QuantLib::Days);
                    grid[i][k][tau] = qe::hestonCallPrice(S[i], K[k], T[tau], r, q, v[i], kappa_Q, theta_Q, xi, rho, evalDate);
                }
            }
        }

        // compute the values of true C
        C.resize(t+1);
        for(int i = 0 ; i <= t ; i++){
            QuantLib::Date evalDate = cal.advance(todayDate, i, QuantLib::Days);
            int nTradingDays = cal.businessDaysBetween(evalDate, maturityDate, false, true);
            C[i] = qe::hestonCallPrice(S[i], strike, nTradingDays/252.0, r, q, v[i], kappa_Q, theta_Q, xi, rho, evalDate);
        }
    }

    // get value of S at a particular day
    double get_S(int i){
        return S[i];
    }

    // Get the value of v at a particular day
    double get_v(int i){
        return v[i];
    }

    vector<double> get_T(){
        return T;
    }

    vector<double> get_K(){
        return K;
    }

    double get_r(){
        return r;
    }

    double get_q(){
        return q;
    }

    // get log return at time i
    double get_log_return(int i){
        return log(S[i] / S[i-1]);
    }

    string get_maturity_date(){
        return maturity;
    }

    int get_time_steps(){
        return t;
    }

    // get the grid for a particular day
    vector<vector<double>> get_grid(int i){
        return grid[i];
    }

    // Get the value of true call option price at a particular day
    double get_C(int i){
        return C[i];
    }

    // today + i days
    QuantLib::Date get_todays_date(int i){
        QuantLib::Date todayDate = QuantLib::DateParser::parseISO(today);
        QuantLib::Calendar cal = QuantLib::TARGET();
        QuantLib::Date evalDate = cal.advance(todayDate, i, QuantLib::Days);
        return evalDate;
    }

    QuantLib::Date get_maturity(){
        QuantLib::Date maturityDate = QuantLib::DateParser::parseISO(maturity);
        return maturityDate;
    }

    // after i days from today
    double get_time_to_maturity(int i){
        QuantLib::Date todayDate = QuantLib::DateParser::parseISO(today);
        QuantLib::Date maturityDate = QuantLib::DateParser::parseISO(maturity);
        QuantLib::Calendar cal = QuantLib::TARGET();
        QuantLib::Date evalDate = cal.advance(todayDate, i, QuantLib::Days);
        int nTradingDays = cal.businessDaysBetween(evalDate, maturityDate, false, true);

        return nTradingDays/252.0;
    }

    double get_strike(){
        return strike;
    }

    vector<double> get_true_params(){
        return {mu, xi, kappa_P, theta_P, rho, lambda, v0};
    }
};