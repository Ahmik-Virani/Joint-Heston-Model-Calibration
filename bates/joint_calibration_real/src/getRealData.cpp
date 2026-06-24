#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <map>
#include <random>
#include <string>
#include <sstream>
#include <ql/time/date.hpp>
#include <ql/time/daycounters/actual365fixed.hpp>

#include "callPrice.cpp"

using namespace std;
using namespace QuantLib;

class getRealData{
private:
    // number of days we will be going ahead
    int t;

    vector<double> S;
    vector<pair<string, int>> dates_and_hour;
    map<string, double> rates;      // a map which maps the date to the rate
    map<int, pair<string, int>> index_to_date_and_hour; // a map to match the index to the date and hour

    // Grid for each time steps
    // time, K, T
    vector<vector<vector<double>>> grid;

    // convert string to QunatLib Date
    QuantLib::Date parseDate(const std::string& s) {
        std::stringstream ss(s);

        std::string dayStr, monthStr, yearStr;

        // here date is yyyy-mm-dd
        getline(ss, yearStr, '-');
        getline(ss, monthStr, '-');
        getline(ss, dayStr, '-');

        int day = std::stoi(dayStr);
        int monthInt = std::stoi(monthStr);
        int year = std::stoi(yearStr);

        Month month;

        if      (monthInt == 1) month = January;
        else if (monthInt == 2) month = February;
        else if (monthInt == 3) month = March;
        else if (monthInt == 4) month = April;
        else if (monthInt == 5) month = May;
        else if (monthInt == 6) month = June;
        else if (monthInt == 7) month = July;
        else if (monthInt == 8) month = August;
        else if (monthInt == 9) month = September;
        else if (monthInt == 10) month = October;
        else if (monthInt == 11) month = November;
        else if (monthInt == 12) month = December;
        else throw std::runtime_error("Invalid month");

        return QuantLib::Date(day, month, year);
    }

public:
    // Here we will initialize all the values
    getRealData(){
        cout<<"I am here"<<endl;
        // get the stock path
        ifstream s_inp("S_path.csv");
        if (!s_inp.is_open()) {
            throw runtime_error("Could not open S_path.csv");
        }
        string line;
        // ignore the header
        getline(s_inp, line);

        while(getline(s_inp, line)){
            stringstream ss(line);
            string date;
            string hour;
            string _val;
            getline(ss, date, ',');
            getline(ss, hour, ',');
            getline(ss, _val, ',');

            double val = stod(_val);
            int hr = stoi(hour);
            S.push_back(val);
            dates_and_hour.push_back({date,hr});
        }

        // now we are done getting the stock path
        s_inp.close();

        // get the risk free rate
        // please note that this is daily
        ifstream r_inp("r_path.csv");
        if (!r_inp.is_open()) {
            throw runtime_error("Could not open r_path.csv");
        }
        // ignore the header
        getline(r_inp, line);

        // ensure initially the rates map is clear
        rates.clear();
        while(getline(r_inp, line)){
            stringstream ss(line);
            string date, _r;
            getline(ss, date, ',');
            getline(ss, _r, ',');

            double rate = stod(_r);

            // since rates are daily, we store it using index of date
            rates[date] = rate;
        }

        // we are done with the rates
        r_inp.close();

        // we need to get the call option grid
        ifstream C_inp("C_grid.csv");
        if (!C_inp.is_open()) {
            throw runtime_error("Could not open C_grid.csv");
        }
        getline(C_inp, line);

        int cur_ind = 0;
        
        grid.resize(S.size());
       
        // need fixed 365 day calendar
        while(getline(C_inp, line)){
            stringstream ss(line);
            string date, hour, Expiry, Strike_Price, Close;
            getline(ss, date, ',');
            getline(ss, hour, ',');
            getline(ss, Expiry, ',');
            getline(ss, Strike_Price, ',');
            getline(ss, Close, ',');

            QuantLib::Date this_date = parseDate(date);
            QuantLib::Date this_expiry = parseDate(Expiry);
            
            QuantLib::Integer calendarDays = this_expiry - this_date;
            // [TODO] check if expiry hour on expiry date is correct or not
            int expiry_hour = 8;
            int hr = stoi(hour);
            int hours = expiry_hour - hr;
            double maturityYears = std::max(1.0 / (365.0 * 24.0), (calendarDays / 365.0 + hours / (24.0 * 365.0)));

            if(index_to_date_and_hour.empty() || index_to_date_and_hour[cur_ind-1]!=make_pair(date, hr)){
                index_to_date_and_hour[cur_ind++] = {date,hr};
            }
            grid[cur_ind-1].push_back({stod(Strike_Price), maturityYears, stod(Close)});
        }
        C_inp.close();
    }

    // get value of S at a particular day
    double get_S(int i){
        return S[i];
    }

    // get rate of a particular index
    double get_r(int i){
        // get which date we need
        string date = index_to_date_and_hour[i].first;

        return (rates[date]/100.0);
    }

    // [TODO] - can do this
    double get_q(int i){
        return 0.0;
    }

    // get log return at time i
    double get_log_return(int i){
        return log(S[i] / S[i-1]);
    }

    int get_time_steps(){
        return index_to_date_and_hour.size();
    }

    // get the grid for a particular day
    vector<vector<double>> get_grid(int i){
        return grid[i];
    }

    // [TODO] - update to include bates terms
    // get some guess for the parameters
    vector<double> get_guess(){
        // 0 - mu
        // 1 - sigma (vol-of-vol)
        // 2 - kappa
        // 3 - theta
        // 4 - rho
        // 5 - lambda
        // 6 - Instantaneous volatility
        // 7 - Jump Intensity
        // 8 - Jump Mean
        // 9 - Jump Volatility
        return {0.1, 0.5, 2, 0.04, -0.7, 0.5, 0.20, 0.5, -0.05, 0.2};
    }

    // get date for index i
    pair<QuantLib::Date,int> get_date_and_hour(int i){
        return {parseDate(index_to_date_and_hour[i].first), index_to_date_and_hour[i].second};
    }

    // a function which takes in heston parameters at time t
    // and checks how will they fit the current grid
    // it returns the error
    double get_penalty(int t, double v0, double kappaQ, double thetaQ, double xi, double rho, double jumpIntensity, double jumpMean, double jumpVolatility){

        double total_error = 0.0;

        // go through all the vectors of the grid at timestep t
        for(int i = 0 ; i < grid[t].size() ; i++){
            double computed_call_price = 0.0;
            computed_call_price = qe::batesCallPrice(
                get_S(t),                          // spot
                grid[t][i][0],                     // fixed strike
                grid[t][i][1],                     // maturity in years (scalar)
                get_r(t),                          // r
                get_q(t),                          // q
                max(1e-8, v0),                     // v0
                max(1e-8, kappaQ),                 // kappaQ
                max(1e-8, thetaQ),                 // thetaQ
                max(1e-8, xi),                     // xi
                clamp(rho, -0.999, 0.999),         // rho
                max(1e-8, jumpIntensity),
                jumpMean,
                max(1e-8, jumpVolatility),
                get_date(t)                        // anchor pricing at this day
            );   

            double true_price = grid[t][i][2];
            total_error += abs(true_price - computed_call_price);
        }
        double average_error = total_error / grid[t].size();
        return exp(-average_error);
    }

    // same as earlier, but we print to a log file
    void get_penalty(int t, double v0, double kappaQ, double thetaQ, double xi, double rho, ofstream &log_file){

        // go through all the vectors of the grid at timestep t
        for(int i = 0 ; i < grid[t].size() ; i++){
            double computed_call_price = 0.0;
            computed_call_price = qe::batesCallPrice(
                get_S(t),                          // spot
                grid[t][i][0],                     // fixed strike
                grid[t][i][1],                     // maturity in years (scalar)
                get_r(t),                          // r
                get_q(t),                          // q
                max(1e-8, v0),                     // v0
                max(1e-8, kappaQ),                 // kappaQ
                max(1e-8, thetaQ),                 // thetaQ
                max(1e-8, xi),                     // xi
                clamp(rho, -0.999, 0.999),         // rho
                max(1e-8, jumpIntensity),
                jumpMean,
                max(1e-8, jumpVolatility),
                get_date(t)                        // anchor pricing at this day
            );   

            double true_price = grid[t][i][2];
            double error = abs(true_price - computed_call_price);
            // date, hour, strike, maturity, true_price, computed_price, abs_error
            log_file << index_to_date_and_hour[t].first << ',' << index_to_date_and_hour[t].second << ',' << grid[t][i][0] << ',' << grid[t][i][1] << ',' << true_price << ',' << computed_call_price << ',' << error << '\n';
        }
    }
};