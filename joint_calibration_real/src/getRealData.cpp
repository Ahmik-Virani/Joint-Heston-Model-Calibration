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
// #include <ql/time/daycounters/actual365fixed.hpp>
#include <ql/time/calendars/india.hpp>

#include "callPrice.cpp"

using namespace std;
using namespace QuantLib;

class getRealData{
private:
    // number of days we will be going ahead
    int t;

    vector<double> S;
    vector<string> dates;
    map<string, double> rates;      // a map which maps the time to the rate
    map<int, string> index_to_date; // a map to match the index to the date

    // Grid for each time steps
    // time, K, T
    vector<vector<vector<double>>> grid;

    // convert string to QunatLib Date
    QuantLib::Date parseDate(const std::string& s) {
        std::stringstream ss(s);

        std::string dayStr, monthStr, yearStr;

        getline(ss, dayStr, '-');
        getline(ss, monthStr, '-');
        getline(ss, yearStr, '-');

        int day = std::stoi(dayStr);
        int year = std::stoi(yearStr);

        Month month;

        if      (monthStr == "Jan") month = January;
        else if (monthStr == "Feb") month = February;
        else if (monthStr == "Mar") month = March;
        else if (monthStr == "Apr") month = April;
        else if (monthStr == "May") month = May;
        else if (monthStr == "Jun") month = June;
        else if (monthStr == "Jul") month = July;
        else if (monthStr == "Aug") month = August;
        else if (monthStr == "Sep") month = September;
        else if (monthStr == "Oct") month = October;
        else if (monthStr == "Nov") month = November;
        else if (monthStr == "Dec") month = December;
        else throw std::runtime_error("Invalid month");

        return QuantLib::Date(day, month, year);
    }

public:
    // Here we will initialize all the values
    getRealData(){
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
            string _val;
            getline(ss, date, ',');
            getline(ss, _val, ',');

            double val = stod(_val);
            S.push_back(val);
            dates.push_back(date);
        }

        // now we are done getting the stock path
        s_inp.close();

        // get the risk free rate
        // please note that this is monthly
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

            // since rates are monthly
            // we only store month and year
            string store_date = date.substr(3);
            rates[store_date] = rate;
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
       
        // get the indian business calendar
        QuantLib::India calendar;
        while(getline(C_inp, line)){
            stringstream ss(line);
            string date, Expiry, Strike_Price, Close;
            getline(ss, date, ',');
            getline(ss, Expiry, ',');
            getline(ss, Strike_Price, ',');
            getline(ss, Close, ',');

            QuantLib::Date this_date = parseDate(date);
            QuantLib::Date this_expiry = parseDate(Expiry);
            Size tradingDays = calendar.businessDaysBetween(this_date, this_expiry);

            if(index_to_date.empty() || index_to_date[cur_ind-1]!=date){
                index_to_date[cur_ind++] = date;
            }
            grid[cur_ind-1].push_back({stod(Strike_Price), tradingDays/252.0, stod(Close)});
        }
        C_inp.close();
    }

    // get value of S at a particular day
    double get_S(int i){
        return S[i];
    }

    // get rate of a particular date
    double get_r(int i){
        // get which date we need
        string date = index_to_date[i];
        
        // we only need month and year
        string fetch_date = date.substr(3);

        return rates[fetch_date];
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
        return index_to_date.size();
    }

    // get the grid for a particular day
    vector<vector<double>> get_grid(int i){
        return grid[i];
    }

    // get some guess for the parameters
    vector<double> get_guess(){
        return {0.1, 0.5, 2, 0.04, -0.7, 0.5, 0.20};
    }

    // get date for index i
    QuantLib::Date get_date(int i){
        return parseDate(index_to_date[i]);
    }

    // a function which takes in heston parameters at time t
    // and checks how will they fit the current grid
    // it returns the error
    double get_penalty(int t, double v0, double kappaQ, double thetaQ, double xi, double rho){

        double total_error = 0.0;

        // go through all the vectors of the grid at timestep t
        for(int i = 0 ; i < grid[t].size() ; i++){
            double computed_call_price = 0.0;
            computed_call_price = qe::hestonCallPrice(
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
            computed_call_price = qe::hestonCallPrice(
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
                get_date(t)                        // anchor pricing at this day
            );   

            double true_price = grid[t][i][2];
            double error = abs(true_price - computed_call_price);
            // date, strike, maturity, true_price, computed_price, abs_error
            log_file << index_to_date[t] << ',' << grid[t][i][0] << ',' << grid[t][i][1] << ',' << true_price << ',' << computed_call_price << ',' << error << '\n';
        }
    }
};