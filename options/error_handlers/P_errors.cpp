#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <stdexcept>

using namespace std;


std::string toISODate(const std::string& input) {
    // Expected: "January 15th, 2025"

    std::stringstream ss(input);

    std::string month, day, year;
    ss >> month >> day >> year;

    // Remove trailing comma from year
    if (!year.empty() && year.back() == ',')
        year.pop_back();

    // Remove "st", "nd", "rd", "th", and comma from day
    while (!day.empty() && !isdigit(day.back()))
        day.pop_back();

    static const std::unordered_map<std::string, int> monthMap = {
        {"January",1}, {"February",2}, {"March",3}, {"April",4},
        {"May",5}, {"June",6}, {"July",7}, {"August",8},
        {"September",9}, {"October",10}, {"November",11}, {"December",12}
    };

    std::ostringstream out;
    out << year << '-'
        << std::setw(2) << std::setfill('0') << monthMap.at(month)
        << '-'
        << std::setw(2) << std::setfill('0') << std::stoi(day);

    return out.str();
}

std::string marketDateToISODate(const std::string& input) {
    // Expected: "01-Jan-2025"
    std::stringstream ss(input);

    std::string day, month, year;
    getline(ss, day, '-');
    getline(ss, month, '-');
    getline(ss, year, '-');

    static const std::unordered_map<std::string, int> monthMap = {
        {"Jan",1}, {"Feb",2}, {"Mar",3}, {"Apr",4},
        {"May",5}, {"Jun",6}, {"Jul",7}, {"Aug",8},
        {"Sep",9}, {"Oct",10}, {"Nov",11}, {"Dec",12}
    };

    std::ostringstream out;
    out << year << '-'
        << std::setw(2) << std::setfill('0') << monthMap.at(month)
        << '-'
        << std::setw(2) << std::setfill('0') << std::stoi(day);

    return out.str();
}

template <typename T>
void printVec(const vector<T>& arr){
    for(size_t i = 0; i < arr.size(); i++){
        cout<<arr[i]<<",";
    }
    cout<<"\n";
}

void getStrikes(string filename,vector<double>& strikes,vector<string>& dates){
    ifstream file(filename);
    string line;
    getline(file,line);
    while (getline(file,line)){
        stringstream ss(line);
        string token;
        vector<string>f;
        while(getline(ss,token,',')){
            if(!token.empty()) f.push_back(token);
        }
        dates.push_back(marketDateToISODate(f[0]));
        strikes.push_back(stod(f[1]));
    }
}

vector<double> getReturns(const vector<double>& strikes){
    vector<double>returns;
    for(size_t i = 1; i < strikes.size(); i++){
        returns.push_back(strikes[i]/strikes[i-1]);
    }
    return returns;
}


void readFile(string filename,string filetype,vector<string>& p_dates,vector<vector<double>>& p_params,
            unordered_map<string,double> &returns_with_dates){
    ifstream file(filename);
    string line;
    getline(file,line);

    while (getline(file,line)){
        stringstream ss(line);
        string token;
        vector<string>f;
        while (getline(ss,token,',')){
            if(!token.empty()) f.push_back(token);
        }
        if (filetype == "garch_mcmc"){
            string date = toISODate(f[0] + "," + f[1]);
            auto it = returns_with_dates.find(date);
            if (it != returns_with_dates.end()){
                vector<double>single_p_params(7,0.0);
                single_p_params[0] = stod(f[2]); //mu_mean
                single_p_params[1] = stod(f[4]); //kappa_mean
                single_p_params[2] = stod(f[6]); //theta_mean
                single_p_params[3] = stod(f[8]); //volvol_mean
                single_p_params[4] = stod(f[10]); //rho_mean
                single_p_params[5] = stod(f[12]); //vproxy
                single_p_params[6] = it->second; 
                p_dates.push_back(date);
                p_params.push_back(single_p_params);

            }
        }
        else if(filetype == "pmcmc"){
            string date = toISODate(f[0] + "," + f[1]);
            auto it = returns_with_dates.find(date);
            if (it != returns_with_dates.end()){
                vector<double>single_p_params(7,0);
                single_p_params[0] = stod(f[2]); //mu_mean
                single_p_params[1] = stod(f[4]); //kappa_mean
                single_p_params[2] = stod(f[6]); //theta_mean
                single_p_params[3] = stod(f[8]); //volvol_mean
                single_p_params[4] = stod(f[10]); //rho_mean
                single_p_params[5] = stod(f[12]); //v_t_mean
                single_p_params[6] = it->second; 
                p_dates.push_back(date);
                p_params.push_back(single_p_params);

            }

        }
        
    }
}

void writeFile(string filename,const string& variance_header,const vector<string>& p_dates,
                const vector<vector<double>>& p_params){
    std::ofstream out(filename);

    if (!out) {
        throw std::runtime_error("Could not open file");
    }

    out << "date,mu_mean,kappa_mean,theta_mean,vol-of-vol_mean,rho_mean,"
        << variance_header << ",return,one_step_nll\n";

    for (size_t row_idx = 0; row_idx < p_params.size(); ++row_idx) {
        const auto& row = p_params[row_idx];
        out << p_dates[row_idx];
        for (size_t i = 0; i < row.size(); ++i) {
            out << "," << std::setprecision(10) << row[i];
        }
        if (row.size() == 7)
            out << ",";
        out << '\n';
    }

    out.close();
}



double nll(double r_next, double mu, double var)
{
    const double dt = 1.0 / 252.0;

    // Numerical safety
    const double eps = 1e-12;
    var = std::max(var, eps);

    double mean = (mu - 0.5 * var) * dt;
    double variance = var * dt;

    return 0.5 * std::log(2.0 * M_PI * variance)
         + 0.5 * (r_next - mean) * (r_next - mean) / variance;
}

void one_step_nll(vector<vector<double>>& p_params){
    if (p_params.size() < 2)
        return;

    for(size_t i = 0; i < p_params.size() - 1; i++){
        double r_next = p_params[i+1][6];
        double mu = p_params[i][0];
        double var = p_params[i][5];
        double one_step_nll = nll(r_next,mu,var);
        p_params[i].push_back(one_step_nll);
    }   
}

int main(){
    //garch starts
    string s_path = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/S_path.csv";
    vector<double>garch_strikes;
    vector<string>garch_market_dates;
    getStrikes(s_path,garch_strikes,garch_market_dates);
    vector<double>garch_returns = getReturns(garch_strikes);
    unordered_map<string,double> garch_returns_with_dates;
    for(size_t i = 1; i < garch_market_dates.size(); i++){
        garch_returns_with_dates[garch_market_dates[i]] = garch_returns[i - 1];
    }
    vector<string>garch_p_dates;
    vector<vector<double>> garch_p_params;
    string garch_params = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/garch_calibration_params.csv";
    readFile(garch_params,"garch_mcmc",garch_p_dates,garch_p_params,garch_returns_with_dates);
    one_step_nll(garch_p_params);
    for(size_t i = 0; i < garch_p_params.size(); i++){
        printVec(garch_p_params[i]);
    }
    string garch_output_path = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/garch_nll.csv";
    writeFile(garch_output_path,"vproxy",garch_p_dates,garch_p_params);
    //garch ends
    
    //pmcmc starts
    vector<double>pmcmc_strikes;
    vector<string>pmcmc_market_dates;
    getStrikes(s_path,pmcmc_strikes,pmcmc_market_dates);
    vector<double>pmcmc_returns = getReturns(pmcmc_strikes);
    unordered_map<string,double> pmcmc_returns_with_dates;
    for(size_t i = 1; i < pmcmc_market_dates.size(); i++){
        pmcmc_returns_with_dates[pmcmc_market_dates[i]] = pmcmc_returns[i - 1];
    }
    vector<string>pmcmc_p_dates;
    vector<vector<double>> pmcmc_p_params;
    string pmcmc_params = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/pmcmc_calibration_params.csv";
    readFile(pmcmc_params,"pmcmc",pmcmc_p_dates,pmcmc_p_params,pmcmc_returns_with_dates);
    one_step_nll(pmcmc_p_params);
    for(size_t i = 0; i < pmcmc_p_params.size(); i++){
        printVec(pmcmc_p_params[i]);
    }
    string pmcmc_output_path = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/options/baseline_real/Output/pmcmc_nll.csv";
    writeFile(pmcmc_output_path,"v_t_mean",pmcmc_p_dates,pmcmc_p_params);
    //pmcmc ends

}
