#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct MarketPoint {
    string date;
    int hour;
    double spot;
};

string toISODate(const string& input) {
    // Expected: "October 23rd, 2023"
    stringstream ss(input);

    string month, day, year;
    ss >> month >> day >> year;

    if (!year.empty() && year.back() == ',')
        year.pop_back();

    while (!day.empty() && !isdigit(day.back()))
        day.pop_back();

    static const unordered_map<string, int> monthMap = {
        {"January",1}, {"February",2}, {"March",3}, {"April",4},
        {"May",5}, {"June",6}, {"July",7}, {"August",8},
        {"September",9}, {"October",10}, {"November",11}, {"December",12}
    };

    ostringstream out;
    out << year << '-'
        << setw(2) << setfill('0') << monthMap.at(month)
        << '-'
        << setw(2) << setfill('0') << stoi(day);

    return out.str();
}

vector<string> splitCSVLine(const string& line) {
    stringstream ss(line);
    string token;
    vector<string> fields;

    while (getline(ss, token, ',')) {
        if (!token.empty())
            fields.push_back(token);
    }

    return fields;
}

vector<MarketPoint> readMarketPath(const string& filename) {
    ifstream file(filename);
    if (!file)
        throw runtime_error("Could not open market path file: " + filename);

    string line;
    getline(file, line);

    vector<MarketPoint> market_points;
    while (getline(file, line)) {
        vector<string> fields = splitCSVLine(line);
        if (fields.size() < 3)
            continue;

        market_points.push_back({fields[0], stoi(fields[1]), stod(fields[2])});
    }

    return market_points;
}

vector<double> getReturns(const vector<MarketPoint>& market_points) {
    vector<double> returns;

    for (size_t i = 1; i < market_points.size(); i++)
        returns.push_back(market_points[i].spot / market_points[i - 1].spot);

    return returns;
}

void readParams(const string& filename, vector<string>& p_dates,
                vector<vector<double>>& p_params) {
    ifstream file(filename);
    if (!file)
        throw runtime_error("Could not open parameter file: " + filename);

    string line;
    getline(file, line);

    while (getline(file, line)) {
        vector<string> fields = splitCSVLine(line);
        if (fields.size() < 13)
            continue;

        string date = toISODate(fields[0] + "," + fields[1]);

        vector<double> single_p_params(7, 0.0);
        single_p_params[0] = stod(fields[2]);  // mu_mean
        single_p_params[1] = stod(fields[4]);  // kappa_mean
        single_p_params[2] = stod(fields[6]);  // theta_mean
        single_p_params[3] = stod(fields[8]);  // volvol_mean
        single_p_params[4] = stod(fields[10]); // rho_mean
        single_p_params[5] = stod(fields[12]); // vproxy or v_t_mean

        p_dates.push_back(date);
        p_params.push_back(single_p_params);
    }
}

void attachHourlyReturns(const vector<MarketPoint>& market_points,
                         const vector<double>& returns,
                         vector<string>& p_dates,
                         vector<int>& p_hours,
                         vector<vector<double>>& p_params) {
    if (p_params.empty())
        return;

    if (market_points.size() < p_params.size() + 1)
        throw runtime_error("Not enough market path rows to align hourly parameter rows");

    const size_t first_market_idx = market_points.size() - p_params.size();

    for (size_t row_idx = 0; row_idx < p_params.size(); row_idx++) {
        const size_t market_idx = first_market_idx + row_idx;
        if (market_idx == 0)
            throw runtime_error("Cannot attach a return for the first market path row");
        if (p_dates[row_idx] != market_points[market_idx].date) {
            throw runtime_error("Calibration date does not match aligned hourly market row");
        }

        p_hours.push_back(market_points[market_idx].hour);
        p_params[row_idx][6] = returns[market_idx - 1];
    }
}

void writeFile(const string& filename, const string& variance_header,
               const vector<string>& p_dates, const vector<int>& p_hours,
               const vector<vector<double>>& p_params) {
    ofstream out(filename);
    if (!out)
        throw runtime_error("Could not open output file: " + filename);

    out << "date,hour,mu_mean,kappa_mean,theta_mean,vol-of-vol_mean,rho_mean,"
        << variance_header << ",return,one_step_nll\n";

    for (size_t row_idx = 0; row_idx < p_params.size(); row_idx++) {
        const auto& row = p_params[row_idx];
        out << p_dates[row_idx] << "," << p_hours[row_idx];

        for (size_t i = 0; i < row.size(); i++)
            out << "," << setprecision(10) << row[i];

        if (row.size() == 7)
            out << ",";

        out << '\n';
    }
}

double nll(double r_next, double mu, double var) {
    const double dt = 1.0 / 252.0;
    const double eps = 1e-12;

    var = max(var, eps);

    double mean = (mu - 0.5 * var) * dt;
    double variance = var * dt;

    return 0.5 * log(2.0 * M_PI * variance)
         + 0.5 * (r_next - mean) * (r_next - mean) / variance;
}

void oneStepNLL(vector<vector<double>>& p_params) {
    if (p_params.size() < 2)
        return;

    for (size_t i = 0; i < p_params.size() - 1; i++) {
        double r_next = p_params[i + 1][6];
        double mu = p_params[i][0];
        double var = p_params[i][5];
        p_params[i].push_back(nll(r_next, mu, var));
    }
}

int main() {
    const string base_path = "/data1/sandesh/arka/Joint-Heston-Model-Calibration/bates/baseline_real";
    const string s_path = base_path + "/S_path.csv";

    vector<MarketPoint> bates_market_points = readMarketPath(s_path);
    vector<double> bates_returns = getReturns(bates_market_points);

    vector<string> bates_garch_dates;
    vector<int> bates_garch_hours;
    vector<vector<double>> bates_garch_p_params;
    string bates_garch_params = base_path + "/Output/garch_calibration_params.csv";
    readParams(bates_garch_params, bates_garch_dates, bates_garch_p_params);
    attachHourlyReturns(bates_market_points, bates_returns, bates_garch_dates,
                        bates_garch_hours, bates_garch_p_params);
    oneStepNLL(bates_garch_p_params);
    writeFile(base_path + "/Output/garch_nll.csv", "vproxy", bates_garch_dates,
              bates_garch_hours, bates_garch_p_params);

    vector<string> bates_pmcmc_dates;
    vector<int> bates_pmcmc_hours;
    vector<vector<double>> bates_pmcmc_p_params;
    string bates_pmcmc_params = base_path + "/Output/pmcmc_calibration_params.csv";
    readParams(bates_pmcmc_params, bates_pmcmc_dates, bates_pmcmc_p_params);
    attachHourlyReturns(bates_market_points, bates_returns, bates_pmcmc_dates,
                        bates_pmcmc_hours, bates_pmcmc_p_params);
    oneStepNLL(bates_pmcmc_p_params);
    writeFile(base_path + "/Output/pmcmc_nll.csv", "v_t_mean", bates_pmcmc_dates,
              bates_pmcmc_hours, bates_pmcmc_p_params);

    return 0;
}
