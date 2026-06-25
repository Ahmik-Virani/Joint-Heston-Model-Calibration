#include <qe/path.hpp>
#include <qe/garch.hpp>

#include <ql/quantlib.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../src/getRealData.cpp"

using namespace qe;
using namespace QuantLib;
using namespace std;

PPath buildExpandingPath(getRealData& data, int end_t) {
    PPath ppath;

    ppath.logS.resize(end_t + 1);
    ppath.returns.resize(end_t);

    for (int i = 0; i <= end_t; i++) {
        ppath.logS[i] = log(data.get_S(i));
    }

    for (int i = 1; i <= end_t; ++i) {
        ppath.returns[i - 1] = ppath.logS[i] - ppath.logS[i - 1];
    }

    return ppath;
}

int main() {
    getRealData data;
    const int time_steps = data.get_time_steps();
    const Size steps = 252;

    cout << "Data time steps: " << time_steps << "\n";
    cout << fixed << setprecision(8);
    int starting_steps = 10;
    if (time_steps <= 1) {
        cout << "Not enough data to build a return path; need at least 2 prices.\n";
        return 0;
    }

    if (starting_steps >= time_steps) {
        const int end_t = time_steps - 1;
        PPath ppath = buildExpandingPath(data, end_t);
        GarchParams gParams = garchPathFit(ppath);
        vector<double> daily_hPath = getGarchPath(gParams, ppath);
        vector<double> annual_hPath(daily_hPath.size(), 0.0);
        for (int i = 0; i < static_cast<int>(daily_hPath.size()); ++i) {
            annual_hPath[i] = daily_hPath[i] * steps;
        }

        cout << "\nOnly " << time_steps << " data points available, which is not more than "
            << starting_steps << "; skipping expanding window and fitting the full path once.\n";
        cout << "Full path through t=" << end_t
            << " S=" << data.get_S(end_t)
            << " returns=" << ppath.returns.size()
            << " garch_path=" << annual_hPath.size() << "\n";
        cout << "Returns: ";
        for (double r : ppath.returns) {
            cout << r << " ";
        }
        cout << "\n" << gParams;
        cout << "Annualized GARCH h: ";
        for (double h : annual_hPath) {
            cout << h << " ";
        }
        cout << "\n";
    }
    else{
        for (int end_t = starting_steps; end_t < time_steps; end_t++) {
            PPath ppath = buildExpandingPath(data, end_t);
            GarchParams gParams = garchPathFit(ppath);
            vector<double> daily_hPath = getGarchPath(gParams, ppath);
            vector<double> annual_hPath(daily_hPath.size(), 0.0);

            for (int i = 0; i < static_cast<int>(daily_hPath.size()); ++i) {
                annual_hPath[i] = daily_hPath[i] * steps;
            }

            cout << "\nExpanding window through t=" << end_t
                << " S=" << data.get_S(end_t)
                << " returns=" << ppath.returns.size()
                << " garch_path=" << annual_hPath.size() << "\n";
            cout << "Returns: ";
            for (double r : ppath.returns) {
                cout << r << " ";
            }
            cout << "\n" << gParams;
            cout << "Annualized GARCH h: ";
            for (double h : annual_hPath) {
                cout << h << " ";
            }
            cout << "\n";
        }
    }

    return 0;
}
