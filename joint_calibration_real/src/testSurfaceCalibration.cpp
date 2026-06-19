#include <ql/quantlib.hpp>
#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "surfaceCalibration.cpp"
#include "getRealData.cpp"

// Build from joint_calibration_real/ so getRealData can find S_path.csv, r_path.csv, and C_grid.csv:
// g++ src/testSurfaceCalibration.cpp -std=c++17 -I$HOME/local/include -I/data1/sandesh/arka/eigen-3.4.0 -L$HOME/local/lib -lQuantLib -Wl,-rpath,$HOME/local/lib -o testSurfaceCalibration
// Run:
// ./testSurfaceCalibration [particle_index] [time_index] -> ./testSurfaceCalibration 0 0

using namespace std;

class testCalibration {
private:
    static void print_particle(const string& label, const vector<double>& particle) {
        static const vector<string> names = {
            "mu", "sigma", "kappa", "theta", "rho", "lambda", "v0"
        };

        cout << label << '\n';
        cout << fixed << setprecision(10);
        for (size_t i = 0; i < particle.size(); ++i) {
            const string& name = i < names.size() ? names[i] : string("param");
            cout << "  " << name << "[" << i << "] = " << particle[i] << '\n';
        }
    }

public:
    getRealData data;

    vector<double> create_guess() {
        return data.get_guess();
    }

    vector<double> initialization_simplex(
        int particle_index,
        int time_index,
        const vector<double>& guess
    ) {
        cout << "Running simplex initialization for particle "
             << particle_index << ", time " << time_index << '\n';

        surfaceCalibration this_particle(
            data.get_grid(time_index),
            guess,
            data.get_date(time_index),
            data.get_S(time_index),
            data.get_r(time_index),
            data.get_q(time_index)
        );

        return this_particle.getCalibration();
    }

    surfaceCalibrationLaplacian::MAP initialization_lms(
        int particle_index,
        int time_index,
        const vector<double>& guess
    ) {
        cout << "Running LMS initialization for particle "
             << particle_index << ", time " << time_index << '\n';

        surfaceCalibrationLaplacian this_particle(
            data.get_grid(time_index),
            guess,
            data.get_date(time_index),
            data.get_S(time_index),
            data.get_r(time_index),
            data.get_q(time_index)
        );

        surfaceCalibrationLaplacian::MAP map_result = this_particle.getCalibration();
        return map_result;
    }

    void run_single_particle_time(int particle_index, int time_index) {
        if (time_index < 0 || time_index >= data.get_time_steps()) {
            throw out_of_range("time_index is outside the available data range");
        }

        const vector<double> guess = create_guess();

        cout << "Grid size at time " << time_index
             << " = " << data.get_grid(time_index).size() << '\n';
        cout << "Spot = " << data.get_S(time_index)
             << ", r = " << data.get_r(time_index)
             << ", q = " << data.get_q(time_index) << '\n';

        print_particle("Initial guess", guess);

        const vector<double> simplex_particle =
            initialization_simplex(particle_index, time_index, guess);
        print_particle("Simplex calibration result", simplex_particle);

        const surfaceCalibrationLaplacian::MAP lms_map =
            initialization_lms(particle_index, time_index, guess);
        print_particle("LMS/Laplacian calibration result", lms_map.result);

        cout << "LMS/Laplacian covariance matrix in P-space:\n"
             << lms_map.covariance_p << '\n';
    }
};

int main(int argc, char** argv) {
    int particle_index = 0;
    int time_index = 0;

    if (argc > 1) {
        particle_index = stoi(argv[1]);
    }
    if (argc > 2) {
        time_index = stoi(argv[2]);
    }

    try {
        testCalibration test;
        test.run_single_particle_time(particle_index, time_index);
    } catch (const exception& e) {
        cerr << "testSurfaceCalibration failed: " << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
