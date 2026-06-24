// #include <ql/quantlib.hpp>
// #include <iostream>

// int main() {
//     std::cout << "QuantLib version: "
//               << QL_VERSION << std::endl;

//     for (int i = 0; i < 5; ++i){
//         std::cout<<i<<",";
//     }
//     return 0;
// }

#include <iostream>

// Armadillo
#include <armadillo>

// Eigen
#include <Eigen/Core>

// (Optional) QuantLib
#include <ql/version.hpp>

int main() {

    std::cout << "===== Library Versions =====\n\n";

    // Armadillo
    std::cout << "Armadillo version: "
              << arma::arma_version::major << "."
              << arma::arma_version::minor << "."
              << arma::arma_version::patch << "\n";

    // Eigen
    std::cout << "Eigen version: "
              << EIGEN_WORLD_VERSION << "."
              << EIGEN_MAJOR_VERSION << "."
              << EIGEN_MINOR_VERSION << "\n";

    // QuantLib
    std::cout << "QuantLib version: "
              << QL_VERSION << "\n";

    std::cout << "\n============================\n";

    return 0;
}
