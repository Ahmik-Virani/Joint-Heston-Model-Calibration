#pragma once
#include <iomanip>
#include "qe/path.hpp"
#include <vector>
#include <Eigen/Core>

namespace qe{

    using Eigen::VectorXd;
    using std::vector;

    
    struct GarchParams{
        double mu;
        double omega;
        double alpha;
        double beta;
    };

    int get_data(PPath& returns_data);
    double getNLL(GarchParams& params,PPath& ppath);
    GarchParams xToParams(const VectorXd& x);
    VectorXd Paramstox(GarchParams& params);
    vector<double> getGarchPath(GarchParams& gParams,PPath ppath);
    GarchParams garchPathFit(PPath& ppath);
    std::ostream& operator<<(std::ostream& os, const GarchParams& p);

}