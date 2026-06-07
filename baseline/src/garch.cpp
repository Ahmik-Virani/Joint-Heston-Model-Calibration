#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <Eigen/Core>
#include <LBFGS.h>
#include <iomanip>
#include "qe/garch.hpp"
#include "qe/path.hpp"


using namespace std;
using namespace Eigen;

namespace qe{
int get_data(PPath& returns_data){
    
    ifstream file("./logReturns.csv");
    if(!file.is_open()){
        cerr<<"Error opening file!"<<endl;
        return -1;
    }
    vector<double>returns;
    vector<double>v;
    
    string line;

    getline(file,line); //this is for skipping the header;

    while(getline(file,line)){
        stringstream ss(line);
        string token;
        getline(ss,token,',');
        double logS_val = stod(token);

        getline(ss,token,',');
        double v_val = stod(token);

        getline(ss,token,',');
        double returns_val = stod(token);

        returns_data.logS.push_back(logS_val);
        returns_data.v.push_back(v_val);
        returns_data.returns.push_back(returns_val);

    }
    return 1;
}

std::ostream& operator<<(std::ostream& os, const GarchParams& p)
{
    os << std::fixed << setprecision(8);
    os << "GARCH Parameters\n";
    os << "----------------\n";
    os << "mu     : " << p.mu    << "\n";
    os << "omega  : " << p.omega << "\n";
    os << "alpha  : " << p.alpha << "\n";
    os << "beta   : " << p.beta  << "\n";
    os << "alpha+beta : " << (p.alpha + p.beta) << "\n";
    return os;
}



double getNLL(GarchParams& params,PPath& ppath){
    
    int n = ppath.returns.size();
    if (ppath.returns[n-1] == 0.0) n = n - 1;
    
    if (params.omega < 0) throw invalid_argument("Omega must be grater than 0");
    if (params.alpha < 0) throw invalid_argument("Alpha must be grater than 0");
    if (params.beta < 0) throw invalid_argument("Beta must be grater than 0");
    if (params.alpha + params.beta > 1.0) throw invalid_argument("Alpha + Beta must be less than 1");
    
    double mu = params.mu;
    double omega = params.omega;
    double alpha = params.alpha;
    double beta = params.beta;
    double h0 = omega/(1 - alpha - beta);
    double nll = 0.0;
    if(!(h0 > 1e-16) || !isfinite(h0)) throw invalid_argument("h0 is infinite.");
    
    double t_h = h0;

    for(int i = 0; i < n;i++){
        double t_err = ppath.returns[i] - mu;
        nll += 0.5*(log(t_h) + t_err * t_err/t_h);
        t_h = omega + alpha * t_err * t_err + beta * t_h;
        if(!(t_h > 1e-16) || !isfinite(t_h)){
            return numeric_limits<double>::infinity();
        } 
    }
    return nll;
}


GarchParams xToParams(const VectorXd& x){
    GarchParams initParams;
    initParams.mu = x[0];
    initParams.omega = exp(x[1]);
    double a = exp(x[2]);
    double b = exp(x[3]);
    double denom = 1 + a + b;
    initParams.alpha = a/denom;
    initParams.beta = b/denom;
    return initParams;
}

VectorXd Paramstox(GarchParams& params){
    VectorXd x(4);

    x[0] = params.mu;
    x[1] = log(params.omega);

    // invert softmax
    // a = alpha / (1 - alpha - beta)
    // b = beta  / (1 - alpha - beta)

    double denom = 1.0 - params.alpha - params.beta;

    x[2] = log(params.alpha / denom);
    x[3] = log(params.beta  / denom);

    return x;
}

vector<double> getGarchPath(GarchParams& gParams,PPath ppath){
    int n = ppath.returns.size();
    if (ppath.returns[n-1] == 0.0) n = n - 1;
    vector<double>hPath(n+1,0.0);
    double mu = gParams.mu;
    double omega = gParams.omega;
    double alpha = gParams.alpha;
    double beta = gParams.beta;
    hPath[0] = omega/(1 - alpha - beta);
    for (int i = 0; i < n;i++){
        double t_err = ppath.returns[i] - mu;
        hPath[i+1] = omega + alpha * t_err * t_err + beta * hPath[i];
    }
    return hPath;

}

struct GarchObjective{
    PPath& data;
    double eps = 1e-6;
    GarchObjective(PPath& path): data(path){}

    double evaluate(const VectorXd& x){
        GarchParams params = xToParams(x);
        return getNLL(params,data);
    }

    double operator()(const VectorXd& x, VectorXd& grad){
        double fx = evaluate(x);
        int dim = x.size();
        grad.resize(dim);
        for (int i = 0; i < dim;i++){
            VectorXd x_plus = x;
            VectorXd x_minus = x;
            x_plus[i] = x[i] + eps;
            x_minus[i] = x[i] - eps;
            double f_plus = evaluate(x_plus);
            double f_minus = evaluate(x_minus);
            grad[i] = (f_plus - f_minus)/(2 * eps);
        }
        return fx;
    }

};


GarchParams garchPathFit(PPath& ppath){
   
    double init_mu = 0;
    double init_var = 0;
    for(int i = 0; i < ppath.returns.size() - 1; i++){
        init_mu += ppath.returns[i];
        init_var += (ppath.returns[i] * ppath.returns[i]);
    }
    init_mu = init_mu/(ppath.returns.size() - 1);
    init_var = init_var/(ppath.returns.size() - 1) - init_mu * init_mu; 

    // cout<<"Initial Mu:"<<init_mu<<endl;
    // cout<<"Initial Var:"<<init_var<<endl;
    
    GarchParams initGuess;
    initGuess.mu = init_mu;
    initGuess.alpha = 0.05;
    initGuess.beta = 0.90;
    initGuess.omega = max(1e-6,init_var *(1 - initGuess.alpha - initGuess.beta));

    VectorXd x = Paramstox(initGuess);
    // cout<< "Initial X array:"<<endl;
    // cout<<"Mu unconstrained:"<<x[0]<<endl;
    // cout<<"Omega unconstrained:"<<x[1]<<endl;
    // cout<<"Alpha unconstrained:"<<x[2]<<endl;
    // cout<<"Beta unconstrained:"<<x[3]<<endl;

    GarchObjective garchobj(ppath);
    LBFGSpp::LBFGSParam<double> optParam;
    optParam.epsilon = 1e-8;
    optParam.max_iterations = 200;
    
    LBFGSpp::LBFGSSolver<double> solver(optParam);
    double fx = 0.0;
    int niter = solver.minimize(garchobj,x,fx);
    // cout<<"After optimization:"<<endl;
    // cout<<"Mu:"<<x[0] << "\n" << "Omega:" << x[1] << "\n" << "Alpha:" << x[2] << "\n" << "Beta:" << x[3] << "\n"; 
    GarchParams gParams = xToParams(x);
    //cout<<gParams<<endl;
    //vector<double>hPath = getGarchPath(gParams,ppath);
    //cout<<"v path length:"<<ppath.v.size()<<endl;
    //cout<<"h path length:"<<hPath.size()<<endl;
    return gParams;
}
} //namespace qe;