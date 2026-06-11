#pragma once

#include <ql/quantlib.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <thread>
#include <atomic>

#include "surfaceCalibration.cpp"
#include "surfaceCalibration_cmaes.cpp"
#include "surfaceCalibration_multistart_lm.cpp"
#include "surfaceCalibration_multistart_simplex.cpp"
#include "getRealData.cpp"

using namespace std;

class smoothStateCalibration {
private:
    getRealData data;

    int no_of_threads;  // the number of threads
    // define some locks
    mutex set_lock_1, set_lock_2;

    // index to next particle
    atomic<int> next_particle{0};

    int N;      // Number of particles

    int number_of_samples; // Number of samples of ancestral sampling

    const double dt = 1.0 / 252.0;

    // Initialize how many time steps to move ahead
    int time_steps;

    // counter for the number of time reinitialization happens
    atomic<int> ctr_reinit{0};

    // A global random generator
    mt19937 rng{random_device{}()};

    // lambda assumed constant

    // Each particle stores 
    // 1. Heston Parameters
    // 2. Lambda (P->Q) mapping
    // 3. Instantaneous Volatility

    // Thus each particle has size 5+1+1 = 7

    // 0 - mu
    // 1 - sigma (vol-of-vol)
    // 2 - kappa
    // 3 - theta
    // 4 - rho
    // 5 - lambda
    // 6 - Instantaneous volatility
    map<int, string> index_to_name;

    vector<vector<double>> particles;

    // Next we need a map defining how much each variable can move
    map<int, double> ou_sigma_map;

    // mapping each particle to its max value of likelihood * penalty
    vector<double> max_likelihood_map;

    // Weights for each of the particle
    vector<double> weights;

    vector<double> P_to_Q(const vector<double> &particle){
        vector<double> Q_space_params = particle;
        Q_space_params[2] = particle[2] + particle[5];
        Q_space_params[3] = (particle[2] * particle[3]) / (particle[2] + particle[5] + 1e-8);
        return Q_space_params;
    }

    void moveAhead(vector<double> &particle, double log_return){
        // Heston Parameters need to be done via an OU process

        thread_local mt19937 local_rng{random_device{}()};
        normal_distribution<double> distribution(0.0, 1.0);
        for(int i = 0 ; i < 5 ; i++){
            double dW = distribution(local_rng) * sqrt(dt);
            particle[i] += ou_sigma_map[i] * dW;
        }

        // theta, kappa, sigma > 0
        particle[3] = max(particle[3], 1e-8);
        particle[2] = max(particle[2], 1e-8);
        particle[1] = max(particle[1], 1e-8);

        // ensure rho is within bounds
        particle[4] = clamp(particle[4], -0.999, 0.999);

        // volatility moves via heston discretization
        const double z = distribution(local_rng);
        const double rho = particle[4];

        const double v_prev = max(1e-8, particle[6]);
        const double mu     = particle[0];
        const double sigma  = particle[1];
        const double kappaP = particle[2];
        const double thetaP = particle[3];

        const double eps1 = (log_return - (mu - 0.5 * v_prev) * dt) / sqrt(v_prev * dt);
        const double eps2 = rho * eps1 + sqrt(max(1.0 - rho * rho, 1e-12)) * z;

        // Variance Euler step under 
        const double v_next = v_prev
                            + kappaP * (thetaP - v_prev) * dt
                            + sigma * sqrt(v_prev * dt) * eps2;

        particle[6] = max(1e-8, v_next);
    }

    // a function which normalizes the weights
    void normalizeWeights(){
        double sum = accumulate(weights.begin(), weights.end(), 0.0);

        if(sum <= 0){
            for(int i = 0 ; i < N ; i++){
                weights[i] = 1.0 / N;
            }
            return ;
        }

        for(int i = 0 ; i < N ; i++){
            weights[i] = weights[i] / sum;
        }
    }

    // function to compute ESS
    double ESS(){
        double sum = 0.0;
        for(int i = 0 ; i < N ; i++){
            sum += weights[i] * weights[i];
        }

        if(sum >= 1e-9){
            return 1/sum;
        }

        return 0.0;
    }

    // resample based on updated weights
    vector<int> systematicResample(){
        // booking for ancestral sampling
        vector<int> idx(N, -1);

        vector<vector<double>> newParticles(N, vector<double> (7, 0.0));
        vector<double> cdf(N, 0.0);
        cdf[0] = weights[0];
        for (int i = 1; i < N; i++) cdf[i] = cdf[i - 1] + weights[i];

        uniform_real_distribution<double> u01(0.0, 1.0 / max(1, N));
        double u0 = u01(rng);

        int j = 0;
        for (int i = 0; i < N; i++) {
            double u = u0 + static_cast<double>(i) / static_cast<double>(N);
            while (j < N - 1 && u > cdf[j]) j++;
            newParticles[i] = particles[j];

            // keep track of parent
            idx[i] = j;
        }

        particles.swap(newParticles);
        fill(weights.begin(), weights.end(), 1.0 / max(1, N));

        return idx;
    }

    // function which computes normal returns given particle state
    double likelihood_log_return_given_v(const vector<double>& particle, double observed_log_return){

        const double mu = particle[0];
        const double v_prev = max(1e-8, particle[6]);

        const double mean = (mu - 0.5 * v_prev) * dt;
        const double var  = max(v_prev * dt, 1e-12);

        const double z = observed_log_return - mean;
        const double logLik = -0.5 * (log(2.0 * M_PI * var) + (z * z) / var);

        return logLik;
    }

    // [TODO] - shall we remove the randomness?
    vector<double> create_guess(){
        vector<double> guess(7);

        vector<double> guess_params = data.get_guess();

        thread_local mt19937 local_rng{random_device{}()};
        for(int i = 0 ; i < 7 ; i++){
            uniform_real_distribution<double> distr(0.5 * guess_params[i], 1.5 * guess_params[i]);
            guess[i] = distr(local_rng);
        }

        return guess;
    }

    // initialze this particles using surface calibration
    void initialization(int particle_index, int time_index){
        vector<double> guess = create_guess();
        surfaceCalibration this_particle(data.get_grid(time_index), guess, data.get_date(time_index), data.get_S(time_index), data.get_r(time_index), data.get_q(time_index));
        // surfaceCalibrationCMAES this_particle(data.get_grid(time_index), guess, data.get_date(time_index), data.get_S(time_index), data.get_r(time_index), data.get_q(time_index));
        // surfaceCalibrationMultiStartLM this_particle(data.get_grid(time_index), guess, data.get_date(time_index), data.get_S(time_index), data.get_r(time_index), data.get_q(time_index));
        // surfaceCalibrationMultiStartSimplex this_particle(data.get_grid(time_index), guess, data.get_date(time_index), data.get_S(time_index), data.get_r(time_index), data.get_q(time_index));

        particles[particle_index] = this_particle.getCalibration();
    }

    // function to compute the mean and variance of each parameter
    // [TODO] - I am doing average of variance just like others, is it fine
    vector<pair<double,double>> compute_mean_and_variance(){
        // first is mean, second is variance
        vector<pair<double, double>> mean_variance(7);

        for(int i = 0 ; i < N ; i++){
            for(int parameter = 0 ; parameter < 7 ; parameter++){
                mean_variance[parameter].first += weights[i] * particles[i][parameter];
                mean_variance[parameter].second += weights[i] * particles[i][parameter] * particles[i][parameter];
            }
        }

        for(int parameter = 0 ; parameter < 7 ; parameter++){
            mean_variance[parameter].second -= (mean_variance[parameter].first * mean_variance[parameter].first);
        }

        return mean_variance;
    }

    // Book keeping vectors for ancestral sampling
    vector<vector<int>> ancestor;
    vector<vector<double>> v_value;

    // Do ancestral sampling to get volatility paths
    vector<double> ancestralSampling(){   
        vector<double> sampled_path(time_steps + 1, 0.0);

        discrete_distribution<>dist(weights.begin(),weights.end());
        int k = dist(rng);

        for(int t = time_steps ; t >= 0 ; t--){ 
            sampled_path[t] = v_value[t][k];
            k = ancestor[t][k];
            if((k<0) && (t>0)) break;
        }
        return sampled_path;
    }

    // a boolean which says to use ESS or not
    bool use_ESS;

public:
    // n is the number of particles
    smoothStateCalibration(getRealData &passed_data, int n=10, int n_samples=50, bool to_use_ESS=true, int p=8){
        // Get the true data
        N=n;
        no_of_threads=p;
        number_of_samples = n_samples;
        data = passed_data;
        
        particles.assign(N, vector<double>(7, 0));       // ensure that size of params is N

        // kappa and xi can take larger diffusion than theta.
        ou_sigma_map[0] = 0.12;   // mu
        ou_sigma_map[1] = 0.20;   // sigma / xi (vol-of-vol)
        ou_sigma_map[2] = 0.35;   // kappa
        ou_sigma_map[3] = 0.015;  // theta (long-run variance level) - theta needs much smaller noise in level-space (it is usually around ~0.01 to 0.10).
        ou_sigma_map[4] = 0.08;   // rho - rho should be moderate to avoid constant boundary clipping.
        ou_sigma_map[5] = 0.18;   // lambda

        index_to_name[0] = "mu";
        index_to_name[1] = "sigma";
        index_to_name[2] = "kappa";
        index_to_name[3] = "theta";
        index_to_name[4] = "rho";
        index_to_name[5] = "lambda";
        index_to_name[6] = "volatility";

        time_steps = data.get_time_steps();
        
        weights.assign(N, 1.0/N);

        // initialize matrices required for ancestral sampling
        ancestor.resize(time_steps + 1, vector<int> (N));
        v_value.resize(time_steps + 1, vector<double> (N));

        // Initialize all the particles with values of heston parameters based on surface calibration
        max_likelihood_map.assign(N,-1);
        for(int i = 0 ; i < N ; i++){
            initialization(i, 0);

            // Do the book keeping
            v_value[0][i] = particles[i][6];
            ancestor[0][i] = i;
        }

        // initialize the boolean use_ESS
        use_ESS = to_use_ESS;
    }

    // for a particular time step, we update particle of index i
    void _update_particle(int t, int i, vector<double> &cur_weights, set<int> &reinitialized_indexes, set<int> &surviving_indexes, const double &observed_log_return){
        moveAhead(particles[i], observed_log_return);

        // Compute the log likelihood
        double this_likelihood = likelihood_log_return_given_v(particles[i], observed_log_return);
        weights[i] *= exp(this_likelihood);
        cur_weights[i] *= exp(this_likelihood);
        
        // Then we update weights
        // Add penatly of real - computed, or could make normal error, either way will work
        
        const vector<double> qParams = P_to_Q(particles[i]);
        
        // we pass the parameters and compute the average error, i.e. how well it fits the current surface
        double penalty = data.get_penalty(t, particles[i][6], qParams[2], qParams[3], particles[i][1], particles[i][4]);

        weights[i] *= penalty;
        cur_weights[i] *= penalty;
        
        // if any of the values goes less than 50% of max value, then re-initialize
        if(max_likelihood_map[i] == -1){
            max_likelihood_map[i] = cur_weights[i];
            set_lock_1.lock();
            surviving_indexes.insert(i);
            set_lock_1.unlock();
        }else{
            max_likelihood_map[i] = max(max_likelihood_map[i], cur_weights[i]);
            if(cur_weights[i] < 0.5 * max_likelihood_map[i]){
                initialization(i, t);
                max_likelihood_map[i] = -1;
                ctr_reinit.fetch_add(1);

                // Check how to update weightsf

                // option 1 - Carry forward the old weights (original code)

                // Option 2 - Set weights to 1/N (but does not make too much sense)
                // weights[i] = 1.0/N

                // Option 3 - Set weights to average before normalizing them
                // weights[i] = accumulate(weights.begin(), weights.end(), 0.0) / N;

                // Option 4 - Rescore from penalty and weights and use some alpha term for weighting
                set_lock_2.lock();
                reinitialized_indexes.insert(i);
                set_lock_2.unlock();
            }else{
                set_lock_1.lock();
                surviving_indexes.insert(i);
                set_lock_1.unlock();
            }
        }
    }

    void main(){

        // Initialization is done by the constructor 
        
        ofstream error_file("Errors.csv");
        error_file << "date,strike,maturity,true_price,computed_price,abs_error" << '\n';
        ofstream log_file("smooth_calibration.log");
        ofstream statistical_file("statistics.log");
        ofstream statistical_csv("parameters.csv");
        statistical_csv << "mu,vol-of-vol,kappa,theta,rho,lambda,v_t" << '\n';
        for(int t = 1 ; t < time_steps ; t++){
            log_file << "Time: " << t << '\n';
    
            // We now have updated values for each particle
            // book keeping for each particle
            vector<double> cur_weights(N, 1.0);

            // sets which keep track of all the indices which need to be reinitialized
            set<int> reinitialized_indexes, surviving_indexes;

            double observed_log_return = data.get_log_return(t);

            // Next we need each particle to move to time step t+1
            // This for loop can be parallelized

            // We have a lambda function called by threads
            auto worker = [&]() {
                while(true){
                    int i = next_particle.fetch_add(1);
                    if(i>=N) break;

                    _update_particle(t, i, cur_weights, reinitialized_indexes, surviving_indexes, observed_log_return);
                }
            };

            // We create threads
            vector<thread> threads;

            for(int k = 0 ; k < no_of_threads ; k++){
                threads.emplace_back(worker);
            }

            for(int k = 0 ; k < no_of_threads ; k++){
                threads[k].join();
            }

            next_particle.store(0);

            // The weights only need to be modified if atleast one of the indexes are re-initalized 
            if(reinitialized_indexes.size() != 0){
                // assign alpha in [0, 1]
                double alpha = 0.05;
                for(int i : reinitialized_indexes){
                    double likelihood_for_i = likelihood_log_return_given_v(particles[i], observed_log_return);

                    const vector<double> qParams = P_to_Q(particles[i]);

                    double penalty_for_i = data.get_penalty(t, particles[i][6], qParams[2], qParams[3], particles[i][1], particles[i][4]);

                    weights[i] = exp(likelihood_for_i) * penalty_for_i * alpha / reinitialized_indexes.size();
                }

                for(int i : surviving_indexes){
                    weights[i] *= (1 - alpha);
                }
            }

            // normalize the weights so they sum up to 1
            normalizeWeights();

            vector<pair<double, double>> mean_variance = compute_mean_and_variance();
            statistical_file << "Time " << t << '\n';
            for(int i = 0 ; i < 7 ; i++){
                statistical_file << index_to_name[i] << ": mean=" << mean_variance[i].first << ", variance=" << mean_variance[i].second << '\n';
                statistical_csv << mean_variance[i].first << ',';
            }
            statistical_csv << '\n';
            statistical_file << "-----------------------------------------------\n";

            for(int i = 0 ; i < N ; i++){
                log_file << "Particle " << i << ":";
                for(int j = 0 ; j < 7 ; j++){
                    log_file << particles[i][j] << ' ';
                } 
                log_file << " - w = " << weights[i];
                log_file << endl;
            }

            // errors to print
            data.get_penalty(
                t, 
                mean_variance[6].first, 
                mean_variance[2].first + mean_variance[5].first, 
                (mean_variance[2].first * mean_variance[3].first)/(mean_variance[2].first + mean_variance[5].first + 1e-8), 
                mean_variance[1].first, 
                mean_variance[4].first, 
                error_file
            );

            // Then we resample

            if(use_ESS){
                // resample based on ESS
                if(ESS() < 0.5 * N){
                    vector<int> idx = systematicResample();
                    for(int j = 0 ; j < N ; j++){
                        v_value[t][j] = particles[j][6];
                        ancestor[t][j] = idx[j];
                        weights[j] = 1.0/N;
                    }
                }else{
                    for(int j = 0 ; j < N ; j++){
                        v_value[t][j] = particles[j][6];
                        ancestor[t][j] = j;
                    }
                }
            }else{
                vector<int> idx = systematicResample();
                for(int j = 0 ; j < N ; j++){
                    v_value[t][j] = particles[j][6];
                    ancestor[t][j] = idx[j];
                    weights[j] = 1.0/N;
                }
            }
        }

        log_file.close();
        statistical_file.close();
        error_file.close();
        statistical_csv.close();

        // At the end just get the final value of the parameters
        auto final_mean_and_variance = compute_mean_and_variance();
        ofstream final_mean("calibrated_params.csv");

        int ctr = 0;
        for(auto [u, v] : final_mean_and_variance){
            ctr++;
            final_mean << u;
            if(ctr != 7){
                final_mean << ',';
            }
        }
        final_mean << '\n';

        ctr=0;
        for(auto [u, v] : final_mean_and_variance){
            ctr++;
            final_mean << v;
            if(ctr != 6){
                final_mean << ',';
            }
        }
        final_mean.close();

        // file to write sampled paths
        ofstream sampled_paths("sampled_paths.csv");

        // Call ancestral sampling
        for(int i = 0 ; i < number_of_samples ; i++){
            vector<double> sample = ancestralSampling();
            int sz = sample.size();
            for(int j = 0 ; j < sz ; j++){
                sampled_paths << sample[j];
                if(j!=sz-1){
                    sampled_paths << ',';
                }
            }
            if(i!=number_of_samples-1){
                sampled_paths << '\n';
            }
        }

        sampled_paths.close();

        cout << "The number of reinitializations are : " << ctr_reinit << '\n';
    }
};