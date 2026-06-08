Here We do joint calibration of P and Q space as mentioned in the paper

## Before Running the code
Ensure you have downloaded QuantLib, Eigen, Armadillo (optional), LBFGSPP.

## Repository Structure
1. app
    - main.cpp: this will be the file which will be run when terminal command given, calling the and linking the other files
2. src
    - callPrice.cpp: Given the spot price, strike price, time to maturity, risk free rate, dividents price, heston parameters, and date of evaluation, it will use QuantLib Library and return the call option.
    - generateData.cpp: Based on the initialization of parameters in the constructor (listed in the later section), it will generate 
        - The stock price path for time t=1...T (where T is the total number of time steps this program will run for).
        - The volatility path for time t=1...T (where T is the total number of time steps this program will run for).
        - the call price grid for each time from t=1...T (where T is the total number of time steps this program will run for).
        - The true value of call prices which will be tested (we are only testing how value of one call option with fixed maturity and strike changes over time)
    - surfaceCalibration.cpp: Given the call option grid for a particular day, the spot price, the risk-free rate, and price of divident, it will return the heston parameters in P-space, and the volatlility risk using surface calibration
    - smoothStateCalibration.cpp: The core function of our paper, this is responsible for the calibration and create the following log/csv files which can later be used for visualization
        - smooth_calibration.log: Shows the current state of all particles at each time
        - statistics.log: Displays the weighted mean and variance of all the parameters at each time step.
        - calibrated_params.csv: Will store the final value of mean and variance after the program has run for all the time steps in the order: mu, xi, kappa, theta, rho, lambda
        - S_path.csv: The generated stock price path
        - v_path.csv: The generated latent volatility path
        - true_params.csv: The parameters used to generate all the data in the order - mu, xi, kappa_P, theta_P, rho, lambda, v0
        - sampled_paths.csv: each row represents one sample (i.e. entire latent volatility path)
        
3. result_visualization.ipynb
    - Used for showing the results and plotting the latent paths

## Parameter Modification

These values can be modified as per need.

1. In src/generateData.cpp file, the constructor has assigned the values of the following parameters (please modify as per need):
    - The number of time steps to evalue this program for
    - The spot price at time t=0
    - The true value of initial volatility v0
    - The true value heston parameter mu
    - The true value heston parameter kappa in the P-space
    - The true value heston parameter theta in the P-space
    - The true value heston parameter xi (vol-of-vol)
    - The true value heston parameter rho
    - The value of risk-free-rate (r)
    - The value of dividends (q)
    - The true value of volatility risk parameter (lambda)
    - The values of strike prices to make the call grid
    - The values of maturity time (in years) to make the call grid
    - The strike price for the data to be evaluated at
    - The date of maturity for the data being evaluated (YYYY-MM-DD format)
    - The date of today (or evaluation day of choice), when the evaluation is being started (YYYY-MM-DD format)

2. In app/main.cpp file (Line Number 15):
    - The second parameter is the number of particles of particle filter
    - The third parameter is the number of samples to be sampled of latent volatility path
    - The fourth parameter is a boolean, indicating weather use ESS for resampling on not
3. In src/smoothStateCalibration.cpp:
    - Line 399, the value of alpha as explained in the paper (for modifying weights of particles in case of re-initialization)
    - Line Number 450, the ESS condition for resampling

## Deployment
To run the file:
```bash
g++ app/main.cpp -std=c++17 -I$HOME/local/include -L$HOME/local/lib -lQuantLib -o calibrate
```

Followed by
```bash
./calibrate
```