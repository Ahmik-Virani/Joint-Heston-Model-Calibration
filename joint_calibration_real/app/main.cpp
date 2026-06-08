#include<iostream>
#include<fstream>
#include<string>
#include<vector>

#include "../src/smoothStateCalibration.cpp"
#include "../src/getRealData.cpp"

using namespace std;

int main(){
    // First we generate the data
    getRealData data;
    cout<<"Data Generation Done."<<endl;
    // Then we need to do calibration
    // We pass data, number of particles, number of volatility path samples, use ESS, number of processes/threads
    // [TODO] - pass alpha from here
    smoothStateCalibration cali(data, 200, 2, true, 1);
    cali.main();

    return 0;
}