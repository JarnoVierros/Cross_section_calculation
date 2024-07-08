
#include <iostream>
#include <vector>
#include <rapidcsv.h>

using namespace std;

const double x_0 = 0.01;

double calc_b(double r, double b_min, double phi) {
    return sqrt(pow(r/2, 2) + b_min*b_min - r*b_min*cos(phi));
}

double calc_x(double Y) {
    return exp(-Y)*x_0;
}

int main() {
    rapidcsv::Document doc("data/dipole_amplitude_with_IP_dependence.csv");

    vector<double> r = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");

    double fr = r[0];
    double fb = calc_b(r[0], b_min[0], phi[0]);
    double fx = calc_x(Y[0]);
    double fN = N[0];

    for (int i=1; i<r.size(); i++) {
        if (r[i]=fr&&calc_b(r[i], b_min[i], phi[i])==fb&&calc_x(Y[i])) {
            cout << "r: " << r[i] << endl;
            cout << "b_min: " << b_min[i] << endl;
            cout << "phi: " << phi[i] << endl;
            cout << "Y: " << Y[i] << endl;
            cout << "N: " << N[i] << endl;
            cout << "fN: " << fN << endl;
            cout << endl;
        }
    }


}