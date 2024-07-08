
#include <iostream>
#include <vector>
#include <rapidcsv.h>

using namespace std;

const double x_0 = 0.01;

const int b_res = 30*30;
static array<array<array<array<double, 4>, 81>, b_res>, 30> table;

double calc_b(double r, double b_min, double phi) {
    return sqrt(pow(r/2, 2) + b_min*b_min - r*b_min*cos(phi));
}

double calc_x(double Y) {
    return exp(-Y)*x_0;
}

int main() {
    cout << "start" << endl;
    rapidcsv::Document doc("data/dipole_amplitude_with_IP_dependence.csv");
    
    vector<double> r = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");
    
    const int icof = 30*30*81;

    for (int i=0; i<30; i++) {
        for (int j=0; j<30*30; j++) {
            for (int k=0; k<81; k++) {
                int index = i*icof + j*81 + k;
                cout << index << endl;
                double b = calc_b(r[index], b_min[index], phi[index]);
                double x = calc_x(Y[index]);
                table[i][j][k][0] = r[index];
                table[i][j][k][1] = b;
                table[i][j][k][2] = x;
                table[i][j][k][3] = N[index];
            }
        }
    }

    cout << endl << endl;

    int i = 1;
    bool ordered = true;
    while (true) {
        cout << i << endl;
        if (table[0][i][0][1] < table[0][i-1][0][1]) {
            table[0][i].swap(table[0][i-1]);
            ordered = false;
            if (i > 1) {
                i -= 1;
                continue;
            }
        }
        i++;
        if (i == table[0][i].size()) {
            if (ordered == true) {
                break;
            } else {
                i = 1;
                ordered = true;
            }
        }
    }
    cout << "sorting finished" << endl;
    for (int i=0; i<1; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<5; k++) {
                int index = i*icof + j*81 + k;
                cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
            }
        }
    }

    
}