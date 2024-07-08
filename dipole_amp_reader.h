
#include <iostream>
#include <vector>
#include <rapidcsv.h>

using namespace std;

const double x_0 = 0.01;

//const int b_res = 30*30;
//static array<array<array<array<double, 4>, 81>, b_res>, 30> table;

double calc_b(double r, double b_min, double phi) {
    return sqrt(pow(r/2, 2) + b_min*b_min - r*b_min*cos(phi));
}

double calc_x(double Y) {
    return exp(-Y)*x_0;
}

double get_dipole_amplitude(array<array<array<array<double, 4>, 81>, 900>, 30> &table, double r, double b, double x) {
    int i = 0;
    while (table[i][0][0][0] < r) {
        i++;
    }
    if (i != 0) {
        if (abs(table[i][0][0][0] - r) > r - table[i-1][0][0][0]) { //absolute value in case r is very large
            i = i - 1;
        }
    }

    int j = 0;
    while (table[i][j][0][1] < b) {
        j++;
    }
    if (j != 0) {
        int lower_j = j - 1;
        while (table[i][j][0][1] == table[i][lower_j][0][1]) {
            lower_j--;
        }
        if (abs(table[i][j][0][1] - b) > b - table[i][lower_j][0][1]) { //absolute value in case b is very large
            j = lower_j;
        }
    }

    int k = 0;
    while (table[i][j][k][2] < x) {
        k++;
    }
    if (k != 0) {
        if (abs(table[i][j][k][2] - x) > x - table[i][j][k-1][2]) { //absolute value in case x is very large
            k = k - 1;
        }
    }

    return table[i][j][k][3];
}

void load_dipole_amplitudes(array<array<array<array<double, 4>, 81>, 900>, 30> &table, string filename) {
    cout << "Reading " << filename << endl;
    rapidcsv::Document doc(filename);
    
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
                //cout << index << endl;
                double b = calc_b(r[index], b_min[index], phi[index]);
                double x = calc_x(Y[index]);
                table[i][j][k][0] = r[index];
                table[i][j][k][1] = b;
                table[i][j][k][2] = x;
                table[i][j][k][3] = N[index];
            }
        }
    }

    //cout << endl << endl;

    for (int j=0; j<30; j++) {
        int i = 1;
        bool ordered = true;
        while (true) {
            //cout << i << ", " << table[j][i][0][1] << endl;
            if (table[j][i][0][1] < table[j][i-1][0][1]) {
                table[j][i].swap(table[j][i-1]);
                ordered = false;
                if (i > 1) {
                    i -= 1;
                    continue;
                }
            }
            i++;
            if (i == table[j].size()) {
                if (ordered == true) {
                    break;
                } else {
                    i = 1;
                    ordered = true;
                }
            }
        }
    }

    cout << "b sorting finished" << endl;
    
    for (int i=0; i<30; i++) {
        for (int j=0; j<30*30; j++) {
            int k = 1;
            bool ordered = true;
            while (true) {
                //cout << k << ", " << table[i][j][k][2] << endl;
                if (table[i][j][k][2] < table[i][j][k-1][2]) {
                    table[i][j][k].swap(table[i][j][k-1]);
                    ordered = false;
                    if (k > 1) {
                        k -= 1;
                        continue;
                    }
                }
                k++;
                if (k == table[i][j].size()) {
                    if (ordered == true) {
                        break;
                    } else {
                        k = 1;
                        ordered = true;
                    }
                }
            }
        }
    }

    cout << "x sorting finished, table ready" << endl;

}
