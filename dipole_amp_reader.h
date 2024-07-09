
#include <iostream>
#include <vector>
#include <rapidcsv.h>
#include <math.h>

using namespace std;

const double x_0 = 0.01;

//const int b_res = 30*30;
//static array<array<array<array<double, 4>, 81>, b_res>, 30> table;

double calc_b(double r, double b_min, double phi) {
    return sqrt(pow(r/2, 2) + b_min*b_min - r*b_min*cos(M_PI - phi));
}

double calc_x(double Y) {
    return exp(-Y)*x_0;
}

double calc_max_phi(double r, double b_min) {
    return M_PI - acos(r/(2*b_min));
}

double get_dipole_amplitude(array<array<array<array<double, 4>, 81>, 30>, 30> &table, double r, double b, double x) {
    //cout << "Getting amplitude, r=" << r << ", b=" << b << ", x=" << x << endl;
    long unsigned int i = 0;
    while (table[i][0][0][0] < r) {
        i++;
        if (i == table.size() - 1) {
            break;
        }
    }
    if (i != 0) {
        if (abs(table[i][0][0][0] - r) > r - table[i-1][0][0][0]) { //absolute value in case r is very large
            i = i - 1;
        }
    }

    long unsigned int j = 0;
    while (table[i][j][0][1] < b) {
        j++;
        if (j == table[i].size() - 1) {
            break;
        }
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

    long unsigned int k = 0;
    while (table[i][j][k][2] < x) {
        k++;
        if (k == table[i][j].size() - 1) {
            break;
        }
    }
    if (k != 0) {
        if (abs(table[i][j][k][2] - x) > x - table[i][j][k-1][2]) { //absolute value in case x is very large
            k = k - 1;
        }
    }

    return table[i][j][k][3];
}

void load_dipole_amplitudes(array<array<array<array<double, 4>, 81>, 30>, 30> &table, string filename) {
    cout << "Reading " << filename << endl;
    rapidcsv::Document doc(filename);
    
    vector<double> r = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");

    cout << "Data read successfully" << endl;
    
    const int icof = 30*30*81;
    const int jcof = 30*81;

    for (int i=0; i<30; i++) {
        for (int j=0; j<30; j++) {
            for (int k=0; k<81; k++) {
                double N_mean = 0;
                double b_mean = 0;
                int mean_denominator = 0;
                int sub_index = i*icof + j*jcof;
                double max_phi = calc_max_phi(r[sub_index], b_min[sub_index]);
                for (int l=0; l<30; l++) {
                    int index = sub_index + l*81 + k;
                    if (phi[index] > max_phi) {
                        break;
                    }
                    b_mean += calc_b(r[index], b_min[index], phi[index]); // factor 2 not needed since the end result is an average, i.e. Nb also doesn't have factor 2
                    N_mean += N[index];
                    mean_denominator++;
                }
                N_mean = N_mean/mean_denominator;
                b_mean = b_mean/mean_denominator;
                double x = calc_x(Y[sub_index + k]); //value of l is set to zero here
                table[i][j][k][0] = r[sub_index];
                table[i][j][k][1] = b_mean;
                table[i][j][k][2] = x;
                table[i][j][k][3] = N_mean;
            }
        }
    }

    //cout << endl << endl;

    for (int j=0; j<30; j++) {
        long unsigned int i = 1;
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
            long unsigned int k = 1;
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
