
#include <iostream>
#include <vector>
#include <rapidcsv.h>
#include <math.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"

using namespace std;

const double x_0 = 0.01;

double calc_max_phi(double r, double b_min) {
    if (r > 2*b_min) {
        return 2*M_PI;
    } else {
        return M_PI - acos(r/(2*b_min));
    }
}

double calc_x(double Y) {
    return exp(-Y)*x_0;
}

double get_dipole_amplitude(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, double r, double b, double phi, double x) {
    //cout << "Getting amplitude, r=" << r << ", b=" << b << ", x=" << x << endl;

    if (phi > calc_max_phi(r, b)) {
        return 0;
    }

    long unsigned int i = 0;
    int i_closer = 0;
    bool exact_r = false;
    if (r <= table[0][0][0][0][0]) {
        i = 0;
        i_closer = i;
        exact_r = true;
    } else if (r >= table[table.size()-1][0][0][0][0]) {
        i = table.size()-1;
        i_closer = i;
        exact_r = true;
    } else {
        while (table[i][0][0][0][0] < r) {
            i++;
            if (i > table.size()-1) {
                throw 1;
            }
        }
        if (table[i][0][0][0][0] == r) {
            i_closer = i;
            exact_r = true;
        }else if (table[i][0][0][0][0] - r < r - table[i-1][0][0][0][0]) {
            i_closer = i;
        } else {
            i_closer = i-1;
        }
    }

    long unsigned int j = 0;
    int j_closer = 0;
    bool exact_b = false;
    if (b <= table[i_closer][0][0][0][1]) {
        j = 0;
        j_closer = j;
        exact_b = true;
    } else if (b >= table[i_closer][table[i_closer].size()-1][0][0][1]) {
        j = table[i_closer].size()-1;
        j_closer = j;
        exact_b = true;
    } else {
        while (table[i_closer][j][0][0][1] < b) {
            j++;
            if (j > table[i_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i_closer][j][0][0][1] == b) {
            j_closer = j;
            exact_b = true;
        }else if (table[i_closer][j][0][0][1] - b < b - table[i_closer][j-1][0][0][1]) {
            j_closer = j;
        } else {
            j_closer = j-1;
        }
    }

    long unsigned int k = 0;
    int k_closer = 0;
    bool exact_phi = false;
    if (phi <= table[i_closer][j_closer][0][0][2]) {
        k = 0;
        k_closer = k;
        exact_phi = true;
    } else if (phi >= table[i_closer][j_closer][table[i_closer][j_closer].size()-1][0][2]) {
        k = table[i_closer][j_closer].size()-1;
        k_closer = k;
        exact_phi = true;
    } else {
        while (table[i_closer][j_closer][k][0][2] < phi) {
            k++;
            if (k > table[i_closer][j_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i_closer][j_closer][k][0][2] == phi) {
            k_closer = k;
            exact_phi = true;
        }else if (table[i_closer][j_closer][k][0][2] - phi < phi - table[i_closer][j_closer][k-1][0][2]) {
            k_closer = k;
        } else {
            k_closer = k-1;
        }
    }

    long unsigned int l = 0;
    int l_closer = 0;
    bool exact_x = false;
    if (x <= table[i_closer][j_closer][k_closer][0][3]) {
        l = 0;
        l_closer = l;
        exact_x = true;
    } else if (x >= table[i_closer][j_closer][k_closer][table[i_closer][j_closer][k_closer].size()-1][3]) {
        l = table[i_closer][j_closer][k_closer].size()-1;
        l_closer = l;
        exact_x = true;
    } else {
        while (table[i_closer][j_closer][k_closer][l][3] < x) {
            l++;
            if (l > table[i_closer][j_closer][k_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i_closer][j_closer][k_closer][l][3] == x) {
            l_closer = l;
            exact_x = true;
        }else if (table[i_closer][j_closer][k_closer][l][3] - x < x - table[i_closer][j_closer][k_closer][l-1][3]) {
            l_closer = l;
        } else {
            l_closer = l-1;
        }
    }

    double r_corr, b_corr, phi_corr, x_corr;

    if (exact_r) {
        r_corr = 0;
    } else {
        if (i_closer == i) {
            r_corr = -1*(table[i][j_closer][k_closer][l_closer][4] - table[i-1][j_closer][k_closer][l_closer][4])*(1-(r - table[i-1][0][0][0][0])/(table[i][0][0][0][0] - table[i-1][0][0][0][0]));
        } else {
            r_corr = (r - table[i-1][0][0][0][0])/(table[i][0][0][0][0] - table[i-1][0][0][0][0])*(table[i][j_closer][k_closer][l_closer][4] - table[i-1][j_closer][k_closer][l_closer][4]);
        }
        //cout << endl;
        //cout << "r=" << r << ", low=" << table[i-1][0][0][0] << ", high=" << table[i][0][0][0] << ", closer=" << r_closer << endl;
        //cout << "location=" << (r - table[i-1][0][0][0])/(table[i][0][0][0] - table[i-1][0][0][0]) << ", higher_N=" << table[i][j+b_closer][k+x_closer][3] << ", lower_N=" << table[i-1][j+b_closer][k+x_closer][3] << ", corr=" << r_corr << endl;
        //cout << endl;
    }

    if (exact_b) {
        b_corr = 0;
    } else {
        if (j_closer == j) {
            b_corr = -1*(table[i_closer][j][k_closer][l_closer][4] - table[i_closer][j-1][k_closer][l_closer][4])*(1-(b - table[i_closer][j-1][0][0][1])/(table[i_closer][j][0][0][1] - table[i_closer][j-1][0][0][1]));
        } else {
            b_corr = (b - table[i_closer][j-1][0][0][1])/(table[i_closer][j][0][0][1] - table[i_closer][j-1][0][0][1])*(table[i_closer][j][k_closer][l_closer][4] - table[i_closer][j-1][k_closer][l_closer][4]);
        }
    }

    if (exact_phi) {
        phi_corr = 0;
    } else {
        if (k_closer == k) {
            phi_corr = -1*(table[i_closer][j_closer][k][l_closer][4] - table[i_closer][j_closer][k-1][l_closer][4])*(1-(phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2]));
        } else {
            phi_corr = (phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2])*(table[i_closer][j_closer][k][l_closer][4] - table[i_closer][j_closer][k-1][l_closer][4]);
        }
    }

    if (exact_x) {
        x_corr = 0;
    } else {
        if (l_closer == l) {
            x_corr = -1*(table[i_closer][j_closer][k_closer][l][4] - table[i_closer][j_closer][k_closer][l-1][4])*(1-(x - table[i_closer][j_closer][k_closer][l-1][3])/(table[i_closer][j_closer][k_closer][l][3] - table[i_closer][j_closer][k_closer][l-1][3]));
        } else {
            x_corr = (x - table[i_closer][j_closer][k_closer][l-1][3])/(table[i_closer][j_closer][k_closer][l][3] - table[i_closer][j_closer][k_closer][l-1][3])*(table[i_closer][j_closer][k_closer][l][4] - table[i_closer][j_closer][k_closer][l-1][4]);
        }
    }
    
    //double return_value = table[i+r_closer][j+b_closer][k+x_closer][3] + r_corr + b_corr + x_corr;
    //cout << r_corr << ", " << b_corr << ", " << x_corr << endl;
    //cout << "return value: " << table[i+r_closer][j+b_closer][k+x_closer][3] << ", " << table[i+r_closer-1][j+b_closer][k+x_closer][3] << ", " << return_value << endl << endl;
    //return return_value;

    return table[i_closer][j_closer][k_closer][l_closer][4] + r_corr + b_corr + phi_corr + x_corr;
}

void load_dipole_amplitudes(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, string filename, double plot=false) {
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
                for (int k=0; k<30; k++) {
                    for (int l=0; l<81; l++) {
                    int index = i*icof + j*jcof + k*81 + l;
                    table[i][j][k][l][0] = r[index];
                    table[i][j][k][l][1] = b_min[index];
                    table[i][j][k][l][2] = phi[index];
                    table[i][j][k][l][3] = calc_x(Y[index]);
                    table[i][j][k][l][4] = N[index];
                }
            }
        }
    }
    cout << "Table ready" << endl;
    return;
}