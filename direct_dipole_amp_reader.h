
#include <iostream>
#include <vector>
#include <rapidcsv.h>
#include <math.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"

#include <ctime>
#include <linterp.h>

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

double fn (double x1, double x2) { return sin(x1 + x2); }

void create_p_interpolator(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, InterpMultilinear<4, double>* &return_interp_ML) {
    std::vector<double> r_vec, b_min_vec, phi_vec, x_vec;

    array<int, 4> grid_sizes;
    grid_sizes[0] = 30;
    grid_sizes[1] = 30;
    grid_sizes[2] = 30;
    grid_sizes[3] = 81;//81
    
    int num_elements = grid_sizes[0]*grid_sizes[1]*grid_sizes[2]*grid_sizes[3];

    std::vector<double> N_vec(num_elements);

    for (int i=0; i<grid_sizes[0]; i++) {
        r_vec.push_back(log(table[i][0][0][0][0]));
    }

    for (int i=0; i<grid_sizes[1]; i++) {
        b_min_vec.push_back(log(table[0][i][0][0][1]));
    }

    for (int i=0; i<grid_sizes[2]; i++) {
        phi_vec.push_back(table[0][0][i][0][2]);
    }

    for (int i=0; i<grid_sizes[3]; i++) {
        x_vec.push_back(log(table[0][0][0][i][3]));
    }

    std::vector<std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(r_vec.begin());
    grid_iter_list.push_back(b_min_vec.begin());
    grid_iter_list.push_back(phi_vec.begin());
    grid_iter_list.push_back(x_vec.begin());

    cout << endl;

    for (int i=0; i<grid_sizes[0]; i++) {
        for (int j=0; j<grid_sizes[1]; j++) {
            for (int k=0; k<grid_sizes[2]; k++) {
                for (int l=0; l<grid_sizes[3]; l++) {
                    N_vec[i*grid_sizes[1]*grid_sizes[2]*grid_sizes[3] + j*grid_sizes[2]*grid_sizes[3] + k*grid_sizes[3] + l] = log(table[i][j][k][l][4]);
                }
            }
        }
    }

    return_interp_ML = new InterpMultilinear<4, double>(grid_iter_list.begin(), grid_sizes.begin(), N_vec.data(), N_vec.data() + num_elements);
}

void create_Pb_interpolator(array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> &table, InterpMultilinear<4, double>* &return_interp_ML) {
    std::vector<double> r_vec, b_min_vec, phi_vec, x_vec;

    array<int, 4> grid_sizes;
    grid_sizes[0] = 40;
    grid_sizes[1] = 40;
    grid_sizes[2] = 40;
    grid_sizes[3] = 81;
    
    int num_elements = grid_sizes[0]*grid_sizes[1]*grid_sizes[2]*grid_sizes[3];

    std::vector<double> N_vec(num_elements);

    for (int i=0; i<grid_sizes[0]; i++) {
        r_vec.push_back(log(table[i][0][0][0][0]));
    }

    for (int i=0; i<grid_sizes[1]; i++) {
        b_min_vec.push_back(log(table[0][i][0][0][1]));
    }

    for (int i=0; i<grid_sizes[2]; i++) {
        phi_vec.push_back(table[0][0][i][0][2]);
    }

    for (int i=0; i<grid_sizes[3]; i++) {
        x_vec.push_back(log(table[0][0][0][i][3]));
    }

    std::vector<std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(r_vec.begin());
    grid_iter_list.push_back(b_min_vec.begin());
    grid_iter_list.push_back(phi_vec.begin());
    grid_iter_list.push_back(x_vec.begin());

    cout << endl;

    for (int i=0; i<grid_sizes[0]; i++) {
        for (int j=0; j<grid_sizes[1]; j++) {
            for (int k=0; k<grid_sizes[2]; k++) {
                for (int l=0; l<grid_sizes[3]; l++) {
                    N_vec[i*grid_sizes[1]*grid_sizes[2]*grid_sizes[3] + j*grid_sizes[2]*grid_sizes[3] + k*grid_sizes[3] + l] = log(table[i][j][k][l][4]);
                }
            }
        }
    }

    return_interp_ML = new InterpMultilinear<4, double>(grid_iter_list.begin(), grid_sizes.begin(), N_vec.data(), N_vec.data() + num_elements);
}

double get_p_dipole_amplitude(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, double r, double b, double phi, double x, bool limit_phi=true) {
    //cout << "Getting amplitude, r=" << r << ", b=" << b << ", phi=" << phi << ", x=" << x << endl;

    if (limit_phi) {
        if (phi > calc_max_phi(r, b)) {
            return 0;
        }
    }


    long unsigned int i = 0;
    long unsigned int i_closer = 0;
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
        }else if (log(table[i][0][0][0][0]) - log(r) < log(r) - log(table[i-1][0][0][0][0])) {
            i_closer = i;
        } else {
            i_closer = i-1;
        }
    }

    long unsigned int j = 0;
    long unsigned int j_closer = 0;
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
        }else if (log(table[i_closer][j][0][0][1]) - log(b) < log(b) - log(table[i_closer][j-1][0][0][1])) {
            j_closer = j;
        } else {
            j_closer = j-1;
        }
    }

    long unsigned int k = 0;
    long unsigned int k_closer = 0;
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
    long unsigned int l_closer = 0;
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
            //cout << table[i_closer][j_closer][k_closer][l][3] << endl;
            l++;
            if (l > table[i_closer][j_closer][k_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i_closer][j_closer][k_closer][l][3] == x) {
            l_closer = l;
            exact_x = true;
        }else if (log(table[i_closer][j_closer][k_closer][l][3]) - log(x) < log(x) - log(table[i_closer][j_closer][k_closer][l-1][3])) {
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
            r_corr = -1*(log(table[i][j_closer][k_closer][l_closer][4]) - log(table[i-1][j_closer][k_closer][l_closer][4]))*(1-(log(r) - log(table[i-1][0][0][0][0]))/(log(table[i][0][0][0][0]) - log(table[i-1][0][0][0][0])));
        } else {
            r_corr = (log(r) - log(table[i-1][0][0][0][0]))/(log(table[i][0][0][0][0]) - log(table[i-1][0][0][0][0]))*(log(table[i][j_closer][k_closer][l_closer][4]) - log(table[i-1][j_closer][k_closer][l_closer][4]));
        }
        /*
        cout << endl;
        cout << "r=" << r << ", low=" << table[i-1][0][0][0][0] << ", high=" << table[i][0][0][0][0] << ", closer=" << i-i_closer << endl;
        cout << "location=" << (r - table[i-1][0][0][0][0])/(table[i][0][0][0][0] - table[i-1][0][0][0][0]) << ", higher_N=" << table[i][j_closer][k_closer][l_closer][4] << ", lower_N=" << table[i-1][j_closer][k_closer][l_closer][4] << ", corr=" << r_corr << endl;
        cout << endl;
        */
    }

    if (exact_b) {
        b_corr = 0;
    } else {
        if (j_closer == j) {
            b_corr = -1*(log(table[i_closer][j][k_closer][l_closer][4]) - log(table[i_closer][j-1][k_closer][l_closer][4]))*(1-(log(b) - log(table[i_closer][j-1][0][0][1]))/(log(table[i_closer][j][0][0][1]) - log(table[i_closer][j-1][0][0][1])));
        } else {
            b_corr = (log(b) - log(table[i_closer][j-1][0][0][1]))/(log(table[i_closer][j][0][0][1]) - log(table[i_closer][j-1][0][0][1]))*(log(table[i_closer][j][k_closer][l_closer][4]) - log(table[i_closer][j-1][k_closer][l_closer][4]));
        }
        /*
        cout << endl;
        cout << "j=" << j << endl;
        cout << "j_closer=" << j_closer << endl;
        cout << "i_closer=" << i_closer << endl;
        cout << "b=" << b << ", low=" << table[i_closer][j-1][0][0][1] << ", high=" << table[i_closer][j][0][0][1] << ", closer=" << j-j_closer << endl;
        cout << "location=" << (b - table[i_closer][j-1][0][0][1])/(table[i_closer][j][0][0][1] - table[i_closer][j-1][0][0][1]) << ", higher_N=" << table[i_closer][j][k_closer][l_closer][4] << ", lower_N=" << table[i_closer][j-1][k_closer][l_closer][4] << ", corr=" << b_corr << endl;
        cout << endl;
        */
    }

    if (exact_phi) {
        phi_corr = 0;
    } else {
        if (k_closer == k) {
            phi_corr = -1*(log(table[i_closer][j_closer][k][l_closer][4]) - log(table[i_closer][j_closer][k-1][l_closer][4]))*(1-(phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2]));
        } else {
            phi_corr = (phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2])*(log(table[i_closer][j_closer][k][l_closer][4]) - log(table[i_closer][j_closer][k-1][l_closer][4]));
        }
    }

    if (exact_x) {
        x_corr = 0;
        //cout << "exact" << endl;
    } else {
        if (l_closer == l) {
            x_corr = -1*(log(table[i_closer][j_closer][k_closer][l][4]) - log(table[i_closer][j_closer][k_closer][l-1][4]))*(1-(log(x) - log(table[i_closer][j_closer][k_closer][l-1][3]))/(log(table[i_closer][j_closer][k_closer][l][3]) - log(table[i_closer][j_closer][k_closer][l-1][3])));
        } else {
            x_corr = (log(x) - log(table[i_closer][j_closer][k_closer][l-1][3]))/(log(table[i_closer][j_closer][k_closer][l][3]) - log(table[i_closer][j_closer][k_closer][l-1][3]))*(log(table[i_closer][j_closer][k_closer][l][4]) - log(table[i_closer][j_closer][k_closer][l-1][4]));
        }
        /*
        cout << endl;
        cout << "l=" << l << endl;
        cout << "l_closer=" << l_closer << endl;
        cout << "x=" << x << ", low=" << table[i_closer][j_closer][k_closer][l-1][3] << ", high=" << table[i_closer][j_closer][k_closer][l][3] << ", closer=" << l-l_closer << endl;
        cout << "location=" << (x - table[i_closer][j_closer][k_closer][l-1][3])/(table[i_closer][j_closer][k_closer][l][3] - table[i_closer][j_closer][k_closer][l-1][3]) << ", higher_N=" << table[i_closer][j_closer][k_closer][l][4] << ", lower_N=" << table[i_closer][j_closer][k_closer][l-1][4] << ", corr=" << x_corr << endl;
        cout << endl;
        */
    }
    /*
    double return_value = table[i_closer][j_closer][k_closer][l_closer][4] + r_corr + b_corr + x_corr;
    cout << r_corr << ", " << b_corr << ", " << phi_corr << ", " << x_corr << endl;
    if (!exact_r) {
        cout << "return value: " << table[i][j_closer][k_closer][l_closer][4] << ", " << table[i-1][j_closer][k_closer][j_closer][4] << ", " << return_value << endl << endl;
    }
    */

    return exp(log(table[i_closer][j_closer][k_closer][l_closer][4]) + r_corr + b_corr + phi_corr + x_corr);
}

double get_Pb_dipole_amplitude(array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> &table, double r, double b, double phi, double x, bool limit_phi=true) {
    //cout << "Getting amplitude, r=" << r << ", b=" << b << ", phi=" << phi << ", x=" << x << endl;

    if (limit_phi) {
        if (phi > calc_max_phi(r, b)) {
            return 0;
        }
    }

    long unsigned int i = 0;
    long unsigned int i_closer = 0;
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
        }else if (log(table[i][0][0][0][0]) - log(r) < log(r) - log(table[i-1][0][0][0][0])) {
            i_closer = i;
        } else {
            i_closer = i-1;
        }
    }

    long unsigned int j = 0;
    long unsigned int j_closer = 0;
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
        }else if (log(table[i_closer][j][0][0][1]) - log(b) < log(b) - log(table[i_closer][j-1][0][0][1])) {
            j_closer = j;
        } else {
            j_closer = j-1;
        }
    }

    long unsigned int k = 0;
    long unsigned int k_closer = 0;
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
    long unsigned int l_closer = 0;
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
            //cout << table[i_closer][j_closer][k_closer][l][3] << endl;
            l++;
            if (l > table[i_closer][j_closer][k_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i_closer][j_closer][k_closer][l][3] == x) {
            l_closer = l;
            exact_x = true;
        }else if (log(table[i_closer][j_closer][k_closer][l][3]) - log(x) < log(x) - log(table[i_closer][j_closer][k_closer][l-1][3])) {
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
            r_corr = -1*(log(table[i][j_closer][k_closer][l_closer][4]) - log(table[i-1][j_closer][k_closer][l_closer][4]))*(1-(log(r) - log(table[i-1][0][0][0][0]))/(log(table[i][0][0][0][0]) - log(table[i-1][0][0][0][0])));
        } else {
            r_corr = (log(r) - log(table[i-1][0][0][0][0]))/(log(table[i][0][0][0][0]) - log(table[i-1][0][0][0][0]))*(log(table[i][j_closer][k_closer][l_closer][4]) - log(table[i-1][j_closer][k_closer][l_closer][4]));
        }
        /*
        cout << endl;
        cout << "r=" << r << ", low=" << table[i-1][0][0][0][0] << ", high=" << table[i][0][0][0][0] << ", closer=" << i-i_closer << endl;
        cout << "location=" << (r - table[i-1][0][0][0][0])/(table[i][0][0][0][0] - table[i-1][0][0][0][0]) << ", higher_N=" << table[i][j_closer][k_closer][l_closer][4] << ", lower_N=" << table[i-1][j_closer][k_closer][l_closer][4] << ", corr=" << r_corr << endl;
        cout << endl;
        */
    }

    if (exact_b) {
        b_corr = 0;
    } else {
        if (j_closer == j) {
            b_corr = -1*(log(table[i_closer][j][k_closer][l_closer][4]) - log(table[i_closer][j-1][k_closer][l_closer][4]))*(1-(log(b) - log(table[i_closer][j-1][0][0][1]))/(log(table[i_closer][j][0][0][1]) - log(table[i_closer][j-1][0][0][1])));
        } else {
            b_corr = (log(b) - log(table[i_closer][j-1][0][0][1]))/(log(table[i_closer][j][0][0][1]) - log(table[i_closer][j-1][0][0][1]))*(log(table[i_closer][j][k_closer][l_closer][4]) - log(table[i_closer][j-1][k_closer][l_closer][4]));
        }
        /*
        cout << endl;
        cout << "j=" << j << endl;
        cout << "j_closer=" << j_closer << endl;
        cout << "i_closer=" << i_closer << endl;
        cout << "b=" << b << ", low=" << table[i_closer][j-1][0][0][1] << ", high=" << table[i_closer][j][0][0][1] << ", closer=" << j-j_closer << endl;
        cout << "location=" << (b - table[i_closer][j-1][0][0][1])/(table[i_closer][j][0][0][1] - table[i_closer][j-1][0][0][1]) << ", higher_N=" << table[i_closer][j][k_closer][l_closer][4] << ", lower_N=" << table[i_closer][j-1][k_closer][l_closer][4] << ", corr=" << b_corr << endl;
        cout << endl;
        */
    }

    if (exact_phi) {
        phi_corr = 0;
    } else {
        if (k_closer == k) {
            phi_corr = -1*(log(table[i_closer][j_closer][k][l_closer][4]) - log(table[i_closer][j_closer][k-1][l_closer][4]))*(1-(phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2]));
        } else {
            phi_corr = (phi - table[i_closer][j_closer][k-1][0][2])/(table[i_closer][j_closer][k][0][2] - table[i_closer][j_closer][k-1][0][2])*(log(table[i_closer][j_closer][k][l_closer][4]) - log(table[i_closer][j_closer][k-1][l_closer][4]));
        }
    }

    if (exact_x) {
        x_corr = 0;
        //cout << "exact" << endl;
    } else {
        if (l_closer == l) {
            x_corr = -1*(log(table[i_closer][j_closer][k_closer][l][4]) - log(table[i_closer][j_closer][k_closer][l-1][4]))*(1-(log(x) - log(table[i_closer][j_closer][k_closer][l-1][3]))/(log(table[i_closer][j_closer][k_closer][l][3]) - log(table[i_closer][j_closer][k_closer][l-1][3])));
        } else {
            x_corr = (log(x) - log(table[i_closer][j_closer][k_closer][l-1][3]))/(log(table[i_closer][j_closer][k_closer][l][3]) - log(table[i_closer][j_closer][k_closer][l-1][3]))*(log(table[i_closer][j_closer][k_closer][l][4]) - log(table[i_closer][j_closer][k_closer][l-1][4]));
        }
        /*
        cout << endl;
        cout << "l=" << l << endl;
        cout << "l_closer=" << l_closer << endl;
        cout << "x=" << x << ", low=" << table[i_closer][j_closer][k_closer][l-1][3] << ", high=" << table[i_closer][j_closer][k_closer][l][3] << ", closer=" << l-l_closer << endl;
        cout << "location=" << (x - table[i_closer][j_closer][k_closer][l-1][3])/(table[i_closer][j_closer][k_closer][l][3] - table[i_closer][j_closer][k_closer][l-1][3]) << ", higher_N=" << table[i_closer][j_closer][k_closer][l][4] << ", lower_N=" << table[i_closer][j_closer][k_closer][l-1][4] << ", corr=" << x_corr << endl;
        cout << endl;
        */
    }
    /*
    double return_value = table[i_closer][j_closer][k_closer][l_closer][4] + r_corr + b_corr + x_corr;
    cout << r_corr << ", " << b_corr << ", " << phi_corr << ", " << x_corr << endl;
    if (!exact_r) {
        cout << "return value: " << table[i][j_closer][k_closer][l_closer][4] << ", " << table[i-1][j_closer][k_closer][j_closer][4] << ", " << return_value << endl << endl;
    }
    */

    return exp(log(table[i_closer][j_closer][k_closer][l_closer][4]) + r_corr + b_corr + phi_corr + x_corr);
}

double get_raw_p_dipole_amplitude(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, double r, double b, double phi, double x) {

    if (phi > calc_max_phi(r, b)) {
        return 0;
    }

    int r_index = 0;
    if (table[29][0][0][0][0] <= r) {
        r_index = 29;
    } else if (r <= table[0][0][0][0][0]) {
        r_index = 0;
    } else {
        while (table[r_index][0][0][0][0] < r) {
            r_index++;
        }
        if (r - table[r_index-1][0][0][0][0] < table[r_index][0][0][0][0] - r) {
            r_index = r_index - 1;
        }
    }

    int b_index = 0;
    if (table[0][29][0][0][1] <= b) {
        b_index = 29;
    } else if (b <= table[0][0][0][0][1]) {
        b_index = 0;
    } else {
        while (table[0][b_index][0][0][1] < b) {
            b_index++;
        }
        if (b - table[0][b_index-1][0][0][1] < table[0][b_index][0][0][1] - b) {
            b_index = b_index - 1;
        }
    }

    int phi_index = 0;
    if (table[0][0][29][0][2] <= phi) {
        phi_index = 29;
    } else if (phi <= table[0][0][0][0][2]) {
        phi_index = 0;
    } else {
        while (table[0][0][phi_index][0][2] < phi) {
            phi_index++;
        }
        if (phi - table[0][0][phi_index-1][0][2] < table[0][0][phi_index][0][2] - phi) {
            phi_index = phi_index - 1;
        }
    }

    int x_index = 0;
    if (table[0][0][0][80][3] <= x) {
        x_index = 80;
    } else if (x <= table[0][0][0][0][3]) {
        x_index = 0;
    } else {
        while (table[0][0][0][x_index][3] < x) {
            x_index++;
        }
        if (x - table[0][0][0][x_index-1][3] < table[0][0][0][x_index][3] - x) {
            x_index = x_index - 1;
        }
    }
    return table[r_index][b_index][phi_index][x_index][4];
}

void load_p_dipole_amplitudes(array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> &table, string filename, double plot=false) {
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
                    table[i][j][k][80-l][0] = r[index];
                    table[i][j][k][80-l][1] = b_min[index];
                    table[i][j][k][80-l][2] = phi[index];
                    table[i][j][k][80-l][3] = calc_x(Y[index]);
                    table[i][j][k][80-l][4] = N[index];
                }
            }
        }
    }
    cout << "Table ready" << endl;
    return;
}

void load_Pb_dipole_amplitudes(array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> &table, string filename, double plot=false) {
    cout << "Reading " << filename << endl;
    rapidcsv::Document doc(filename);
    
    vector<double> r = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");

    cout << "Data read successfully" << endl;
    
    const int icof = 40*40*81;
    const int jcof = 40*81;

    for (int i=0; i<40; i++) {
            for (int j=0; j<40; j++) {
                for (int k=0; k<40; k++) {
                    for (int l=0; l<81; l++) {
                    int index = i*icof + j*jcof + k*81 + l;
                    table[i][j][k][80-l][0] = r[index];
                    table[i][j][k][80-l][1] = b_min[index];
                    table[i][j][k][80-l][2] = phi[index];
                    table[i][j][k][80-l][3] = calc_x(Y[index]);
                    table[i][j][k][80-l][4] = N[index];
                }
            }
        }
    }
    cout << "Table ready" << endl;
    return;
}
