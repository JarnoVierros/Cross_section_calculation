
#include <iostream>
#include <vector>
#include <rapidcsv.h>
#include <math.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

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


double get_dipole_amplitude(array<array<array<array<double, 4>, 81>, 30>, 30> &table, double r, double b, double x, bool restrict_b=false) {
    //cout << "Getting amplitude, r=" << r << ", b=" << b << ", x=" << x << endl;
    long unsigned int i = 0;
    int r_closer = 0; //shifts i to the cell closer to the desired value
    bool exact_r = false;
    if (r <= table[0][0][0][0]) {
        i = 0;
        exact_r = true;
    } else if (r >= table[table.size()-1][0][0][0]) {
        i = table.size()-1;
        exact_r = true;
    } else {
        while (table[i][0][0][0] < r) {
            i++;
            if (i > table.size()-1) {
                throw 1;
            }
        }
        if (table[i][0][0][0] == r) {
            r_closer = 0;
        }else if (table[i][0][0][0] - r < r - table[i-1][0][0][0]) {
            r_closer = 0;
        } else {
            r_closer = -1;
        }
    }

    long unsigned int j = 0;
    int b_closer = 0;
    bool exact_b = false;
    if (b <= table[i+r_closer][0][0][1]) {
        if (restrict_b) {
            return 0;
        }
        j = 0;
        exact_b = true;
        //cout << "Getting amplitude, r=" << r << ", b=" << b << ", x=" << x << endl;
        //cout << "b less than " << table[i+r_closer][0][0][1] << endl;
    } else if (b >= table[i+r_closer][table[i+r_closer].size()-1][0][1]) {
        if (restrict_b) {
            return 0;
        }
        j = table[i+r_closer].size()-1;
        exact_b = true;
        //cout << "Getting amplitude, r=" << r << ", b=" << b << ", x=" << x << endl;
        //cout << "b more than " << table[i+r_closer][table[i+r_closer].size()-1][0][1] << endl;
    } else {
        while (table[i+r_closer][j][0][1] < b) {
            j++;
            if (j > table[i+r_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i+r_closer][j][0][1] == b) {
            b_closer = 0;
        }else if (table[i+r_closer][j][0][1] - b < b - table[i+r_closer][j-1][0][1]) {
            b_closer = 0;
        } else {
            b_closer = -1;
        }
    }

    long unsigned int k = 0;
    int x_closer = 0; // 2 meas exact match
    bool exact_x = false;
    if (x <= table[i+r_closer][j+b_closer][0][2]) {
        k = 0;
        exact_x = true;
    } else if (x >= table[i+r_closer][j+b_closer][table[i+r_closer][j+b_closer].size()-1][2]) {
        k = table[i+r_closer][j+b_closer].size()-1;
        exact_x = true;
    } else {
        while (table[i+r_closer][j+b_closer][k][2] < x) {
            k++;
            if (k > table[i+r_closer][j+b_closer].size()-1) {
                throw 1;
            }
        }
        if (table[i+r_closer][j+b_closer][k][2] == x) {
            x_closer = 0;
        }else if (table[i+r_closer][j+b_closer][k][2] - x < x - table[i+r_closer][j+b_closer][k-1][2]) {
            x_closer = 0;
        } else {
            x_closer = -1;
        }
    }

    double r_corr, b_corr, x_corr;

    if (exact_r) {
        r_corr = 0;
    } else {
        if (r_closer == 0) {
            r_corr = -1*(table[i][j+b_closer][k+x_closer][3] - table[i-1][j+b_closer][k+x_closer][3])*(1-(r - table[i-1][0][0][0])/(table[i][0][0][0] - table[i-1][0][0][0]));
        } else {
            r_corr = (r - table[i-1][0][0][0])/(table[i][0][0][0] - table[i-1][0][0][0])*(table[i][j+b_closer][k+x_closer][3] - table[i-1][j+b_closer][k+x_closer][3]);
        }
        //cout << endl;
        //cout << "r=" << r << ", low=" << table[i-1][0][0][0] << ", high=" << table[i][0][0][0] << ", closer=" << r_closer << endl;
        //cout << "location=" << (r - table[i-1][0][0][0])/(table[i][0][0][0] - table[i-1][0][0][0]) << ", higher_N=" << table[i][j+b_closer][k+x_closer][3] << ", lower_N=" << table[i-1][j+b_closer][k+x_closer][3] << ", corr=" << r_corr << endl;
        //cout << endl;
    }

    if (exact_b) {
        b_corr = 0;
    } else {
        if (b_closer == 0) {
            b_corr = -1*(table[i+r_closer][j][k+x_closer][3] - table[i+r_closer][j-1][k+x_closer][3])*(1 - (b - table[i+r_closer][j-1][0][1])/(table[i+r_closer][j][0][1] - table[i+r_closer][j-1][0][1]));
        } else {
            b_corr = (b - table[i+r_closer][j-1][0][1])/(table[i+r_closer][j][0][1] - table[i+r_closer][j-1][0][1])*(table[i+r_closer][j][k+x_closer][3] - table[i+r_closer][j-1][k+x_closer][3]);
        }
        //cout << endl;
        //cout << "b=" << b << ", low=" << table[i+r_closer][j-1][0][1] << ", high=" << table[i+r_closer][j][0][1] << ", closer=" << b_closer << endl;
        //cout << "location=" << (b - table[i+r_closer][j-1][0][1])/(table[i+r_closer][j][0][1] - table[i+r_closer][j-1][0][1]) << ", higher_N=" << table[i+r_closer][j][k+x_closer][3] << ", lower_N=" << table[i+r_closer][j-1][k+x_closer][3] << ", corr=" << b_corr << endl;
        //cout << endl;
    }

    if (exact_x) {
        x_corr = 0;
    } else {
        if (x_closer == 0) {
            x_corr = -1*(table[i+r_closer][j+b_closer][k][3] - table[i+r_closer][j+b_closer][k-1][3])*(1 - (x - table[i+r_closer][j+b_closer][k-1][2])/(table[i+r_closer][j+b_closer][k][2] - table[i+r_closer][j+b_closer][k-1][2]));
        } else {
            x_corr = (x - table[i+r_closer][j+b_closer][k-1][2])/(table[i+r_closer][j+b_closer][k][2] - table[i+r_closer][j+b_closer][k-1][2])*(table[i+r_closer][j+b_closer][k][3] - table[i+r_closer][j+b_closer][k-1][3]);
        }
    }
    
    double return_value = table[i+r_closer][j+b_closer][k+x_closer][3] + r_corr + b_corr + x_corr;
    //cout << r_corr << ", " << b_corr << ", " << x_corr << endl;
    //cout << "return value: " << table[i+r_closer][j+b_closer][k+x_closer][3] << ", " << table[i+r_closer-1][j+b_closer][k+x_closer][3] << ", " << return_value << endl << endl;
    //return return_value;
    return table[i+r_closer][j+b_closer][k+x_closer][3] + r_corr + b_corr + x_corr;
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

    int count = 0;
    double peak_b_min, peak_phi, peak_r, peak_x;
    double peak_N = 0;
    for (int i=15; i<30; i++) {
        int past_count = 630;
        double sub_N[past_count], sub_b[past_count];
        for (int l=40; l<81; l++) {
            for (int j=0; j<30; j++) {
                for (int k=0; k<30; k++) {
                    int index = i*icof + j*jcof + k*81 + l;
                    
                    if (N[index] > 1e-9) {
                        //continue; // 4.32171e-09
                    }
                    if (calc_b(r[index], b_min[index], phi[index]) < 0.5) {
                        //cout << N[index] << endl;
                        count++;
                        sub_N[j*30 + k] = N[index];
                        sub_b[j*30 + k] = calc_b(r[index], b_min[index], phi[index]);
                        if (peak_N < N[index]) {
                            peak_N = N[index];
                            peak_x = Y[index];
                            peak_b_min = b_min[index];
                            peak_phi = phi[index];
                            peak_r = r[index];
                        }
                    }
                }
            }
            cout << peak_r << ", " << peak_b_min << ", " << peak_phi << ", " << peak_x << ", " << peak_N << ", " << endl;
            cout << "index: " << count << endl;
            TGraph* graph = new TGraph(past_count, sub_b, sub_N);
            TString title = "Dipole amplitude at r="+to_string(r[i*icof])+", x="+to_string(exp(-Y[i*icof+l])*0.01);
            graph->SetTitle(title);
            graph->GetXaxis()->SetTitle("b");
            graph->GetYaxis()->SetTitle("N");
            TCanvas* canvas = new TCanvas();
            graph->Draw("AP");
            gPad->SetLogx();
            //gPad->SetLogy();
            canvas->Print("test.pdf");
            return;
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
        for (int j=0; j<30; j++) {
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