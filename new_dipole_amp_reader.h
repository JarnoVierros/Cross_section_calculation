
#include <iostream>
#include <vector>
#include <rapidcsv.h>
#include <math.h>
/*
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLine.h"
*/
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
    if (r > 2*b_min) {
        return 2*M_PI;
    } else {
        return M_PI - acos(r/(2*b_min));
    }
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
    
    //double return_value = table[i+r_closer][j+b_closer][k+x_closer][3] + r_corr + b_corr + x_corr;
    //cout << r_corr << ", " << b_corr << ", " << x_corr << endl;
    //cout << "return value: " << table[i+r_closer][j+b_closer][k+x_closer][3] << ", " << table[i+r_closer-1][j+b_closer][k+x_closer][3] << ", " << return_value << endl << endl;
    //return return_value;

    return table[i+r_closer][j+b_closer][k+x_closer][3] + r_corr + b_corr + x_corr;
}

void generate_dipole_amplitudes(string in_filename, string out_filename) {
    cout << "Reading " << in_filename << endl;
    rapidcsv::Document doc(in_filename);
    
    vector<double> r_vec = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min_vec = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi_vec = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");

    cout << "Dipole amplitude read successfully" << endl;
    
    const int icof = 30*30*81;
    const int jcof = 30*81;

    ofstream outfile1("data/variable_change_table.txt");
    outfile1 << "r1, r2, b1, b2, z, x, N\n";

    int index;
    double r, b_min, phi, theta, z;
    double r1, r2, b1, b2, x, N_val;
    int i, j, k, l, m, n;
    array<double, 7> cell;

    for (i=0; i<30; i++) {
        for (j=0; j<30; j++) {
            for (k=0; k<30; k++) {
                for (l=0; l<81; l++) {
                    index = i*icof + j*jcof + k*81 + l;

                    r = r_vec[index];
                    b_min = b_min_vec[index];
                    phi = phi_vec[index];

                    if (phi > calc_max_phi(r, b_min)) {
                        //cout << "skipping" << endl;
                        continue; //Skip if phi is in forbidden region
                    }

                    for (m=0; m<100; m++) {
                        theta = m/100.0*2*M_PI;
                        for (n=0; n<100; n++) {
                            z = n/100.0;

                            r1 = -r*cos(theta+phi);
                            r2 = -r*sin(theta+phi);
                            b1 = b_min*cos(theta) + (1-z)*r*cos(theta+phi);
                            b2 = b_min*sin(theta) + (1-z)*r*sin(theta+phi);
                            x = calc_x(Y[index]);
                            N_val = N[index];

                            cell = {r1, r2, b1, b2, z, x, N_val};
                            //raw_table.push_back(cell);
                            outfile1 << r1 << "," << r2 << "," << b1 << "," << b2 << "," << z << "," << x << "," << N_val << "\n";

                            r1 = r*cos(theta+phi);
                            r2 = r*sin(theta+phi);

                            cell = {r1, r2, b1, b2, z, x, N_val};
                            //raw_table.push_back(cell);
                            outfile1 << r1 << "," << r2 << "," << b1 << "," << b2 << "," << z << "," << x << "," << N_val << "\n";

                            r1 = -r*cos(theta-phi);
                            r2 = -r*sin(theta-phi);
                            b1 = b_min*cos(theta) + (1-z)*r*cos(theta-phi);
                            b2 = b_min*sin(theta) + (1-z)*r*sin(theta-phi);

                            cell = {r1, r2, b1, b2, z, x, N_val};
                            //raw_table.push_back(cell);
                            outfile1 << r1 << "," << r2 << "," << b1 << "," << b2 << "," << z << "," << x << "," << N_val << "\n";

                            r1 = r*cos(theta-phi);
                            r2 = r*sin(theta-phi);

                            cell = {r1, r2, b1, b2, z, x, N_val};
                            //raw_table.push_back(cell);
                            outfile1 << r1 << "," << r2 << "," << b1 << "," << b2 << "," << z << "," << x << "," << N_val << "\n";
                        }
                    }
                }
                //cout << raw_table.size() << endl;
            }
            cout << "b_min" << endl;
        }
        cout << "r" << endl;
    }
    outfile1.close();
    cout << "Variable change completed" << endl;

    //array<array<array<array<array<array<array<double, 7>, 81>, 100>, 100>, 100>, 100>, 100> table;

    ofstream outfile(out_filename);
    outfile << "r1, r2, b1, b2, z, x, N\n";

    const int x_steps = 100;
    const double x_start = 1e-9; //2.37024e-05
    const double x_stop = 0.01;
    const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

    const int table_size = 100;
    const double r_limit = 40;
    const double b_limit = 20;
    for (int r1i=0; r1i<table_size; r1i++) {
        for (int r2i=0; r2i<table_size; r2i++) {
            for (int b1i=0; b1i<table_size; b1i++) {
                for (int b2i=0; b2i<table_size; b2i++) {
                    for (int zi=0; zi<table_size; zi++) {
                        for (int xi=0; xi<81; xi++) {
                            double r1 = r1i/table_size*2.0*r_limit - r_limit;
                            double r2 = r2i/table_size*2.0*r_limit - r_limit;
                            double b1 = b1i/table_size*2.0*b_limit - b_limit;
                            double b2 = b2i/table_size*2.0*b_limit - b_limit;
                            double z = zi/table_size*1.0;
                            double x = pow(10, log10(x_start) + xi*x_step);

                            int best_index = 0;
                            double best_distance = abs(r1-raw_table[0][0])/34.0 + abs(r2-raw_table[0][1])/34.0
                                 + abs(b1-raw_table[0][3])/17.0 + abs(b1-raw_table[0][4])/17.0
                                 + abs(z-raw_table[0][5]) + abs(x - raw_table[0][6]);;

                            for (int i=0; i<raw_table.size(); i++) {
                                double distance = abs(r1-raw_table[i][0])/34.0 + abs(r2-raw_table[i][1])/34.0
                                 + abs(b1-raw_table[i][3])/17.0 + abs(b1-raw_table[i][4])/17.0
                                 + abs(z-raw_table[i][5]) + abs(x - raw_table[i][6]); // might need to boost x distance
                                if (distance < best_distance) {
                                best_index = i;
                                best_distance = distance;
                                }
                            }
                            /*
                            table[r1i][r2i][b1i][b2i][zi][xi][0] = raw_table[best_index][0];
                            table[r1i][r2i][b1i][b2i][zi][xi][1] = raw_table[best_index][1];
                            table[r1i][r2i][b1i][b2i][zi][xi][2] = raw_table[best_index][2];
                            table[r1i][r2i][b1i][b2i][zi][xi][3] = raw_table[best_index][3];
                            table[r1i][r2i][b1i][b2i][zi][xi][4] = raw_table[best_index][4];
                            table[r1i][r2i][b1i][b2i][zi][xi][5] = raw_table[best_index][5];
                            table[r1i][r2i][b1i][b2i][zi][xi][6] = raw_table[best_index][6];       
                            */
                            outfile << raw_table[best_index][0] << "," << raw_table[best_index][1] << "," << raw_table[best_index][3] << "," << raw_table[best_index][4] << "," << raw_table[best_index][5] << "," << raw_table[best_index][6] << "\n";
                        }
                    }
                }
            }
        }
    }
    cout << "Program finished successfully" << endl;
}


/*
void load_dipole_amplitudes(array<array<array<array<double, 4>, 81>, 30>, 30> &table, string filename, double plot=false) {
    cout << "Reading " << filename << endl;
    rapidcsv::Document doc(filename);
    
    vector<double> r = doc.GetColumn<double>("r [GeV^-1]");
    vector<double> b_min = doc.GetColumn<double>("b_min [GeV^-1]");
    vector<double> phi = doc.GetColumn<double>("phi");
    vector<double> Y = doc.GetColumn<double>("Y");
    vector<double> N = doc.GetColumn<double>("N");

    cout << "Dipole amplitude read successfully" << endl;
    
    const int icof = 30*30*81;
    const int jcof = 30*81;

    const int b_steps = 30;
    const double b_start = 4e-5; //2.37024e-05
    const double b_stop = 17.3;
    const double b_step = 1.0/(b_steps-1)*log10(b_stop/b_start);

    for (int i=0; i<30; i++) {
        for (int l=0; l<81; l++) {
            vector<double> sub_b, sub_x, sub_N;
            for (int j=0; j<30; j++) {
                for (int k=0; k<30; k++) {
                    int index = i*icof + j*jcof + k*81 + l;
                    //cout << "phi=" << phi[index] << ", max=" << calc_max_phi(r[index], b_min[index]) << endl;
                    if (phi[index] > calc_max_phi(r[index], b_min[index])) {
                        //cout << "skipping" << endl;
                        continue; //Skip if phi is in forbidden region
                    }
                    double sub_b_value = calc_b(r[index], b_min[index], phi[index]);
                    if (sub_b_value == 0) {
                        //cout << "zero found" << endl;
                        sub_b_value = 2.5e-5;
                    }
                    //if (sub_b_value>2) {continue;}
                    sub_b.push_back(sub_b_value);
                    sub_x.push_back(calc_x(Y[index]));
                    sub_N.push_back(N[index]);
                }
            }
            double centers[30], averages[30], borders[31];
            double max_N = 0;
            for (int n=0; n<30; n++) {
                int hits = 0;
                double total_N = 0;
                double b_center = pow(10, log10(b_start) + n*b_step);
                double low_limit;
                if (n == 0) {
                    low_limit = 2e-5;
                } else {
                    low_limit = sqrt(pow(10, log10(b_start) + (n-1)*b_step)*pow(10, log10(b_start) + n*b_step));
                }
                double high_limit;
                if (n == 29) {
                    high_limit = 1.5*b_center;
                } else {
                    high_limit = sqrt(pow(10, log10(b_start) + n*b_step)*pow(10, log10(b_start) + (n+1)*b_step));
                }
                //cout << "options: " << sub_b.size() << endl;
                for (long unsigned int m=0; m<sub_b.size(); m++) {
                    if (low_limit < sub_b[m] && sub_b[m] < high_limit) {
                        hits++;
                        total_N += sub_N[m];
                    }
                }
                //cout << "hits: " << hits << endl;
                if (hits < 1) {
                    if (n == 0) {
                        hits = 1;
                        total_N = 0;
                    } else {
                        hits = 1;
                        total_N = table[i][n-1][l][3];
                    }
                }
                table[i][n][l][0] = r[i*icof];
                table[i][n][l][1] = b_center;
                table[i][n][l][2] = sub_x[n];
                table[i][n][l][3] = total_N/hits;
                if (plot) {
                    if (max_N < total_N/hits) {max_N = total_N/hits;}
                    centers[n] = b_center;
                    averages[n] = total_N/hits;
                    borders[n] = low_limit;
                    if (n==29) {borders[30] = high_limit;}
                    //cout << table[i][n][l][0] << ", " << table[i][n][l][1] << ", " << table[i][n][l][2] << ", " << table[i][n][l][3] << endl;
                }
            }
            if (plot) {
                double* sub_b_arr = &sub_b[0];
                double* sub_N_arr = &sub_N[0];
                int array_size = sub_b.size();
                TGraph* datapoints = new TGraph(array_size, sub_b_arr, sub_N_arr);
                datapoints->SetMarkerStyle(6);
                datapoints->SetMarkerColor(2);

                TGraph* averages_graph = new TGraph(30, centers, averages);

                vector<TLine*> lines(31);
                for (int n=0; n<31; n++) {
                    lines[n] = new TLine();
                    lines[n]->SetX1(borders[n]);
                    lines[n]->SetX2(borders[n]);
                    lines[n]->SetY1(0);
                    lines[n]->SetY2(1.1*max_N);
                    //cout << borders[n] << ", " << max_N << endl;
                }

                TCanvas* canvas = new TCanvas();

                averages_graph->GetXaxis()->SetLimits(borders[0]/1.5, borders[30]*1.5);
                //averages_graph->GetXaxis()->SetLimits(5e-7, 2);
                averages_graph->Draw("A*");
                TString title = "N as a function of b at r=" + to_string(r[i*icof]) + ", x="+to_string(calc_x(Y[i*icof+l]));
                averages_graph->SetTitle(title);
                averages_graph->GetXaxis()->SetTitle("b");
                averages_graph->GetYaxis()->SetTitle("N");
                datapoints->Draw("P");

                for (int n=0; n<31; n++) {
                    lines[n]->Draw();
                }

                gPad->SetLogx();
                //gPad->SetLogy();

                title = "b_figs/r_" + to_string(r[i*icof]) + "_x_" + to_string(calc_x(Y[i*icof+l])) + ".pdf";
                canvas->Print(title);
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
*/