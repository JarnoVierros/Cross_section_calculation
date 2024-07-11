
#include "TGraph.h"
#include "TCanvas.h"

#include "new_dipole_amp_reader.h"
#include <chrono>
#include <iostream>

int main() {

    static array<array<array<array<double, 4>, 81>, 30>, 30> table;
    string filename = "data/dipole_amplitude_with_IP_dependence.csv";
    
    load_dipole_amplitudes(table, filename);
    cout << "Amplitude loaded" << endl;
    
    /*
    double min_bs[30*81];
    double x[30*81];
    int count = 0;
    for (int i=0; i<30; i++) {
        for (int k=0; k<81; k++) {
            //cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
            min_bs[count] = table[i][0][k][1];
            x[count] = count;
            count++;
        }
    }
    */
    /*
    double b[30];
    double N[30];
    for (int i=0; i<30; i++) {
        //cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
        b[i] = table[0][i][0][1];
        N[i] = table[0][i][0][3];
    }

    double edge[31], edgeN[31];
    const int b_steps = 30;
    const double b_start = 1e-05; //2.37024e-05
    const double b_stop = 17.3;
    const double b_step = 1.0/(b_steps-1)*log10(b_stop/b_start);
    const int edge_b_steps = b_steps + 2;
    const double edge_b_start = pow(10, log10(b_start) - 1/2*b_step);
    const double edge_b_stop = pow(10, log10(b_stop) + 1/2*b_step);;
    const double edge_b_step = 1.0/(edge_b_steps-1)*log10(edge_b_stop/edge_b_start);
    cout << edge_b_start << endl;
    cout << edge_b_step << endl;
    cout << endl;
    for (int i=0; i<29; i++) {
        double local_start = pow(10, log10(b_start) + i*b_step);
        double local_stop = pow(10, log10(b_start) + (i+1)*b_step);
        //edge[i] = pow(10, 1/2*log10(local_start*local_stop));
        edge[i+1] = sqrt(local_start*local_stop);
        edgeN[i+1] = 1.1;
        cout << edge[i+1] << endl;
    }
    edge[0] = edge[1]/2;
    edgeN[0] = 1.1;
    edge[30] = edge[29]*2;
    edgeN[30] = 1.1;

    TGraph* edge_graph = new TGraph(31, edge, edgeN);
    TGraph* graph = new TGraph(30, b, N);

    TCanvas* canvas = new TCanvas();
    edge_graph->Draw("A*");
    graph->Draw("*");
    

    gPad->SetLogx();
    //gPad->SetLogy();

    canvas->Print("test_b.pdf");
    */

    /*
    for (int i=0; i<30; i++) {
        for (int j=0; j<30; j++) {
            for (int k=0; k<81; k++) {
                if (i == 29) {
                    cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
                }
            }
            if (i==29) {
                getchar();
            }
        }
    }
    */
    /*
    double b[30], r[30];
    
    for (int i=0; i<30; i++) {
        for (int j=0; j<1; j++) {
            b[i] = table[i][29][0][1];
            r[i] = table[i][29][0][0];
            cout << i << ": " << r[i] << ", " << b[i] << endl;
        }
    }

    TGraph* graph = new TGraph(30, r, b);
    graph->SetLineColor(2);

    TCanvas* canvas = new TCanvas();

    graph->Draw("AL");

    //gPad->SetLogx();
    //gPad->SetLogy();

    canvas->Print("test_b.pdf");
    */
    /*
    auto start = chrono::high_resolution_clock::now();
    double res = get_dipole_amplitude(table, 121.77, 11.84, 1e-05);
    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::nanoseconds>(stop-start);
    cout << "duration: " << duration.count() << endl;

    cout << "amplitude: " << res << endl;
    */
    //0.000155246, 2.43256e-05, 1.12535e-09, 4.95126e-06
}
