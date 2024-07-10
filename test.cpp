
#include "TGraph.h"
#include "TCanvas.h"

#include "dipole_amp_reader.h"
#include <chrono>
#include <iostream>

int main() {

    static array<array<array<array<double, 4>, 81>, 30>, 30> table;
    string filename = "data/dipole_amplitude_with_IP_dependence.csv";
    
    load_dipole_amplitudes(table, filename);
    cout << "Amplitude loaded" << endl;
    string nothing;
    /*
    double min = 9999;
    double max = 0;
    for (int i=0; i<30; i++) {
        for (int j=0; j<30; j++) {
            for (int k=0; k<81; k++) {
                cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
                if (table[i][j][k][1] > max) {
                    max = table[i][j][k][1];
                }
                if (table[i][j][k][1] < min) {
                    min = table[i][j][k][1];
                }
            }
            //getchar();
        }
    }
    cout << "min: " << min << endl;
    cout << "max: " << max << endl;
    */
    
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
