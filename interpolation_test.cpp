#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <thread>
using namespace std;

#include "direct_dipole_amp_reader.h"

int main() {
    
    string dipole_amp_type = "bk";
    string nucleus_type = "p";
    string filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
    static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;
    
    load_dipole_amplitudes(table, filename);
    int data_size = table.size();

    TMultiGraph* multigraph = new TMultiGraph();
   

    double r = table[20][0][0][0][0];
    double b_min = table[0][25][0][0][1];
    //double phi = table[0][0][6][0][2];
    double x_pom = table[0][0][0][50][3];

    cout << "r=" << r << endl;
    cout << "b_min=" << b_min << endl;
    //cout << "phi=" << phi << endl;
    cout << "x_pom=" << x_pom << endl;

    TString title = "Dipole amplitude with r="+to_string(r)+", b_min="+to_string(b_min)+", x="+to_string(x_pom)+";phi (rad);N";
    multigraph->SetTitle(title);

    const int x_steps = 10000;

    //r, b_min
    /*
    const double x_start = 1e-5;
    const double x_stop = 100;
    const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);
    */

    //x
    const double x_start = 1e-10;
    const double x_stop = 0.1;
    const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);
    

    double raw_x[x_steps], raw_y[x_steps], interpolated_x[x_steps], interpolated_y[x_steps];

    for (int i=0; i<x_steps; i++) {
        //double x = pow(10, log10(x_start) + i*x_step);
        double x = 1.0*i/x_steps*M_PI;
        double phi = x;
        raw_x[i] = x;
        raw_y[i] = get_raw_dipole_amplitude(table, r, b_min, phi, x_pom);
        cout << "x=" << x << ", raw=" << raw_y[i];
        interpolated_x[i] = x;
        interpolated_y[i] = get_dipole_amplitude(table, r, b_min, phi, x_pom);
        cout << ", interpolated=" << interpolated_y[i] << endl;
    }

    TGraph* raw_graph = new TGraph(x_steps, raw_x, raw_y);
    raw_graph->SetTitle("Raw");
    multigraph->Add(raw_graph, "L");

    TGraph* interpolated_graph = new TGraph(x_steps, interpolated_x, interpolated_y);
    interpolated_graph->SetTitle("Interpolated");
    multigraph->Add(interpolated_graph, "L");

    TCanvas* canvas = new TCanvas("canvas1", "", 1000, 800);
    //gPad->SetLogx();
    //gPad->SetLogy();
    multigraph->GetYaxis()->SetRangeUser(0.024, 0.036);
    multigraph->Draw("A PMC PLC");
    canvas->BuildLegend(0.7, 0.7, 0.9, 0.9);
    canvas->Print("figures/interpolation_test.pdf");
    
    return 0;
}