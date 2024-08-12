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
    array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> table;

    load_dipole_amplitudes(table, filename);

    int data_size = table.size();

    TMultiGraph* multigraph = new TMultiGraph();

    //const double r = 0;
    const double b_min = 0;
    const double phi = 0;
    const double x = 0;

    const int x_steps = 50;
    const double x_start = 2e1;
    const double x_stop = 2e4;
    const double x_step = 1.0/(x_steps-1)*log10(x_stop/x_start);

    double raw_x[x_steps], raw_y[x_steps], interpolated_x[x_steps], interpolated_y[x_steps];

    for (int i=0; i<x_steps; i++) {
        double x = pow(10, log10(x_start) + i*x_step);
        raw_x[i] = x;
        raw_y[i] = get_raw_dipole_amplitude(table, x, b_min, phi, x);
        interpolated_x[i] = x;
        interpolated_y[i] = get_dipole_amplitude(table, x, b_min, phi, x);
    }

    TGraph* raw_graph = new TGraph(x_steps, raw_x, interpolated_x);
    raw_graph->SetTitle("Raw");
    multigraph->Add(raw_graph);

    TGraph* interpolated_graph = new TGraph(x_steps, interpolated_x, interpolated_y);
    interpolated_graph->SetTitle("Interpolated");
    multigraph->Add(interpolated_graph);

    TCanvas* canvas = new TCanvas("canvas1", "", 1000, 800);

    multigraph->Draw("A PMC PLC");
    canvas->BuildLegend();
    canvas->Print("interpoaltion_test.pdf");
}