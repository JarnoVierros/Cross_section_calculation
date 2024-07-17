
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "cross_section_file_reader.h"

int main() {

    string filenames[] = {
        "data/J_T_inclusive_r_34.txt",
        "data/J_T_inclusive_r_30.txt",
        "data/J_T_inclusive_r_20.txt",
        "data/J_T_inclusive_r_10.txt",
        "data/J_T_inclusive_r_5.txt",
        "data/J_T_inclusive_r_4.txt",
        "data/J_T_inclusive_r_3.txt",
        "data/J_T_inclusive_r_2.txt",
        "data/J_T_inclusive_r_1.txt",
        "data/J_T_inclusive_r_0_50.txt"
    };

    string titles[] {
        "r limit=34",
        "r limit=30",
        "r limit=20",
        "r limit=10",
        "r limit=5",
        "r limit=4",
        "r limit=3",
        "r limit=2",
        "r limit=1",
        "r limit=0.5"
    };

    TMultiGraph* multigraph = new TMultiGraph();
    multigraph->SetTitle("Transverse cross sections with different dipole size integration limits at Q=2.5 GeV;x;sigma (mb)");

    for (int i=0; i<size(filenames); i++) {
        vector<double> Q2, x, sigma, sigma_error;
        cout << "reading: " << filenames[i] << endl;
        read_sigma_file(filenames[i], Q2, x, sigma, sigma_error);

        double* x_arr = &x[0];
        double* sigma_arr = &sigma[0];
        double* sigma_error_arr = &sigma_error[0];
        double x_error[x.size()];
        for (int j=0; j<x.size(); j++) {
            x_error[j] = 0;
        }

        TGraphErrors* subgraph = new TGraphErrors(x.size(), x_arr, sigma_arr, x_error);
        TString subgraph_name = titles[i];
        subgraph->SetTitle(subgraph_name);
        multigraph->Add(subgraph);
    }

    TCanvas* multigraph_canvas = new TCanvas("multigraph_canvas", "", 1100, 600);
    gPad->SetLogx();
    multigraph->Draw("A PMC PLC");
    multigraph_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
    multigraph_canvas->Print("figures/transverse_r_limit_comparison.pdf");
}