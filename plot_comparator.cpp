
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "cross_section_file_reader.h"

int main() {

    string filenames[] = {
        "output/rapidity_LHC_inclusive_D0_c_bfkl_Pb.txt",
        "output/rapidity_LHC_inclusive_D0_c_bk_Pb.txt"
    };

    string titles[] {
        "bfkl",
        "bk"
    };

    Color_t colors[] {
        kBlue,
        kRed
    };

    TMultiGraph* multigraph = new TMultiGraph();
    multigraph->SetTitle("inclusive charm cross section at pT=0 and Q^2=0;y;sigma");

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
        subgraph->SetLineColor(colors[i]);
        multigraph->Add(subgraph);
    }

    TCanvas* multigraph_canvas = new TCanvas("multigraph_canvas", "", 1100, 600);
    //gPad->SetLogx();
    multigraph->Draw("A"); //"A PMC PLC"
    multigraph_canvas->BuildLegend(0.55, 0.55, 0.7, 0.9);
    multigraph_canvas->Print("figures/inclusive_charm_cross_section_comparison.pdf");
}