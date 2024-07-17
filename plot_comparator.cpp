
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "cross_section_file_reader.h"

int main() {

    string filenames[] = {
        "asd",
        "dsa",
        "wsad"
    };

    string titles[] {
        "a",
        "b",
        "c"
    };

    TMultiGraph* multigraph = new TMultiGraph();
    multigraph->SetTitle("Longitudinal cross sections with different dipole size integration limits at Q=2.5 GeV;x;sigma (mb)");

    for (int i=0; i<size(filenames); i++) {
        vector<double> Q2, x, sigma, sigma_error;
        read_sigma_file(filenames[i], Q2, x, sigma, sigma_error);

        double* x_arr = &x[0];
        double* sigma_arr = &sigma[0];
        double* sigma_error_arr = &sigma_error[0];
        double x_error[x.size()];

        TGraphErrors* subgraph = new TGraphErrors(x.size(), x_arr, sigma_arr, x_error);
        TString subgraph_name = titles[i];
        subgraph->SetTitle(subgraph_name);
        multigraph->Add(subgraph);
    }

    TCanvas* multigraph_canvas = new TCanvas("multigraph_canvas");
    multigraph_canvas->Print("longitudinal_r_limit_comparison.pdf");
}