
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "cross_section_file_reader.h"

#include <iostream>

using namespace std;

int main() {


    string filename = "output/Phi_table_p_bk_2A.txt";

    vector<double> P2, y, xpom, result, error;
    cout << "reading: " << filename << endl;
    read_Phi_file(filename, P2, y, xpom, result, error);

    double P2_choice = 13.8017;
    double y_choice = 0.317241;
    double xpom_choice = 0.00126896;

    vector<double> X, Y, X_error, Y_error;

    for (int i=0; i<P2.size(); i++) {
        if (P2[i] == P2_choice && y[i] == y_choice && xpom[i] < 0.01) {
            X.push_back(xpom[i]);
            Y.push_back(result[i]);
            X_error.push_back(0);
            Y_error.push_back(error[i]);
        }
    }

    double* X_arr = &X[0];
    double* Y_arr = &Y[0];
    double* X_error_arr = &X_error[0];
    double* Y_error_arr = &Y_error[0];

    TCanvas* canvas = new TCanvas("asd", "asd");

    TGraphErrors* graph = new TGraphErrors(X.size(), X_arr, Y_arr, X_error_arr, Y_error_arr);

    //gPad->SetLogy();
    gPad->SetLogx();
    graph->Draw("APL");

    canvas->Print("figures/Phi_plot.pdf");
}