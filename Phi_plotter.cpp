
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


    string filename = "output/Phi_archive/Phi_table_p_bk_2A.txt";

    vector<double> P2, y, xpom, result, error;
    cout << "reading: " << filename << endl;
    read_Phi_file(filename, P2, y, xpom, result, error);

    double P2_choice_up = 0.227603; //0.805927
    double P2_choice_down = P2_choice_up;

    double y_choice_up = 0.13895; //0.0316228
    double y_choice_down = y_choice_up;

    double xpom_choice_up = 0.000372759; //0.000316228
    double xpom_choice_down = xpom_choice_up;

    vector<double> X, Y, X_error, Y_error;

    string show = "y";

    if (show == "P2") {
        for (int i=0; i<P2.size(); i++) {
            if ((y_choice_down <= y[i] && y[i] <= y_choice_up) && (xpom_choice_down <= xpom[i] && xpom[i] <= xpom_choice_up)) {
                X.push_back(P2[i]);
                Y.push_back(result[i]);
                X_error.push_back(0);
                Y_error.push_back(error[i]);
            }
        }
    } else if (show == "y") {
        for (int i=0; i<P2.size(); i++) {
            if ((P2_choice_down <= P2[i] && P2[i] <= P2_choice_up) && (xpom_choice_down <= xpom[i] && xpom[i] <= xpom_choice_up)) {
                X.push_back(y[i]);
                Y.push_back(result[i]);
                X_error.push_back(0);
                Y_error.push_back(error[i]);
            }
        }
    } else if (show == "xpom") {
        for (int i=0; i<P2.size(); i++) {
            if ((P2_choice_down <= P2[i] && P2[i] <= P2_choice_up) && (y_choice_down <= y[i] && y[i] <= y_choice_up)) {
                X.push_back(xpom[i]);
                Y.push_back(result[i]);
                X_error.push_back(0);
                Y_error.push_back(error[i]);
            }
        }
    } else {
        return 0;
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