
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include "cross_section_file_reader.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <thread>

int main() {

    string filenames[] = {
        "output/rapidity_LHC_inclusive_D0_c_bfkl_Pb.txt",
        "output/rapidity_LHC_inclusive_D0_c_bk_Pb.txt"
    };

    string titles[] {
        "bfkl",
        "bk"
    };

    int linestyles[] {
        2, //9
        1
    };

    Color_t colors[] {
        kBlue,
        kRed,
        kGreen
    };

    TMultiGraph* multigraph = new TMultiGraph();
    multigraph->SetTitle("inclusive charm cross section;rapidity;cross section (GeV^{-2})");

    for (int i=0; i<size(filenames); i++) {
        vector<double> pT, y, sigma, sigma_error;
        cout << "reading: " << filenames[i] << endl;
        read_sigma_file(filenames[i], pT, y, sigma, sigma_error);

        vector<vector<double>> y_columns, sigma_columns, sigma_err_columns;
        vector<double> pT_values;
        double current_pT = pT[0];
        vector<double> current_y_column, current_sigma_column, current_sigma_err_column;
        for (int j=0; j<size(pT); j++) {
            if (pT[j] != current_pT) {
                y_columns.push_back(current_y_column);
                sigma_columns.push_back(current_sigma_column);
                sigma_err_columns.push_back(current_sigma_err_column);
                pT_values.push_back(current_pT);
                current_y_column = {};
                current_sigma_column = {};
                current_pT = pT[j];
            }
            current_y_column.push_back(y[j]);
            current_sigma_column.push_back(sigma[j]);
            current_sigma_err_column.push_back(sigma_error[j]);
        }
        y_columns.push_back(current_y_column);
        sigma_columns.push_back(current_sigma_column);
        sigma_err_columns.push_back(current_sigma_err_column);
        pT_values.push_back(current_pT);

        for (int j=0; j<size(pT_values); j++) {
            double* y_arr = &y_columns[j][0];
            double* sigma_arr = &sigma_columns[j][0];
            double* sigma_error_arr = &sigma_err_columns[j][0];
            double y_error[y.size()];
            for (int k=0; k<y.size(); k++) {
                y_error[k] = 0;
            }

            TGraphErrors* subgraph = new TGraphErrors(size(y_columns[j]), y_arr, sigma_arr, y_error);

            stringstream pT_stream;
            pT_stream << fixed << setprecision(0) << pT_values[j];
            string pT_string = pT_stream.str();
            string name_str = titles[i] + " pT=" + pT_string;
            TString subgraph_name = name_str;

            subgraph->SetTitle(subgraph_name);
            subgraph->SetLineColor(colors[j]);
            subgraph->SetLineStyle(linestyles[i]);
            multigraph->Add(subgraph);
        }
    }

    TCanvas* multigraph_canvas = new TCanvas("multigraph_canvas", "", 1100, 600);
    //gPad->SetLogx();
    multigraph->Draw("A"); //"A PMC PLC"
    multigraph_canvas->BuildLegend(0.25, 0.55, 0.4, 0.9);
    multigraph_canvas->Print("figures/inclusive_charm_cross_section_rapidity.pdf");
}