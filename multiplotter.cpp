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

#include "cross_section_file_reader.h"

void vector_to_array(double array[], vector<double> &vector) {
  for (int i=0; i<vector.size(); i++) {
    array[i] = vector[i];
  }
}

void zero_array(double array[], int size) {
  for (int i=0; i<size; i++) {
    array[i] = 0;
  }
}

int main() {

  string filenames[] = {
    "data/J_L_inclusive_bk.txt",
    "data/J_L_inclusive_bfkl.txt"
  };
  int filecount = size(filenames);

  vector<double> Q2[filecount];
  vector<vector<double>> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    vector<double> raw_Q2, raw_x, raw_sigma, raw_sigma_error;
    read_sigma_file(filenames[i], raw_Q2, raw_x, raw_sigma, raw_sigma_error);
    //vector<double> new_Q2, 
    split_by_Q2(Q2[i], x[i], sigma[i], sigma_error[i], raw_Q2, raw_x, raw_sigma, raw_sigma_error);
  }

  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("Longitudinal inclusive cross section comparison between BK and BFKL");

  for (int i=0; i<filecount; i++) {
    for (int j=0; j<Q2[i].size(); j++) {

      double x_arr[x[i][j].size()];
      vector_to_array(x_arr, x[i][j]);

      double sigma_arr[sigma[i][j].size()];
      vector_to_array(sigma_arr, sigma[i][j]);

      double sigma_error_arr[sigma_error[i][j].size()];
      vector_to_array(sigma_error_arr, sigma_error[i][j]);

      double x_error_arr[x[i][j].size()];
      zero_array(x_error_arr, sigma_error[i][j].size());

      TGraphErrors* subgraph = new TGraphErrors(x[i][j].size(), x_arr, sigma_arr, x_error_arr, sigma_error_arr);
      stringstream Q2_stream;
      Q2_stream << fixed << setprecision(0) << Q2[i][j];
      TString subgraph_name;

      if (i == 0) {
        subgraph_name = "BK, Q^{2}=" + Q2_stream.str();
      } else if(i == 1) {
        subgraph_name = "BFKL, Q^{2}=" + Q2_stream.str();
        subgraph->SetLineStyle(2);
      }
      
      subgraph->SetTitle(subgraph_name);
      comparison_graphs->Add(subgraph);
    }
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  comparison_graphs->GetXaxis()->SetTitle("x");
  comparison_graphs->GetYaxis()->SetTitle("#sigma (GeV^-2)");

  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/L_BFKL_BK_comp.pdf");
  
  return 0;
}
