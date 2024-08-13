#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLatex.h"

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
    "data/GBW_diff_LHC_T_sigma_W_bk_Pb.txt",
    "data/GBW_var_change_diff_LHC_T_sigma_W_bk_Pb.txt"
  };
  TString titles[] = {
    "Original",
    "Variable change"
  };
  int line_styles[] = {
    1,
    2
  };

  int filecount = size(filenames);

  vector<double> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    cout << "reading " << filenames[i] << endl;
    read_LHC_sigma_file(filenames[i], x[i], sigma[i], sigma_error[i]);
  }

  TMultiGraph* multigraph = new TMultiGraph();
  multigraph->SetTitle("Multigraph");


  for (int j=0; j<size(filenames); j++) {
    double x_array[x[j].size()];
    double sigma_array[x[j].size()];
    double x_error_array[x[j].size()];
    double sigma_error_array[x[j].size()];
    for (int i=0; i<x[j].size(); i++) {
      x_array[i] = x[j][i];
      sigma_array[i] = sigma[j][i];
      x_error_array[i] = 0;
      sigma_error_array[i] = sigma_error[j][i];
    }
    TGraphErrors* subgraph = new TGraphErrors(x[j].size(), x_array, sigma_array, x_error_array, sigma_error_array);
    subgraph->SetTitle(titles[j]);
    subgraph->SetLineStyle(line_styles[j]);
    multigraph->Add(subgraph);
  }

  TCanvas* canvas1 = new TCanvas("canvas1", "", 1000, 600);
  multigraph->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();

  if (true) {
    canvas1->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    canvas1->BuildLegend(0.2, 0.7, 0.35, 0.9);
  }

  canvas1->Print("figures/LHC_multiplot.pdf");
  
  return 0;
}
