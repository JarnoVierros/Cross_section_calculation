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

  vector<double> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    read_LHC_sigma_file(filenames[i], x[i], sigma[i], sigma_error[i]);
  }

  double BK_sigma_tot[x[0].size()], BFKL_sigma_tot[x[0].size()];

  for (int i=0; i<x[0].size(); i++) {
    BK_sigma_tot[i] = sigma[0][i] + sigma[1][i];
    BFKL_sigma_tot[i] = sigma[2][i] + sigma[3][i];
  }

  TMultiGraph* comparison_graph = new TMultiGraph();
  comparison_graph->SetTitle("Total inclusive cross section comparison between BK and BFKL");

  double* BK_x_arr;
  vector_to_array(BK_x_arr, x[0]);

  TGraph* BK_graph = new TGraph(x[0].size(), BK_x_arr, BK_sigma_tot);
  BK_graph->SetTitle("BK");
  comparison_graph->Add(BK_graph);

  double* BFKL_x_arr;
  vector_to_array(BFKL_x_arr, x[0]);

  TGraph* BFKL_graph = new TGraph(x[0].size(), BFKL_x_arr, BFKL_sigma_tot);
  BFKL_graph->SetTitle("BFKL");
  comparison_graph->Add(BFKL_graph);

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graph->Draw("A PMC PLC");

  gPad->SetLogx();
  gPad->SetLogy();
  comparison_graph->GetXaxis()->SetTitle("W (GeV)");
  comparison_graph->GetYaxis()->SetTitle("#sigma_{tot} (GeV^-2)");

  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/LHC_BFKL_BK_comp.pdf");
  
  return 0;
}
