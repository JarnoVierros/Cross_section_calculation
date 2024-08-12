#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
using namespace std;

#include "cross_section_file_reader.h"


int main() {

  string numerator_filename("/home/jarno/Cross_section_calculation/archive/data/LHC/diff_LHC_T_sigma_W_bk_Pb.txt");
  string denominator_filename("/home/jarno/Cross_section_calculation/archive/data/LHC/diff_LHC_T_sigma_W_bk_p.txt");

  vector<double> initial_numerator_x, initial_numerator_sigma, initial_numerator_sigma_error, initial_denominator_x, initial_denominator_sigma, initial_denominator_sigma_error;

  read_LHC_sigma_file(numerator_filename, initial_numerator_x, initial_numerator_sigma, initial_numerator_sigma_error);
  read_LHC_sigma_file(denominator_filename, initial_denominator_x, initial_denominator_sigma, initial_denominator_sigma_error);

  vector<double> numerator_Q2, denominator_Q2;
  vector<vector<double>> numerator_x, numerator_sigma, numerator_sigma_error, 

  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("bk Pb p ratio");
  for (long unsigned int i=0; i < Q2_values.size(); i++) {
    double ratio[x_values[i].size()];
    double x[x_values[i].size()];
    for (long unsigned int j=0; j<x_values[i].size(); j++) {
      ratio[j] = numerator[i][j]/denominator[i][j];
      x[j] = x_values[i][j];
    }
    TGraph* subgraph = new TGraph(x_values[i].size(), x, ratio);
    TString subgraph_name = "Q^{2}=" + to_string(Q2_values[i]);
    subgraph->SetTitle(subgraph_name);
    comparison_graphs->Add(subgraph);
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/kb_Pb_p_ratio.pdf");
  
  return 0;
}
