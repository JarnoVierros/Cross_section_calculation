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

const double c = 1.0;
const double A = pow(208, 4.0/3); 

int main() {

  string numerator_filename("/home/jarno/Cross_section_calculation/data/J_LHC_T_inclusive_bk_p.txt");
  string denominator_filename("/home/jarno/Cross_section_calculation/archive/data/LHC/inclusive/J_LHC_T_inclusive_bk_p.txt");

  vector<double> initial_numerator_Q2, initial_denominator_Q2, initial_numerator_W, initial_numerator_sigma, initial_numerator_sigma_error, initial_denominator_W, initial_denominator_sigma, initial_denominator_sigma_error;

  read_sigma_file(numerator_filename, initial_numerator_W, initial_numerator_Q2, initial_numerator_sigma, initial_numerator_sigma_error);
  read_sigma_file(denominator_filename, initial_denominator_W, initial_denominator_Q2, initial_denominator_sigma, initial_denominator_sigma_error);

  vector<double> numerator_W, denominator_W;
  vector<vector<double>> numerator_Q2, numerator_sigma, numerator_sigma_error, denominator_Q2, denominator_sigma, denominator_sigma_error;

  split_by_Q2(numerator_W, numerator_Q2, numerator_sigma, numerator_sigma_error, initial_numerator_Q2, initial_numerator_W, initial_numerator_sigma, initial_numerator_sigma_error);
  split_by_Q2(denominator_W, denominator_Q2, denominator_sigma, denominator_sigma_error, initial_denominator_Q2, initial_denominator_W, initial_denominator_sigma, initial_denominator_sigma_error);

  TMultiGraph* comparison_graphs = new TMultiGraph();
  comparison_graphs->SetTitle("Ratio between different integration methods in longitudinal case");
  for (long unsigned int i=0; i < numerator_W.size(); i++) {
    double ratio[numerator_Q2[i].size()];
    double x[numerator_Q2[i].size()];
    for (long unsigned int j=0; j<numerator_Q2[i].size(); j++) {
      ratio[j] = numerator_sigma[i][j]/(c*A*denominator_sigma[i][j]);
      x[j] = numerator_Q2[i][j];
    }
    TGraph* subgraph = new TGraph(numerator_Q2[i].size(), x, ratio);
    TString subgraph_name = "W=" + to_string(numerator_W[i]) + " GeV";
    subgraph->SetTitle(subgraph_name);
    comparison_graphs->Add(subgraph, "C");
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  comparison_canvas->Print("figures/nuclear_suppression_ratio.pdf");
  
  return 0;
}
