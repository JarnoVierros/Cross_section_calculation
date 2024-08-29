#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <sstream>
using namespace std;

#include "cross_section_file_reader.h"

const double c = 1;//0.18;
//const double A = pow(208, 4.0/3);
static double A;

int main() {

  string polarization = "L";
  bool inclusive = true;
  
  string numerator_filename, denominator_filename;
  
  if (inclusive) {
    numerator_filename = "data/J_"+polarization+"_inclusive_Q2_bk_Pb.txt";
    denominator_filename = "data/J_"+polarization+"_inclusive_Q2_bk_p.txt";
    A = pow(208, 1);
  } else {
    numerator_filename = "data/diff_Q2_"+polarization+"_sigma_bk_Pb.txt";
    denominator_filename = "data/diff_Q2_"+polarization+"_sigma_bk_p.txt";
    A = pow(208, 4.0/3);
  }


  vector<double> initial_numerator_Q2, initial_denominator_Q2, initial_numerator_W, initial_numerator_sigma, initial_numerator_sigma_error, initial_denominator_W, initial_denominator_sigma, initial_denominator_sigma_error;

  read_sigma_file(numerator_filename, initial_numerator_W, initial_numerator_Q2, initial_numerator_sigma, initial_numerator_sigma_error);
  read_sigma_file(denominator_filename, initial_denominator_W, initial_denominator_Q2, initial_denominator_sigma, initial_denominator_sigma_error);

  vector<double> numerator_W, denominator_W;
  vector<vector<double>> numerator_Q2, numerator_sigma, numerator_sigma_error, denominator_Q2, denominator_sigma, denominator_sigma_error;

  split_by_Q2(numerator_W, numerator_Q2, numerator_sigma, numerator_sigma_error, initial_numerator_W, initial_numerator_Q2, initial_numerator_sigma, initial_numerator_sigma_error);
  split_by_Q2(denominator_W, denominator_Q2, denominator_sigma, denominator_sigma_error, initial_denominator_W, initial_denominator_Q2, initial_denominator_sigma, initial_denominator_sigma_error);

  TMultiGraph* comparison_graphs = new TMultiGraph();
  TString title;
  if (inclusive) {
    if (polarization == "L") {
      title = "Nuclear suppression ratio for longitudinal inclusive cross section; Q^{2} (GeV^{2})";
    } else {
      title = "Nuclear suppression ratio for transverse inclusive cross section; Q^{2} (GeV^{2})";
    }
  } else {
    if (polarization == "L") {
      title = "Nuclear suppression ratio for longitudinal diffractive cross section; Q^{2} (GeV^{2})";
    } else {
      title = "Nuclear suppression ratio for transverse diffractive cross section; Q^{2} (GeV^{2})";
    }
  }
  comparison_graphs->SetTitle(title);
  for (long unsigned int i=0; i < numerator_W.size(); i++) {
    double ratio[numerator_Q2[i].size()];
    double x[numerator_Q2[i].size()];
    for (long unsigned int j=0; j<numerator_Q2[i].size(); j++) {
      ratio[j] = numerator_sigma[i][j]/(c*A*denominator_sigma[i][j]);
      x[j] = numerator_Q2[i][j];
    }
    TGraph* subgraph = new TGraph(numerator_Q2[i].size(), x, ratio);
    stringstream W_stream;
    W_stream << fixed << setprecision(0) << numerator_W[i];
    TString subgraph_name = "W=" + W_stream.str() + " GeV";
    subgraph->SetTitle(subgraph_name);
    comparison_graphs->Add(subgraph, "C");
  }

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graphs->Draw("A PMC PLC");

  gPad->SetLogx();
  if (true) {
    comparison_canvas->BuildLegend(0.75, 0.2, 0.9, 0.55);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  if (inclusive) {
    title = "figures/"+polarization+"_inclusive_nuclear_suppression_ratio.pdf";
  } else {
    title = "figures/"+polarization+"_diffractive_nuclear_suppression_ratio.pdf";
  }
  comparison_canvas->Print(title);
  
  return 0;
}
