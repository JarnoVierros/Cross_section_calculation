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

  string numerator_filename, denominator_filename;
  TString title, outfile_name;

  const string sigma_type = "Pb BK";
  numerator_filename = "data/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt";
  denominator_filename = "data/J_LHC_T_inclusive_c_bk_Pb.txt";
  title = "Diffractive "+sigma_type+" cross section divided by inclusive "+sigma_type+" cross section";
  outfile_name = "figures/Pb_BK_mixed_dipamp_diff_inc_ratio.pdf";
  
  vector<double> numerator_W, numerator_sigma, numerator_sigma_error, denominator_W, denominator_sigma, denominator_sigma_error;

  read_LHC_sigma_file(numerator_filename, numerator_W, numerator_sigma, numerator_sigma_error);
  read_LHC_sigma_file(denominator_filename, denominator_W, denominator_sigma, denominator_sigma_error);

  
  double ratio[numerator_W.size()];
  double W[numerator_W.size()];
  for (long unsigned int j=0; j<numerator_W.size(); j++) {
    ratio[j] = numerator_sigma[j]/denominator_sigma[j];
    W[j] = numerator_W[j];
  }

  TGraph* comparison_graph = new TGraph(numerator_W.size(), W, ratio);
  comparison_graph->SetTitle(title);
  comparison_graph->GetXaxis()->SetTitle("W (GeV)");
  //comparison_graph->GetYaxis()->SetTitle("");

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graph->Draw("AC");

  gPad->SetLogx();

  comparison_canvas->Print(outfile_name);
  
  return 0;
}
