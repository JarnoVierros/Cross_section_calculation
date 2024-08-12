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

  bool diffractive = true;

  if (diffractive) {
    numerator_filename = "data/diff_LHC_T_sigma_W_bk_Pb.txt";
    denominator_filename = "data/diff_LHC_T_sigma_W_bk_p.txt";
    title = "Diffractive BK Pb p cross section ration; W (GeV); #sigma_{Pb}/#sigma_{p}";
    outfile_name = "figures/diffractive_bk_Pb_p_ratio.pdf";
  } else {
    numerator_filename = "data/J_LHC_T_inclusive_bk_Pb.txt";
    denominator_filename = "data/J_LHC_T_inclusive_bk_p.txt";
    title = "Inclusive BK Pb p cross section ration; W (GeV); #sigma_{Pb}/#sigma_{p}";
    outfile_name = "figures/inclusive_bk_Pb_p_ratio.pdf";
  }
  
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

  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graph->Draw("AL");

  gPad->SetLogx();

  comparison_canvas->Print(outfile_name);
  
  return 0;
}
