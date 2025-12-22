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

  vector<string> numerator_filenames = {"bk_version", "bfkl_version"};
  vector<string> denominator_filenames = {"bk_version", "bfkl_version"};

  TString title, outfile_name;


  //const string sigma_type = "Pb BK";
  //numerator_filename = "data/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt";
  //denominator_filename = "data/J_LHC_T_inclusive_c_bk_Pb.txt";
  title = "Test";
  outfile_name = "figures/Test.pdf";

  TMultiGraph* ratios_graph = new TMultiGraph();
  ratios_graph->SetTitle(title);
  ratios_graph->GetXaxis()->SetTitle("W (GeV)");
  
  for (int i=0;i<numerator_filenames.size();i++) {
    vector<double> numerator_W, numerator_sigma, numerator_sigma_error, denominator_W, denominator_sigma, denominator_sigma_error;
    read_LHC_sigma_file(numerator_filenames[i], numerator_W, numerator_sigma, numerator_sigma_error);
    read_LHC_sigma_file(denominator_filenames[i], denominator_W, denominator_sigma, denominator_sigma_error);
    double ratio[numerator_W.size()];
    double W[numerator_W.size()];
    for (long unsigned int j=0; j<numerator_W.size(); j++) {
      ratio[j] = numerator_sigma[j]/denominator_sigma[j];
      W[j] = numerator_W[j];
    }
    TGraph* ratio_graph = new TGraph(numerator_W.size(), W, ratio);
    ratio_graph->SetTitle("asd");
    //ratio_graph->SetLineColor(2);
    ratios_graph->Add(ratio_graph);
  }

  TCanvas* ratios_canvas = new TCanvas("ratios_canvas", "", 1000, 600);
  ratios_graph->Draw("AC");

  gPad->SetLogx();

  ratios_canvas->Print(outfile_name);

  return 0;
}
