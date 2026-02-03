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

  vector<string> numerator_filenames = {
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_c_bfkl_Pb_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_c_bfkl_p_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_c_bk_p_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_b_bfkl_Pb_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_b_bk_Pb_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_b_bfkl_p_diffraction.txt",
    "/home/jarno/Cross_section_calculation/output/diff_LHC_T_sigma_W_b_bk_p_diffraction.txt",
  };
  vector<string> denominator_filenames = {
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_c_bfkl_Pb.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_c_bk_Pb.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_c_bfkl_p.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_c_bk_p.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_b_bfkl_Pb.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_b_bk_Pb.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_b_bfkl_p.txt",
    "/home/jarno/Cross_section_calculation/output/J_LHC_T_inclusive_b_bk_p.txt",
  };
  vector<TString> graph_titles = {"PbPb#rightarrowc#bar{c} bfkl", "PbPb#rightarrowc#bar{c} bk", "pp#rightarrowc#bar{c} bfkl", "pp#rightarrowc#bar{c} bk", "PbPb#rightarrowb#bar{b} bfkl", "PbPb#rightarrowb#bar{b} bk", "pp#rightarrowb#bar{b} bfkl", "pp#rightarrowb#bar{b} bk"};
  vector<int> line_colors = {2, 2, 4, 4, 2, 2, 4, 4};
  vector<int> line_styles = {2, 1, 5, 7, 2, 1, 5, 7};

  TString title, outfile_name;


  //const string sigma_type = "Pb BK";
  //numerator_filename = "data/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt";
  //denominator_filename = "data/J_LHC_T_inclusive_c_bk_Pb.txt";
  //title = "Diffractive c#bar{c} cross section divided by the inclusive cross section";
  outfile_name = "figures/diffractive_inclusive_ccbar_ratio.pdf";

  TCanvas* ratios_canvas = new TCanvas("ratios_canvas", "", 1000, 2*600);
  ratios_canvas->Divide(1, 2);

  ratios_canvas->cd(1);



  TMultiGraph* ratios_graph = new TMultiGraph();
  //ratios_graph->SetTitle(title);
  
  for (int i=0;i<4;i++) {
    vector<double> numerator_W, numerator_sigma, numerator_sigma_error, denominator_W, denominator_sigma, denominator_sigma_error;
    read_LHC_sigma_file(numerator_filenames[i], numerator_W, numerator_sigma, numerator_sigma_error);
    read_LHC_sigma_file(denominator_filenames[i], denominator_W, denominator_sigma, denominator_sigma_error);
    double ratio[numerator_W.size()];
    double W[numerator_W.size()];
    int data_size = 0;
    for (long unsigned int j=0; j<numerator_W.size(); j++) {
      double new_ratio = numerator_sigma[j]/denominator_sigma[j];
      if (new_ratio > 1) {
        break;
      }
      if (numerator_W[j] > 2.1e3) {
        continue;
      }
      ratio[j] = numerator_sigma[j]/denominator_sigma[j];
      W[j] = numerator_W[j];
      data_size++;
      //cout << W[j] << ", " << ratio[j] << endl;
    }
    TGraph* ratio_graph = new TGraph(data_size, W, ratio);
    ratio_graph->SetTitle(graph_titles[i]);
    ratio_graph->SetLineColor(line_colors[i]);
    ratio_graph->SetLineStyle(line_styles[i]);
    //ratio_graph->GetXaxis()->SetRangeUser(9e1, 2e4);
    ratios_graph->Add(ratio_graph);
  }
  //ratios_graph->GetXaxis()->SetRangeUser(9e1, 2e4);

  ratios_graph->GetYaxis()->SetTitle("#frac{diffractive}{inclusive}");
  gPad->SetLogx();

  ratios_graph->Draw("AC");
  ratios_graph->GetXaxis()->SetTitle("W (GeV)");
  ratios_graph->GetXaxis()->SetLimits(2.6e1, 2.1e3);

  ratios_canvas->cd(1)->BuildLegend(0.15, 0.6, 0.5, 0.9);


  ratios_canvas->cd(2);
  TMultiGraph* ratios_graph_2 = new TMultiGraph();
  //ratios_graph->SetTitle(title);
  
  for (int i=4;i<8;i++) {
    vector<double> numerator_W, numerator_sigma, numerator_sigma_error, denominator_W, denominator_sigma, denominator_sigma_error;
    read_LHC_sigma_file(numerator_filenames[i], numerator_W, numerator_sigma, numerator_sigma_error);
    read_LHC_sigma_file(denominator_filenames[i], denominator_W, denominator_sigma, denominator_sigma_error);
    double ratio[numerator_W.size()];
    double W[numerator_W.size()];
    int data_size = 0;
    for (long unsigned int j=0; j<numerator_W.size(); j++) {
      double new_ratio = numerator_sigma[j]/denominator_sigma[j];
      if (new_ratio > 1) {
        break;
      }
      if (numerator_W[j] > 2.1e3) {
        continue;
      }
      ratio[j] = numerator_sigma[j]/denominator_sigma[j];
      W[j] = numerator_W[j];
      data_size++;
      //cout << W[j] << ", " << ratio[j] << endl;
    }
    TGraph* ratio_graph = new TGraph(data_size, W, ratio);
    ratio_graph->SetTitle(graph_titles[i]);
    ratio_graph->SetLineColor(line_colors[i]);
    ratio_graph->SetLineStyle(line_styles[i]);
    ratios_graph_2->Add(ratio_graph);
  }

  gPad->SetLogx();
  //ratios_graph_2->GetYaxis()->SetRangeUser(0, 01);
  ratios_graph_2->Draw("AC");
  ratios_graph_2->GetXaxis()->SetTitle("W (GeV)");
  ratios_graph_2->GetYaxis()->SetTitle("#frac{diffractive}{inclusive}");
  ratios_graph_2->GetXaxis()->SetLimits(2.6e1, 2.1e3);
  
  ratios_canvas->cd(2)->BuildLegend(0.15, 0.6, 0.5, 0.9);

  ratios_canvas->Print(outfile_name);

  return 0;
}
