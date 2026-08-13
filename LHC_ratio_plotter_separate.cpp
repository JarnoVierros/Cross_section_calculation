#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TLatex.h"

#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <sstream>
using namespace std;

#include "cross_section_file_reader.h"

const string main_dir = "/home/jarno/Cross_section_calculation";

int main() {

  double linewidth = 3;

  vector<string> numerator_filenames = {
    main_dir + "/output/diff_LHC_T_sigma_W_c_bfkl_Pb_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_c_bfkl_p_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_c_bk_p_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_b_bfkl_Pb_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_b_bk_Pb_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_b_bfkl_p_diffraction.txt",
    main_dir + "/output/diff_LHC_T_sigma_W_b_bk_p_diffraction.txt",
  };
  vector<string> denominator_filenames = {
    main_dir + "/output/J_LHC_T_inclusive_c_bfkl_Pb.txt",
    main_dir + "/output/J_LHC_T_inclusive_c_bk_Pb.txt",
    main_dir + "/output/J_LHC_T_inclusive_c_bfkl_p.txt",
    main_dir + "/output/J_LHC_T_inclusive_c_bk_p.txt",
    main_dir + "/output/J_LHC_T_inclusive_b_bfkl_Pb.txt",
    main_dir + "/output/J_LHC_T_inclusive_b_bk_Pb.txt",
    main_dir + "/output/J_LHC_T_inclusive_b_bfkl_p.txt",
    main_dir + "/output/J_LHC_T_inclusive_b_bk_p.txt",
  };
  vector<TString> graph_titles = {"Pb bfkl", "Pb bk", "p bfkl", "p bk", "Pb bfkl", "Pb bk", "p bfkl", "p bk"};
  vector<int> line_colors = {2, 2, 4, 4, 2, 2, 4, 4};
  vector<int> line_styles = {2, 1, 5, 7, 2, 1, 5, 7};

  TString title, outfile_name_1, outfile_name_2;


  //const string sigma_type = "Pb BK";
  //numerator_filename = "data/diff_LHC_T_sigma_W_c_bk_Pb_diffraction.txt";
  //denominator_filename = "data/J_LHC_T_inclusive_c_bk_Pb.txt";
  //title = "Diffractive c#bar{c} cross section divided by the inclusive cross section";
  outfile_name_1 = "figures/diffractive_inclusive_ccbar_ratio.ps";
  outfile_name_2 = "figures/diffractive_inclusive_bbbar_ratio.ps";


  TCanvas* top_canvas = new TCanvas("top_canvas", "", 0.5*1.5*2000, 1.5*800);
  top_canvas->SetBottomMargin(0.14);
  top_canvas->SetLeftMargin(0.2);

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
    ratio_graph->SetLineWidth(linewidth);
    ratio_graph->SetMarkerSize(0);
    //ratio_graph->GetXaxis()->SetRangeUser(9e1, 2e4);
    ratios_graph->Add(ratio_graph);
  }
  //ratios_graph->GetXaxis()->SetRangeUser(9e1, 2e4);

  ratios_graph->GetYaxis()->SetTitle("#frac{#it{#sigma}_{diffractive}}{#it{#sigma}_{inclusive}}");
  gPad->SetLogx();

  ratios_graph->Draw("AC");
  ratios_graph->GetXaxis()->SetTitle("#it{W} [GeV]");
  ratios_graph->GetXaxis()->SetLimits(2.6e1, 2.1e3);
  ratios_graph->GetYaxis()->SetTitleSize(0.05);
  ratios_graph->GetXaxis()->SetTitleSize(0.05);
  ratios_graph->GetYaxis()->SetLabelSize(0.05);
  ratios_graph->GetXaxis()->SetLabelSize(0.05);

  ratios_graph->GetXaxis()->SetTitleOffset(1.3);

  top_canvas->BuildLegend(0.25, 0.6, 0.45, 0.9);

  TLatex* Q2_text = new TLatex(3.5e1, 0.4, "Q^{2}=0\\text{ GeV}^{2}");
  Q2_text->SetTextSize(0.06);
  Q2_text->Draw("same");
  TLatex* gamma = new TLatex(3.5e1, 0.5, "#gamma");
  gamma->SetTextSize(0.06);
  gamma->Draw("same");
  TLatex* process_text = new TLatex(3.5e1, 0.5, "\\ \\, +\\text{Pb}/\\text{p} \\rightarrow c\\bar{c}+X");
  process_text->SetTextSize(0.06);
  process_text->Draw("same");

  top_canvas->Print(outfile_name_1);

  
  TCanvas* bottom_canvas = new TCanvas("bottom_canvas", "", 0.5*1.5*2000, 1.5*800);
  bottom_canvas->SetBottomMargin(0.14);
  bottom_canvas->SetLeftMargin(0.2);

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
    ratio_graph->SetLineWidth(linewidth);
    ratio_graph->SetMarkerSize(0);
    ratios_graph_2->Add(ratio_graph);
  }

  gPad->SetLogx();
  //ratios_graph_2->GetYaxis()->SetRangeUser(0, 01);
  ratios_graph_2->Draw("AC");
  ratios_graph_2->GetXaxis()->SetTitle("#it{W} [GeV]");
  ratios_graph_2->GetYaxis()->SetTitle("#frac{#it{#sigma}_{diffractive}}{#it{#sigma}_{inclusive}}");
  ratios_graph_2->GetXaxis()->SetLimits(2.6e1, 2.1e3);
  ratios_graph_2->GetYaxis()->SetTitleSize(0.05);
  ratios_graph_2->GetXaxis()->SetTitleSize(0.05);
  ratios_graph_2->GetYaxis()->SetLabelSize(0.05);
  ratios_graph_2->GetXaxis()->SetLabelSize(0.05);

  ratios_graph_2->GetXaxis()->SetTitleOffset(1.3);

  bottom_canvas->BuildLegend(0.25, 0.6, 0.45, 0.9);

  Q2_text = new TLatex(3.5e1, 0.075, "Q^{2}=0\\text{ GeV}^{2}");
  Q2_text->SetTextSize(0.06);
  Q2_text->Draw("same");
  gamma = new TLatex(3.5e1, 0.09, "#gamma");
  gamma->SetTextSize(0.06);
  gamma->Draw("same");
  process_text = new TLatex(3.5e1, 0.09, "\\ \\, +\\text{Pb}/\\text{p} \\rightarrow b\\bar{b}+X");
  process_text->SetTextSize(0.06);
  process_text->Draw("same");

  bottom_canvas->Print(outfile_name_2);

  return 0;
}
