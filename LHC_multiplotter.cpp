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
    "data/J_LHC_T_inclusive_bk_Pb.txt",
    "data/J_LHC_L_inclusive_bk_Pb.txt",
    "data/J_LHC_T_inclusive_bfkl_Pb.txt",
    "data/J_LHC_L_inclusive_bfkl_Pb.txt",
    "data/J_LHC_T_inclusive_bk_Pb_xpom.txt",
    "data/J_LHC_L_inclusive_bk_Pb_xpom.txt",
    "data/J_LHC_T_inclusive_bfkl_Pb_xpom.txt",
    "data/J_LHC_L_inclusive_bfkl_Pb_xpom.txt",
    "data/J_LHC_T_inclusive_bk.txt",
    "data/J_LHC_L_inclusive_bk.txt",
    "data/J_LHC_T_inclusive_bfkl.txt",
    "data/J_LHC_L_inclusive_bfkl.txt",
    "data/J_LHC_T_inclusive_bk_Pb_xpom_Q2=0.txt",
    "data/J_LHC_L_inclusive_bk_Pb_xpom_Q2=0.txt",
    "data/J_LHC_T_inclusive_bfkl_Pb_xpom_Q2=0.txt",
    "data/J_LHC_L_inclusive_bfkl_Pb_xpom_Q2=0.txt",
    "data/J_LHC_T_inclusive_bk_p_xpom_Q2=0.txt",
    "data/J_LHC_L_inclusive_bk_p_xpom_Q2=0.txt",
    "data/J_LHC_T_inclusive_bfkl_p_xpom_Q2=0.txt",
    "data/J_LHC_L_inclusive_bfkl_p_xpom_Q2=0.txt",
  };
  int filecount = size(filenames);

  vector<double> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    cout << "reading " << filenames[i] << endl;
    read_LHC_sigma_file(filenames[i], x[i], sigma[i], sigma_error[i]);
  }

  double BK_sigma_tot_Pb[x[0].size()], BFKL_sigma_tot_Pb[x[0].size()], BK_sigma_tot_Pb_xpom[x[0].size()], BFKL_sigma_tot_Pb_xpom[x[0].size()], BK_sigma_tot_p[x[0].size()], BFKL_sigma_tot_p[x[0].size()], BK_sigma_tot_Pb_Q20[x[0].size()], BFKL_sigma_tot_Pb_Q20[x[0].size()], BK_sigma_tot_p_Q20[x[0].size()], BFKL_sigma_tot_p_Q20[x[0].size()];

  for (int i=0; i<x[0].size(); i++) {
    BK_sigma_tot_Pb[i] = sigma[0][i] + sigma[1][i];
    BFKL_sigma_tot_Pb[i] = sigma[2][i] + sigma[3][i];
    BK_sigma_tot_Pb_xpom[i] = sigma[4][i] + sigma[5][i];
    BFKL_sigma_tot_Pb_xpom[i] = sigma[6][i] + sigma[7][i];
    BK_sigma_tot_p[i] = sigma[8][i] + sigma[9][i];
    BFKL_sigma_tot_p[i] = sigma[10][i] + sigma[11][i];
    BK_sigma_tot_Pb_Q20[i] = sigma[12][i] + sigma[13][i];
    BFKL_sigma_tot_Pb_Q20[i] = sigma[14][i] + sigma[15][i];
    BK_sigma_tot_p_Q20[i] = sigma[16][i] + sigma[17][i];
    BFKL_sigma_tot_p_Q20[i] = sigma[18][i] + sigma[19][i];
  }

  TMultiGraph* comparison_graph = new TMultiGraph();
  comparison_graph->SetTitle("Inclusive cross section comparison between BK and BFKL");

  double BK_Pb_x_arr[x[0].size()];
  vector_to_array(BK_Pb_x_arr, x[0]);

  TGraph* BK_graph = new TGraph(x[0].size(), BK_Pb_x_arr, BK_sigma_tot_Pb);
  BK_graph->SetTitle("Pb BK");
  BK_graph->SetLineColor(1);
  //comparison_graph->Add(BK_graph);

  double BFKL_Pb_x_arr[x[2].size()];
  vector_to_array(BFKL_Pb_x_arr, x[2]);

  TGraph* BFKL_graph = new TGraph(x[2].size(), BFKL_Pb_x_arr, BFKL_sigma_tot_Pb);
  BFKL_graph->SetTitle("Pb BFKL");
  BFKL_graph->SetLineColor(1);
  BFKL_graph->SetLineStyle(2);
  //comparison_graph->Add(BFKL_graph);

  double BK_Pb_xpom_x_arr[x[0].size()];
  vector_to_array(BK_Pb_xpom_x_arr, x[0]);

  TGraph* BK_xpom_graph = new TGraph(x[0].size(), BK_Pb_xpom_x_arr, BK_sigma_tot_Pb_xpom);
  BK_xpom_graph->SetTitle("Pb BK xpom");
  BK_xpom_graph->SetLineColor(2);
  //comparison_graph->Add(BK_xpom_graph);

  double BFKL_Pb_xpom_x_arr[x[2].size()];
  vector_to_array(BFKL_Pb_xpom_x_arr, x[2]);

  TGraph* BFKL_xpom_graph = new TGraph(x[2].size(), BFKL_Pb_xpom_x_arr, BFKL_sigma_tot_Pb_xpom);
  BFKL_xpom_graph->SetTitle("Pb BFKL xpom");
  BFKL_xpom_graph->SetLineColor(2);
  BFKL_xpom_graph->SetLineStyle(2);
  //comparison_graph->Add(BFKL_xpom_graph);

  double BK_p_x_arr[x[0].size()];
  vector_to_array(BK_p_x_arr, x[0]);

  TGraph* BK_p_graph = new TGraph(x[0].size(), BK_p_x_arr, BK_sigma_tot_p);
  BK_p_graph->SetTitle("p BK");
  BK_p_graph->SetLineColor(3);
  //comparison_graph->Add(BK_p_graph);

  double BFKL_p_x_arr[x[2].size()];
  vector_to_array(BFKL_p_x_arr, x[2]);

  TGraph* BFKL_p_graph = new TGraph(x[2].size(), BFKL_p_x_arr, BFKL_sigma_tot_p);
  BFKL_p_graph->SetTitle("p BFKL");
  BFKL_p_graph->SetLineColor(3);
  BFKL_p_graph->SetLineStyle(2);
  //comparison_graph->Add(BFKL_p_graph);

  double BK_Pb_Q20_x_arr[x[0].size()];
  vector_to_array(BK_Pb_Q20_x_arr, x[0]);

  TGraph* BK_Q20_graph = new TGraph(x[0].size(), BK_Pb_Q20_x_arr, BK_sigma_tot_Pb_Q20);
  BK_Q20_graph->SetTitle("Pb BK");
  BK_Q20_graph->SetLineColor(4);
  comparison_graph->Add(BK_Q20_graph);

  double BFKL_Pb_Q20_x_arr[x[2].size()];
  vector_to_array(BFKL_Pb_Q20_x_arr, x[2]);

  TGraph* BFKL_Q20_graph = new TGraph(x[2].size(), BFKL_Pb_Q20_x_arr, BFKL_sigma_tot_Pb_Q20);
  BFKL_Q20_graph->SetTitle("Pb BFKL");
  BFKL_Q20_graph->SetLineColor(4);
  BFKL_Q20_graph->SetLineStyle(2);
  comparison_graph->Add(BFKL_Q20_graph);

  double BK_p_Q20_x_arr[x[0].size()];
  vector_to_array(BK_p_Q20_x_arr, x[0]);

  TGraph* BK_p_Q20_graph = new TGraph(x[0].size(), BK_p_Q20_x_arr, BK_sigma_tot_p_Q20);
  BK_p_Q20_graph->SetTitle("p BK");
  BK_p_Q20_graph->SetLineColor(2);
  comparison_graph->Add(BK_p_Q20_graph);

  double BFKL_p_Q20_x_arr[x[2].size()];
  vector_to_array(BFKL_p_Q20_x_arr, x[2]);

  TGraph* BFKL_p_Q20_graph = new TGraph(x[2].size(), BFKL_p_Q20_x_arr, BFKL_sigma_tot_p_Q20);
  BFKL_p_Q20_graph->SetTitle("p BFKL");
  BFKL_p_Q20_graph->SetLineColor(2);
  BFKL_p_Q20_graph->SetLineStyle(2);
  comparison_graph->Add(BFKL_p_Q20_graph);


  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graph->Draw("A");

  gPad->SetLogx();
  gPad->SetLogy();
  comparison_graph->GetXaxis()->SetTitle("W (GeV)");
  comparison_graph->GetYaxis()->SetTitle("#sigma_{L+T} (GeV^{-2})");

  if (false) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.55, 0.35, 0.9);
  }

  TLatex* Q2_text = new TLatex(2e3, 1e-2, "Q^{2} = 0 GeV^{2}");
  Q2_text->Draw("same");

  comparison_canvas->Print("figures/LHC_BFKL_BK_comp_all.pdf");
  
  return 0;
}
