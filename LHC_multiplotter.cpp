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
    "data/J_LHC_T_exclusive_bk_Pb.txt",
    "data/J_LHC_T_exclusive_bfkl_Pb.txt",
    "data/J_LHC_T_exclusive_bk_p.txt",
    "data/J_LHC_T_exclusive_bfkl_p.txt",
  };
  int filecount = size(filenames);

  vector<double> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    cout << "reading " << filenames[i] << endl;
    read_LHC_sigma_file(filenames[i], x[i], sigma[i], sigma_error[i]);
  }

  double BK_sigma_tot_Pb_Q20[x[0].size()], BFKL_sigma_tot_Pb_Q20[x[0].size()], BK_sigma_tot_p_Q20[x[0].size()], BFKL_sigma_tot_p_Q20[x[0].size()];

  for (int i=0; i<x[0].size(); i++) {
    BK_sigma_tot_Pb_Q20[i] = sigma[0][i];
    BFKL_sigma_tot_Pb_Q20[i] = sigma[1][i];
    BK_sigma_tot_p_Q20[i] = sigma[2][i];
    BFKL_sigma_tot_p_Q20[i] = sigma[3][i];
  }

  TMultiGraph* comparison_graph = new TMultiGraph();
  comparison_graph->SetTitle("Exclusive cross section comparison between BK and BFKL");

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
  //Q2_text->Draw("same");

  comparison_canvas->Print("figures/LHC_exclusive_prediction.pdf");
  
  return 0;
}
