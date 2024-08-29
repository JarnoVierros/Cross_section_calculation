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

const string nucleus = "Pb";
const bool diffractive = false;

int main() {

  string filenames[4];
  /*
  string filenames[] = {
    "archive/data/LHC/exclusive/J_LHC_T_exclusive_bk_Pb.txt",
    "archive/data/LHC/exclusive/J_LHC_T_exclusive_bfkl_Pb.txt",
    "archive/data/LHC/exclusive/J_LHC_T_exclusive_bk_p.txt",
    "archive/data/LHC/exclusive/J_LHC_T_exclusive_bfkl_p.txt",
  };
  */

  if (diffractive) {
    filenames[0] = "data/diff_LHC_T_sigma_W_bk_Pb.txt";
    filenames[1] = "data/diff_LHC_T_sigma_W_bfkl_Pb.txt";
    filenames[2] = "data/diff_LHC_T_sigma_W_bk_p.txt";
    filenames[3] = "data/diff_LHC_T_sigma_W_bfkl_p.txt";
  } else {
    /*
    filenames[0] = "archive/data/LHC/inclusive/J_LHC_T_inclusive_bk_Pb.txt";
    filenames[1] = "archive/data/LHC/inclusive/J_LHC_T_inclusive_bfkl_Pb.txt";
    filenames[2] = "archive/data/LHC/inclusive/J_LHC_T_inclusive_bk_p.txt";
    filenames[3] = "archive/data/LHC/inclusive/J_LHC_T_inclusive_bfkl_p.txt";
    */
    filenames[0] = "data/J_LHC_T_inclusive_bk_Pb.txt";
    filenames[1] = "data/J_LHC_T_inclusive_bfkl_Pb.txt";
    filenames[2] = "data/J_LHC_T_inclusive_bk_p.txt";
    filenames[3] = "data/J_LHC_T_inclusive_bfkl_p.txt";
  }

  int filecount = size(filenames);

  vector<double> x[filecount], sigma[filecount], sigma_error[filecount];

  for (int i=0; i<size(filenames); i++) {
    cout << "reading " << filenames[i] << endl;
    read_LHC_sigma_file(filenames[i], x[i], sigma[i], sigma_error[i]);
  }

  double BK_sigma_tot_Pb_Q20[x[0].size()], BFKL_sigma_tot_Pb_Q20[x[0].size()], BK_sigma_tot_p_Q20[x[0].size()], BFKL_sigma_tot_p_Q20[x[0].size()];

  const double GeV_to_nb_conversion = 1e6*0.389379;

  for (int i=0; i<x[0].size(); i++) {
    BK_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*sigma[0][i];
    BFKL_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*sigma[1][i];
    BK_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*sigma[2][i];
    BFKL_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*sigma[3][i];
  }

  TMultiGraph* comparison_graph = new TMultiGraph();
  TString title;
  if (diffractive) {
    if (nucleus == "Pb") {
      title = "c#bar{c} production cross section in diffractive #gamma lead scattering";
    } else {
      title = "c#bar{c} production cross section in diffractive #gamma proton scattering";
    }
  } else {
    if (nucleus == "Pb") {
      title = "c#bar{c} production cross section in inclusive #gamma lead scattering";
    } else {
      title = "c#bar{c} production cross section in inclusive #gamma proton scattering";
    }
  }
  comparison_graph->SetTitle(title);


  if (nucleus == "Pb") {
    double BK_Pb_Q20_x_arr[x[0].size()];
    vector_to_array(BK_Pb_Q20_x_arr, x[0]);

    TGraph* BK_Q20_graph = new TGraph(x[0].size(), BK_Pb_Q20_x_arr, BK_sigma_tot_Pb_Q20);
    BK_Q20_graph->SetTitle("BK"); //Pb BK
    BK_Q20_graph->SetLineColor(4);
    comparison_graph->Add(BK_Q20_graph);

    double BFKL_Pb_Q20_x_arr[x[2].size()];
    vector_to_array(BFKL_Pb_Q20_x_arr, x[2]);

    TGraph* BFKL_Q20_graph = new TGraph(x[2].size(), BFKL_Pb_Q20_x_arr, BFKL_sigma_tot_Pb_Q20);
    BFKL_Q20_graph->SetTitle("BFKL"); //Pb BFKL
    BFKL_Q20_graph->SetLineColor(4);
    BFKL_Q20_graph->SetLineStyle(2);
    comparison_graph->Add(BFKL_Q20_graph);
  }

  if (nucleus == "p") {
    double BK_p_Q20_x_arr[x[0].size()];
    vector_to_array(BK_p_Q20_x_arr, x[0]);

    TGraph* BK_p_Q20_graph = new TGraph(x[0].size(), BK_p_Q20_x_arr, BK_sigma_tot_p_Q20);
    BK_p_Q20_graph->SetTitle("BK"); //p BK
    BK_p_Q20_graph->SetLineColor(2);
    comparison_graph->Add(BK_p_Q20_graph);

    double BFKL_p_Q20_x_arr[x[2].size()];
    vector_to_array(BFKL_p_Q20_x_arr, x[2]);

    TGraph* BFKL_p_Q20_graph = new TGraph(x[2].size(), BFKL_p_Q20_x_arr, BFKL_sigma_tot_p_Q20);
    BFKL_p_Q20_graph->SetTitle("BFKL"); //p BFKL
    BFKL_p_Q20_graph->SetLineColor(2);
    BFKL_p_Q20_graph->SetLineStyle(2);
    comparison_graph->Add(BFKL_p_Q20_graph);
  }



  TCanvas* comparison_canvas = new TCanvas("comparison_canvas", "", 1000, 600);
  comparison_graph->Draw("A");

  gPad->SetLogx();
  gPad->SetLogy();
  comparison_graph->GetXaxis()->SetTitle("W (GeV)");
  if (diffractive) {
    title = "#sigma^{#gamma"+nucleus+"#rightarrowc#bar{c}"+nucleus+"X} (nb)";
  } else {
    title = "#sigma^{#gamma"+nucleus+"#rightarrowc#bar{c}X} (nb)";
  }
  comparison_graph->GetYaxis()->SetTitle(title);

  if (false) {
    comparison_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
  } else {
    comparison_canvas->BuildLegend(0.2, 0.7, 0.35, 0.9);
  }


  if (diffractive) {
    if (nucleus == "Pb") {
      TLatex* Q2_text = new TLatex(3e3, 2e4, "Q^{2} = 0 GeV^{2}");
      Q2_text->Draw("same");
    } else {
      TLatex* Q2_text = new TLatex(3e3, 100, "Q^{2} = 0 GeV^{2}");
      Q2_text->Draw("same");
    }
  } else {
    if (nucleus == "Pb") {
      TLatex* Q2_text = new TLatex(3e3, 6e5, "Q^{2} = 0 GeV^{2}");
      Q2_text->Draw("same");
    } else {
      TLatex* Q2_text = new TLatex(3e3, 4e3, "Q^{2} = 0 GeV^{2}");
      Q2_text->Draw("same");
    }
  }
  /*
  if (true) {
    TLatex* Q2_text = new TLatex(3e3, 20, "Q^{2} = 0 GeV^{2}");
    Q2_text->Draw("same");
  } else if (true) {
    TLatex* Q2_text = new TLatex(3e3, 5e-2, "Q^{2} = 0 GeV^{2}");
    Q2_text->Draw("same");
  } else if (false) {
    TLatex* Q2_text = new TLatex(2e3, 7e-2, "Q^{2} = 0 GeV^{2}");
    Q2_text->Draw("same");
  } else {
    TLatex* Q2_text = new TLatex(3e3, 2, "Q^{2} = 0 GeV^{2}");
    Q2_text->Draw("same");
  }
  */
  
  if (diffractive) {
    title = "figures/LHC_diffractive_"+nucleus+"_prediction.pdf";
  } else {
    title = "figures/LHC_inclusive_"+nucleus+"_prediction.pdf";
  }

  comparison_canvas->Print(title);
  
  return 0;
}
