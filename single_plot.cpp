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

const double min_b_W = 83.6; //83.6

string nucleus = "p";
bool diffractive = false;
string diff_dipole = ""; //_diffraction
string particle_name = "b";

int main() {


  vector<TString> graph_titles;
  if (nucleus == "Pb" && diffractive) {
    graph_titles = {"#gammaPb#rightarrowc#bar{c}PbX bk", "#gammaPb#rightarrowc#bar{c}PbX bfkl", "#gammaPb#rightarrowb#bar{b}PbX bk", "#gammaPb#rightarrowb#bar{b}PbX bfkl"};
  } else if (nucleus == "Pb" && !diffractive) {
    graph_titles = {"#gammaPb#rightarrowc#bar{c}X bk", "#gammaPb#rightarrowc#bar{c}X bfkl", "#gammaPb#rightarrowb#bar{b}X bk", "#gammaPb#rightarrowb#bar{b}X bfkl"};
  } else if (nucleus == "p" && diffractive) {
    graph_titles = {"#gammap#rightarrowc#bar{c}pX bk", "#gammap#rightarrowc#bar{c}pX bfkl", "#gammap#rightarrowb#bar{b}pX bk", "#gammap#rightarrowb#bar{b}pX bfkl"};
  } else if (nucleus == "p" && !diffractive) {
    graph_titles = {"#gammap#rightarrowc#bar{c}X bk", "#gammap#rightarrowc#bar{c}X bfkl", "#gammap#rightarrowb#bar{b}X bk", "#gammap#rightarrowb#bar{b}X bfkl"};
  }
  

  string c_filenames[4];
  string b_filenames[4];

  if (diffractive) {
    c_filenames[0] = "output/diff_LHC_T_sigma_W_c_bk_Pb"+diff_dipole+".txt";
    c_filenames[1] = "output/diff_LHC_T_sigma_W_c_bfkl_Pb"+diff_dipole+".txt";
    c_filenames[2] = "output/diff_LHC_T_sigma_W_c_bk_p"+diff_dipole+".txt";
    c_filenames[3] = "output/diff_LHC_T_sigma_W_c_bfkl_p"+diff_dipole+".txt";

    b_filenames[0] = "output/diff_LHC_T_sigma_W_b_bk_Pb"+diff_dipole+".txt";
    b_filenames[1] = "output/diff_LHC_T_sigma_W_b_bfkl_Pb"+diff_dipole+".txt";
    b_filenames[2] = "output/diff_LHC_T_sigma_W_b_bk_p"+diff_dipole+".txt";
    b_filenames[3] = "output/diff_LHC_T_sigma_W_b_bfkl_p"+diff_dipole+".txt";
  } else {

    c_filenames[0] = "output/J_LHC_T_inclusive_c_bk_Pb"+diff_dipole+".txt";
    c_filenames[1] = "output/J_LHC_T_inclusive_c_bfkl_Pb"+diff_dipole+".txt";
    c_filenames[2] = "output/J_LHC_T_inclusive_c_bk_p"+diff_dipole+".txt";
    c_filenames[3] = "output/J_LHC_T_inclusive_c_bfkl_p"+diff_dipole+".txt";

    b_filenames[0] = "output/J_LHC_T_inclusive_b_bk_Pb"+diff_dipole+".txt";
    b_filenames[1] = "output/J_LHC_T_inclusive_b_bfkl_Pb"+diff_dipole+".txt";
    b_filenames[2] = "output/J_LHC_T_inclusive_b_bk_p"+diff_dipole+".txt";
    b_filenames[3] = "output/J_LHC_T_inclusive_b_bfkl_p"+diff_dipole+".txt";
  }

  int filecount = size(c_filenames);

  vector<double> c_x[filecount], c_sigma[filecount], c_sigma_error[filecount];
  vector<double> raw_b_x[filecount], b_sigma[filecount], b_sigma_error[filecount];

  for (int i=0; i<size(c_filenames); i++) {
    cout << "reading " << c_filenames[i] << endl;
    read_LHC_sigma_file(c_filenames[i], c_x[i], c_sigma[i], c_sigma_error[i]);
  }

  for (int i=0; i<size(b_filenames); i++) {
    cout << "reading " << b_filenames[i] << endl;
    read_LHC_sigma_file(b_filenames[i], raw_b_x[i], b_sigma[i], b_sigma_error[i]);
  }

  int b_skip = 0;
  vector<double> b_x[1];
  for (int i=0; i<raw_b_x[0].size(); i++) {
    if (raw_b_x[0][i] < min_b_W) {
      b_skip++;
      continue;
    }
    b_x[0].push_back(raw_b_x[0][i]);
  }

  double c_BK_sigma_tot_Pb_Q20[c_x[0].size()], c_BFKL_sigma_tot_Pb_Q20[c_x[0].size()], c_BK_sigma_tot_p_Q20[c_x[0].size()], c_BFKL_sigma_tot_p_Q20[c_x[0].size()];
  double b_BK_sigma_tot_Pb_Q20[raw_b_x[0].size()], b_BFKL_sigma_tot_Pb_Q20[raw_b_x[0].size()], b_BK_sigma_tot_p_Q20[raw_b_x[0].size()], b_BFKL_sigma_tot_p_Q20[raw_b_x[0].size()];

  const double GeV_to_nb_conversion = 1e6*0.389379;

  for (int i=0; i<c_x[0].size(); i++) {
    c_BK_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*c_sigma[0][i];
    c_BFKL_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*c_sigma[1][i];
    c_BK_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*c_sigma[2][i];
    c_BFKL_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*c_sigma[3][i];
  }
  
  for (int i=0; i<b_x[0].size(); i++) {
    b_BK_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*b_sigma[0][i+b_skip];
    b_BFKL_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*b_sigma[1][i+b_skip];
    b_BK_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*b_sigma[2][i+b_skip];
    b_BFKL_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*b_sigma[3][i+b_skip];
  }

  TCanvas *double_canvas = new TCanvas("double_canvas", "", 0.5*1.5*2000, 1.5*800);

  //double_canvas->Draw();
  //double_canvas->cd();
  /*
  TPad *left_pad = new TPad("left_pad", "left_pad", 0, 0, 0.5, 1.0);
  left_pad->SetRightMargin(0.045);
  left_pad->Draw();
  left_pad->cd();
  */

  if (!diffractive) {
    //cout << "TEST" << endl;

    TMultiGraph* left_comparison_graph = new TMultiGraph();
    TString left_title;
    if (diffractive) {
      if (nucleus == "Pb") {
        //left_title = "Diffractive c#bar{c} and b#bar{b} cross sections in #gamma lead scattering";
      } else {
        //left_title = "Diffractive c#bar{c} and b#bar{b} cross sections in #gamma proton scattering";
      }
    } else {
      if (nucleus == "Pb") {
        //left_title = "Inclusive c#bar{c} and b#bar{b} cross sections in #gamma lead scattering";
      } else {
        //left_title = "Inclusive c#bar{c} and b#bar{b} cross sections in #gamma proton scattering";
      }
    }
    left_comparison_graph->SetTitle(left_title);

    if (nucleus == "Pb") {
      double c_BK_Pb_Q20_x_arr[c_x[0].size()];
      vector_to_array(c_BK_Pb_Q20_x_arr, c_x[0]);
      TGraph* c_BK_Q20_graph = new TGraph(c_x[0].size(), c_BK_Pb_Q20_x_arr, c_BK_sigma_tot_Pb_Q20);
      c_BK_Q20_graph->SetTitle(graph_titles[0]); //Pb BK
      c_BK_Q20_graph->SetLineColor(2);
      left_comparison_graph->Add(c_BK_Q20_graph);
      double c_BFKL_Pb_Q20_x_arr[c_x[2].size()];
      vector_to_array(c_BFKL_Pb_Q20_x_arr, c_x[2]);

      TGraph* c_BFKL_Q20_graph = new TGraph(c_x[2].size(), c_BFKL_Pb_Q20_x_arr, c_BFKL_sigma_tot_Pb_Q20);
      c_BFKL_Q20_graph->SetTitle(graph_titles[1]); //Pb BFKL
      c_BFKL_Q20_graph->SetLineColor(2);
      c_BFKL_Q20_graph->SetLineStyle(2);
      left_comparison_graph->Add(c_BFKL_Q20_graph);

      double b_BK_Pb_Q20_x_arr[b_x[0].size()];
      vector_to_array(b_BK_Pb_Q20_x_arr, b_x[0]);

      TGraph* b_BK_Q20_graph = new TGraph(b_x[0].size(), b_BK_Pb_Q20_x_arr, b_BK_sigma_tot_Pb_Q20);
      b_BK_Q20_graph->SetTitle(graph_titles[2]); //Pb BK
      b_BK_Q20_graph->SetLineColor(4);
      b_BK_Q20_graph->SetLineStyle(7);
      left_comparison_graph->Add(b_BK_Q20_graph);

      double b_BFKL_Pb_Q20_x_arr[b_x[0].size()];
      vector_to_array(b_BFKL_Pb_Q20_x_arr, b_x[0]);

      TGraph* b_BFKL_Q20_graph = new TGraph(b_x[0].size(), b_BFKL_Pb_Q20_x_arr, b_BFKL_sigma_tot_Pb_Q20);
      b_BFKL_Q20_graph->SetTitle(graph_titles[3]); //Pb BFKL
      b_BFKL_Q20_graph->SetLineColor(4);
      b_BFKL_Q20_graph->SetLineStyle(3);
      left_comparison_graph->Add(b_BFKL_Q20_graph);
    }

    if (nucleus == "p") {
      double c_BK_p_Q20_x_arr[c_x[0].size()];
      vector_to_array(c_BK_p_Q20_x_arr, c_x[0]);

      TGraph* c_BK_p_Q20_graph = new TGraph(c_x[0].size(), c_BK_p_Q20_x_arr, c_BK_sigma_tot_p_Q20);
      c_BK_p_Q20_graph->SetTitle(graph_titles[0]); //p BK
      c_BK_p_Q20_graph->SetLineColor(2);
      left_comparison_graph->Add(c_BK_p_Q20_graph);

      double c_BFKL_p_Q20_x_arr[c_x[2].size()];
      vector_to_array(c_BFKL_p_Q20_x_arr, c_x[2]);

      TGraph* c_BFKL_p_Q20_graph = new TGraph(c_x[2].size(), c_BFKL_p_Q20_x_arr, c_BFKL_sigma_tot_p_Q20);
      c_BFKL_p_Q20_graph->SetTitle(graph_titles[1]); //p BFKL
      c_BFKL_p_Q20_graph->SetLineColor(2);
      c_BFKL_p_Q20_graph->SetLineStyle(2);
      left_comparison_graph->Add(c_BFKL_p_Q20_graph);


      double b_BK_p_Q20_x_arr[b_x[0].size()];
      vector_to_array(b_BK_p_Q20_x_arr, b_x[0]);
      TGraph* b_BK_p_Q20_graph = new TGraph(b_x[0].size(), b_BK_p_Q20_x_arr, b_BK_sigma_tot_p_Q20);
      b_BK_p_Q20_graph->SetTitle(graph_titles[2]); //p BK
      b_BK_p_Q20_graph->SetLineColor(4);
      b_BK_p_Q20_graph->SetLineStyle(7);
      left_comparison_graph->Add(b_BK_p_Q20_graph);

      double b_BFKL_p_Q20_x_arr[b_x[0].size()];
      vector_to_array(b_BFKL_p_Q20_x_arr, b_x[0]);

      TGraph* b_BFKL_p_Q20_graph = new TGraph(b_x[0].size(), b_BFKL_p_Q20_x_arr, b_BFKL_sigma_tot_p_Q20);
      b_BFKL_p_Q20_graph->SetTitle(graph_titles[3]); //p BFKL
      b_BFKL_p_Q20_graph->SetLineColor(4);
      b_BFKL_p_Q20_graph->SetLineStyle(3);
      left_comparison_graph->Add(b_BFKL_p_Q20_graph);
    }

    left_comparison_graph->Draw("A");

    gPad->SetLogx();
    gPad->SetLogy();
    left_comparison_graph->GetXaxis()->SetTitle("W (GeV)");
    if (diffractive) {
      left_title = "#sigma^{#gamma"+nucleus+"#rightarrowq#bar{q}"+nucleus+"X} (nb)";
    } else {
      left_title = "#sigma^{#gamma"+nucleus+"#rightarrowq#bar{q}X} (nb)";
    }
    left_comparison_graph->GetYaxis()->SetTitle(left_title);
    left_comparison_graph->GetYaxis()->SetTitleSize(0.045);
    left_comparison_graph->GetYaxis()->SetTitleOffset(0.95);

    if (false) {
      double_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
    } else {
      double_canvas->BuildLegend(0.14, 0.7, 0.415, 0.9);
    }

    //left_comparison_graph->GetXaxis()->SetLabelOffset(-0.05);
    left_comparison_graph->GetXaxis()->SetTitleOffset(1.2);

    if (particle_name == "c") {
      if (diff_dipole=="_diffraction") {
        if (diffractive) {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 8e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 2e2, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        } else {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e6, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      } else {
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
      }
    } else if (particle_name == "b") {
      if (diff_dipole=="_diffraction") {
        if (diffractive) {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e3, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 3e0, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        } else {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e6, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      } else {
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
            TLatex* Q2_text = new TLatex(3e3, 2e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e2, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      }
    }
  } else {
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //double_canvas->cd();

    //string nucleus = "p";
    string diff_dipole = "_diffraction"; //_diffraction
    string particle_name = "b";



    if (diffractive) {
      c_filenames[0] = "output/diff_LHC_T_sigma_W_c_bk_Pb"+diff_dipole+".txt";
      c_filenames[1] = "output/diff_LHC_T_sigma_W_c_bfkl_Pb"+diff_dipole+".txt";
      c_filenames[2] = "output/diff_LHC_T_sigma_W_c_bk_p"+diff_dipole+".txt";
      c_filenames[3] = "output/diff_LHC_T_sigma_W_c_bfkl_p"+diff_dipole+".txt";

      b_filenames[0] = "output/diff_LHC_T_sigma_W_b_bk_Pb"+diff_dipole+".txt";
      b_filenames[1] = "output/diff_LHC_T_sigma_W_b_bfkl_Pb"+diff_dipole+".txt";
      b_filenames[2] = "output/diff_LHC_T_sigma_W_b_bk_p"+diff_dipole+".txt";
      b_filenames[3] = "output/diff_LHC_T_sigma_W_b_bfkl_p"+diff_dipole+".txt";
    } else {

      c_filenames[0] = "output/J_LHC_T_inclusive_c_bk_Pb"+diff_dipole+".txt";
      c_filenames[1] = "output/J_LHC_T_inclusive_c_bfkl_Pb"+diff_dipole+".txt";
      c_filenames[2] = "output/J_LHC_T_inclusive_c_bk_p"+diff_dipole+".txt";
      c_filenames[3] = "output/J_LHC_T_inclusive_c_bfkl_p"+diff_dipole+".txt";

      b_filenames[0] = "output/J_LHC_T_inclusive_b_bk_Pb"+diff_dipole+".txt";
      b_filenames[1] = "output/J_LHC_T_inclusive_b_bfkl_Pb"+diff_dipole+".txt";
      b_filenames[2] = "output/J_LHC_T_inclusive_b_bk_p"+diff_dipole+".txt";
      b_filenames[3] = "output/J_LHC_T_inclusive_b_bfkl_p"+diff_dipole+".txt";
    }

    vector<double> right_c_x[filecount], right_c_sigma[filecount], right_c_sigma_error[filecount];
    vector<double> raw_right_b_x[filecount], right_b_sigma[filecount], right_b_sigma_error[filecount];

    for (int i=0; i<size(c_filenames); i++) {
      cout << "reading " << c_filenames[i] << endl;
      read_LHC_sigma_file(c_filenames[i], right_c_x[i], right_c_sigma[i], right_c_sigma_error[i]);
    }

    for (int i=0; i<size(b_filenames); i++) {
      cout << "reading " << b_filenames[i] << endl;
      read_LHC_sigma_file(b_filenames[i], raw_right_b_x[i], right_b_sigma[i], right_b_sigma_error[i]);
    }

    b_skip = 0;
    vector<double> right_b_x[1];
    for (int i=0; i<raw_right_b_x[0].size(); i++) {
      if (raw_right_b_x[0][i] < min_b_W) {
        b_skip++;
        continue;
      }
      right_b_x[0].push_back(raw_right_b_x[0][i]);
    }

    double right_c_BK_sigma_tot_Pb_Q20[right_c_x[0].size()], right_c_BFKL_sigma_tot_Pb_Q20[right_c_x[0].size()], right_c_BK_sigma_tot_p_Q20[right_c_x[0].size()], right_c_BFKL_sigma_tot_p_Q20[right_c_x[0].size()];
    double right_b_BK_sigma_tot_Pb_Q20[right_b_x[0].size()], right_b_BFKL_sigma_tot_Pb_Q20[right_b_x[0].size()], right_b_BK_sigma_tot_p_Q20[right_b_x[0].size()], right_b_BFKL_sigma_tot_p_Q20[right_b_x[0].size()];

    for (int i=0; i<right_c_x[0].size(); i++) {
      right_c_BK_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*right_c_sigma[0][i];
      right_c_BFKL_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*right_c_sigma[1][i];
      right_c_BK_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*right_c_sigma[2][i];
      right_c_BFKL_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*right_c_sigma[3][i];
    }

    for (int i=0; i<right_b_x[0].size(); i++) {
      right_b_BK_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*right_b_sigma[0][i+b_skip];
      right_b_BFKL_sigma_tot_Pb_Q20[i] = GeV_to_nb_conversion*right_b_sigma[1][i+b_skip];
      right_b_BK_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*right_b_sigma[2][i+b_skip];
      right_b_BFKL_sigma_tot_p_Q20[i] = GeV_to_nb_conversion*right_b_sigma[3][i+b_skip];
    }


    TMultiGraph* right_comparison_graph = new TMultiGraph();
    TString right_title;
    if (diffractive) {
      if (nucleus == "Pb") {
        //right_title = "Diffractive c#bar{c} and b#bar{b} cross sections in #gamma lead scattering";
      } else {
        //right_title = "Diffractive c#bar{c} and b#bar{b} cross sections in #gamma proton scattering";
      }
    } else {
      if (nucleus == "Pb") {
        //right_title = "Inclusive c#bar{c} and b#bar{b} cross sections in #gamma lead scattering";
      } else {
        //right_title = "Inclusive c#bar{c} and b#bar{b} cross sections in #gamma proton scattering";
      }
    }
    right_comparison_graph->SetTitle(right_title);


    if (nucleus == "Pb") {
      double c_BK_Pb_Q20_x_arr[right_c_x[0].size()];
      vector_to_array(c_BK_Pb_Q20_x_arr, right_c_x[0]);

      TGraph* c_BK_Q20_graph = new TGraph(right_c_x[0].size(), c_BK_Pb_Q20_x_arr, right_c_BK_sigma_tot_Pb_Q20);
      c_BK_Q20_graph->SetTitle(graph_titles[0]); //Pb BK
      c_BK_Q20_graph->SetLineColor(2);
      right_comparison_graph->Add(c_BK_Q20_graph);

      double c_BFKL_Pb_Q20_x_arr[right_c_x[2].size()];
      vector_to_array(c_BFKL_Pb_Q20_x_arr, right_c_x[2]);

      TGraph* c_BFKL_Q20_graph = new TGraph(right_c_x[2].size(), c_BFKL_Pb_Q20_x_arr, right_c_BFKL_sigma_tot_Pb_Q20);
      c_BFKL_Q20_graph->SetTitle(graph_titles[1]); //Pb BFKL
      c_BFKL_Q20_graph->SetLineColor(2);
      c_BFKL_Q20_graph->SetLineStyle(2);
      right_comparison_graph->Add(c_BFKL_Q20_graph);


      double b_BK_Pb_Q20_x_arr[right_b_x[0].size()];
      vector_to_array(b_BK_Pb_Q20_x_arr, right_b_x[0]);

      TGraph* b_BK_Q20_graph = new TGraph(right_b_x[0].size(), b_BK_Pb_Q20_x_arr, right_b_BK_sigma_tot_Pb_Q20);
      b_BK_Q20_graph->SetTitle(graph_titles[2]); //Pb BK
      b_BK_Q20_graph->SetLineColor(4);
      b_BK_Q20_graph->SetLineStyle(7);
      right_comparison_graph->Add(b_BK_Q20_graph);

      double b_BFKL_Pb_Q20_x_arr[right_b_x[0].size()];
      vector_to_array(b_BFKL_Pb_Q20_x_arr, right_b_x[0]);

      TGraph* b_BFKL_Q20_graph = new TGraph(right_b_x[0].size(), b_BFKL_Pb_Q20_x_arr, right_b_BFKL_sigma_tot_Pb_Q20);
      b_BFKL_Q20_graph->SetTitle(graph_titles[3]); //Pb BFKL
      b_BFKL_Q20_graph->SetLineColor(4);
      b_BFKL_Q20_graph->SetLineStyle(3);
      right_comparison_graph->Add(b_BFKL_Q20_graph);
    }

    if (nucleus == "p") {
      double c_BK_p_Q20_x_arr[right_c_x[0].size()];
      vector_to_array(c_BK_p_Q20_x_arr, right_c_x[0]);

      TGraph* c_BK_p_Q20_graph = new TGraph(right_c_x[0].size(), c_BK_p_Q20_x_arr, right_c_BK_sigma_tot_p_Q20);
      c_BK_p_Q20_graph->SetTitle(graph_titles[0]); //p BK
      c_BK_p_Q20_graph->SetLineColor(2);
      right_comparison_graph->Add(c_BK_p_Q20_graph);

      double c_BFKL_p_Q20_x_arr[right_c_x[2].size()];
      vector_to_array(c_BFKL_p_Q20_x_arr, right_c_x[2]);

      TGraph* c_BFKL_p_Q20_graph = new TGraph(right_c_x[2].size(), c_BFKL_p_Q20_x_arr, right_c_BFKL_sigma_tot_p_Q20);
      c_BFKL_p_Q20_graph->SetTitle(graph_titles[1]); //p BFKL
      c_BFKL_p_Q20_graph->SetLineColor(2);
      c_BFKL_p_Q20_graph->SetLineStyle(2);
      right_comparison_graph->Add(c_BFKL_p_Q20_graph);



      double b_BK_p_Q20_x_arr[right_b_x[0].size()];
      vector_to_array(b_BK_p_Q20_x_arr, right_b_x[0]);

      TGraph* b_BK_p_Q20_graph = new TGraph(right_b_x[0].size(), b_BK_p_Q20_x_arr, right_b_BK_sigma_tot_p_Q20);
      b_BK_p_Q20_graph->SetTitle(graph_titles[2]); //p BK
      b_BK_p_Q20_graph->SetLineColor(4);
      b_BK_p_Q20_graph->SetLineStyle(7);
      right_comparison_graph->Add(b_BK_p_Q20_graph);

      double b_BFKL_p_Q20_x_arr[right_b_x[0].size()];
      vector_to_array(b_BFKL_p_Q20_x_arr, right_b_x[0]);

      TGraph* b_BFKL_p_Q20_graph = new TGraph(right_b_x[0].size(), b_BFKL_p_Q20_x_arr, right_b_BFKL_sigma_tot_p_Q20);
      b_BFKL_p_Q20_graph->SetTitle(graph_titles[3]); //p BFKL
      b_BFKL_p_Q20_graph->SetLineColor(4);
      b_BFKL_p_Q20_graph->SetLineStyle(3);
      right_comparison_graph->Add(b_BFKL_p_Q20_graph);
    }
    //right_comparison_graph->GetYaxis()->SetLimits(3e-1, 5e4);
    right_comparison_graph->Draw("A");

    gPad->SetLogx();
    gPad->SetLogy();
    right_comparison_graph->GetXaxis()->SetTitle("W (GeV)");
    if (diffractive) {
      right_title = "#sigma^{#gamma"+nucleus+"#rightarrowq#bar{q}"+nucleus+"X} (nb)";
    } else {
      right_title = "#sigma^{#gamma"+nucleus+"#rightarrowq#bar{q}X} (nb)";
    }
    right_comparison_graph->GetYaxis()->SetTitle(right_title);
    right_comparison_graph->GetYaxis()->SetTitleSize(0.045);
    right_comparison_graph->GetYaxis()->SetTitleOffset(0.95);

    if (false) {
      double_canvas->BuildLegend(0.75, 0.55, 0.9, 0.9);
    } else {
      double_canvas->BuildLegend(0.14, 0.7, 0.45, 0.9);
    }

    right_comparison_graph->GetXaxis()->SetTitleOffset(1.2);

    if (particle_name == "c") {
      if (diff_dipole=="_diffraction") {
        if (diffractive) {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 8e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 2e2, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        } else {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e6, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      } else {
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
      }
    } else if (particle_name == "b") {
      if (diff_dipole=="_diffraction") {
        if (diffractive) {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e3, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 3e0, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        } else {
          if (nucleus == "Pb") {
            TLatex* Q2_text = new TLatex(3e3, 1e6, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      } else {
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
            TLatex* Q2_text = new TLatex(3e3, 2e4, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          } else {
            TLatex* Q2_text = new TLatex(3e3, 1e2, "Q^{2} = 0 GeV^{2}");
            Q2_text->Draw("same");
          }
        }
      }
    }
    
  }

  TString filename = "figures/single_plot";
  if (nucleus == "Pb") {
    filename += "_Pb";
  } else {
    filename += "_p";
  }
  if (diffractive) {
    filename += "_diffractive";
  } else {
    filename += "_inclusive";
  }
  filename += ".pdf";
  double_canvas->Print(filename);
  return 0;
}
