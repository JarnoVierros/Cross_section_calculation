
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TText.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "cross_section_file_reader.h"

const double alpha_em = 1.0/137;

void read_data_file(string filename, vector<double> &Q2_values, vector<double> &beta_values, vector<double> &x_values, vector<double> &x_pom_F2_values, vector<double> &delta_stat_values, vector<double> &delta_sys_values) {
  ifstream data_file(filename);

  cout << "Reading: " << filename << endl;
  string line;
  while(getline (data_file, line)) {

    long unsigned int i = 0;
    string value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    Q2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    beta_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    x_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    x_pom_F2_values.push_back(stod(value));
    i++;

    value = "";
    while(line[i] != ' ') {
      value += line[i];
      i++;
    }
    delta_stat_values.push_back(stod(value));
    i++;

    value = "";
    while(i < line.size()) {
      value += line[i];
      i++;
    }
    delta_sys_values.push_back(stod(value));
    i++;

  }
  cout << "Finished reading file" << endl;
}

struct plot {
  TMultiGraph* comparison_graph;
  TGraphErrors* measurement_data;
  TGraphErrors* prediction;
  TGraph* FL_prediction;
  TGraph* FT_prediction;
  double Q2;
  double beta;
};

const double x_limits[2] = {1e-4, 0.1};
const double y_limits[2] = {0.0, 0.1};

int main() {

  vector<plot> plots;

  vector<double> L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit;
  vector<double> charm_L_prediction_Q2, charm_L_prediction_beta, charm_L_prediction_x, charm_L_prediction_sigma, charm_L_prediction_sigma_error, charm_L_prediction_fit;
  vector<double> T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit;
  vector<double> charm_T_prediction_Q2, charm_T_prediction_beta, charm_T_prediction_x, charm_T_prediction_sigma, charm_T_prediction_sigma_error, charm_T_prediction_fit;

  bool vector_dipamp = true;
  if (vector_dipamp) {
    string L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/data/differential_diffractive_L_vector_all.txt"
      };

    for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
      read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
    }

    string charm_L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/data/differential_diffractive_L_vector_charm_all.txt"
      };

    for (long unsigned int i=0; i<size(charm_L_prediction_filenames); i++) {
      read_differential_sigma_file(charm_L_prediction_filenames[i], charm_L_prediction_Q2, charm_L_prediction_beta, charm_L_prediction_x, charm_L_prediction_sigma, charm_L_prediction_sigma_error, charm_L_prediction_fit);
    }

    string T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/data/differential_diffractive_T_vector_all.txt"
    };

    for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
      read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
    }

    string charm_T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/data/differential_diffractive_T_vector_charm_all.txt"
    };

    for (long unsigned int i=0; i<size(charm_T_prediction_filenames); i++) {
      read_differential_sigma_file(charm_T_prediction_filenames[i], charm_T_prediction_Q2, charm_T_prediction_beta, charm_T_prediction_x, charm_T_prediction_sigma, charm_T_prediction_sigma_error, charm_T_prediction_fit);
    }

  } else {
    string L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/data/differential_diffractive_L_all.txt"
      };

    for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
      read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
    }

    string charm_L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/data/differential_diffractive_L_charm_all.txt"
      };

    for (long unsigned int i=0; i<size(charm_L_prediction_filenames); i++) {
      read_differential_sigma_file(charm_L_prediction_filenames[i], charm_L_prediction_Q2, charm_L_prediction_beta, charm_L_prediction_x, charm_L_prediction_sigma, charm_L_prediction_sigma_error, charm_L_prediction_fit);
    }

    string T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/data/differential_diffractive_T_all.txt"
    };

    for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
      read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
    }

    string charm_T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/data/differential_diffractive_T_charm_all.txt"
    };

    for (long unsigned int i=0; i<size(charm_T_prediction_filenames); i++) {
      read_differential_sigma_file(charm_T_prediction_filenames[i], charm_T_prediction_Q2, charm_T_prediction_beta, charm_T_prediction_x, charm_T_prediction_sigma, charm_T_prediction_sigma_error, charm_T_prediction_fit);
    }
  }

  vector<double> measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys;
  string measurement_filename = "data/differential_HERA_data.dat";

  read_data_file(measurement_filename, measurement_Q2, measurement_beta, measurement_x, measurement_xpomF2, measurement_delta_stat, measurement_delta_sys);

/*
  double Q2_selections[] = {
    4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 
    7.5, 7.5, 7.5, 7.5, 7.5, 7.5, 
    9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 
    12, 12, 12, 12, 12, 12, 
    18, 18, 18, 18, 18, 18, 
    28, 28, 28, 28, 28, 28, 
    45, 45, 45, 45, 45, 45, 
        75, 75, 75, 75, 75
    };

  double beta_selections[] = {
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9,
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
          0.1, 0.2, 0.4, 0.65, 0.9
    };
*/

  double Q2_selections[] = {
    4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 
    };

  double beta_selections[] = {
    0.04, 0.1, 0.2, 0.4, 0.65, 0.9, 
    };

  if (size(Q2_selections) != size(beta_selections)) {
    cout << "Selections not the same size" << endl;
    throw 1;
  }


  for (int k=0; k<size(Q2_selections); k++) {

    vector<double> x_selection, chosen_measurement_xpomF2, chosen_delta;

    for (int i=0; i<measurement_Q2.size(); i++) {
      if (measurement_Q2[i] != Q2_selections[k]) {
        continue;
      }
      if (measurement_beta[i] != beta_selections[k]) {
        continue;
      }
      x_selection.push_back(measurement_x[i]);
      chosen_measurement_xpomF2.push_back(measurement_xpomF2[i]);
      chosen_delta.push_back(measurement_delta_stat[i]+measurement_delta_sys[i]);
    }


    vector<double> chosen_prediction_xpomF2, chosen_prediction_xpomFL, chosen_prediction_xpomFT, chosen_prediction_error;

    for (int i=0; i<x_selection.size(); i++) {
      double L_sigma, L_error, T_sigma, T_error;
      double L_found = false, T_found = false;
      for (int j=0; j<L_prediction_Q2.size(); j++) {

        if (L_prediction_Q2[j] != Q2_selections[k]) {
          continue;
        }
        if (L_prediction_beta[j] != beta_selections[k]) {
          continue;
        }
        if (L_prediction_x[j] != x_selection[i]) {
          continue;
        }
        L_sigma = L_prediction_sigma[j] + charm_L_prediction_sigma[j];
        L_error = sqrt(L_prediction_sigma_error[j]*L_prediction_sigma_error[j] + charm_L_prediction_sigma_error[j]*charm_L_prediction_sigma_error[j]);
        L_found = true;
        if (L_sigma < 0) {
          L_sigma = 0;
        }
        break;
      }
      for (int j=0; j<T_prediction_Q2.size(); j++) {
        if (T_prediction_Q2[j] != Q2_selections[k]) {
          continue;
        }
        if (T_prediction_beta[j] != beta_selections[k]) {
          continue;
        }
        if (T_prediction_x[j] != x_selection[i]) {
          continue;
        }
        T_sigma = T_prediction_sigma[j] + charm_T_prediction_sigma[j];
        T_error = sqrt(T_prediction_sigma_error[j]*T_prediction_sigma_error[j] + T_prediction_sigma_error[j]*T_prediction_sigma_error[j]);
        T_found = true;
        if (T_sigma < 0) {
          T_sigma = 0;
        }
        break;
      }
      if (!L_found) {
        cout << "Warning: L prediction not found" << endl;
        L_sigma = 0;
      }
      if (!T_found) {
        cout << "Warning: T prediction not found" << endl;
        T_sigma = 0;
      }
      const double correction = 1;//12
      chosen_prediction_xpomF2.push_back(correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*(L_sigma + T_sigma));
      chosen_prediction_xpomFL.push_back(correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*L_sigma);
      chosen_prediction_xpomFT.push_back(correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*T_sigma);
      chosen_prediction_error.push_back(correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*beta_selections[k]*alpha_em)*sqrt(L_error*L_error + T_error*T_error));
    }
/*
struct plot {
  TString title;
  TGraphErrors* measurement_data;
  TGraphErrors* prediction;
  TGraph* FL_prediction;
  TGraph* FT_prediction;
};
*/

    double zeroes[x_selection.size()];
    for (int i=0; i<x_selection.size(); i++) {
      zeroes[i] = 0;
    }
    double* x_selection_arr = &x_selection[0];
    double* chosen_measurement_xpomF2_arr = &chosen_measurement_xpomF2[0];
    double* chosen_delta_arr = &chosen_delta[0];

    TMultiGraph* comparison_graph = new TMultiGraph();
    comparison_graph->GetXaxis()->SetLimits(x_limits[0], x_limits[1]);
    comparison_graph->GetYaxis()->SetRangeUser(y_limits[0], y_limits[1]);
    
    comparison_graph->GetXaxis()->SetLabelSize(0.05);
    comparison_graph->GetYaxis()->SetLabelSize(0.05);

    comparison_graph->GetXaxis()->SetTitle("x");
    comparison_graph->GetYaxis()->SetTitle("x_{pom}F_{2}^{D(3)}");

    /*
    stringstream Q2_stream;
    Q2_stream << fixed << setprecision(1) << Q2_selections[k];
    stringstream beta_stream;
    beta_stream << fixed << setprecision(2) << beta_selections[k];
    TString title = "x_{pom}F_{2}^{D(3)} measurement prediction comparison at Q^{2}="+Q2_stream.str()+", #beta="+beta_stream.str()+";x;x_{pom}F_{2}^{D(3)}";
    comparison_graph->SetTitle(title);
    */

    TGraphErrors* measurement_data = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_measurement_xpomF2_arr, zeroes, chosen_delta_arr);
    measurement_data->SetMarkerStyle(7);
    measurement_data->SetTitle("Measurement data");
    comparison_graph->Add(measurement_data, "P");

    double* chosen_prediction_xpomF2_arr = &chosen_prediction_xpomF2[0];
    double* chosen_prediction_error_arr = &chosen_prediction_error[0];

    TGraphErrors* prediction = new TGraphErrors(x_selection.size(), x_selection_arr, chosen_prediction_xpomF2_arr, zeroes, chosen_prediction_error_arr);
    prediction->SetTitle("Prediction");
    prediction->SetMarkerStyle(8);
    prediction->SetMarkerColor(2);
    prediction->SetLineColor(2);
    comparison_graph->Add(prediction, "C");

    double* chosen_prediction_xpomFL_arr = &chosen_prediction_xpomFL[0];

    TGraph* FL_prediction = new TGraph(x_selection.size(), x_selection_arr, chosen_prediction_xpomFL_arr);
    FL_prediction->SetLineColor(3);
    FL_prediction->SetLineStyle(3);
    comparison_graph->Add(FL_prediction, "C");

    double* chosen_prediction_xpomFT_arr = &chosen_prediction_xpomFT[0];

    TGraph* FT_prediction = new TGraph(x_selection.size(), x_selection_arr, chosen_prediction_xpomFT_arr);
    FT_prediction->SetLineColor(4);
    FT_prediction->SetLineStyle(2);
    comparison_graph->Add(FT_prediction, "C");

    plot new_plot = {comparison_graph, measurement_data, prediction, FL_prediction, FT_prediction, Q2_selections[k], beta_selections[k]};
    plots.push_back(new_plot);

  }

  int figure_width = 6;
  int figure_height = 1; //8
  TCanvas* multicanvas = new TCanvas("multicanvas", "multipads", figure_width*10000, figure_height*10000);
  multicanvas->Divide(figure_width, figure_height, 0, 0);

  int offset = 0;
  for (int i=0; i<plots.size(); i++) {

    if (i == 42) {
      offset = 1;
    }
    //gPad->SetTopMargin(0.1);
    multicanvas->cd(i+1+offset);
    //gPad->SetTickx(2);

    plots[i].comparison_graph->Draw("A");

    gPad->SetLogx();

    stringstream Q2_stream;
    Q2_stream << setprecision(2) << plots[i].Q2;
    TString Q2_string = "Q2=" + Q2_stream.str();
    TText* Q2_text = new TText(1e-2, 0.09, Q2_string);
    Q2_text->Draw("Same");

    stringstream beta_stream;
    beta_stream << setprecision(2) << plots[i].beta;
    TString beta_string = "beta=" + beta_stream.str();
    TText* beta_text = new TText(1e-2, 0.08, beta_string);
    beta_text->Draw("Same");

    /*
    float location[4];
    location[0] = 0.5;
    location[1] = 0.7;
    location[2] = 0.75;
    location[3] = 0.9;
    
    TLegend* legend = new TLegend(location[0], location[1], location[2], location[3]);
    legend->AddEntry(plots[i].measurement_data,"Measurement data");
    legend->AddEntry(plots[i].prediction,"Prediction");
    legend->AddEntry(plots[i].FL_prediction,"F_{L}");
    legend->AddEntry(plots[i].FT_prediction,"F_{T}");
    //legend->SetTextSize(0.04);
    legend->Draw();
    */
  }
  multicanvas->cd(0);

  

  //TText* title = new TText(0.48, 0.9, "Title");
  //title->Draw("Same");

  //TPad *top_pad = new TPad("top_pad", "top", 0, 0.45, 1, 0.9);
  //top_pad->Draw();

  TString figure_filename = "figures/F2D_data_comparison.pdf";
  multicanvas->Print(figure_filename);

  return 0;
}
