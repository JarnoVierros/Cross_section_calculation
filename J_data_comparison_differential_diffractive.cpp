
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TText.h"
#include "TLatex.h"
#include "TColor.h"

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
  TGraph* qqg_correction;
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
  vector<double> qqg_low_beta_correction_Q2, qqg_low_beta_correction_beta, qqg_low_beta_correction_x, qqg_low_beta_correction_sigma, qqg_low_beta_correction_error, qqg_low_beta_correction_fit;
  vector<double> charm_qqg_low_beta_correction_Q2, charm_qqg_low_beta_correction_beta, charm_qqg_low_beta_correction_x, charm_qqg_low_beta_correction_sigma, charm_qqg_low_beta_correction_error, charm_qqg_low_beta_correction_fit;

  double unit_scaler = 1;
  bool Jani_qqg_correction = false;
  bool vector_dipamp = true;
  if (vector_dipamp) {
    string L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_L_vector_all.txt"
      };

    for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
      read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
    }

    string charm_L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_L_vector_charm_all.txt"
      };

    for (long unsigned int i=0; i<size(charm_L_prediction_filenames); i++) {
      read_differential_sigma_file(charm_L_prediction_filenames[i], charm_L_prediction_Q2, charm_L_prediction_beta, charm_L_prediction_x, charm_L_prediction_sigma, charm_L_prediction_sigma_error, charm_L_prediction_fit);
    }

    string T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_T_vector_all.txt"
    };

    for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
      read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
    }

    string charm_T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_T_vector_charm_all.txt"
    };

    for (long unsigned int i=0; i<size(charm_T_prediction_filenames); i++) {
      read_differential_sigma_file(charm_T_prediction_filenames[i], charm_T_prediction_Q2, charm_T_prediction_beta, charm_T_prediction_x, charm_T_prediction_sigma, charm_T_prediction_sigma_error, charm_T_prediction_fit);
    }

    if (Jani_qqg_correction) {
      string qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/Jani_low_beta_corrections_all.txt"
      };

      for (long unsigned int i=0; i<size(qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(qqg_low_beta_correction_filenames[i], qqg_low_beta_correction_Q2, qqg_low_beta_correction_beta, qqg_low_beta_correction_x, qqg_low_beta_correction_sigma, qqg_low_beta_correction_error, qqg_low_beta_correction_fit);
      }

      string charm_qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/Jani_low_beta_corrections_all_charm.txt"
      };

      for (long unsigned int i=0; i<size(charm_qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(charm_qqg_low_beta_correction_filenames[i], charm_qqg_low_beta_correction_Q2, charm_qqg_low_beta_correction_beta, charm_qqg_low_beta_correction_x, charm_qqg_low_beta_correction_sigma, charm_qqg_low_beta_correction_error, charm_qqg_low_beta_correction_fit);
      }
    } else {
      string qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/low_beta_corrections_all.txt"
      };

      for (long unsigned int i=0; i<size(qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(qqg_low_beta_correction_filenames[i], qqg_low_beta_correction_Q2, qqg_low_beta_correction_beta, qqg_low_beta_correction_x, qqg_low_beta_correction_sigma, qqg_low_beta_correction_error, qqg_low_beta_correction_fit);
      }

      string charm_qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/low_beta_corrections_all_charm.txt"
      };

      for (long unsigned int i=0; i<size(charm_qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(charm_qqg_low_beta_correction_filenames[i], charm_qqg_low_beta_correction_Q2, charm_qqg_low_beta_correction_beta, charm_qqg_low_beta_correction_x, charm_qqg_low_beta_correction_sigma, charm_qqg_low_beta_correction_error, charm_qqg_low_beta_correction_fit);
      }
    }



  } else {
    string L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_L_10k_all_diffdiff.txt"
      };

    for (long unsigned int i=0; i<size(L_prediction_filenames); i++) {
      read_differential_sigma_file(L_prediction_filenames[i], L_prediction_Q2, L_prediction_beta, L_prediction_x, L_prediction_sigma, L_prediction_sigma_error, L_prediction_fit);
    }

    string charm_L_prediction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_L_10k_all_diffdiff_charm.txt"
      };

    for (long unsigned int i=0; i<size(charm_L_prediction_filenames); i++) {
      read_differential_sigma_file(charm_L_prediction_filenames[i], charm_L_prediction_Q2, charm_L_prediction_beta, charm_L_prediction_x, charm_L_prediction_sigma, charm_L_prediction_sigma_error, charm_L_prediction_fit);
    }

    string T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_T_10k_all_diffdiff.txt"
    };

    for (long unsigned int i=0; i<size(T_prediction_filenames); i++) {
      read_differential_sigma_file(T_prediction_filenames[i], T_prediction_Q2, T_prediction_beta, T_prediction_x, T_prediction_sigma, T_prediction_sigma_error, T_prediction_fit);
    }

    string charm_T_prediction_filenames[] = {
      "/home/jarno/Cross_section_calculation/output/old/differential_diffractive_T_10k_all_diffdiff_charm.txt"
    };

    for (long unsigned int i=0; i<size(charm_T_prediction_filenames); i++) {
      read_differential_sigma_file(charm_T_prediction_filenames[i], charm_T_prediction_Q2, charm_T_prediction_beta, charm_T_prediction_x, charm_T_prediction_sigma, charm_T_prediction_sigma_error, charm_T_prediction_fit);
    }

    if (Jani_qqg_correction) {
      string qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/Jani_low_beta_corrections_all.txt"
      };

      for (long unsigned int i=0; i<size(qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(qqg_low_beta_correction_filenames[i], qqg_low_beta_correction_Q2, qqg_low_beta_correction_beta, qqg_low_beta_correction_x, qqg_low_beta_correction_sigma, qqg_low_beta_correction_error, qqg_low_beta_correction_fit);
      }

      string charm_qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/Jani_low_beta_corrections_all_charm.txt"
      };

      for (long unsigned int i=0; i<size(charm_qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(charm_qqg_low_beta_correction_filenames[i], charm_qqg_low_beta_correction_Q2, charm_qqg_low_beta_correction_beta, charm_qqg_low_beta_correction_x, charm_qqg_low_beta_correction_sigma, charm_qqg_low_beta_correction_error, charm_qqg_low_beta_correction_fit);
      }
    } else {
      string qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/low_beta_corrections_all.txt"
      };

      for (long unsigned int i=0; i<size(qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(qqg_low_beta_correction_filenames[i], qqg_low_beta_correction_Q2, qqg_low_beta_correction_beta, qqg_low_beta_correction_x, qqg_low_beta_correction_sigma, qqg_low_beta_correction_error, qqg_low_beta_correction_fit);
      }

      string charm_qqg_low_beta_correction_filenames[] = {
        "/home/jarno/Cross_section_calculation/output/old/low_beta_corrections_all_charm.txt"
      };

      for (long unsigned int i=0; i<size(charm_qqg_low_beta_correction_filenames); i++) {
        read_differential_sigma_file(charm_qqg_low_beta_correction_filenames[i], charm_qqg_low_beta_correction_Q2, charm_qqg_low_beta_correction_beta, charm_qqg_low_beta_correction_x, charm_qqg_low_beta_correction_sigma, charm_qqg_low_beta_correction_error, charm_qqg_low_beta_correction_fit);
      }
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

  if (size(Q2_selections) != size(beta_selections)) {
    cout << "Selections not the same size" << endl;
    throw 1;
  }

  int figure_width = 6;
  int figure_height = 8;
  double margin_fraction = 0.1;

  TLegend* legend = new TLegend(0.5*margin_fraction, 0.5*margin_fraction, margin_fraction+0.95*(1-2*margin_fraction)/figure_width, margin_fraction+0.95*(1-2*margin_fraction)/figure_height);

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
      chosen_measurement_xpomF2.push_back(unit_scaler*measurement_xpomF2[i]);
      chosen_delta.push_back(unit_scaler*(measurement_delta_stat[i]+measurement_delta_sys[i]));
    }


    vector<double> chosen_prediction_xpomF2, chosen_prediction_xpomFL, chosen_prediction_xpomFT, chosen_prediction_error, chosen_qqg_correction;

    for (int i=0; i<x_selection.size(); i++) {
      double L_sigma, L_error, T_sigma, T_error, qqg_correction, qqg_correction_error;
      double L_found = false, T_found = false, qqg_found = false;
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
        T_error = sqrt(T_prediction_sigma_error[j]*T_prediction_sigma_error[j] + charm_T_prediction_sigma_error[j]*charm_T_prediction_sigma_error[j]);
        T_found = true;
        if (T_sigma < 0) {
          T_sigma = 0;
        }
        break;
      }
      for (int j=0; j<qqg_low_beta_correction_Q2.size(); j++) {
        if (qqg_low_beta_correction_Q2[j] != Q2_selections[k]) {
          continue;
        }
        if (qqg_low_beta_correction_beta[j] != beta_selections[k]) {
          continue;
        }
        if (qqg_low_beta_correction_x[j] != x_selection[i]) {
          continue;
        }
        qqg_correction = qqg_low_beta_correction_sigma[j] + charm_qqg_low_beta_correction_sigma[j];
        qqg_correction_error = sqrt(qqg_low_beta_correction_error[j]*qqg_low_beta_correction_error[j] + charm_qqg_low_beta_correction_error[j]*charm_qqg_low_beta_correction_error[j]);
        qqg_found = true;
        if (qqg_correction < 0) {
          qqg_correction = 0;
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
      if (!qqg_found) {
        cout << "Warning: qqg correction not found" << endl;
        qqg_correction = 0;
      }

      const double Jani_dipamp_normalization = 1.0/20; //1.0/40
      if (Jani_qqg_correction) {
        qqg_correction = Jani_dipamp_normalization*qqg_correction;
        qqg_correction_error = Jani_dipamp_normalization*qqg_correction_error;
      }

      const double correction = 1;//12
      chosen_prediction_xpomFL.push_back(unit_scaler*correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*L_sigma);
      chosen_prediction_xpomFT.push_back(unit_scaler*correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*T_sigma);
      double combined_error = correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*beta_selections[k]*alpha_em)*sqrt(L_error*L_error + T_error*T_error);
      combined_error = sqrt(combined_error*combined_error + qqg_correction_error*qqg_correction_error);
      chosen_prediction_error.push_back(unit_scaler*combined_error);
      chosen_qqg_correction.push_back(unit_scaler*qqg_correction);
      chosen_prediction_xpomF2.push_back(unit_scaler*(correction*Q2_selections[k]*Q2_selections[k]/(pow(2*M_PI, 2)*alpha_em*beta_selections[k])*(L_sigma + T_sigma) + qqg_correction));
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

    vector<double> valid_prediction_x;
    for (int m=0;m<x_selection.size();m++) {
      if (x_selection[m] < 0.01) {
        valid_prediction_x.push_back(x_selection[m]);
      }
    }
    double* prediction_x_selection_arr = &valid_prediction_x[0];

    int valid_prediction_size = valid_prediction_x.size();

    vector<double> valid_chosen_prediction_xpomF2, valid_chosen_prediction_error, valid_chosen_prediction_xpomFL, valid_chosen_prediction_xpomFT, valid_chosen_qqg_correction;
    for (int m=0;m<valid_prediction_size;m++) {
      valid_chosen_prediction_xpomF2.push_back(chosen_prediction_xpomF2[m]);
      valid_chosen_prediction_error.push_back(chosen_prediction_error[m]);
      valid_chosen_prediction_xpomFL.push_back(chosen_prediction_xpomFL[m]);
      valid_chosen_prediction_xpomFT.push_back(chosen_prediction_xpomFT[m]);
      valid_chosen_qqg_correction.push_back(chosen_qqg_correction[m]);
    }
    double* chosen_prediction_xpomF2_arr = &valid_chosen_prediction_xpomF2[0];
    double* chosen_prediction_error_arr = &valid_chosen_prediction_error[0];
    double* chosen_prediction_xpomFL_arr = &valid_chosen_prediction_xpomFL[0];
    double* chosen_prediction_xpomFT_arr = &valid_chosen_prediction_xpomFT[0];
    double* chosen_qqg_correction_arr = &valid_chosen_qqg_correction[0];


    TMultiGraph* comparison_graph = new TMultiGraph();
    comparison_graph->GetXaxis()->SetLimits(x_limits[0], x_limits[1]);
    comparison_graph->GetYaxis()->SetRangeUser(unit_scaler*y_limits[0], unit_scaler*y_limits[1]);
    
    comparison_graph->GetXaxis()->SetLabelSize(0.15);
    comparison_graph->GetYaxis()->SetLabelSize(0.15);

    comparison_graph->GetXaxis()->SetNdivisions(3);
    comparison_graph->GetYaxis()->SetNdivisions(6);

    //comparison_graph->GetXaxis()->SetLabelFont(21);
    comparison_graph->GetXaxis()->SetLabelOffset(-0.05);
    if (k>41) {
      comparison_graph->GetXaxis()->SetTickSize(0.04);
    } else {
      comparison_graph->GetXaxis()->SetTickSize(0.075);
    }

    //comparison_graph->GetXaxis()->SetTitle("x");
    //comparison_graph->GetYaxis()->SetTitle("x_{pom}F_{2}^{D(3)}");

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

    TGraphErrors* prediction = new TGraphErrors(valid_prediction_size, x_selection_arr, chosen_prediction_xpomF2_arr, zeroes, chosen_prediction_error_arr);
    prediction->SetTitle("Prediction");
    prediction->SetMarkerStyle(8);
    prediction->SetMarkerColor(2);
    prediction->SetLineWidth(2);

    Int_t cred = TColor::GetFreeColorIndex();
    TColor* color_red = new TColor(cred, 228.0/255, 26.0/255, 28.0/255);
    prediction->SetLineColor(cred);

    comparison_graph->Add(prediction, "C");

    TGraph* FL_prediction = new TGraph(valid_prediction_size, x_selection_arr, chosen_prediction_xpomFL_arr);

    Int_t cblue = TColor::GetFreeColorIndex();
    TColor* color_blue = new TColor(cblue, 55.0/255, 126.0/255, 184.0/255);
    FL_prediction->SetLineColor(cblue);

    FL_prediction->SetLineStyle(3);
    FL_prediction->SetLineWidth(3);
    comparison_graph->Add(FL_prediction, "C");

    TGraph* FT_prediction = new TGraph(valid_prediction_size, x_selection_arr, chosen_prediction_xpomFT_arr);

    Int_t cgreen = TColor::GetFreeColorIndex();
    TColor* color_green = new TColor(cgreen, 77.0/255, 175.0/255, 74.0/255);
    FT_prediction->SetLineColor(cgreen);

    FT_prediction->SetLineStyle(2);
    FT_prediction->SetLineWidth(3);
    comparison_graph->Add(FT_prediction, "C");

    TGraph* qqg_prediction = new TGraph(valid_prediction_size, x_selection_arr, chosen_qqg_correction_arr);

    Int_t cpurple = TColor::GetFreeColorIndex();
    TColor* color_purple = new TColor(cpurple, 152.0/255, 78.0/255, 163.0/255);
    qqg_prediction->SetLineColor(cpurple);

    qqg_prediction->SetLineStyle(8);
    qqg_prediction->SetLineWidth(3);
    comparison_graph->Add(qqg_prediction, "C");

    plot new_plot = {comparison_graph, measurement_data, prediction, FL_prediction, FT_prediction, qqg_prediction, Q2_selections[k], beta_selections[k]};
    plots.push_back(new_plot);

    if (k==0) {
      legend->AddEntry(measurement_data, "q#bar{q} measurement");
      legend->AddEntry(prediction, "total", "L");
      legend->AddEntry(FL_prediction, "longitudinal", "L");
      legend->AddEntry(FT_prediction, "transverse", "L");
      legend->AddEntry(qqg_prediction, "gluon emission", "L");
    }

  }



  double fig_size_x = 100;
  double fig_size_y = 100;
  TCanvas* multicanvas = new TCanvas("multicanvas", "multipads", figure_width*fig_size_x/(1-2*margin_fraction), figure_height*fig_size_y/(1-2*margin_fraction));
  //TCanvas* multicanvas = new TCanvas("multicanvas", "multipads", 769, 1115);
  
  multicanvas->Draw();
  //TPad* multipad = new TPad("multipad", "multipad", margin_fraction, margin_fraction, 1-margin_fraction, 1-margin_fraction);
  //multipad->Draw();
  multicanvas->cd(0);
  TPad* subpads[48];
  /*
  subpads[1] = new TPad("asd", "dsa", 0.1, 0.5, 0.3, 0.7);
  subpads[1]->SetFillColor(1);
  subpads[1]->Draw();

  subpads[2] = new TPad("asd", "dsa", 0.6, 0.1, 0.9, 0.2);
  subpads[2]->SetFillColor(1);
  subpads[2]->Draw();
  */
  int offset = 0;
  //for (int i=0; i<plots.size(); i++) {
  for (int i=0; i<48; i++) {
    if (i==42) {continue;}
    
    multicanvas->cd(0);
    //cout << i << endl;
    double x1 = margin_fraction+(i%figure_width)*1.0/figure_width*(1-2*margin_fraction);
    if (i%6==0) {x1=0;}
    double x2 = margin_fraction+(i%figure_width+1)*1.0/figure_width*(1-2*margin_fraction);
    double y1 = 1-margin_fraction-(i/figure_width+1)*1.0/figure_height*(1-2*margin_fraction);
    if (i/6==7) {y1=0;}
    double y2 = 1-margin_fraction-(i/figure_width)*1.0/figure_height*(1-2*margin_fraction);
    //cout << x1 << ", " << y1 << ", " << x2 << ", " << y2 << endl;
    subpads[i] = new TPad("subpad", "subpad", x1, y1, x2, y2);
    subpads[i]->SetMargin(0, 0, 0, 0);
    if (i%6==0) {
      subpads[i]->SetLeftMargin(margin_fraction/(margin_fraction+(1-2*margin_fraction)/figure_width));
    }
    if (i/6==7) {
      subpads[i]->SetBottomMargin(margin_fraction/(margin_fraction+(1-2*margin_fraction)/figure_width));
    }
    subpads[i]->Draw();
    subpads[i]->cd(0);

    int offset;
    if (i>42) {
      offset = -1;
    } else {
      offset = 0;
    }

    plots[i+offset].comparison_graph->Draw("A");

    gPad->SetLogx();

    stringstream Q2_stream;
    Q2_stream << setprecision(2) << plots[i+offset].Q2;
    TString Q2_string = "Q^{2}=" + Q2_stream.str();
    TLatex* Q2_text = new TLatex(4e-3, unit_scaler*0.083, Q2_string);
    Q2_text->SetTextSize(0.15);
    Q2_text->Draw("Same");

    stringstream beta_stream;
    beta_stream << setprecision(2) << plots[i+offset].beta;
    TString beta_string = "#beta=" + beta_stream.str();
    TLatex* beta_text = new TLatex(4e-3, unit_scaler*0.070, beta_string);
    beta_text->SetTextSize(0.15);
    beta_text->Draw("Same");

    
    float location[4];
    location[0] = 0.5;
    location[1] = 0.7;
    location[2] = 0.75;
    location[3] = 0.9;
    
    /*
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
  TString y_unit_string = "x_{P}F_{2}^{D(3)}";
  TLatex* y_unit_text = new TLatex(0.1*margin_fraction, 1-1.05*margin_fraction, y_unit_string);
  y_unit_text->SetTextSize(0.03);
  y_unit_text->Draw("Same");


  multicanvas->cd(0);
  TString x_unit_string = "x";
  TLatex* x_unit_text = new TLatex(1-1.05*margin_fraction, 0.4*margin_fraction, x_unit_string);
  x_unit_text->SetTextSize(0.03);
  x_unit_text->Draw("Same");
  

  multicanvas->cd(0);

  

  //TText* title = new TText(0.48, 0.9, "Title");
  //title->Draw("Same");

  //TPad *top_pad = new TPad("top_pad", "top", 0, 0.45, 1, 0.9);
  //top_pad->Draw();

  legend->SetTextSize(0.02);
  legend->Draw("Same");

  TString figure_filename = "figures/F2D_data_comparison.pdf";
  multicanvas->Print(figure_filename);

  return 0;
}
