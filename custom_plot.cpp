#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"

#include <iostream>

using namespace std;

const double alpha_em = 1.0/137;
const double s = 320;

int main() {
    int array_size = 3;
    int prediction_array_size = 5;
    double Q2 =35;
    double x_pom[array_size] = {0.004, 0.010, 0.018};
    double beta[array_size] = {0.25, 0.10, 0.04};
    double prediction_x_pom[prediction_array_size] = {0.004, 0.007, 0.010, 0.014, 0.018};
    double prediction_beta[prediction_array_size] = {0.25, 0.175, 0.10, 0.07, 0.04};
    double measurement_values[array_size] = {1.5, 0.63, 0.62};
    double relative_measurement_error[array_size] = {25+27, 23+29, 18+47};
    double measurement_error[array_size] = {relative_measurement_error[0]/100*measurement_values[0], relative_measurement_error[1]/100*measurement_values[1], relative_measurement_error[2]/100*measurement_values[2]};
    double predicted_T[prediction_array_size] = {3.16e-8, 1.3e-8, 3.87e-9, 1.94e-9, 6.4e-10};
    double predicted_L[prediction_array_size] = {5.3e-10, 2e-10, 6e-11, 3e-11, 1.95696e-11};
    double x_errors[array_size] = {0, 0, 0};

    double prediction[prediction_array_size];
    for (int i=0; i<prediction_array_size; i++) {
        double FL = 1/prediction_x_pom[i]*Q2/(pow(2*M_PI, 2)*alpha_em)*Q2/prediction_beta[i]*predicted_L[i];
        double FT = 1/prediction_x_pom[i]*Q2/(pow(2*M_PI, 2)*alpha_em)*Q2/prediction_beta[i]*predicted_T[i];
        double F2 = FL + FT;
        double M_X2 = Q2*(1/prediction_beta[i]-1);
        double y = (Q2 + M_X2)/(s*prediction_x_pom[i]);
        double sigma_r = F2 - y*y/(1+pow(1-y, 2))*FL;
        double correction = 12;
        prediction[i] = correction*sigma_r;
        cout << sigma_r << endl;
    }

    TMultiGraph* multigraph = new TMultiGraph();

    TGraphErrors* mesurement_data = new TGraphErrors(3, x_pom, measurement_values, x_errors, measurement_error);
    mesurement_data->SetMarkerStyle(7);

    TGraph* prediction_graph = new TGraph(prediction_array_size, prediction_x_pom, prediction);
    
    multigraph->Add(mesurement_data, "P");
    multigraph->Add(prediction_graph, "C");

    multigraph->SetTitle("Diffractive reduced c#bar{c} cross section #tilde{#sigma}_{D}^{c#bar{c}}(x_{pom}, #beta, Q^{2})");  

    TCanvas* MC = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);

    multigraph->GetYaxis()->SetTitle("#tilde{#sigma}_{D}^{c#bar{c}}");
    multigraph->GetXaxis()->SetTitle("x_{pom}");
    multigraph->GetYaxis()->SetRangeUser(0, 2.4);

    multigraph->Draw("A");

    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->SetTextSize(0.045);
    legend->AddEntry(mesurement_data,"H1 data");
    legend->AddEntry(prediction_graph,"BK prediction");
    legend->Draw();

    TLatex* beta_text_1 = new TLatex(0.0045, 1.6, "#beta=0.25");
    beta_text_1->Draw("Same");
    TLatex* Q2_text_1 = new TLatex(0.0045, 1.8, "Q^{2}=35 GeV^{2}");
    Q2_text_1->Draw("Same");

    TLatex* beta_text_2 = new TLatex(0.0105, 0.6, "#beta=0.10");
    beta_text_2->Draw("Same");
    TLatex* Q2_text_2 = new TLatex(0.0105, 0.8, "Q^{2}=35 GeV^{2}");
    Q2_text_2->Draw("Same");

    TLatex* beta_text_3 = new TLatex(0.0155, 0.9, "#beta=0.04");
    beta_text_3->Draw("Same");
    TLatex* Q2_text_3 = new TLatex(0.0155, 1.1, "Q^{2}=35 GeV^{2}");
    Q2_text_3->Draw("Same");

    MC->Print("figures/corrected_F2Dcharm_data_comparison.pdf");

}