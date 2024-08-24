#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"

int main() {
    double Q2 =35;
    double x_pom[] = {0.004, 0.010, 0.018};
    double beta[] = {0.25, 0.10, 0.04};
    double x[] = {x_pom[0]*beta[0], x_pom[1]*beta[1], x_pom[2]*beta[2]};
    double measurement_values[] = {1.5, 0.63, 0.62};
    double relative_measurement_error[] = {25+27, 23+29, 18+47};
    double measurement_error[] = {relative_measurement_error[0]/100*measurement_values[0], relative_measurement_error[1]/100*measurement_values[1], relative_measurement_error[2]/100*measurement_values[2]};
    
    double predicted_T = {};
    double predicted_L = {};
    double prediction[] = {1.441945751106215, 0.234253640879954, 0.05283514645762909};
    double x_errors[] = {0, 0, 0};

    TMultiGraph* multigraph = new TMultiGraph();

    TGraphErrors* mesurement_data = new TGraphErrors(3, x_pom, measurement_values, x_errors, measurement_error);
    mesurement_data->SetMarkerStyle(7);

    TGraph* prediction_graph = new TGraph(3, x_pom, prediction);
    
    multigraph->Add(mesurement_data, "P");
    multigraph->Add(prediction_graph, "PL");

    multigraph->SetTitle("Diffractive reduced c#bar{c} cross section #tilde{#sigma}_{D}^{c#bar{c}}(x_{pom}, #beta, Q^{2})");  

    TCanvas* MC = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);

    multigraph->GetYaxis()->SetTitle("#tilde{#sigma}_{D}^{c#bar{c}}");
    multigraph->GetXaxis()->SetTitle("x_{pom}");
    multigraph->GetYaxis()->SetRangeUser(0, 2.4);

    multigraph->Draw("A");

    TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
    legend->SetTextSize(0.045);
    legend->AddEntry(mesurement_data,"Measurement data");
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

    MC->Print("F2Dcharm_data_comparison.pdf");

}