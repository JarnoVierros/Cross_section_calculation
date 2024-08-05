#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

int main() {
    double Q2 =35;
    double x_pom[] = {0.004, 0.010, 0.018};
    double beta[] = {0.25, 0.10, 0.04};
    double x[] = {x_pom[0]*beta[0], x_pom[1]*beta[1], x_pom[2]*beta[2]};
    double measurement_values[] = {1.5, 0.63, 0.62};
    double relative_measurement_error[] = {25+27, 23+29, 18+47};
    double measurement_error[] = {relative_measurement_error[0]/100*measurement_values[0], relative_measurement_error[1]/100*measurement_values[1], relative_measurement_error[2]/100*measurement_values[2]};
    double prediction[] = {1.441945751106215, 0.234253640879954, 0.051945612798961774};
    double x_errors[] = {0, 0, 0};

    TMultiGraph* multigraph = new TMultiGraph();

    TGraphErrors* mesurement_data = new TGraphErrors(3, x, measurement_values, x_errors, measurement_error);
    mesurement_data->SetMarkerStyle(7);

    TGraph* prediction_graph = new TGraph(3, x, prediction);
    
    multigraph->Add(mesurement_data, "P");
    multigraph->Add(prediction_graph, "PC");

    TCanvas* MC = new TCanvas("diff_T_sigma_canvas", "", 1000, 600);

    multigraph->Draw("A");

    MC->Print("F2Dcharm_data_comparison.pdf");

}