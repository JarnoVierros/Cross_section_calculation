
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

void read_sigma_file(string filename, vector<double> &Q2, vector<double> &x, vector<double> &sigma, vector<double> &sigma_error) {

    ifstream data_file(filename);
    string line;

    bool first_line = true;
    while (getline (data_file, line)) {
        if (first_line) {
            first_line = false;
            continue;
        }
        string Q2_string = "";
        long unsigned int i = 0;
        while (line[i] != ';') {
            Q2_string += line[i];
            i++;
        }
        i++;
        double Q2_value = stod(Q2_string);

        string x_value_string = "";
        while (line[i] != ';') {
            x_value_string += line[i];
            i++;
        }
        i++;
        double x_value = stod(x_value_string);

        string sigma_value_string = "";
        while (line[i] != ';') {
            sigma_value_string += line[i];
            i++;
        }
        i++;
        double sigma_value = stod(sigma_value_string);

        string sigma_error_string = "";
        while (i < line.size()) {
            sigma_error_string += line[i];
            i++;
        }
        double sigma_error_value = stod(sigma_error_string);

        Q2.push_back(Q2_value);
        x.push_back(x_value);
        sigma.push_back(sigma_value);
        sigma_error.push_back(sigma_error_value);
    }
}

void read_LHC_sigma_file(string filename, vector<double> &x, vector<double> &sigma, vector<double> &sigma_error) {

    ifstream data_file(filename);
    string line;

    bool first_line = true;
    while (getline (data_file, line)) {
        if (first_line) {
            first_line = false;
            continue;
        }
        long unsigned int i = 0;
        string x_value_string = "";
        while (line[i] != ';') {
            x_value_string += line[i];
            i++;
        }
        i++;
        double x_value = stod(x_value_string);

        string sigma_value_string = "";
        while (line[i] != ';') {
            sigma_value_string += line[i];
            i++;
        }
        i++;
        double sigma_value = stod(sigma_value_string);

        string sigma_error_string = "";
        while (i < line.size()) {
            sigma_error_string += line[i];
            i++;
        }
        double sigma_error_value = stod(sigma_error_string);

        x.push_back(x_value);
        sigma.push_back(sigma_value);
        sigma_error.push_back(sigma_error_value);
    }
}

void read_differential_sigma_file(string filename, vector<double> &Q2, vector<double> &beta, vector<double> &x, vector<double> &sigma, vector<double> &sigma_error, vector<double> &fit) {

    ifstream data_file(filename);
    string line;

    bool first_line = true;
    while (getline (data_file, line)) {
        if (first_line) {
            first_line = false;
            continue;
        }

        string Q2_string = "";
        long unsigned int i = 0;
        while (line[i] != ';') {
            Q2_string += line[i];
            i++;
        }
        i++;
        double Q2_value = stod(Q2_string);

        string beta_string = "";
        while (line[i] != ';') {
            beta_string += line[i];
            i++;
        }
        i++;
        double beta_value = stod(beta_string);

        string x_value_string = "";
        while (line[i] != ';') {
            x_value_string += line[i];
            i++;
        }
        i++;
        double x_value = stod(x_value_string);

        string sigma_value_string = "";
        while (line[i] != ';') {
            sigma_value_string += line[i];
            i++;
        }
        i++;
        double sigma_value = stod(sigma_value_string);

        string sigma_error_string = "";
        while (line[i] != ';') {
            sigma_error_string += line[i];
            i++;
        }
        i++;
        double sigma_error_value = stod(sigma_error_string);

        string fit_string = "";
        while (i < line.size()) {
            fit_string += line[i];
            i++;
        }
        double fit_value = stod(fit_string);

        Q2.push_back(Q2_value);
        beta.push_back(beta_value);
        x.push_back(x_value);
        sigma.push_back(sigma_value);
        sigma_error.push_back(sigma_error_value);
        fit.push_back(fit_value);
    }
}

void split_by_Q2(vector<double> &new_Q2, vector<vector<double>> &new_x, vector<vector<double>> &new_sigma, vector<vector<double>> &new_sigma_error, vector<double> &Q2, vector<double> &x, vector<double> &sigma, vector<double> &sigma_error) {
    vector<double> current_x, current_sigma, current_sigma_error;
    double current_Q2 = Q2[0];
    for (long unsigned int i=0; i<size(Q2); i++) {
        if (Q2[i] != current_Q2) {
            new_x.push_back(current_x);
            new_sigma.push_back(current_sigma);
            new_sigma_error.push_back(current_sigma_error);
            new_Q2.push_back(current_Q2);
            current_x.clear();
            current_sigma.clear();
            current_sigma_error.clear();
            current_Q2 = Q2[i];
        }
        current_x.push_back(x[i]);
        current_sigma.push_back(sigma[i]);
        current_sigma_error.push_back(sigma_error[i]);
    }
    new_x.push_back(current_x);
    new_sigma.push_back(current_sigma);
    new_sigma_error.push_back(current_sigma_error);
    new_Q2.push_back(current_Q2);
}