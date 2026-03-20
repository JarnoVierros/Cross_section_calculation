

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include "gsl/gsl_pow_int.h"

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include <thread>

#include "direct_dipole_amp_reader.h"

using namespace std;


const double alpha_em = 1.0/137;
const int N_c = 3;
const double C_A = 3.0;
const double C_F = 4.0/3;

const bool is_charm = false;

static double e_f;
static double m_f;

const int warmup_calls = 10000;
const int integration_calls = 100000;//100000000
const int integration_iterations = 1;

const double absolute_error = 0;
const double relative_error = 0.001;
const int max_iterations = 1000000;

string filename_end = "";
bool diffraction_dipamp = true;
const bool detailed = true;

const int i_start = 0; // number of data points to skip
const int data_inclusion_count = 287;

static array<array<array<array<double, 4>, Phi_xpom_size>, Phi_y_size>, Phi_P2_size> Phi_table;
static array<array<array<array<double, 4>, Phi2A_xpom_size>, Phi2A_y_size>, Phi2A_P2_size> Phi2A_table;

static InterpMultilinear<3, double>* Phi_interpolator;
static InterpMultilinear<3, double>* Phi2A_interpolator;

const string Phi_filename = "output/Phi_table_p_bk.txt";
const string Phi2A_filename = "output/Phi_table_p_bk_2A.txt";

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

double Phi(double P2, double y, double xpom) {
    array<double, 3> args = {log(P2), log(y), log(xpom)};
    double interp = Phi_interpolator->interp(args.begin());
    if (interp == 0) {
        return 0;
    } else {
        return exp(interp);
    };
    
}

double Phi2A(double P2, double y, double xpom) {
    array<double, 3> args = {log(P2), log(y), log(xpom)};
    double interp = Phi2A_interpolator->interp(args.begin());
    if (interp == 0) {
        return 0;
    } else {
        return exp(interp);
    };
}

//////////////////////////////////////////

struct subintegration_parameters {
    double K2;
    double beta;
    double xpom;
    double Q2;
    double sigma;
    double sigma_error;
    double sigma_fit;
};

struct integration_parameters {
    double beta;
    double xpom;
    double Q2;
    double sigma;
    double sigma_error;
    double sigma_fit;
};

double integrand2(double z, double K2, double beta, double xpom, double Q2) {
    //cout << "part1=" << 1/(1-z)*(1+z*z)*Phi(K2, beta/z, xpom) << endl;
    //cout << "part2=" << 1/(1-z)*(-2*Phi(K2, beta, xpom)) << endl;
    return 1/(1-z)*((1+z*z)*Phi(K2, beta/z, xpom) - 2*Phi(K2, beta, xpom));
}

double wrong_integrand2(double z, double K2, double beta, double xpom, double Q2) {
    //cout << "part1=" << 1/(1-z)*(1+z*z)*Phi(K2, beta/z, xpom) << endl;
    //cout << "part2=" << 1/(1-z)*(-2*Phi(K2, beta, xpom)) << endl;
    return (1+z*z)*Phi(K2, beta/z, xpom);
}

double f2(double z, void * params) {
    struct subintegration_parameters *par = (struct subintegration_parameters *)params;
    return integrand2(z, par->K2, par->beta, par->xpom, par->Q2);
}

double integrand1(double K2, double beta, double xpom, double Q2) {

    double res, err;
    int n = max_iterations;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);

    subintegration_parameters params2;
    params2.K2 = K2;
    params2.beta = beta;
    params2.xpom = xpom;
    params2.Q2 = Q2;

    gsl_function F2;
    F2.function = &f2;
    F2.params = &params2;

    gsl_integration_qags(&F2, beta, 1, absolute_error, relative_error, n, w, &res, &err);
    //cout << res << " pm " << err << ", rel=" << err/res << endl;
    //cout << "part3=" << (2*log(1-beta)+3.0/2)*Phi(K2, beta, xpom) << endl;
    gsl_integration_workspace_free(w);


    double factor = log(Q2/(beta*K2));
    double part_1 = res;
    double part_2 = (2*log(1-beta)+3.0/2)*Phi(K2, beta, xpom);

    return factor*(part_1 + part_2);
}

double wrong_integrand1(double K2, double beta, double xpom, double Q2) {

    double res, err;
    int n = max_iterations;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);

    subintegration_parameters params2;
    params2.K2 = K2;
    params2.beta = beta;
    params2.xpom = xpom;
    params2.Q2 = Q2;

    gsl_function F2;
    F2.function = &f2;
    F2.params = &params2;

    gsl_integration_qaws_table* qaws_table = gsl_integration_qaws_table_alloc(0, 1, 0, 0);

    gsl_integration_qaws(&F2, beta, 1, qaws_table, absolute_error, relative_error, n, w, &res, &err);
    //cout << res << " pm " << err << ", rel=" << err/res << endl;
    //cout << "part3=" << (2*log(1-beta)+3.0/2)*Phi(K2, beta, xpom) << endl;
    gsl_integration_workspace_free(w);


    double factor = log(Q2/(beta*K2));
    double part_1 = res;
    double part_2 = 3.0/2*(1.0/beta*Phi(K2, 1, xpom)-Phi(K2, beta, xpom));

    return factor*(part_1 + part_2);
}

double f1(double K2, void * params) {
    struct integration_parameters *par = (struct integration_parameters *)params;
    return integrand1(K2, par->beta, par->xpom, par->Q2);
}

double xpomFqq(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

    double normalization = 2*N_c*C_F*alpha_em/pow(2*M_PI, 6)*e_f*e_f;

    double res, err;

    int n = max_iterations;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (n);

    integration_parameters params1;
    params1.beta = beta;
    params1.xpom = xpom;
    params1.Q2 = Q2;

    gsl_function F1;
    F1.function = &f1;
    F1.params = &params1;

    gsl_integration_qags(&F1, 0, Q2, absolute_error, relative_error, n, w, &res, &err);

    result = normalization*res;
    error = normalization*err;
    //cout << "xpomFqq: Q²=" << params1.Q2 << ", xpom=" << params1.xpom << ", beta=" << params1.beta << ", res: " << result << ", err: " << error << endl;

    gsl_integration_workspace_free(w);

    return 0;
}

/////////////////


double xpomFgq_integrand(double beta, double xpom, double Q2, double z, double K2) {
    double normalization = 2*N_c*C_F*alpha_em/pow(2*M_PI, 6)*e_f*e_f;
    double integrand = normalization*(z*z+gsl_pow_2(1-z))*log(Q2/(beta*K2))*Phi2A(K2, beta/z, xpom);
    return integrand;
}

double integration_function_gq(double *k, size_t dim, void * params) {
  double z = k[0];
  double K2 = k[1];
  struct integration_parameters *par = (struct integration_parameters *)params;

  return xpomFgq_integrand(par->beta, par->xpom, par->Q2, z, K2);
}

double xpomFgq(double beta, double xpom, double Q2, double &result, double &error, double &fit) {

  const int dim = 2;
  double res, err;

  double xl[dim] = {beta, 0};
  double xu[dim] = {1, Q2};

  struct integration_parameters params = {1, 1, 1};
  params.beta = beta;
  params.xpom = xpom;
  params.Q2 = Q2;

  const gsl_rng_type *qg_rng;
  gsl_rng *rng;

  gsl_monte_function qg_rng_G = {&integration_function_gq, dim, &params};

  gsl_rng_env_setup();
  int status = 0;

  qg_rng = gsl_rng_default;
  rng = gsl_rng_alloc(qg_rng);
  gsl_rng_set(rng, 1);

  gsl_monte_vegas_state *T_s = gsl_monte_vegas_alloc(dim);
  status = gsl_monte_vegas_integrate(&qg_rng_G, xl, xu, dim, warmup_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qg_warmup_error: " << status << endl; throw (status);}
  status = gsl_monte_vegas_integrate(&qg_rng_G, xl, xu, dim, integration_calls, rng, T_s, &res, &err);
  if (status != 0) {cout << "qg_integration_error: " << status << endl; throw (status);}

  if (gsl_isnan(res)) {
    res = 0;
    cout << "nan found at xpom=" << params.xpom << endl;
  }
  result = res;
  error = err;
  fit = gsl_monte_vegas_chisq(T_s);
  //cout << "xpomFgq: Q²=" << params.Q2 << ", xpom=" << params.xpom << ", beta=" << params.beta << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(T_s) << endl;

  gsl_monte_vegas_free(T_s);

  return 0;
}

int calc_total_xpomF_Tqqg_contribution(double beta, double xpom, double Q2, double &result, double &error, double &fit, double &gq_result, double &qq_result) {
    double fits[2];

    double gq_error, gq_fit;
    xpomFgq(beta, xpom, Q2, gq_result, gq_error, gq_fit);

    double qq_error, qq_fit;
    xpomFqq(beta, xpom, Q2, qq_result, qq_error, qq_fit);
    qq_fit = 0;

    result  = gq_result + 2*qq_result;
    error = sqrt(gsl_pow_2(gq_error) + gsl_pow_2(2*qq_error));
    fits[0] = gq_fit;
    fits[1] = qq_fit;
    if (gq_fit > qq_fit) {
        fit = gq_fit;
    } else {
        fit = qq_fit;
    }

    cout << "gq=" << gq_result << "±" << gq_error << ", qq=" << qq_result << "±" << qq_error << ", total=" << result << "±" << error << endl;

    return 0;
}

struct thread_par_struct
{
  double Q2;
  double xpom;
  double beta;
  double &result;
  double &error;
  double &fit;
  double &gq_result;
  double &qq_result;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6, double &a7, double &a8) : Q2(a1), xpom(a2), beta(a3), result(a4), error(a5), fit(a6), gq_result(a7), qq_result(a8) {}
};

int low_beta_thread_func(thread_par_struct par) {
  calc_total_xpomF_Tqqg_contribution(par.beta, par.xpom, par.Q2, par.result, par.error, par.fit, par.gq_result, par.qq_result);
  return 0;
}

int main() {

  
  if (is_charm) {
    e_f = 2.0/3;
    m_f = 1.27;
    filename_end += "_charm";
  } else {
    e_f = sqrt(2.0/3*2.0/3+2*1.0/3*1.0/3);
    m_f = 0;
  }
  

  //e_f = sqrt(2*2.0/3*2.0/3+2*1.0/3*1.0/3);

  gsl_set_error_handler_off();


  vector<double> Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values;

  read_data_file("data/expanded_differential_HERA_data.dat", Q2_values, beta_values, x_values, x_pom_F2_values, delta_stat_values, delta_sys_values);

  cout << "loading Phi table" << endl;
  load_Phi_table(Phi_table, Phi_filename);
  cout << "creating Phi interpolator" << endl;
  create_Phi_interpolator(Phi_table, Phi_interpolator);
  cout << "loading Phi2A table" << endl;
  load_Phi2A_table(Phi2A_table, Phi2A_filename);
  cout << "creating Phi2A interpolator" << endl;
  create_Phi2A_interpolator(Phi2A_table, Phi2A_interpolator);
  cout << "Phi tables loaded and interpolators created" << endl;

  /*
  double Q2 = 4.5;
  double beta = 0.04;
  double x = 0.00012;
  double resulta, errora, fita;
  calc_total_xpomF_Tqqg_contribution(beta, x, Q2, resulta, errora, fita);
  cout << "x=" << x << ", Phi=" << resulta << endl;

  for (int i=0; i<100; i++) {
    double xpom = 1e-5 + i/99.0*0.01;
    double result, error, fit;
    calc_total_xpomF_Tqqg_contribution(beta, xpom*beta, Q2, result, error, fit);
    cout << "xpom=" << xpom << ", Phi=" << result << endl;
  }
  */

  thread low_beta_integration_threads[data_inclusion_count];
  double result[data_inclusion_count], error[data_inclusion_count], fit[data_inclusion_count], gq_result[data_inclusion_count], qq_result[data_inclusion_count];

  static auto t1 = chrono::high_resolution_clock::now();

  for (int i=0; i<data_inclusion_count; i++) {

    thread_par_struct low_beta_par(Q2_values[i+i_start], x_values[i+i_start]/beta_values[i+i_start], beta_values[i+i_start], result[i], error[i], fit[i], gq_result[i], qq_result[i]);
    low_beta_integration_threads[i] = thread(low_beta_thread_func, low_beta_par);

  }

  for (int i=0; i<data_inclusion_count; i++) {
    low_beta_integration_threads[i].join();
  }

  static auto t2 = chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
  cout << "Calculation finished in " << duration.count() << " seconds" << endl;

  if (!detailed) {

    ofstream L_output_file("output/J_low_beta_corrections"+filename_end+".txt");
    L_output_file << "Q2 (GeV);beta;x;sigma (mb);sigma error (mb);fit" << endl;

    for (int i=0; i<data_inclusion_count; i++) {
      ostringstream Q2;
      Q2 << Q2_values[i+i_start];
      ostringstream beta;
      beta << beta_values[i+i_start];
      ostringstream x;
      x << x_values[i+i_start];
      ostringstream result_stream;
      result_stream << result[i];
      ostringstream error_stream;
      error_stream << error[i];
      ostringstream fit_stream;
      fit_stream << fit[i];
      string line = Q2.str() + ";" + beta.str() + ";" + x.str() + ";" + result_stream.str() + ";" + error_stream.str() + ";" + fit_stream.str();
      L_output_file << line << endl;
    }
    L_output_file.close();
  } else {
    ofstream L_output_file("output/J_low_beta_corrections"+filename_end+"_detailed.txt");
    L_output_file << "Q2 (GeV);beta;x;gq;qq;total;error;fit" << endl;

    for (int i=0; i<data_inclusion_count; i++) {
      ostringstream Q2;
      Q2 << Q2_values[i+i_start];
      ostringstream beta;
      beta << beta_values[i+i_start];
      ostringstream x;
      x << x_values[i+i_start];
      ostringstream gq;
      gq << gq_result[i+i_start];
      ostringstream qq;
      qq << qq_result[i+i_start];
      ostringstream result_stream;
      result_stream << result[i];
      ostringstream error_stream;
      error_stream << error[i];
      ostringstream fit_stream;
      fit_stream << fit[i];
      string line = Q2.str() + ";" + beta.str() + ";" + x.str() + ";" + gq.str() + ";" + qq.str() + ";" + result_stream.str() + ";" + error_stream.str() + ";" + fit_stream.str();
      L_output_file << line << endl;
    }
    L_output_file.close();
  }
  return 0;
}

