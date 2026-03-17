

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include <thread>

#include "direct_dipole_amp_reader.h"
#include "gsl/gsl_pow_int.h"

using namespace std;

static array<array<array<array<array<double, 5>, 81>, 30>, 30>, 30> p_table;
static array<array<array<array<array<double, 5>, 81>, 40>, 40>, 40> Pb_table;

static InterpMultilinear<4, double>* interpolator;

const string dipole_amp_type = "bk";
const string nucleus_type = "p";
const bool diffraction_dipamp = true;
const int Phi_n = 2;
const bool A_dipole = true;
const string filename_end = "_p_bk_2A";


// 50000000*50*20*20 too much

const int warmup_calls = 100; //10000
const int integration_calls = 1000; //100000
const int integration_iterations = 1;

const int P2_resolution = 30; //100
const double min_P2 = 0.01;
const double max_P2 = 100;

const int y_resolution = 30; //50
const double y_min = 0.01;
const double y_max = 1;

const int xpom_resolution = 30; //50
const double xpom_start = 1e-4;
const double xpom_stop = 0.01;

static double r_limit; // 34.64
static double b_min_limit; // 17.32

//const double r_radius = 100;
//const double b_radius = 100;

const double C_A = 3;
const double C_F = 4.0/3;

struct thread_par_struct
{
  double P2;
  double y;
  double xpom;
  double &result;
  double &error;
  double &fit;
  thread_par_struct(double a1, double a2, double a3, double &a4, double &a5, double &a6) : P2(a1), y(a2), xpom(a3), result(a4), error(a5), fit(a6) {}
};

struct parameters {double P2; double y; double xpom;};

double dipole_amplitude(double r, double b_min, double phi, double xpom) {
    if (calc_max_phi(r, b_min) < phi) {
        return 0;
    } else {
        array<double, 4> args = {log(r), log(b_min), phi, log(xpom)};
        return exp(interpolator->interp(args.begin()));
    }
}

double sub_integrand(double j, double r, double R, double v, double b, double theta, double P2, double y, double xpom) {
    double P = sqrt(P2);

    double integrand = gsl_sf_bessel_J0(sqrt(gsl_pow_2(r*cos(j+v)-R*cos(v-j))+gsl_pow_2(r*sin(j+v)-R*sin(v-j)))*P*sqrt(1-y));
    //cout << "1: " << gsl_sf_bessel_J0(r_minus_R*P*sqrt(1-y)) << endl;
    integrand *= gsl_sf_bessel_Kn(Phi_n, r*P*sqrt(y))*gsl_sf_bessel_Kn(Phi_n, R*P*sqrt(y));
    //cout << "r=" << r << ", R=" << R << ", P=" << P << ", y=" << y << endl;
    //cout << "2: " << gsl_sf_bessel_Kn(Phi_n, r*P*sqrt(y))*gsl_sf_bessel_Kn(Phi_n, R*P*sqrt(y)) << endl;
    //integrand *= cos(Phi_n*j);
    //cout << "3: " << cos(Phi_n*(phi_r-phi_R)) << endl;

    double eff_r = r;
    double eff_bmin = sqrt(b*b+r*b*cos(j+v-theta)+1.0/4*r*r);
    double eff_phi = acos(-(b*cos(j+v-theta)+1.0/2*r)/eff_bmin);

    double N_r = dipole_amplitude(eff_r, eff_bmin, eff_phi, xpom);


    eff_r = R;
    eff_bmin = sqrt(b*b+R*b*cos(v-j-theta)+1.0/4*R*R);
    eff_phi = acos(-(b*cos(v-j-theta)+1.0/2*R)/eff_bmin);

    double N_R = dipole_amplitude(eff_r, eff_bmin, eff_phi, xpom);

    if (A_dipole) {
        N_r = 1 - pow(1-N_r, C_A/C_F);
        N_R = 1 - pow(1-N_R, C_A/C_F);
    }

    integrand *= N_r*N_R;
    //cout << "4: " << N_r*N_R << endl;

    return integrand;
}

struct subint_parameters {double r; double R; double v; double b; double theta; double P2; double y; double xpom;};

double f(double j, void * params) {
    struct subint_parameters *par = (struct subint_parameters *)params;
    return sub_integrand(j, par->r, par->R, par->v, par->b, par->theta, par->P2, par->y, par->xpom);
}
//static int phi_counter = 0;
double Phi_integrand(double r, double R, double v, double b, double theta, double P2, double y, double xpom) {
    
    int n = 1000;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(n);
    gsl_integration_qawo_table* qawo_table = gsl_integration_qawo_table_alloc(Phi_n, 2*M_PI, GSL_INTEG_COSINE, n);

    subint_parameters subint_params;
    subint_params.r = r;
    subint_params.R = R;
    subint_params.v = v;
    subint_params.b = b;
    subint_params.theta=theta;
    subint_params.P2 = P2;
    subint_params.y=y;
    subint_params.xpom=xpom;

    gsl_function F;
    F.function = &f;
    F.params = &subint_params;

    double result, error;
    gsl_integration_qawo(&F, 0, 0, 0.01, n, workspace, qawo_table, &result, &error);

    //cout << "sub_result=" << result << ", sub_error=" << error << ", rer_error=" << error/result << endl;

    gsl_integration_workspace_free(workspace);
    gsl_integration_qawo_table_free(qawo_table);

    //phi_counter += 1;
    //cout << phi_counter << endl;

    return result;

}

double g(double *k, size_t dim, void * params) {
    struct parameters *par = (struct parameters *)params;
    return Phi_integrand(k[0], k[1], k[2], k[3], k[4], par->P2, par->y, par->xpom);
}

void integrate(thread_par_struct par) {

    const int dim = 5;
    double res, err;

    double xl[dim] = {0, 0, 0, 0, 0};
    double xu[dim] = {r_limit, r_limit, 2*M_PI, b_min_limit, 2*M_PI};

    struct parameters params = {1, 1};
    params.P2 = par.P2;
    params.y = par.y;
    params.xpom = par.xpom;
    double &result = par.result;
    double &error = par.error;
    double &fit = par.fit;

    const gsl_rng_type *T;
    gsl_rng *rng;
    gsl_monte_function G = {&g, dim, &params};

    gsl_rng_env_setup();
    int status = 0;

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, 1);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
    status = gsl_monte_vegas_integrate(&G, xl, xu, dim, warmup_calls, rng, s, &res, &err);
    double duration_double;
    if (status != 0) {cout << "integrate error: " << status << endl; throw (status);}
    for (int i=0; i<integration_iterations; i++) {
        static auto t1 = chrono::high_resolution_clock::now();
        status = gsl_monte_vegas_integrate(&G, xl, xu, dim, integration_calls, rng, s, &res, &err);
        if (status != 0) {cout << "integrate error: " << status << endl; throw (status);}
        static auto t2 = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
        duration_double = duration.count();
        cout << "result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(s) << ", duration: " << duration.count() << " seconds" << endl;
    }
    //phi_counter = 0;
    if (gsl_isnan(res)) {
        res = 0;
        cout << "nan found at xpom=" << params.xpom << endl;
    }
    result = res;
    error = err;
    fit = gsl_monte_vegas_chisq(s);
    //cout << "P2=" << params.P2 << ", y=" << params.y << ", xpom=" << params.xpom << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(s) << ", duration: " << duration_double << " seconds" << endl;

    gsl_monte_vegas_free(s);

}

int main() {
    
    gsl_set_error_handler_off();

    string filename;
    if (diffraction_dipamp) {
    filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+"_diffraction.csv";
    } else {
    filename = "data/dipole_amplitude_with_IP_dependence_"+dipole_amp_type+"_"+nucleus_type+".csv";
    }
    if (nucleus_type == "p") {
    load_p_dipole_amplitudes(p_table, filename);
    create_p_interpolator(p_table, interpolator);
    } else if (nucleus_type == "Pb") {
    load_Pb_dipole_amplitudes(Pb_table, filename);
    create_Pb_interpolator(Pb_table, interpolator);
    } else {
    throw 1;
    }

    if (nucleus_type == "Pb") {
        r_limit = 657; // 34.64, 657
        b_min_limit = 328; // 17.32, 328
    } else if (nucleus_type == "p") {
        r_limit = 34.64;
        b_min_limit = 17.32;
    } else {
        cout << "invalid nucleus type" << endl;
        throw 1;
    }

    vector<double> P2_steps, y_steps;

    for (int i=0; i<P2_resolution; i++) {
        P2_steps.push_back(min_P2 + 1.0*i/(P2_resolution-1)*(max_P2-min_P2));
    }
    for (int i=0; i<y_resolution; i++) {
        y_steps.push_back(y_min + 1.0*i/(y_resolution-1)*(y_max-y_min));
    }

    const double xpom_step = 1.0/(xpom_resolution-1)*log10(xpom_stop/xpom_start);

    vector<double> xpom_steps;

    for (int i=0; i<xpom_resolution; i++) {
        double xpom = pow(10, log10(xpom_start) + i*xpom_step);
        xpom_steps.push_back(xpom);
    }


    //thread integration_threads[P2_resolution*y_resolution*xpom_resolution];
    //double results[P2_resolution*y_resolution*xpom_resolution], errors[P2_resolution*y_resolution*xpom_resolution], fits[P2_resolution*y_resolution*xpom_resolution];

    TString out_filename = "output/Phi_table" + filename_end + ".txt";
    ofstream output_file(out_filename);
    output_file << "P2;y;xpom;result;error;fit" << endl;
    output_file.close();


    static auto t1 = chrono::high_resolution_clock::now();

    vector<thread> threads;
    int max_threads = 15;
    int current_threads = 0;

    double results[max_threads], errors[max_threads], fits[max_threads];
    double rP2[max_threads], ry[max_threads], rxpom[max_threads];

    for (int i=0; i<P2_resolution; i++) {
        for (int j=0; j<y_resolution; j++) {
            for (int k=0; k<xpom_resolution; k++) {
                if (current_threads >= max_threads) {
                    for (int l=0; l<max_threads; l++) {
                        threads[l].join();
                    }
                    //cout << "clearing threads, threads: " << current_threads << ", i=" << i << ", j=" << j << ", k=" << k << endl;
                    int index = i*y_resolution*xpom_resolution + j*xpom_resolution + k;
                    auto now = chrono::high_resolution_clock::now();
                    auto duration = chrono::duration_cast<chrono::seconds>(now-t1);
                    cout << index << "/" << P2_resolution*y_resolution*xpom_resolution << " integrals finished in " << duration.count()/60.0 << " minutes" << endl;
                    threads.clear();
                    current_threads = 0;

                    ofstream output_file;
                    output_file.open(out_filename, ios::app);
                    //output_file << "P2;y;xpom;result;error;fit" << endl;

                    for (int l=0; l<max_threads; l++) {

                        ostringstream P2;
                        P2 << rP2[l];
                        ostringstream y;
                        y << ry[l];
                        ostringstream xpom;
                        xpom << rxpom[l];
                        ostringstream result;
                        result << results[l];
                        ostringstream error;
                        error << errors[l];
                        ostringstream fit;
                        fit << fits[l];
                        string line = P2.str() + ";" + y.str() + ";" + xpom.str() + ";" + result.str() + ";" + error.str() + ";" + fit.str();
                        output_file << line << endl;

                    }
                    output_file.close();
                }

                //int index = i*y_resolution*xpom_resolution + j*xpom_resolution + k;
                thread_par_struct par(P2_steps[i], y_steps[j], xpom_steps[k], results[current_threads], errors[current_threads], fits[current_threads]);
                rP2[current_threads] = P2_steps[i];
                ry[current_threads] = y_steps[j];
                rxpom[current_threads] = xpom_steps[k];
                threads.push_back(thread(integrate, par));
                current_threads += 1;
                //cout << "added thread, threads: " << current_threads << endl;

            }
        }
    }


    for (int l=0; l<current_threads; l++) {
        threads[l].join();
    }
    //cout << "final clearing of threads, threads: " << current_threads << endl;
    threads.clear();

    output_file.open(out_filename, ios::app);
    //output_file << "P2;y;xpom;result;error;fit" << endl;

    for (int l=0; l<current_threads; l++) {

        ostringstream P2;
        P2 << rP2[l];
        ostringstream y;
        y << ry[l];
        ostringstream xpom;
        xpom << rxpom[l];
        ostringstream result;
        result << results[l];
        ostringstream error;
        error << errors[l];
        ostringstream fit;
        fit << fits[l];
        string line = P2.str() + ";" + y.str() + ";" + xpom.str() + ";" + result.str() + ";" + error.str() + ";" + fit.str();
        output_file << line << endl;

    }
    output_file.close();

    /*
                TString out_filename = "output/Phi_table" + filename_end + ".txt";
                ofstream output_file;
                output_file.open(out_filename, ios::app);
                //output_file << "P2;y;xpom;result;error;fit" << endl;

                for (int i=0; i<P2_resolution; i++) {
                    for (int j=0; j<y_resolution; j++) {
                        for (int k=0; k<xpom_resolution; k++) {
                            int index = i*y_resolution*xpom_resolution + j*xpom_resolution + k;
                            ostringstream P2;
                            P2 << P2_steps[i];
                            ostringstream y;
                            y << y_steps[j];
                            ostringstream xpom;
                            xpom << xpom_steps[k];
                            ostringstream result;
                            result << results[index];
                            ostringstream error;
                            error << errors[index];
                            ostringstream fit;
                            fit << fits[index];
                            string line = P2.str() + ";" + y.str() + ";" + xpom.str() + ";" + result.str() + ";" + error.str() + ";" + fit.str();
                            output_file << line << endl;
                        }
                    }
                }
                output_file.close();
                }
    */

/*
    for (int i=0; i<P2_resolution*y_resolution*xpom_resolution; i++) {
        integration_threads[i].join();
    }
*/
    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    cout << "Calculation finished in " << duration.count()/60.0 << " minutes" << endl;
    return 0;
}

