

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
const int Phi_n = 1;
const bool A_dipole = false;


// 50000000*50*20*20 too much

const int warmup_calls = 100000;
const int integration_calls = 50000000;
const int integration_iterations = 1;

const int P2_resolution = 50;
const double min_P2 = 0.01;
const double max_P2 = 100;

const int y_resolution = 20;
const double y_min = 0.01;
const double y_max = 1;

const int xpom_resolution = 20;
const double xpom_start = 1e-4;
const double xpom_stop = 1;

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

double Phi_integrand(double r1, double r2, double R1, double R2, double b1, double b2, double P2, double y, double xpom) {
    double r = sqrt(r1*r1+r2*r2);
    double R = sqrt(R1*R1+R2*R2);
    double r_minus_R = sqrt(gsl_pow_2(r1-R1)+gsl_pow_2(r2-R2));
    double P = sqrt(P2);
    double phi_r = atan2(r2, r1);
    double phi_R = atan2(R2, R1);

    double integrand = gsl_sf_bessel_J0(r_minus_R*P*sqrt(1-y));
    //cout << "1: " << gsl_sf_bessel_J0(r_minus_R*P*sqrt(1-y)) << endl;
    integrand *= gsl_sf_bessel_Kn(Phi_n, r*P*sqrt(y))*gsl_sf_bessel_Kn(Phi_n, R*P*sqrt(y));
    //cout << "r=" << r << ", R=" << R << ", P=" << P << ", y=" << y << endl;
    //cout << "2: " << gsl_sf_bessel_Kn(Phi_n, r*P*sqrt(y))*gsl_sf_bessel_Kn(Phi_n, R*P*sqrt(y)) << endl;
    integrand *= cos(Phi_n*(phi_r-phi_R));
    //cout << "3: " << cos(Phi_n*(phi_r-phi_R)) << endl;

    double eff_r1 = r1;
    double eff_r2 = r2;
    double eff_b1 = b1+0.5*r1;
    double eff_b2 = b2+0.5*r2;

    double eff_x1 = eff_b1+0.5*eff_r1;
    double eff_x2 = eff_b2+0.5*eff_r2;
    double eff_y1 = eff_b1-0.5*eff_r1;
    double eff_y2 = eff_b2-0.5*eff_r2;

    double eff_x = sqrt(eff_x1*eff_x1+eff_x2+eff_x2);
    double eff_y = sqrt(eff_y1*eff_y1+eff_y2+eff_y2);

    if (eff_x < eff_y) {
        eff_r1 = -1*eff_r1;
        eff_r2 = -1*eff_r2;
    }

    double eff_r = sqrt(eff_r1*eff_r1+eff_r2*eff_r2);
    double eff_bmin = sqrt(gsl_pow_2(b1+0.5*r1)+gsl_pow_2(b2+0.5*r2));
    double eff_phi = acos(-(eff_r1*(eff_b1+0.5*eff_r1)+eff_r2*(eff_b2+0.5*eff_r2))/(eff_bmin*eff_r));

    double N_r = dipole_amplitude(eff_r, eff_bmin, eff_phi, xpom);



    eff_r1 = R1;
    eff_r2 = R2;
    eff_b1 = b1+0.5*R1;
    eff_b2 = b2+0.5*R2;

    eff_x1 = eff_b1+0.5*eff_r1;
    eff_x2 = eff_b2+0.5*eff_r2;
    eff_y1 = eff_b1-0.5*eff_r1;
    eff_y2 = eff_b2-0.5*eff_r2;

    eff_x = sqrt(eff_x1*eff_x1+eff_x2+eff_x2);
    eff_y = sqrt(eff_y1*eff_y1+eff_y2+eff_y2);

    if (eff_x < eff_y) {
        eff_r1 = -1*eff_r1;
        eff_r2 = -1*eff_r2;
    }

    eff_r = sqrt(eff_r1*eff_r1+eff_r2*eff_r2);
    eff_bmin = sqrt(gsl_pow_2(b1+0.5*r1)+gsl_pow_2(b2+0.5*r2));
    eff_phi = acos(-(eff_r1*(eff_b1+0.5*eff_r1)+eff_r2*(eff_b2+0.5*eff_r2))/(eff_bmin*eff_r));

    double N_R = dipole_amplitude(eff_r, eff_bmin, eff_phi, xpom);

    if (A_dipole) {
        N_r = 1 - pow(1-N_r, C_A/C_F);
        N_R = 1 - pow(1-N_R, C_A/C_F);
    }

    integrand *= N_r*N_R;
    //cout << "4: " << N_r*N_R << endl;

    return integrand;

}

double g(double *k, size_t dim, void * params) {
    struct parameters *par = (struct parameters *)params;
    return Phi_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->P2, par->y, par->xpom);
}

void integrate(thread_par_struct par) {

    const int dim = 6;
    double res, err;

    double xl[dim] = {-r_limit, -r_limit, -r_limit, -r_limit, -b_min_limit, -b_min_limit};
    double xu[dim] = {r_limit, r_limit, r_limit, r_limit, b_min_limit, b_min_limit};

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
        //cout << "Integration iteration " << i << " result: " << res << ", err: " << err  << ", fit: " << gsl_monte_vegas_chisq(s) << ", duration: " << duration.count() << endl;
    }
    if (gsl_isnan(res)) {
        res = 0;
        cout << "nan found at xpom=" << params.xpom << endl;
    }
    result = res;
    error = err;
    fit = gsl_monte_vegas_chisq(s);
    cout << "L, P2=" << params.P2 << ", y=" << params.y << ", xpom=" << params.xpom << ", res: " << result << ", err: " << error << ", fit: " << gsl_monte_vegas_chisq(s) << ", duration: " << duration_double << " seconds" << endl;

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


    thread integration_threads[P2_resolution*y_resolution*xpom_resolution];
    double results[P2_resolution*y_resolution*xpom_resolution], errors[P2_resolution*y_resolution*xpom_resolution], fits[P2_resolution*y_resolution*xpom_resolution];

    static auto t1 = chrono::high_resolution_clock::now();

    for (int i=0; i<P2_resolution; i++) {
        for (int j=0; j<y_resolution; j++) {
            for (int k=0; k<xpom_resolution; k++) {
                int index = i*y_resolution*xpom_resolution + j*xpom_resolution + k;
                thread_par_struct par(P2_steps[i], y_steps[j], xpom_steps[k], results[index], errors[index], fits[index]);
                integration_threads[index] = thread(integrate, par);
            }
        }
    }

    for (int i=0; i<P2_resolution*y_resolution*xpom_resolution; i++) {
        integration_threads[i].join();
    }

    static auto t2 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(t2-t1);
    cout << "Calculation finished in " << duration.count()/60.0 << " minutes" << endl;


    ofstream output_file("output/Phi_table.txt");
    output_file << "P2;y;xpom;result;error;fit" << endl;

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

