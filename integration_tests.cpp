#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

#include <string>
#include <iostream>
#include <thread>
#include <chrono>
#include <fstream>
#include <ostream>
#include <sstream>
#include <iomanip>
using namespace std;

const double radius = 100;
const double width = 200;

struct parameters {double the_parameter;};

double calc_max_phi(double r, double b_min) {
    if (r > 2*b_min) {
        return 2*M_PI;
    } else {
        return M_PI - acos(r/(2*b_min));
    }
}

double trans_integral_1(double r, double b_min, double phi) {
    double z = 1.0/2;
    if (calc_max_phi(r, b_min) < phi) {
        return 0;
    }
    return 4*1/r*1/sqrt(b_min*b_min+2*(1-z)*r*b_min*cos(phi)+gsl_pow_2(1-z)*r*r);
}

double trans_integral_1_wrap(double *k, size_t dim, void * params) {
    return trans_integral_1(k[0], k[1], k[2]);
}

double orig_integral_1(double z, double r, double b_min, double phi, double theta) {
    if (calc_max_phi(r, b_min) < phi) {
        return 0;
    }
    return 2*M_PI*4*1/r*1/sqrt(b_min*b_min+2*(1-z)*r*b_min*cos(phi)+gsl_pow_2(1-z)*r*r);
}

double orig_integral_1_wrap(double *k, size_t dim, void * params) {
    return trans_integral_1(k[0], k[1], k[2]);
}

void trans_integrate(double &return_res, double &return_err, double &return_fit, double (*func)(double *x_array, size_t, void *params), int calls) {

    const int dim = 3;
    double res, err;

    double xl[dim] = {0, 0, 0};
    double xu[dim] = {radius, radius, M_PI};
    struct parameters params = {1};

    const gsl_rng_type *T;
    gsl_rng *rng;

    gsl_monte_function T_G = {func, dim, &params};

    gsl_rng_env_setup ();
    int status = 0;

    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    gsl_rng_set(rng, 1);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);
    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, 10000, rng, s, &res, &err);
    if (status != 0) {cout << "trans_integral_1_wrap: " << status << endl; throw (status);}

    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, calls, rng, s, &res, &err);
    if (status != 0) {cout << "trans_integral_1_wrap: " << status << endl; throw (status);}

    return_res = res;
    return_err = err;
    return_fit = gsl_monte_vegas_chisq(s);

    gsl_monte_vegas_free(s);
}

int main() {

    double res, err, fit;

    trans_integrate(res, err, fit, trans_integral_1_wrap, 100000000);

    cout << "trans integral 1: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    
}