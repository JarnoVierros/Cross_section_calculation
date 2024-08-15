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
const double width = radius;

struct parameters {double the_parameter;};

double calc_max_phi(double r, double b_min) {
    if (r > 2*b_min) {
        return 2*M_PI;
    } else {
        return M_PI - acos(r/(2*b_min));
    }
}

double trans_integral_0(double r, double theta) {
    return r*1/(r*r*r+1);
}

double trans_integral_0_wrap(double *k, size_t dim, void * params) {
    return trans_integral_0(k[0], k[1]);
}

double trans_integral_1(double r, double b_min, double phi, double z) {
    if (calc_max_phi(r, b_min) < phi) {
        return 0;
    }
    return r*b_min*4*2*M_PI*1/(r*r*r+1)*1/(gsl_pow_3(sqrt(b_min*b_min+2*(1-z)*r*b_min*cos(phi)+gsl_pow_2(1-z)*r*r))+1);
}

double trans_integral_1_wrap(double *k, size_t dim, void * params) {
    return trans_integral_1(k[0], k[1], k[2], k[3]);
}

double orig_integral_0(double r1, double r2) {
    return 1/(gsl_pow_3(sqrt(r1*r1+r2*r2))+1);
}

double orig_integral_0_wrap(double *k, size_t dim, void * params) {
    return orig_integral_0(k[0], k[1]);
}

double orig_integral_1(double r1, double r2, double b1, double b2, double z) {
    return 1/(gsl_pow_3(sqrt(r1*r1+r2*r2))+1)*1/(gsl_pow_3(sqrt(b1*b1+b2*b2))+1);
}

double orig_integral_1_wrap(double *k, size_t dim, void * params) {
    return orig_integral_1(k[0], k[1], k[2],  k[3],  k[4]);
}

void trans_integrate(double &return_res, double &return_err, double &return_fit, double (*func)(double *x_array, size_t, void *params), int calls, size_t dim) {

    double res, err;

    double xl[dim];
    double xu[dim];

    if (dim==2) {
        xl[0] = 0;
        xu[0] = radius;
        xl[1] = 0;
        xu[1] = 2*M_PI;
    } else if (dim==3) {
        xl[0] = 0;
        xu[0] = radius;
        xl[1] = 0;
        xu[1] = radius;
        xl[2] = 0;
        xu[2] = M_PI;
    } else if (dim==4) {
        xl[0] = 0;
        xu[0] = radius;
        xl[1] = 0;
        xu[1] = radius;
        xl[2] = 0;
        xu[2] = M_PI;
        xl[3] = 0;
        xu[3] = 1;
    }

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
    if (status != 0) {cout << "trans_integral: " << status << endl; throw (status);}

    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, calls, rng, s, &res, &err);
    if (status != 0) {cout << "trans_integral: " << status << endl; throw (status);}

    return_res = res;
    return_err = err;
    return_fit = gsl_monte_vegas_chisq(s);

    gsl_monte_vegas_free(s);
}

void orig_integrate(double &return_res, double &return_err, double &return_fit, double (*func)(double *x_array, size_t, void *params), int calls, size_t dim) {

    double res, err;

    double xl[dim];
    double xu[dim];

    if (dim==2) {
        xl[0] = -width;
        xu[0] = width;
        xl[1] = -width;
        xu[1] = width;
    } else if (dim==4) {
        xl[0] = -width;
        xu[0] = width;
        xl[1] = -width;
        xu[1] = width;
        xl[2] = -width;
        xu[2] = width;
        xl[3] = -width;
        xu[3] = width;
    } else if (dim==5) {
        xl[0] = -width;
        xu[0] = width;
        xl[1] = -width;
        xu[1] = width;
        xl[2] = -width;
        xu[2] = width;
        xl[3] = -width;
        xu[3] = width;
        xl[4] = 0;
        xu[4] = 1;
    }

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
    if (status != 0) {cout << "orig_integral: " << status << endl; throw (status);}

    status = gsl_monte_vegas_integrate(&T_G, xl, xu, dim, calls, rng, s, &res, &err);
    if (status != 0) {cout << "orig_integral: " << status << endl; throw (status);}

    return_res = res;
    return_err = err;
    return_fit = gsl_monte_vegas_chisq(s);

    gsl_monte_vegas_free(s);
}

int main() {

    double res, err, fit;

    /*
    trans_integrate(res, err, fit, trans_integral_0_wrap, 1e4, 2);
    cout << "trans integral 0: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    //trans integral 1: res=2819.87, err=78.2318, fit=6.22363

    orig_integrate(res, err, fit, orig_integral_0_wrap, 1e4, 2);
    cout << "orig integral 0: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    */
    
    trans_integrate(res, err, fit, trans_integral_1_wrap, 1e7, 3);
    cout << "trans integral 1: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    //trans integral 1: res=2819.87, err=78.2318, fit=6.22363

    orig_integrate(res, err, fit, orig_integral_1_wrap, 1e7, 4);
    cout << "orig integral 1: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    
}