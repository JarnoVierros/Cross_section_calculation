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

const double alpha_em = 1.0/137;
const int N_c = 3;
const double e_f = 2.0/3;
const double m_f = 1.27; //GeV 1.27

const double radius = 100;
const double width = radius;

const double Q2 = 1;
const double M_X = 3;

const double max_theta_root_excess = 1e-6;
const int debug_precision = 10;

struct parameters {double the_parameter;};

double epsilon2(double z, double Q2) {
  return m_f*m_f + z*(1-z)*Q2;
}

double epsilon(double z, double Q2) {
  return sqrt(epsilon2(z, Q2));
}

double calc_h(double r, double b_min, double phi, double z) {
  return b_min*b_min + (1-z)*2*b_min*r*cos(phi) + gsl_pow_2(1-z)*r*r;
}

double calc_b1(double r, double b_min, double phi, double z) {
  return b_min + r*(1-z)*cos(phi);
}

double calc_b2(double r, double b_min, double phi, double z) {
  return r*(1-z)*sin(phi);
}

double calc_j(double b2, double r_bar, double phi_bar, double z) {
  return (1-z)*2*b2*r_bar*sin(phi_bar);
}

double calc_A(double j, double h, double b1, double b2) {
  return sqrt(j*j + 4*h*(b1*b1-gsl_pow_2(j/(2*b2))));
}

bool theta_root_invalid(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  double excess = abs(tan(theta_bar)*(b_min+r*(1-z)*cos(phi)-r_bar*(1-z)*cos(phi_bar+theta_bar)) + r_bar*(1-z)*sin(phi_bar+theta_bar) - r*(1-z)*sin(phi));
  if (excess > max_theta_root_excess) {
    return 1;
  } else {
    return 0;
  }
}

int calc_theta_bar(double return_values[4], double r, double b_min, double phi, double r_bar, double phi_bar, double z) {
  double b1 = calc_b1(r, b_min, phi, z);
  double b2 = calc_b2(r, b_min, phi, z);
  double h = calc_h(r, b_min, phi, z);
  if (r_bar*r_bar > (4*h*b1*b1)/(gsl_pow_2((1-z)*2*sin(phi_bar))*(h-b2*b2))) {
    //r_bar is too large
    return 1;
  }
  double j = calc_j(b2, r_bar, phi_bar, z);
  double A = calc_A(j, h, b1, b2);
  if (gsl_isnan(A)) {
    cout << "A is nan!" << endl;
    cout << "j=" << j << endl;
    cout << "h=" << h << endl;
    cout << "b1=" << b1 << endl;
    cout << "b2=" << b2 << endl;
    cout << "A'=" << j*j + 4*h*(b1*b1-gsl_pow_2(j/(2*b2))) << endl;
    cout << endl;

  }
  double acos_arg_plus_type = (A+j)/(2*h);
  if (acos_arg_plus_type > 1) {
    cout << "acos arg of type plus changed: " << setprecision(debug_precision) << acos_arg_plus_type << "->1" << endl;
    acos_arg_plus_type = 1;
  } else if (acos_arg_plus_type < -1) {
    cout << "acos arg of type plus changed: " << setprecision(debug_precision) << acos_arg_plus_type << "->-1" << endl;
    acos_arg_plus_type = -1;
  }
  double acos_arg_minus_type = (-A+j)/(2*h);
  if (acos_arg_minus_type > 1) {
    cout << "acos arg of type minus changed: " << setprecision(debug_precision) << acos_arg_minus_type << "->1" << endl;
    acos_arg_minus_type = 1;
  } else if (acos_arg_minus_type < -1) {
    cout << "acos arg of type minus changed: " << setprecision(debug_precision) << acos_arg_minus_type << "->-1" << endl;
    acos_arg_minus_type = -1;
  }
  return_values[0] = acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[0], z)) {
    return_values[0] = 10; // 10 means that the theta_bar value is invalid and should be skipped
  }
  return_values[1] = acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[1], z)) {
    return_values[1] = 10;
  }
  return_values[2] = -acos(acos_arg_plus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[2], z)) {
    return_values[2] = 10;
  }
  return_values[3] = -acos(acos_arg_minus_type);
  if (theta_root_invalid(r, b_min, phi, r_bar, phi_bar, return_values[3], z)) {
    return_values[3] = 10;
  }
  return 0;
}

double calc_b_bar(double r, double b_min, double phi, double r_bar, double phi_bar, double theta_bar, double z) {
  return 1/cos(theta_bar)*(b_min + (1-z)*r*cos(phi) - (1-z)*r_bar*cos(theta_bar+phi_bar));
}

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

double trans_integral_2(double r, double b_min, double phi, double r_bar, double phi_bar) {
    double z = 1.0/2;
    if (calc_max_phi(r, b_min) < phi) {
        return 0;
    }
    double total_integrand = 0;
    double theta_bar[4];
    int r_bar_too_large = calc_theta_bar(theta_bar, r, b_min, phi, r_bar, phi_bar, z);
    if (r_bar_too_large == 1) {
        return 0;
    }
    for (int i=0; i<4; i++) {
        if (theta_bar[i] == 10) {
            continue;
        }
        double b_min_bar = calc_b_bar(r, b_min, phi, r_bar, phi_bar, theta_bar[i], z);

        double sub_integrand = 1/abs(cos(theta_bar[i]))*1/abs(b_min_bar*cos(theta_bar[i])+(1-z)*r_bar*cos(theta_bar[i]+phi_bar))
        *r*b_min*r_bar*b_min_bar*4*2*M_PI*exp(-1*(r + sqrt(b_min*b_min+2*(1-z)*r*b_min*cos(phi)+gsl_pow_2(1-z)*r*r) + r_bar));
        if (gsl_isnan(sub_integrand)) {
            cout << "T sub_integrand " << i << " is nan" << endl;
            sub_integrand = 0;
        }

        total_integrand += sub_integrand;
    }
    return total_integrand;
}

double trans_integral_2_wrap(double *k, size_t dim, void * params) {
    return trans_integral_2(k[0], k[1], k[2], k[3], k[4]);
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

double orig_integral_2(double r1, double r2, double b1, double b2, double r1_bar, double r2_bar) {
    return exp(-1*(sqrt(r1*r1+r2*r2) + sqrt(b1*b1+b2*b2) + sqrt(r1_bar*r1_bar+r2_bar*r2_bar)));
}

double orig_integral_2_wrap(double *k, size_t dim, void * params) {
    return orig_integral_2(k[0], k[1], k[2],  k[3],  k[4], k[5]);
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
    } else if (dim==5) {
        xl[0] = 0;
        xu[0] = radius;
        xl[1] = 0;
        xu[1] = radius;
        xl[2] = 0;
        xu[2] = M_PI;
        xl[3] = 0;
        xu[3] = radius;
        xl[4] = 0;
        xu[4] = M_PI;
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
    } else if (dim==6) {
        xl[0] = -width;
        xu[0] = width;
        xl[1] = -width;
        xu[1] = width;
        xl[2] = -width;
        xu[2] = width;
        xl[3] = -width;
        xu[3] = width;
        xl[4] = -width;
        xu[4] = width;
        xl[5] = -width;
        xu[5] = width;
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
    /*
    trans_integrate(res, err, fit, trans_integral_1_wrap, 1e9, 4);
    cout << "trans integral 1: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    //trans integral 1: res=2819.87, err=78.2318, fit=6.22363

    orig_integrate(res, err, fit, orig_integral_1_wrap, 1e9, 5);
    cout << "orig integral 1: res=" << res << ", err=" << err << ", fit=" << fit << endl;
    */
    trans_integrate(res, err, fit, trans_integral_2_wrap, 1e7, 5);
    cout << "trans integral 2: res=" << res << ", err=" << err << ", fit=" << fit << endl;

    orig_integrate(res, err, fit, orig_integral_2_wrap, 1e6, 6);
    cout << "orig integral 2: res=" << res << ", err=" << err << ", fit=" << fit << endl;
}