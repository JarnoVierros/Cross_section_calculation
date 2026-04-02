

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


const int warmup_calls = 100000;
const int integration_calls = 10000000;
const int integration_iterations = 1;

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


double Phi_integrand(double r1, double r2, double R1, double R2, double b1, double b2, double P2, double y, double xpom) {

    double radius = r1*r1+r2*r2+R1*R1+R2*R2+b1*b1+b2*b2;
    //cout << radius << ", " << exp(-sqrt(radius)/10) << endl;
    return exp(-sqrt(radius));

}

double g(double *k, size_t dim, void * params) {
    struct parameters *par = (struct parameters *)params;
    return Phi_integrand(k[0], k[1], k[2], k[3], k[4], k[5], par->P2, par->y, par->xpom);
}

void integrate(thread_par_struct par) {

    const int dim = 6;
    double res, err;

    double limit = 100.0;

    double xl[dim] = {-limit, -limit, -limit, -limit, -limit, -limit};
    double xu[dim] = {limit, limit, limit, limit, limit, limit};

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

    static auto t1 = chrono::high_resolution_clock::now();

    double result, error, fit;

    thread_par_struct par(1, 1, 1, result, error, fit);
    integrate(par);

    cout << "result=" << result << ", error=" << error << ", fit=" << fit << endl;
}

