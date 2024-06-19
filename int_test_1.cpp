#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <string>
#include <iostream>
using namespace std;




double g (double *k, size_t dim, void *params) {
  double A = 1.0 / (M_PI * M_PI * M_PI);
  return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}


int main() {

    const int dim = 5;
    double res, err;

    double xl[3] = { 0, 0, 0 };
    double xu[3] = { M_PI, M_PI, M_PI };

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &g, dim, 0 };

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);


    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);

    gsl_monte_vegas_integrate (&G, xl, xu, dim, 10000, r, s, &res, &err);

    cout << "res: " << res << endl;
    cout << "err: " << err << endl;
    cout << "chisq: " << gsl_monte_vegas_chisq(s) << endl;


    gsl_monte_vegas_free (s);
    gsl_rng_free (r);

    return 0;
}
