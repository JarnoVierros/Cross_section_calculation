#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>

#include <string>
#include <iostream>
using namespace std;


const double alpha_em = 1;
const int N_c = 3;
const double e_f = 2.0/3;
const double m = 1.27; //GeV

const double sigma_0 = 29; //mb
const double Q_0 = sqrt(0.15);
const double x_0 = 0.000041;
const double lambda_star = 0.288;

double Q = 1;
double x = 1;

double normalization = 4*alpha_em*N_c*e_f*e_f/(2*M_PI*2*M_PI);

double epsilon(double z) {
  return sqrt(m*m + z*(1-z)*Q*Q);
}

double dipole_amplitude(double r) {
  return sigma_0*(1 - exp(-1*gsl_pow_2((Q_0*r)/(2*pow(x/x_0, lambda_star/2)))));
}

double L_integrand(double r_x, double r_y, double z) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return 4*Q*Q*z*z*(1-z)*(1-z)*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z)*r))*dipole_amplitude(r);
}

double T_integrand(double r_x, double r_y, double z) {
  double r = sqrt(r_x*r_x + r_y*r_y);
  return (m*m*gsl_pow_2(gsl_sf_bessel_K0(epsilon(z)*r)) + gsl_pow_2(epsilon(z))*(z*z + gsl_pow_2(1-z))*gsl_pow_2(gsl_sf_bessel_K1(epsilon(z)*r)))*dipole_amplitude(r);
}

double g(double *k, size_t dim, void *params) {
  return normalization*L_integrand(k[0], k[1], k[2]);
}

double f(double *k, size_t dim, void *params) {
  return normalization*L_integrand(k[0], k[1], k[2]);
}

int main() {

  const int dim = 3;
  double res, err;

  double xl[3] = {-100, -100, 0};
  double xu[3] = {100, 100, 1};

  cout << "normalization: " << normalization << endl;
  cout << "L_integrand: " << L_integrand(-100, -100, 0.5) << endl;

  const gsl_rng_type *T;
  gsl_rng *rng;

  gsl_monte_function G = {&g, dim, 0};

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  rng = gsl_rng_alloc(T);


  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dim);

  gsl_monte_vegas_integrate (&G, xl, xu, dim, 10000, rng, s, &res, &err);

  cout << "res: " << res << endl;
  cout << "err: " << err << endl;
  cout << "chisq: " << gsl_monte_vegas_chisq(s) << endl;
  cout << endl;

  gsl_monte_vegas_integrate (&G, xl, xu, dim, 100000, rng, s, &res, &err);

  cout << "res: " << res << endl;
  cout << "err: " << err << endl;
  cout << "chisq: " << gsl_monte_vegas_chisq(s) << endl;
  cout << endl;

  gsl_monte_vegas_integrate (&G, xl, xu, dim, 100000, rng, s, &res, &err);

  cout << "res: " << res << endl;
  cout << "err: " << err << endl;
  cout << "chisq: " << gsl_monte_vegas_chisq(s) << endl;
  cout << endl;

  gsl_monte_vegas_integrate (&G, xl, xu, dim, 100000, rng, s, &res, &err);

  cout << "res: " << res << endl;
  cout << "err: " << err << endl;
  cout << "chisq: " << gsl_monte_vegas_chisq(s) << endl;


  gsl_monte_vegas_free(s);
  gsl_rng_free(rng);

  return 0;
}