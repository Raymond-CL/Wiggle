#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <stdlib.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

struct parameters {
  double rad;
};

double integrand(double *dx, size_t ndim, void *params) {
  (void)(ndim);  // unused
  // static cast params to p-pointer
  auto *p = static_cast<parameters *>(params);
  double rsq = 0.0;
  for (size_t i = 0; i < ndim; ++i) rsq += dx[i] * dx[i];
  return (rsq <= p->rad * p->rad) ? 1.0 : 0.0;
}

int main(int argc, char *argv[]) {
  // command-line arguments
  (void)(argc);
  (void)(argv);

  // define number of dimensions and radius
  size_t ndim = 10;
  double radius = 1.0;
  parameters p = {radius};

  // define integration limits
  std::vector<double> dx_lower(ndim, -radius);
  std::vector<double> dx_upper(ndim, +radius);

  // setup integrand
  gsl_monte_function gmf = {&integrand, ndim, &p};

  // setup Monte-Carlo integration environment
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
  gsl_monte_vegas_params vp;

  // define base call and iteration number
  size_t n = 10000000;
  size_t total = 10;
  size_t w = 5, f = total - w;

  // define output format and analytic result
  std::cout << std::setprecision(15) << std::fixed;
  constexpr double pi = 3.14159265358979323846;
  double exact = std::pow(pi * radius * radius, ndim / 2.0) /
                 std::tgamma(ndim / 2.0 + 1.0);
  double result, error;

  // warmup run: iteration=w, calls=n
  gsl_monte_vegas_params_get(s, &vp);
  vp.stage = 0;
  vp.iterations = w;
  gsl_monte_vegas_params_set(s, &vp);
  gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                            n, r, s, &result, &error);

  // final run: iteration=1, calls=f*n
  gsl_monte_vegas_params_get(s, &vp);
  vp.stage = 2;
  vp.iterations = 1;
  gsl_monte_vegas_params_set(s, &vp);
  gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                            f * n, r, s, &result, &error);

  // total calls: w*n + f*n = total*n
  std::cout << "vegas result = " << result << ", error = " << error
            << std::endl;
  std::cout << "exact result = " << exact
            << ", error = " << std::abs(result - exact) << std::endl;

  // free resources
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);

  return 0;
}
