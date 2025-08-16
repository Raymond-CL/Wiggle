#include "vegas.h"

#include <stdlib.h>

#include <iostream>

vegas::vegas(size_t dim) : dim_(dim), xlower_(dim), xupper_(dim) {}

vegas::~vegas() {}

void vegas::set_limits(const std::vector<double> &lower,
                       const std::vector<double> &upper) {
  if (lower.size() != dim_ || upper.size() != dim_) {
    throw std::invalid_argument("Invalid limits");
  }
  xlower_ = lower;
  xupper_ = upper;
}

void vegas::set_limits(const double *lower, const double *upper) {
  for (size_t i = 0; i < dim_; ++i) {
    xlower_[i] = lower[i];
    xupper_[i] = upper[i];
  }
}

void vegas::set_integrand(double (*f)(double *, size_t, void *), void *params) {
  integrand_ = {f, dim_, params};
}

void vegas::integrate(size_t calls, int iterations, double &result,
                      double &error) {
  gsl_rng_env_setup();
  generator_ = gsl_rng_alloc(gsl_rng_default);
  state_ = gsl_monte_vegas_alloc(dim_);
  gsl_monte_vegas_init(state_);
  do {
    gsl_monte_vegas_integrate(&integrand_, xlower_.data(), xupper_.data(), dim_,
                              calls, generator_, state_, &result, &error);
  } while (fabs(gsl_monte_vegas_chisq(state_) - 1.0) > 0.5 && iterations-- > 0);

  gsl_monte_vegas_integrate(&integrand_, xlower_.data(), xupper_.data(), dim_,
                            calls, generator_, state_, &result, &error);
  gsl_monte_vegas_free(state_);
  gsl_rng_free(generator_);
  // integration workflow is initialize state, integrate, and wipe state
  // because there is not state reset, so we need to allocate a new state for each integration
  // otherwise the next integration will remember the previous state, and error propagates
}