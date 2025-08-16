#ifndef VEGAS_H
#define VEGAS_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

#include <vector>

class vegas {
 private:
  size_t dim_;
  std::vector<double> xlower_, xupper_;
  gsl_monte_function integrand_;
  gsl_rng *generator_;
  gsl_monte_vegas_state *state_;

 public:
  vegas(size_t dim);
  ~vegas();
  void set_limits(const std::vector<double> &lower,
                  const std::vector<double> &upper);
  void set_limits(const double *lower, const double *upper);
  void set_integrand(double (*f)(double *, size_t, void *), void *params);
  void integrate(size_t calls, int iterations, double &result, double &error);
};

#endif