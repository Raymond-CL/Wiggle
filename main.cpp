#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdlib.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "physconst.h"

struct parameters {
  double CME;
  double qtmin, qtmax;
  int qtn;
  double ylmin, ylmax;
  double pltmin, pltmax;
  double Mllmin, Mllmax;
  double ktmin, ktmax;
  int atomA, atomZ;
  double m_lep;
  double RA;
  double qt_imb;
};

double integrand(double *dx, size_t ndim, void *params) {
  (void)(ndim);  // unused
  // static cast params to p-pointer
  auto *p = static_cast<parameters *>(params);

  double k1t = dx[0];
  double phik1 = dx[1];
  double kat = dx[2];
  double phika = dx[3];
  double y1 = dx[4];
  double y2 = dx[5];
  double pl_avg = dx[6];
  double phi_pl = dx[7];
  double phi_qt = dx[8];

  // k1t = 2.0;
  // phik1 = 1.0;
  // kat = 2.0;
  // phika = 1.0;
  // y1 = 0.0;
  // y2 = 0.0;
  // pl_avg = 10.0;
  // phi_pl = 1.0;
  // phi_qt = 2.0;
  // p->qt_imb = 0.05;

  double bp = 10.0 / phys::GeVfm;

  // kinematics definition
  double plx = pl_avg * cos(phi_pl);
  double ply = pl_avg * sin(phi_pl);
  double qtx = p->qt_imb * cos(phi_qt);
  double qty = p->qt_imb * sin(phi_qt);

  double p1x = 0.5 * qtx + plx;
  double p1y = 0.5 * qty + ply;
  double p2x = 0.5 * qtx - plx;
  double p2y = 0.5 * qty - ply;
  double p1t = std::sqrt(p1x * p1x + p1y * p1y);
  double p2t = std::sqrt(p2x * p2x + p2y * p2y);

  double ptmax = std::max(p1t, p2t);
  double x1 = ptmax * (std::exp(+y1) + std::exp(+y2)) / p->CME;
  double x2 = ptmax * (std::exp(-y1) + std::exp(-y2)) / p->CME;

  if (x1 <= 0.0 || x1 >= 1.0) return 0.0;
  if (x2 <= 0.0 || x2 >= 1.0) return 0.0;

  double mans = +x1 * x2 * p->CME * p->CME;
  double mant = -x1 * p->CME * ptmax * std::exp(-y1);
  double manu = -x1 * p->CME * ptmax * std::exp(-y2);

  double Mll = std::sqrt(mans);
  if (Mll <= p->Mllmin || Mll >= p->Mllmax) return 0.0;

  double k1x = k1t * cos(phik1);
  double k1y = k1t * sin(phik1);
  double kax = kat * cos(phika);
  double kay = kat * sin(phika);

  double k2x = qtx - k1x;
  double k2y = qty - k1y;
  double kbx = qtx - kax;
  double kby = qty - kay;
  double k2t = std::sqrt(k2x * k2x + k2y * k2y);
  double kbt = std::sqrt(kbx * kbx + kby * kby);
  double phik2 = std::atan2(k2y, k2x);
  double phikb = std::atan2(kby, kbx);

  double mp2 = phys::M_pro * phys::M_pro;
  double k1sq = x1 * x1 * mp2 + k1t * k1t;
  double kasq = x1 * x1 * mp2 + kat * kat;
  double k2sq = x2 * x2 * mp2 + k2t * k2t;
  double kbsq = x2 * x2 * mp2 + kbt * kbt;

  double k1ka =
      std::sqrt((k1x - kax) * (k1x - kax) + (k1y - kay) * (k1y - kay));

  double sigma0 =
      2.0 * phys::alphae * phys::alphae / mans / mans * k1t * kat * k2t * kbt *
      ((mant / manu + manu / mant) * std::cos(phik1 - phika + phik2 - phikb) -
       2.0 * std::cos(phik1 + phika + phik2 + phikb - 4.0 * phi_pl));

  auto ffnuc = [](double k, double r, double a0) {
    double kr = k * r;
    return (std::sin(kr) - kr * std::cos(kr)) * 3.0 / std::pow(kr, 3) /
           (1.0 + a0 * a0 * k * k);
  };

  auto xf = [&ffnuc](double k, double kp, int Z, double r, double a0) {
    if (k <= 0.0 || kp <= 0.0) return 0.0;
    double fac = Z * Z * phys::alphae / phys::PI / phys::PI;
    double nff = ffnuc(k, r, a0) / k / k * ffnuc(kp, r, a0) / kp / kp;
    return fac * nff;
  };

  double x1f1 = xf(std::sqrt(k1sq), std::sqrt(kasq), p->atomZ, p->RA, phys::a0);
  double x2f2 = xf(std::sqrt(k2sq), std::sqrt(kbsq), p->atomZ, p->RA, phys::a0);

  double bessel = bp / phys::twoPI * gsl_sf_bessel_J0(bp * k1ka);

  double costerm = -2.0 * std::cos(4.0 * phi_pl - 4.0 * phi_qt);

  double res =
      pl_avg * p->qt_imb * k1t * kat * bessel * x1f1 * x2f2 * sigma0 * costerm;

  return res;
}

int main(int argc, char *argv[]) {
  // check command line arguments
  (void)(argc);  // unused
  (void)(argv);  // unused

  // define input parameters
  double CME = 5020.0;
  double qtmin = 0.0, qtmax = 0.1;
  int qtn = 50;
  double ylmin = -1.0, ylmax = 1.0;
  double pltmin = 4.0, pltmax = 20.0;
  double Mllmin = 4.0, Mllmax = 45.0;
  double ktmin = 0.0, ktmax = 10.0;
  int atomA = 207, atomZ = 82;
  double m_lep = phys::M_mu;
  double RA = 6.62 / phys::GeVfm;

  // define integration limits
  size_t ndim = 9;
  std::vector<double> dx_lower(ndim);
  std::vector<double> dx_upper(ndim);
  dx_lower[0] = ktmin;
  dx_upper[0] = ktmax;
  dx_lower[1] = 0.0;
  dx_upper[1] = phys::twoPI;
  dx_lower[2] = ktmin;
  dx_upper[2] = ktmax;
  dx_lower[3] = 0.0;
  dx_upper[3] = phys::twoPI;
  dx_lower[4] = ylmin;
  dx_upper[4] = ylmax;
  dx_lower[5] = ylmin;
  dx_upper[5] = ylmax;
  dx_lower[6] = pltmin;
  dx_upper[6] = pltmax;
  dx_lower[7] = 0.0;
  dx_upper[7] = phys::twoPI;
  dx_lower[8] = 0.0;
  dx_upper[8] = phys::twoPI;

  double pro_test = 1.0;
  for (size_t i = 0; i < ndim; ++i) {
    pro_test *= dx_upper[i] - dx_lower[i];
  }

  // setup Monte-Carlo integration environment
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
  gsl_monte_vegas_params vp;

  // define base call and iteration number
  // size_t n = 1000000;
  // size_t total = 10;
  // size_t w = 5, f = total - w;
  // size_t w = 10, f = 10;
  size_t ncall1 = 10000000, itm1 = 10;
  size_t ncall2 = 100000000, itm2 = 1;

  // perform integration loop
  int nbin = qtn;
  double hmin = qtmin;
  double hmax = qtmax;
  double bin = (hmax - hmin) / nbin;
  double result, error;
  for (int i = 0; i < nbin; ++i) {
    double binL = hmin + i * bin;
    double binR = hmin + (i + 1) * bin;
    double binM = (binL + binR) / 2.0;
    double qt_imb = binM;

    parameters p = {CME,    qtmin,  qtmax,  qtn,    ylmin, ylmax,
                    pltmin, pltmax, Mllmin, Mllmax, ktmin, ktmax,
                    atomA,  atomZ,  m_lep,  RA,     qt_imb};

    gsl_monte_function gmf = {&integrand, ndim, &p};

    // warmup run: iteration=w, calls=n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 0;
    vp.iterations = itm1;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim, ncall1,
                              r, s, &result, &error);

    // final run: iteration=1, calls=f*n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 2;
    vp.iterations = itm2;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim, ncall2,
                              r, s, &result, &error);

    std::cout << qt_imb << '\t' << result << '\t' << error << std::endl;
  }

  // free resources
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);

  return 0;
}
