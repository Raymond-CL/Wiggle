/*
Program: Wiggle
Description: Wigner function description of photon-photon to dilepton production
Process: \gamma + \gamma -> l^+ + l^-
*/

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdlib.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>

#include "physconst.h"

/*
Kinematic naming convention:
* center of mass energy: CME
* incoming photon amplitude: k1, k2, conjugate: ka, kb
  - transverse momentum: k1t, kat, k2t, kbt
  - azimuth: phi_k1, phi_ka, phi_k2, phi_kb
  - virtuality: k1sq, kasq, k2sq, kbsq
  - Delta: k1ka=|k1-ka|, k2kb=|k2-kb|
* photon momentum imbalance: qt_imb = k1t + k2t = kat + kbt
  - azimuth: phi_qt
* Sudakov momentum imbalance: lt_imb
  - azimuth: phi_lt
* outgoing lepton: p1, p2
  - transverse momentum: p1t, p2t
  - azimuth: phi_p1, phi_p2
  - rapidity: y1, y2
* lepton momentum imbalance: pt_imb = p1t + p2t = qt_imb + lt_imb
  - azimuth: phi_pt
* lepton relative momentum: Pt_rel = 0.5 * (p1t - p2t)
  - azimuth: phi_Pt
* impact parameter: bp
* Mandelstam variable: mans, mant, manu
* momentum fraction: x1, x2
*/
class kinematics {
 public:
  double CME;
  double Ptmin, Pt_rel, Ptmax, phi_Pt;
  int Ptn;
  double ptmin, pt_imb, ptmax, phi_pt;
  int ptn;
  double bpmin, bp, bpmax;
  int bpn;
  double almin, alpha, almax;
  int aln;
  double ktmin, ktmax;
  double k1t, kat, k2t, kbt;
  double phi_k1, phi_ka, phi_k2, phi_kb;
  double k1sq, kasq, k2sq, kbsq;
  double k1ka, k2kb;
  double qt_imb, phi_qt;
  double lt_imb, phi_lt;
  double p1t, p2t, phi_p1, phi_p2;
  double ylmin, ylmax;
  double y1, y2;
  double m_lep;
  double Mllmin, Mll, Mllmax;
  double x1, x2;
  double mans, mant, manu;
};

double integrand(double *dx, size_t ndim, void *params) {
  (void)(ndim);  // unused
  // static cast params to kinematics pointer
  auto *k = static_cast<kinematics *>(params);

  k->k1t = dx[0];
  k->phi_k1 = dx[1];
  k->kat = dx[2];
  k->phi_ka = dx[3];
  k->y1 = dx[4];
  k->y2 = dx[5];
  k->Pt_rel = dx[6];
  k->phi_Pt = dx[7];
  k->phi_pt = dx[8];
  k->bp = dx[9];

  // dilepton kinematics definition
  double Ptx = k->Pt_rel * cos(k->phi_Pt);
  double Pty = k->Pt_rel * sin(k->phi_Pt);
  double ptx = k->pt_imb * cos(k->phi_pt);
  double pty = k->pt_imb * sin(k->phi_pt);

  // p1t = 0.5pt + PT, p2t = 0.5pt - Pt
  double p1x = 0.5 * ptx + Ptx;
  double p1y = 0.5 * pty + Pty;
  double p2x = 0.5 * ptx - Ptx;
  double p2y = 0.5 * pty - Pty;
  k->p1t = std::sqrt(p1x * p1x + p1y * p1y);
  k->p2t = std::sqrt(p2x * p2x + p2y * p2y);

  // set momentum fraction
  double ptmax = std::max(k->p1t, k->p2t);
  k->x1 = ptmax * (std::exp(+k->y1) + std::exp(+k->y2)) / k->CME;
  k->x2 = ptmax * (std::exp(-k->y1) + std::exp(-k->y2)) / k->CME;
  if (k->x1 <= 0.0 || k->x1 >= 1.0) return 0.0;
  if (k->x2 <= 0.0 || k->x2 >= 1.0) return 0.0;

  // set Mandelstam variables
  k->mans = k->x1 * k->x2 * k->CME * k->CME;
  k->mant = -k->x1 * k->CME * ptmax * std::exp(-k->y1);
  k->manu = -k->x1 * k->CME * ptmax * std::exp(-k->y2);

  // perform dilepton invariant mass cut
  k->Mll = std::sqrt(k->mans);
  if (k->Mll <= k->Mllmin || k->Mll >= k->Mllmax) return 0.0;

  // Sudakov kinematics
  double ltx = 0.0;
  double lty = 0.0;

  // momentum imbalances pt = qt + lt
  double qtx = ptx - ltx;
  double qty = pty - lty;

  // set photon kinematics
  double k1x = k->k1t * cos(k->phi_k1);
  double k1y = k->k1t * sin(k->phi_k1);
  double kax = k->kat * cos(k->phi_ka);
  double kay = k->kat * sin(k->phi_ka);
  double k2x = qtx - k1x;
  double k2y = qty - k1y;
  double kbx = qtx - kax;
  double kby = qty - kay;
  k->k2t = std::sqrt(k2x * k2x + k2y * k2y);
  k->kbt = std::sqrt(kbx * kbx + kby * kby);
  k->phi_k2 = std::atan2(k2y, k2x);
  k->phi_kb = std::atan2(kby, kbx);

  double mp2 = phys::M_pro * phys::M_pro;
  k->k1sq = k->x1 * k->x1 * mp2 + k->k1t * k->k1t;
  k->kasq = k->x1 * k->x1 * mp2 + k->kat * k->kat;
  k->k2sq = k->x2 * k->x2 * mp2 + k->k2t * k->k2t;
  k->kbsq = k->x2 * k->x2 * mp2 + k->kbt * k->kbt;

  k->k1ka = std::sqrt((k1x - kax) * (k1x - kax) + (k1y - kay) * (k1y - kay));

  double sigma0 = 2.0 * phys::alphae * phys::alphae / k->mans / k->mans *
                  k->k1t * k->kat * k->k2t * k->kbt *
                  ((k->mant / k->manu + k->manu / k->mant) *
                       std::cos(k->phi_k1 - k->phi_ka + k->phi_k2 - k->phi_kb) -
                   2.0 * std::cos(k->phi_k1 + k->phi_ka + k->phi_k2 +
                                  k->phi_kb - 4.0 * k->phi_Pt));

  auto ffnuc = [](double kin) {
    constexpr double R = 6.62 / phys::GeVfm;
    double kr = kin * R;
    return (std::sin(kr) - kr * std::cos(kr)) * 3.0 / std::pow(kr, 3) /
           (1.0 + phys::a0 * phys::a0 * kin * kin);
  };

  auto xf = [&ffnuc](double kamp, double kcon) {
    if (kamp <= 0.0 || kcon <= 0.0) return 0.0;
    constexpr int atomZ = 82;
    double fac = atomZ * atomZ * phys::alphae / phys::PI / phys::PI;
    double nff = ffnuc(kamp) / kamp / kamp * ffnuc(kcon) / kcon / kcon;
    return fac * nff;
  };

  double x1f1 = xf(std::sqrt(k->k1sq), std::sqrt(k->kasq));
  double x2f2 = xf(std::sqrt(k->k2sq), std::sqrt(k->kbsq));

  double bessel = k->bp / phys::twoPI * gsl_sf_bessel_J0(k->bp * k->k1ka);

  double costerm = -2.0 * std::cos(4.0 * k->phi_Pt - 4.0 * k->phi_pt);

  double res = k->Pt_rel * k->pt_imb * k->k1t * k->kat * bessel * x1f1 * x2f2 *
               sigma0 * costerm;

  return res;
}

// main program
int main(int argc, char *argv[]) {
  // check command line arguments
  (void)(argc);  // unused
  (void)(argv);  // unused
  // start program timer
  auto start = std::chrono::high_resolution_clock::now();

  // define kinematics input parameters
  kinematics k;
  k.CME = 5020.0;
  k.Ptmin = 4.0, k.Ptmax = 20.0, k.Ptn = 0;
  k.ptmin = 0.0, k.ptmax = 0.1, k.ptn = 50;
  k.bpmin = 9.5, k.bpmax = 10.5, k.bpn = 0;
  k.almin = 0.0, k.almax = 0.0, k.aln = 0;
  k.ktmin = 0.0, k.ktmax = 10.0;
  k.ylmin = -1.0, k.ylmax = +1.0;
  // k.atomA = 207, k.atomZ = 82, k.RA = 6.62 / phys::GeVfm;
  k.m_lep = phys::M_mu;
  k.Mllmin = 4.0, k.Mllmax = 45.0;

  // define integration limits
  size_t ndim = 10;
  std::vector<double> dx_lower(ndim);
  std::vector<double> dx_upper(ndim);
  // dk1t, dphi_k1, dkat, dphi_ka, dy1, dy2
  dx_lower[0] = k.ktmin, dx_upper[0] = k.ktmax;
  dx_lower[1] = 0.0, dx_upper[1] = phys::twoPI;
  dx_lower[2] = k.ktmin, dx_upper[2] = k.ktmax;
  dx_lower[3] = 0.0, dx_upper[3] = phys::twoPI;
  dx_lower[4] = k.ylmin, dx_upper[4] = k.ylmax;
  dx_lower[5] = k.ylmin, dx_upper[5] = k.ylmax;
  // dPt, dphi_Pt, dqt, dphi_qt, dbp
  dx_lower[6] = k.Ptmin, dx_upper[6] = k.Ptmax;
  dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
  dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
  dx_lower[9] = k.bpmin / phys::GeVfm, dx_upper[9] = k.bpmax / phys::GeVfm;

  // setup Monte-Carlo integration environment
  gsl_rng_env_setup();
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
  gsl_monte_vegas_params vp;
  gsl_monte_function gmf = {&integrand, ndim, &k};

  // define call and iteration number
  size_t ncall1 = 10000, itm1 = 10;
  size_t ncall2 = 100000, itm2 = 1;

  // perform integration loop
  int nbin = k.ptn;
  double hmin = k.ptmin;
  double hmax = k.ptmax;
  double bin = (hmax - hmin) / nbin;
  for (int i = 0; i < nbin; ++i) {
    double result, error;
    double binL = hmin + i * bin;
    double binR = hmin + (i + 1) * bin;
    double binM = (binL + binR) / 2.0;
    k.pt_imb = binM;

    // warmup run: iteration=w, calls=n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 0;
    vp.iterations = itm1;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              ncall1, r, s, &result, &error);

    // final run: iteration=1, calls=f*n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 2;
    vp.iterations = itm2;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              ncall2, r, s, &result, &error);

    // print result
    std::cout << binM << '\t' << result << '\t' << error << std::endl;
  }

  // display elapsed time
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

  // free resources
  gsl_monte_vegas_free(s);
  gsl_rng_free(r);

  return 0;
}
