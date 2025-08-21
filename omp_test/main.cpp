#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_sf_bessel.h>
#include <stdlib.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "physconst.h"
#include "vec2d.h"

// Need to cancel this kinematic class and replace it with a simpler struct
// move all variables locally to the integrand
class kinematics {
 public:
  double CME;
  // dilepton relative and total momentum
  double Ptmin, Ptmax;
  double ptmin, ptmax;
  size_t Ptn, ptn;
  vec2d Pt, pt;
  double PT_rel, phi_PT;
  double pt_imb, phi_pt;
  // lepton momentum
  vec2d p1t, p2t;
  // lepton rapidity
  double ylmin, ylmax;
  double y1, y2;
  // lepton mass
  double m_lep;
  // Sudakov momentum
  double ltmin, ltmax;
  vec2d lt;
  // photon momentum imbalance
  vec2d qt;
  double qt_imb, phi_qt;
  // photon momentum
  double ktmin, ktmax;
  vec2d k1t, kat, k2t, kbt;
  double k1sq, kasq, k2sq, kbsq;
  double k1ka, k2kb;
  // impact parameter
  double bpmin, bp, bpmax;
  size_t bpn;
  // azimuthal angle alpha
  double almin, al, almax;
  size_t aln;
  // dilepton mass (cut)
  double Mllmin, Mll, Mllmax;
  // momentum fraction and Mandelstam variables
  double pthard;
  double x1, x2;
  double mans, mant, manu;
  // flags
  int diff;
  int lep;
  // temporary variables
  double sigma0;

  void set_lepton_mom() {
    Pt.setPolar(PT_rel, phi_PT);
    pt.setPolar(pt_imb, phi_pt);
    p1t = 0.5 * pt + Pt;
    p2t = 0.5 * pt - Pt;
  }

  void set_mom_frac() {
    pthard = std::max(p1t.mag(), p2t.mag());
    x1 = pthard * (std::exp(+y1) + std::exp(+y2)) / CME;
    x2 = pthard * (std::exp(-y1) + std::exp(-y2)) / CME;
  }

  void set_mandelstam() {
    mans = x1 * x2 * CME * CME;
    mant = -x1 * CME * pthard * std::exp(-y1);
    manu = -x1 * CME * pthard * std::exp(-y2);
    Mll = std::sqrt(mans);
  }

  void set_photon_mom() {
    qt = pt - lt;
    k2t = qt - k1t;
    kbt = qt - kat;
    k1sq = x1 * x1 * phys::Mp2 + k1t.mag2();
    kasq = x1 * x1 * phys::Mp2 + kat.mag2();
    k2sq = x2 * x2 * phys::Mp2 + k2t.mag2();
    kbsq = x2 * x2 * phys::Mp2 + kbt.mag2();
    k1ka = (k1t - kat).mag();
    k2kb = (k2t - kbt).mag();
  }

  void set_sigma0() {
    double iso = std::cos(k1t.phi() - kat.phi() + k2t.phi() - kbt.phi()) *
                 (mant / manu + manu / mant);
    double aniso = 2.0 * std::cos(k1t.phi() + kat.phi() + k2t.phi() +
                                  kbt.phi() - 4.0 * phi_PT);
    sigma0 = k1t.mag() * kat.mag() * k2t.mag() * kbt.mag() * (iso - aniso) *
             2.0 * phys::alphae * phys::alphae / mans / mans;
  }
};

double integrand(double *dx, size_t ndim, void *params) {
  (void)(ndim);  // unused
  // static cast params to kinematics pointer
  auto *k = static_cast<kinematics *>(params);

  k->k1t.setPolar(dx[0], dx[1]);
  k->kat.setPolar(dx[2], dx[3]);
  k->y1 = dx[4];
  k->y2 = dx[5];

  switch (k->diff) {
    case 0:  // total cross-section
      k->PT_rel = dx[6];
      k->phi_PT = dx[7];
      k->pt_imb = dx[8];
      k->phi_pt = dx[9];
      k->bp = dx[10];
      break;
    case 1:  // dsigma/dPT_rel
      k->phi_PT = dx[6];
      k->pt_imb = dx[7];
      k->phi_pt = dx[8];
      k->bp = dx[9];
      break;
    case 2:  // dsigma/dpt_imb
      k->PT_rel = dx[6];
      k->phi_PT = dx[7];
      k->phi_pt = dx[8];
      k->bp = dx[9];
      break;
    case 3:  // dsigma/dbp
      k->PT_rel = dx[6];
      k->phi_PT = dx[7];
      k->pt_imb = dx[8];
      k->phi_pt = dx[9];
      break;
    case 4:  // dsigma/dal
      k->PT_rel = dx[6];
      k->pt_imb = dx[7];
      k->phi_pt = dx[8];
      k->bp = dx[9];
      break;
    default:
      std::cerr << "Invalid diff value: " << k->diff << std::endl;
      return 0.0;
  }

  bool sud = false;
  if (sud) {
    k->lt.setPolar(dx[10], dx[11]);
  } else {
    k->lt.setPolar(0.0, 0.0);
  }

  k->set_lepton_mom();

  k->set_mom_frac();
  if (k->x1 <= 0.0 || k->x1 >= 1.0) return 0.0;
  if (k->x2 <= 0.0 || k->x2 >= 1.0) return 0.0;

  k->set_mandelstam();
  if (k->Mll <= k->Mllmin || k->Mll >= k->Mllmax) return 0.0;

  k->set_photon_mom();

  k->set_sigma0();

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

  // double costerm = -2.0 * std::cos(4.0 * k->phi_PT - 4.0 * k->phi_pt);
  double costerm = 1.0;

  double res = k->PT_rel * k->pt_imb * k->k1t.mag() * k->kat.mag() * bessel *
               x1f1 * x2f2 * k->sigma0 * costerm;

  return res;
}

// main program
int main(int argc, char *argv[]) {
  // check command line arguments
  (void)(argc);  // unused
  (void)(argv);  // unused

  // start program timer
  auto start = std::chrono::high_resolution_clock::now();

  // setup Monte-Carlo integration environment
  gsl_rng_env_setup();

  // integration calls and iterations
  size_t ncall1, itm1;
  size_t ncall2, itm2;

  // retrieve kinematics parameters from input file
  kinematics k;
  std::ifstream fin("input.dat", std::ios::in);
  if (!fin) return 1;
  fin >> k.CME;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.diff;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.Ptmin >> k.Ptmax >> k.Ptn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.ptmin >> k.ptmax >> k.ptn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.bpmin >> k.bpmax >> k.bpn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.almin >> k.almax >> k.aln;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.ktmin >> k.ktmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.ylmin >> k.ylmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.Mllmin >> k.Mllmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> k.lep;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ncall1 >> itm1;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ncall2 >> itm2;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin.close();

  /*
  We need to set the differential variable for different observables
  Besides from the must-have integrations: "dk1t, dphi_k1, dkat, dphi_ka, dy1,
  dy2" We will select which of the following: "dPT_rel, dphi_PT, dpt_imb,
  dphi_pt, dbp" is required to integrate out, "dphi_bp" was integrated out
  analytically. This is determined by selecting the value of "diff". 0: total
  cross-section, 11D integration. 1: dsigma/dPT_rel = dphi_PT, dpt_imb, dphi_pt,
  dbp 2: dsigma/dpt_imb = dPT_rel, dphi_PT, dphi_PT, dbp 3: dsigma/dbp =
  dPT_rel, dphi_PT, dpt_imb, dphi_pt 4: dsigma/dal = dPT_rel, dpt_imb, dphi_Pp,
  dbp 5: user define
  */
  size_t ndim = 10;
  std::vector<double> dx_lower(ndim);
  std::vector<double> dx_upper(ndim);

  dx_lower[0] = k.ktmin, dx_upper[0] = k.ktmax;
  dx_lower[1] = 0.0, dx_upper[1] = phys::twoPI;
  dx_lower[2] = k.ktmin, dx_upper[2] = k.ktmax;
  dx_lower[3] = 0.0, dx_upper[3] = phys::twoPI;
  dx_lower[4] = k.ylmin, dx_upper[4] = k.ylmax;
  dx_lower[5] = k.ylmin, dx_upper[5] = k.ylmax;

  size_t nbin;
  double hmin;
  double hmax;
  switch (k.diff) {
    case 0:       // sigma total
      ndim += 1;  // add one more dimension
      dx_lower[6] = k.Ptmin, dx_upper[6] = k.Ptmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = k.ptmin, dx_upper[8] = k.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      dx_lower.push_back(k.bpmin / phys::GeVfm);
      dx_upper.push_back(k.bpmax / phys::GeVfm);
      nbin = 1;
      hmin = 0.0;
      hmax = 0.0;
      break;
    case 1:  // dsigma/dPT_rel = dphi_PT, dpt_imb, dphi_pt, dbp
      dx_lower[6] = 0.0, dx_upper[6] = phys::twoPI;
      dx_lower[7] = k.ptmin, dx_upper[7] = k.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = k.bpmin / phys::GeVfm, dx_upper[9] = k.bpmax / phys::GeVfm;
      nbin = k.Ptn;
      hmin = k.Ptmin;
      hmax = k.Ptmax;
      break;
    case 2:  // dsigma/dpt_imb = dPT_rel, dphi_PT, dphi_PT, dbp
      dx_lower[6] = k.Ptmin, dx_upper[6] = k.Ptmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = k.bpmin / phys::GeVfm, dx_upper[9] = k.bpmax / phys::GeVfm;
      nbin = k.ptn;
      hmin = k.ptmin;
      hmax = k.ptmax;
      break;
    case 3:  // dsigma/dbp = dPT_rel, dphi_PT, dpt_imb, dphi_pt
      dx_lower[6] = k.Ptmin, dx_upper[6] = k.Ptmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = k.ptmin, dx_upper[8] = k.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      nbin = k.bpn;
      hmin = k.bpmin / phys::GeVfm;
      hmax = k.bpmax / phys::GeVfm;
      break;
    case 4:  // dsigma/dal = dPT_rel, dpt_imb, dphi, dbp
      dx_lower[6] = k.Ptmin, dx_upper[6] = k.Ptmax;
      dx_lower[7] = k.ptmin, dx_upper[7] = k.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = k.bpmin / phys::GeVfm, dx_upper[9] = k.bpmax / phys::GeVfm;
      nbin = k.aln;
      hmin = k.almin;
      hmax = k.almax;
      break;
    default:
      std::cerr << "Invalid diff value: " << k.diff << std::endl;
      return 1;
  }
  double bin = (hmax - hmin) / static_cast<double>(nbin);
  std::vector<double> bin_mid(nbin), results(nbin), errors(nbin);

  // add two more dimension of integration if Sudakov is active
  bool sud = false;
  if (sud) {
    ndim += 2;
    dx_lower.push_back(k.ltmin), dx_upper.push_back(k.ltmax);
    dx_lower.push_back(0.0), dx_upper.push_back(phys::twoPI);
  }

  // set lepton mass
  switch (k.lep) {
    case 1:
      k.m_lep = phys::M_ele;
      break;
    case 2:
      k.m_lep = phys::M_mu;
      break;
    case 3:
      k.m_lep = phys::M_tau;
      break;
    default:
      break;
  }

// define OpenMP parallel region
#ifdef _OPENMP
#pragma omp parallel for
#endif

  // perform integration loop
  for (size_t i = 0; i < nbin; ++i) {
    // define bin parameters
    double binL = hmin + static_cast<double>(i) * bin;
    double binR = hmin + static_cast<double>(i + 1) * bin;
    bin_mid[i] = (binL + binR) / 2.0;
    double res, err;

    // local kinematics object (copied values from global k)
    kinematics k_local = k;
    switch (k.diff) {
      case 0:
        break;
      case 1:
        k_local.PT_rel = bin_mid[i];
        break;
      case 2:
        k_local.pt_imb = bin_mid[i];
        break;
      case 3:
        k_local.bp = bin_mid[i];
        break;
      case 4:
        k_local.al = bin_mid[i];
        break;
      default:
        break;
    }

    // local GSL monte rng and state
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
    gsl_monte_function gmf = {&integrand, ndim, &k_local};
    gsl_monte_vegas_params vp;

    // warmup run: iteration=w, calls=n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 0;
    vp.iterations = itm1;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              ncall1, r, s, &res, &err);

    // final run: iteration=1, calls=f*n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 2;
    vp.iterations = itm2;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              ncall2, r, s, &res, &err);

    // free resources
    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    // #ifdef _OPENMP
    // std::cout << "Thread " << omp_get_thread_num() << ", index " << i << "
    // finished" << std::endl; #endif

    // store result and error in array
    results[i] = res;
    errors[i] = err;
  }

  // define print function
  auto print = [](std::ostream &out, const std::vector<double> &x,
                  const std::vector<double> &y, const std::vector<double> &e) {
    for (size_t i = 0; i < x.size(); ++i)
      out << x[i] << '\t' << y[i] << '\t' << e[i] << std::endl;
  };

  // output results at the end
  print(std::cout, bin_mid, results, errors);
  std::ofstream fout("results.txt", std::ios::out);
  print(fout, bin_mid, results, errors);
  fout.close();

  // display elapsed time
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;

  return 0;
}
