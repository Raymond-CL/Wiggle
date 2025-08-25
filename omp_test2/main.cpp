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

struct parameters {
  double CME;
  int diff;
  double PT_rel, pt_imb, bp, al;
  bool do_sudakov;
  bool do_Mllcut;
  double Mllmin, Mllmax;
  bool do_Ajcut;
  double Ajmin, Ajmax;
  bool do_decoherent;
  bool do_aniso;
  int lep;
  double m_lep;
};

double integrand(double *dx, size_t ndim, void *params) {
  (void)(ndim);  // unused
  auto *p = static_cast<parameters *>(params);

  vec2d k1t, kat;
  k1t.setPolar(dx[0], dx[1]);
  kat.setPolar(dx[2], dx[3]);
  double y1 = dx[4];
  double y2 = dx[5];

  vec2d PT, pt;
  switch (p->diff) {
    case 0:  // total cross-section
      PT.setPolar(dx[6], dx[7]);
      pt.setPolar(dx[8], dx[9]);
      p->bp = dx[10];
      break;
    case 1:  // dsigma/dPT_rel
      PT.setPolar(p->PT_rel, dx[6]);
      pt.setPolar(dx[7], dx[8]);
      p->bp = dx[9];
      break;
    case 2:  // dsigma/dpt_imb
      PT.setPolar(dx[6], dx[7]);
      pt.setPolar(p->pt_imb, dx[8]);
      p->bp = dx[9];
      break;
    case 3:  // dsigma/dbp
      PT.setPolar(dx[6], dx[7]);
      pt.setPolar(dx[8], dx[9]);
      break;
    case 4:  // dsigma/dal
      p->PT_rel = dx[6];
      p->pt_imb = dx[7];
      // p->phi_pt = dx[8];
      p->bp = dx[9];
      break;
    default:
      std::cerr << "Invalid diff value: " << p->diff << std::endl;
      return 0.0;
  }

  vec2d p1t, p2t;
  p1t = 0.5 * pt + PT;
  p2t = 0.5 * pt - PT;

  vec2d lt;
  if (p->do_sudakov) {
    lt.setPolar(dx[10], dx[11]);  // not correct for total cross-section
  } else {
    lt.setPolar(0.0, 0.0);
  }

  double pthard = std::max(p1t.mag(), p2t.mag());
  // double pthard = PT.mag();
  double x1 = pthard * (std::exp(+y1) + std::exp(+y2)) / p->CME;
  double x2 = pthard * (std::exp(-y1) + std::exp(-y2)) / p->CME;
  if (x1 <= 0.0 || x1 >= 1.0) return 0.0;
  if (x2 <= 0.0 || x2 >= 1.0) return 0.0;

  double mans = +x1 * x2 * p->CME * p->CME;
  double mant = -x1 * p->CME * pthard * std::exp(-y1);
  double manu = -x1 * p->CME * pthard * std::exp(-y2);
  double Mll = std::sqrt(mans);
  if (Mll <= p->Mllmin || Mll >= p->Mllmax) return 0.0;

  vec2d qt;
  qt = pt - lt;
  vec2d k2t, kbt;
  k2t = qt - k1t;
  kbt = qt - kat;
  double k1sq = x1 * x1 * phys::Mp2 + k1t.mag2();
  double kasq = x1 * x1 * phys::Mp2 + kat.mag2();
  double k2sq = x2 * x2 * phys::Mp2 + k2t.mag2();
  double kbsq = x2 * x2 * phys::Mp2 + kbt.mag2();
  double k1ka = (k1t - kat).mag();
  // double k2kb = (k2t - kbt).mag();

  double iso = std::cos(k1t.phi() - kat.phi() + k2t.phi() - kbt.phi()) *
               (mant * mant + manu * manu) / (mant * manu);

  double aniso = 0.0;
  if (p->do_aniso) {
    aniso = 2.0 * std::cos(k1t.phi() + kat.phi() + k2t.phi() + kbt.phi() -
                           4.0 * PT.phi());
  }
  double sigma0 = k1t.mag() * kat.mag() * k2t.mag() * kbt.mag() *
                  (iso - aniso) * 2.0 * phys::alphaesq / (mans * mans);

  auto ffnuc = [](double k) {
    constexpr double R = 6.62 / phys::GeVfm;
    double kR = k * R;
    return (std::sin(kR) - kR * std::cos(kR)) * 3.0 /
           ((kR * kR * kR) * (1.0 + phys::a02 * k * k));
  };

  auto xf = [&ffnuc](double k, double kp, parameters *pp) {
    if (k <= 0.0 || kp <= 0.0) return 0.0;
    constexpr int atomZ2 = 82 * 82;
    constexpr double fac = atomZ2 * phys::alphae / phys::PIsq;
    double nff = ffnuc(k) * ffnuc(kp) / (k * k * kp * kp);
    if (pp->do_decoherent) {
      ;
    }
    return fac * nff;
  };

  double x1f1 = xf(std::sqrt(k1sq), std::sqrt(kasq), p);
  double x2f2 = xf(std::sqrt(k2sq), std::sqrt(kbsq), p);

  double bessel = p->bp / phys::twoPI * gsl_sf_bessel_J0(p->bp * k1ka);

  // double costerm = -2.0 * std::cos(4.0 * k->phi_PT - 4.0 * k->phi_pt);
  double costerm = 1.0;

  double res = PT.mag() * pt.mag() * k1t.mag() * kat.mag() * bessel * x1f1 *
               x2f2 * sigma0 * costerm;

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

  // retrieve kinematics parameters from input file
  parameters p;

  double ktmin, ktmax;
  double ylmin, ylmax;
  double ltmin, ltmax;
  double PTmin, PTmax;
  double ptmin, ptmax;
  double bpmin, bpmax;
  double almin, almax;
  size_t PTn, ptn, bpn, aln;

  // integration calls and iterations
  size_t ncall1, itm1;
  size_t ncall2, itm2;

  std::ifstream fin("input.dat", std::ios::in);
  if (!fin) return 1;
  fin >> p.CME;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ktmin >> ktmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ylmin >> ylmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> p.diff;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> PTmin >> PTmax >> PTn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ptmin >> ptmax >> ptn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> bpmin >> bpmax >> bpn;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> almin >> almax >> aln;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> std::boolalpha >> p.do_sudakov;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ltmin >> ltmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> std::boolalpha >> p.do_Mllcut;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> p.Mllmin >> p.Mllmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> std::boolalpha >> p.do_Ajcut;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> p.Ajmin >> p.Ajmax;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> std::boolalpha >> p.do_decoherent;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> std::boolalpha >> p.do_aniso;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> p.lep;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ncall1 >> itm1;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin >> ncall2 >> itm2;
  fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  fin.close();

  size_t ndim = 10;
  std::vector<double> dx_lower(ndim);
  std::vector<double> dx_upper(ndim);

  dx_lower[0] = ktmin, dx_upper[0] = ktmax;
  dx_lower[1] = 0.0, dx_upper[1] = phys::twoPI;
  dx_lower[2] = ktmin, dx_upper[2] = ktmax;
  dx_lower[3] = 0.0, dx_upper[3] = phys::twoPI;
  dx_lower[4] = ylmin, dx_upper[4] = ylmax;
  dx_lower[5] = ylmin, dx_upper[5] = ylmax;

  size_t nbin;
  double hmin;
  double hmax;
  switch (p.diff) {
    case 0:       // sigma total
      ndim += 1;  // add one more dimension
      dx_lower[6] = PTmin, dx_upper[6] = PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = ptmin, dx_upper[8] = ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      dx_lower.push_back(bpmin / phys::GeVfm);
      dx_upper.push_back(bpmax / phys::GeVfm);
      nbin = 1;
      hmin = 0.0;
      hmax = 0.0;
      break;
    case 1:  // dsigma/dPT_rel = dphi_PT, dpt_imb, dphi_pt, dbp
      dx_lower[6] = 0.0, dx_upper[6] = phys::twoPI;
      dx_lower[7] = ptmin, dx_upper[7] = ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = bpmin / phys::GeVfm, dx_upper[9] = bpmax / phys::GeVfm;
      nbin = PTn;
      hmin = PTmin;
      hmax = PTmax;
      break;
    case 2:  // dsigma/dpt_imb = dPT_rel, dphi_PT, dphi_PT, dbp
      dx_lower[6] = PTmin, dx_upper[6] = PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = bpmin / phys::GeVfm, dx_upper[9] = bpmax / phys::GeVfm;
      nbin = ptn;
      hmin = ptmin;
      hmax = ptmax;
      break;
    case 3:  // dsigma/dbp = dPT_rel, dphi_PT, dpt_imb, dphi_pt
      dx_lower[6] = PTmin, dx_upper[6] = PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = ptmin, dx_upper[8] = ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      nbin = bpn;
      hmin = bpmin / phys::GeVfm;
      hmax = bpmax / phys::GeVfm;
      break;
    case 4:  // dsigma/dal = dPT_rel, dpt_imb, dphi, dbp
      dx_lower[6] = PTmin, dx_upper[6] = PTmax;
      dx_lower[7] = ptmin, dx_upper[7] = ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = bpmin / phys::GeVfm, dx_upper[9] = bpmax / phys::GeVfm;
      nbin = aln;
      hmin = almin;
      hmax = almax;
      break;
    default:
      std::cerr << "Invalid diff value: " << p.diff << std::endl;
      return 1;
  }
  double bin = (hmax - hmin) / static_cast<double>(nbin);
  std::vector<double> bin_mid(nbin), results(nbin), errors(nbin);

  // add two more dimension of integration if Sudakov is active
  if (p.do_sudakov) {
    ndim += 2;
    dx_lower.push_back(ltmin), dx_upper.push_back(ltmax);
    dx_lower.push_back(0.0), dx_upper.push_back(phys::twoPI);
  }

  // set lepton mass
  switch (p.lep) {
    case 1:
      p.m_lep = phys::M_ele;
      break;
    case 2:
      p.m_lep = phys::M_mu;
      break;
    case 3:
      p.m_lep = phys::M_tau;
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
    parameters p_local = p;
    switch (p_local.diff) {
      case 0:
        break;
      case 1:
        p_local.PT_rel = bin_mid[i];
        break;
      case 2:
        p_local.pt_imb = bin_mid[i];
        break;
      case 3:
        p_local.bp = bin_mid[i];
        break;
      case 4:
        p_local.al = bin_mid[i];
        break;
      default:
        break;
    }

    // local GSL monte rng and state
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
    gsl_monte_function gmf = {&integrand, ndim, &p_local};
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
