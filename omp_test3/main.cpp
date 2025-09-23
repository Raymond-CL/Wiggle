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
  // parameters from input file
  double CME;
  double ktmin, ktmax;
  double ylmin, ylmax;
  int diff;
  double PTmin, PTmax;
  double ptmin, ptmax;
  double bpmin, bpmax;
  double almin, almax;
  size_t PTn, ptn, bpn, aln;
  bool do_sudakov;
  double ltmin, ltmax;
  bool do_Mllcut;
  double Mllmin, Mllmax;
  bool do_Ajcut;
  double Ajmin, Ajmax;
  bool do_decoherent;
  bool do_aniso;
  int lep;
  size_t ncall1, itm1;
  size_t ncall2, itm2;
  // parameters in main needed in integrand
  double PT_rel, pt_imb, bp, al;
  double m_lep;
};

template <typename... Args>
bool read_line(std::ifstream &fin, Args &...args);
void read_file(const std::string &filename, parameters &p);
void print_info(const parameters &p);
void init_limits(size_t &ndim, std::vector<double> &dx_lower,
                 std::vector<double> &dx_upper, const parameters &p);

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

  vec2d Dt;
  Dt = k1t - kat;
  // double k1ka = Dt.mag();
  // double phi_Dt = Dt.phi();
  // double k2kb = (k2t - kbt).mag();

  double bessel = p->bp / phys::twoPI * gsl_sf_bessel_J0(p->bp * Dt.mag());
  // double bessel = p->bp / phys::twoPI * gsl_sf_bessel_Jn(4, p->bp *
  // Dt.mag());

  // double costerm = -2.0 * std::cos(4.0 * PT.phi() - 4.0 * pt.phi());
  // double costerm = -2.0 * std::cos(4.0 * PT.phi() - 2.0 * pt.phi() - 2.0 *
  // Dt.phi()); double costerm = -2.0 * std::cos(4.0 * pt.phi() - 4.0 *
  // Dt.phi());
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

  // read input parameters from file
  read_file("input.dat", p);

  // print program details
  print_info(p);
  
  size_t ndim;
  std::vector<double> dx_lower;
  std::vector<double> dx_upper;
  init_limits(ndim, dx_lower, dx_upper, p);

  dx_lower[0] = p.ktmin, dx_upper[0] = p.ktmax;
  dx_lower[1] = 0.0, dx_upper[1] = phys::twoPI;
  dx_lower[2] = p.ktmin, dx_upper[2] = p.ktmax;
  dx_lower[3] = 0.0, dx_upper[3] = phys::twoPI;
  dx_lower[4] = p.ylmin, dx_upper[4] = p.ylmax;
  dx_lower[5] = p.ylmin, dx_upper[5] = p.ylmax;

  size_t nbin;
  double hmin;
  double hmax;
  switch (p.diff) {
    case 0:       // sigma total
      ndim += 1;  // add one more dimension
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = p.ptmin, dx_upper[8] = p.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      dx_lower.push_back(p.bpmin / phys::GeVfm);
      dx_upper.push_back(p.bpmax / phys::GeVfm);
      nbin = 1;
      hmin = 0.0;
      hmax = 0.0;
      break;
    case 1:  // dsigma/dPT_rel = dphi_PT, dpt_imb, dphi_pt, dbp
      dx_lower[6] = 0.0, dx_upper[6] = phys::twoPI;
      dx_lower[7] = p.ptmin, dx_upper[7] = p.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      nbin = p.PTn;
      hmin = p.PTmin;
      hmax = p.PTmax;
      break;
    case 2:  // dsigma/dpt_imb = dPT_rel, dphi_PT, dphi_PT, dbp
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      nbin = p.ptn;
      hmin = p.ptmin;
      hmax = p.ptmax;
      break;
    case 3:  // dsigma/dbp = dPT_rel, dphi_PT, dpt_imb, dphi_pt
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = p.ptmin, dx_upper[8] = p.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      nbin = p.bpn;
      hmin = p.bpmin / phys::GeVfm;
      hmax = p.bpmax / phys::GeVfm;
      break;
    case 4:  // dsigma/dal = dPT_rel, dpt_imb, dphi, dbp
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = p.ptmin, dx_upper[7] = p.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      nbin = p.aln;
      hmin = p.almin;
      hmax = p.almax;
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
    dx_lower.push_back(p.ltmin), dx_upper.push_back(p.ltmax);
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
#ifdef _OPENMP
#pragma omp critical
    {
      std::cout << "Working on bin: " << i
                << "\t under thread: " << omp_get_thread_num() << std::endl;
    }
#else
    std::cout << "Working on bin: " << i << std::endl;
#endif

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
    vp.iterations = p.itm1;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              p.ncall1, r, s, &res, &err);

    // final run: iteration=1, calls=f*n
    gsl_monte_vegas_params_get(s, &vp);
    vp.stage = 2;
    vp.iterations = p.itm2;
    gsl_monte_vegas_params_set(s, &vp);
    gsl_monte_vegas_integrate(&gmf, dx_lower.data(), dx_upper.data(), ndim,
                              p.ncall2, r, s, &res, &err);

    // free resources
    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

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

  // print to console
  std::cout << "----------------------------------------" << std::endl
            << std::scientific << std::setprecision(6)
            << "#   x    \t    y    \t   error  " << std::endl;
  print(std::cout, bin_mid, results, errors);
  // print to file
  std::ofstream fout("results.txt", std::ios::out);
  print(fout, bin_mid, results, errors);
  fout.close();

  // display elapsed time
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "----------------------------------------" << std::endl
            << "Elapsed time: " << std::defaultfloat << elapsed.count()
            << " seconds\n"
            << "----------------------------------------" << std::endl;

  return 0;
}

// read helper function
template <typename... Args>
bool read_line(std::ifstream &fin, Args &...args) {
  std::string line;
  while (std::getline(fin, line)) {
    auto pos = line.find('#');
    if (pos != std::string::npos) line = line.substr(0, pos);
    if (line.empty()) continue;
    std::istringstream iss(line);
    if ((iss >> std::boolalpha >> ... >> args)) return true;
  }
  return false;
}

// read from file
void read_file(const std::string &filename, parameters &p) {
  std::ifstream fin(filename, std::ios::in);
  // if (!fin) return 1;
  read_line(fin, p.CME);
  read_line(fin, p.ktmin, p.ktmax);
  read_line(fin, p.ylmin, p.ylmax);
  read_line(fin, p.diff);
  read_line(fin, p.PTmin, p.PTmax, p.PTn);
  read_line(fin, p.ptmin, p.ptmax, p.ptn);
  read_line(fin, p.bpmin, p.bpmax, p.bpn);
  read_line(fin, p.almin, p.almax, p.aln);
  read_line(fin, p.do_sudakov);
  read_line(fin, p.ltmin, p.ltmax);
  read_line(fin, p.do_Mllcut);
  read_line(fin, p.Mllmin, p.Mllmax);
  read_line(fin, p.do_Ajcut);
  read_line(fin, p.Ajmin, p.Ajmax);
  read_line(fin, p.do_decoherent);
  read_line(fin, p.do_aniso);
  read_line(fin, p.lep);
  read_line(fin, p.ncall1, p.itm1);
  read_line(fin, p.ncall2, p.itm2);
  fin.close();
}

// print input parameters
void print_info(const parameters &p) {
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Wiggle: Dilepton production in UPC" << std::endl;
#if defined(__GNUC__) || defined(__GNUG__)
  std::cout << "Using GCC version " << __VERSION__ << std::endl;
#elif defined(__clang__)
  std::cout << "Using Clang version " << __clang_version__ << std::endl;
#elif defined(_MSC_VER)
  std::cout << "Using MSVC version " << _MSC_VER << std::endl;
#elif defined(__INTEL_COMPILER)
  std::cout << "Using Intel C++ version " << __INTEL_COMPILER << std::endl;
#else
  std::cout << "Using unknown compiler" << std::endl;
#endif
  std::cout << "Using C++ standard version: " << __cplusplus << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Center-of-mass energy: " << p.CME << " GeV" << std::endl;
  std::cout << "k1t: [" << p.ktmin << ", " << p.ktmax << "] GeV" << std::endl;
  std::cout << "y1, y2: [" << p.ylmin << ", " << p.ylmax << "]" << std::endl;
  std::cout << "Differential type: " << p.diff << std::endl;
  switch (p.diff) {
    case 0:
      std::cout << "Total cross-section" << std::endl;
      break;
    case 1:
      std::cout << "dSigma/dPT_rel with PT_rel in [" << p.PTmin << ", "
                << p.PTmax << "] GeV"
                << " with " << p.PTn << " bins" << std::endl;
      break;
    case 2:
      std::cout << "dSigma/dpt_imb with pt_imb in [" << p.ptmin << ", "
                << p.ptmax << "] GeV"
                << " with " << p.ptn << " bins" << std::endl;
      break;
    case 3:
      std::cout << "dSigma/dbp with bp in [" << p.bpmin << ", " << p.bpmax
                << "] fm"
                << " with " << p.bpn << " bins" << std::endl;
      break;
    case 4:
      std::cout << "dSigma/dal with al in [" << p.almin << ", " << p.almax
                << "]"
                << " with " << p.aln << " bins" << std::endl;
      break;
    default:
      std::cerr << "Invalid diff value: " << p.diff << std::endl;
      // return 1;
  }
  if (p.do_sudakov)
    std::cout << "Including Sudakov factor with lt in [" << p.ltmin << ", "
              << p.ltmax << "] GeV" << std::endl;
  if (p.do_Mllcut)
    std::cout << "Including dilepton mass cut with Mll in [" << p.Mllmin << ", "
              << p.Mllmax << "] GeV" << std::endl;
  if (p.do_Ajcut)
    std::cout << "Including Aj cut with Aj in [" << p.Ajmin << ", " << p.Ajmax
              << "]" << std::endl;
  if (p.do_decoherent)
    std::cout << "Including decoherent contribution" << std::endl;
  if (p.do_aniso)
    std::cout << "Including anisotropic contribution" << std::endl;
  switch (p.lep) {
    case 1:
      std::cout << "Lepton type: electron" << std::endl;
      break;
    case 2:
      std::cout << "Lepton type: muon" << std::endl;
      break;
    case 3:
      std::cout << "Lepton type: tau" << std::endl;
      break;
    default:
      std::cout << "Lepton type: all" << std::endl;
      break;
  }
  std::cout << "Warm-up: " << p.ncall1 << " calls for " << p.itm1
            << " iterations" << std::endl;
  std::cout << "Final:   " << p.ncall2 << " calls for " << p.itm2
            << " iterations" << std::endl;
#ifdef _OPENMP
  std::cout << "Using OpenMP with " << omp_get_max_threads() << " threads"
            << std::endl;
#else
  std::cout << "Not using OpenMP" << std::endl;
#endif
  std::cout << "----------------------------------------" << std::endl;
}

void init_limits(size_t &ndim, std::vector<double> &dx_lower,
                 std::vector<double> &dx_upper, const parameters &p) {
  dx_upper = {p.ktmax, phys::twoPI, p.ktmax, phys::twoPI, p.ylmax, p.ylmax};
  dx_lower = {p.ktmin, 0.0, p.ktmin, 0.0, p.ylmin, p.ylmin};

  switch (p.diff) {
    case 0:       // sigma total
      ndim = 11;
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = p.ptmin, dx_upper[8] = p.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      dx_lower.push_back(p.bpmin / phys::GeVfm);
      dx_upper.push_back(p.bpmax / phys::GeVfm);
      break;
    case 1:  // dsigma/dPT_rel = dphi_PT, dpt_imb, dphi_pt, dbp
      ndim = 10;
      dx_lower[6] = 0.0, dx_upper[6] = phys::twoPI;
      dx_lower[7] = p.ptmin, dx_upper[7] = p.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      break;
    case 2:  // dsigma/dpt_imb = dPT_rel, dphi_PT, dphi_PT, dbp
      ndim = 10;
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      break;
    case 3:  // dsigma/dbp = dPT_rel, dphi_PT, dpt_imb, dphi_pt
      ndim = 10;
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = 0.0, dx_upper[7] = phys::twoPI;
      dx_lower[8] = p.ptmin, dx_upper[8] = p.ptmax;
      dx_lower[9] = 0.0, dx_upper[9] = phys::twoPI;
      break;
    case 4:  // dsigma/dal = dPT_rel, dpt_imb, dphi, dbp
      ndim = 10;
      dx_lower[6] = p.PTmin, dx_upper[6] = p.PTmax;
      dx_lower[7] = p.ptmin, dx_upper[7] = p.ptmax;
      dx_lower[8] = 0.0, dx_upper[8] = phys::twoPI;
      dx_lower[9] = p.bpmin / phys::GeVfm, dx_upper[9] = p.bpmax / phys::GeVfm;
      break;
    default:
      std::cerr << "Invalid diff value: " << p.diff << std::endl;
      // return 1;
  }

  if (p.do_sudakov) {
    ndim += 2;
    dx_lower.push_back(p.ltmin), dx_upper.push_back(p.ltmax);
    dx_lower.push_back(0.0), dx_upper.push_back(phys::twoPI);
  }
}