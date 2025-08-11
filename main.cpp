#include <iostream>

#include "kinematics.h"
#include "photonff.h"
#include "vegas.h"

// parameters passed to the integrand function
struct parameters {
  kinematics *kin;
  photonff *pho;
};

// integrand function for VEGAS integration
double integrand(double *dx, size_t dim, void *params_void) {
  auto *params = static_cast<parameters *>(params_void);
  return 1.0;  // params->kin->qt_imb;
  // double rsq = 0.0;
  // for (size_t i = 0; i < dim; ++i) rsq += dx[i] * dx[i];
  // return (rsq <= params->R * params->R) ? 1.0 : 0.0;
}

int main() {
  kinematics *kin = new kinematics();
  photonff *pho = new photonff();

  kin->CME = 5020.0;
  kin->qtmin = 0.0;
  kin->qtmax = 0.1;
  kin->qtn = 50;
  kin->ylmin = -1.0;
  kin->ylmax = +1.0;
  kin->pltmin = 4.0;
  kin->pltmax = 20.0;
  kin->Mllmin = 4.0;
  kin->Mllmax = 45.0;
  kin->ktmin = 0.0;
  kin->ktmax = 10.0;
  kin->m_lep = phys::M_MU;

  pho->atom_A = 207;
  pho->atom_Z = 82;
  pho->radius = 6.62 / phys::GEVfm;

  const size_t DIM = 9;

  double xlow[DIM], xup[DIM];
  xlow[0] = kin->ktmin, xup[0] = kin->ktmax;
  xlow[1] = 0.0, xup[1] = phys::twoPI;
  xlow[2] = kin->ktmin, xup[2] = kin->ktmax;
  xlow[3] = 0.0, xup[3] = phys::twoPI;
  xlow[4] = kin->ylmin, xup[4] = kin->ylmax;
  xlow[5] = kin->ylmin, xup[5] = kin->ylmax;
  xlow[6] = kin->pltmin, xup[6] = kin->pltmax;
  xlow[7] = 0.0, xup[7] = phys::twoPI;
  xlow[8] = 0.0, xup[8] = phys::twoPI;

  vegas *v = new vegas(DIM);
  v->set_limits(xlow, xup);
  parameters params = {kin, pho};
  v->set_integrand(&integrand, &params);

  int nbin = kin->qtn;
  double hmin = kin->qtmin;
  double hmax = kin->qtmax;
  bool islog = false;
  double bin = islog ? log(hmax / hmin) / nbin : (hmax - hmin) / nbin;
  double binL, binM, binR;
  double res, err;
  for (int ind = 1; ind <= nbin; ++ind) {
    binL = islog ? hmin * std::exp(bin * (ind - 1)) : hmin + bin * (ind - 1);
    binR = islog ? hmin * std::exp(bin * ind) : hmin + bin * ind;
    binM = 0.5 * (binL + binR);
    kin->qt_imb = binM;

    v->integrate(1E6, 10, res, err);
    std::cout << kin->qt_imb << '\t' << res << '\t' << err / res << std::endl;
  }

  // clean up
  delete kin;
  delete pho;
  delete v;

  // end program
  return 0;
}