#ifndef KINEMATICS_H
#define KINEMATICS_H

#include "physconst.h"
#include "vec2d.h"

class kinematics {
 public:

  // kinematics ranges
  double CME;
  double pt_imb, ptmin, ptmax, phi_pt;
  double pl_avg, plmin, plmax, phi_pl;
  double qt_imb, qtmin, qtmax, phi_qt;
  double alpha, almin, almax;
  double ktmin, ktmax;
  double phi_Dt;
  double ylmin, ylmax;
  double pltmin, pltmax;
  double Mllmin, Mllmax;

  // bins
  int qtn;

  // photo and lepton momentum
  vec2d lep1, lep2;
  vec2d phoa, phob, pho1, pho2;
  vec2d pt, qt;
  double rap1, rap2, rapmin, rapmax;

  // transverse mass
  double mtsq, mt;
  double m1tsq, m1t;
  double m2tsq, m2t;
  double m_lep;

  // extra kinematic cut
  double yllmin, yll, yllmax;
  double aj, ajcut;

  // momentum fraction
  double x1,x2;

  // impact parameter
  double bpmin, bp, bpmax;
};

#endif  // KINEMATICS_H