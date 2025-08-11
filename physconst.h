#ifndef PHYSCONST_H
#define PHYSCONST_H

namespace phys {
// mathematical constants
constexpr double PI = 3.14159265358979323846;
constexpr double twoPI = 2.0 * PI;

// quantum constants
constexpr double PLANCK = 6.62607015e-34;
constexpr double PLANCK_BAR = PLANCK / twoPI;
constexpr double CHARGE = 1.602176634e-19;
constexpr double HBAR = PLANCK_BAR / CHARGE;
constexpr double LIGHT_SPEED = 299792458.0;
constexpr double HBARC = HBAR * LIGHT_SPEED;

// unit conversion constants
constexpr double GIGA = 1.0e9;
constexpr double femto = 1.0e-15;
constexpr double GEVfm = HBARC / GIGA / femto;

// particle masses
constexpr double M_P = 0.93827208816;
constexpr double M_N = 0.93956542052;
constexpr double M_E = 0.00051;
constexpr double M_MU = 0.10566;
constexpr double M_TAU = 1.77686;
}

#endif  // PHYSCONST_H