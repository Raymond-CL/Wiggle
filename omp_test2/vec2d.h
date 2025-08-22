#pragma once
#include <cmath>

class vec2d {
 private:
  mutable double x_{0.0}, y_{0.0};
  mutable double r_{0.0}, phi_{0.0};

  mutable bool cart_dirty{false};
  mutable bool polar_dirty{false};

  static double normalize_angle(double a) noexcept {
    while (a <= -M_PI) a += 2 * M_PI;
    while (a >  M_PI)  a -= 2 * M_PI;
    return a;
  }

  void updateCartesian() const noexcept {
    if (cart_dirty) {
      x_ = r_ * std::cos(phi_);
      y_ = r_ * std::sin(phi_);
      cart_dirty = false;
    }
  }

  void updatePolar() const noexcept {
    if (polar_dirty) {
      r_   = std::sqrt(x_ * x_ + y_ * y_);
      phi_ = std::atan2(y_, x_);
      polar_dirty = false;
    }
  }

 public:
  // Constructors
  vec2d() noexcept = default;
  // constructor in cartesian
  vec2d(double x, double y) noexcept : x_(x), y_(y) {
    polar_dirty = true;  // polar not initialized
  }
  // constructor in polar (return)
  static vec2d fromPolar(double r, double phi) noexcept {
    vec2d v;
    v.r_   = r;
    v.phi_ = normalize_angle(phi);
    v.cart_dirty = true;   // Cartesian not initialized
    return v;
  }

  // --- Setters ---
  inline void setCartesian(double x, double y) noexcept {
    x_ = x; y_ = y;
    polar_dirty = true;
    cart_dirty  = false;
  }

  inline void setPolar(double r, double phi) noexcept {
    r_ = r; phi_ = normalize_angle(phi);
    cart_dirty  = true;
    polar_dirty = false;
  }

  // --- Getters ---
  inline double x() const noexcept {
    updateCartesian();
    return x_;
  }
  inline double y() const noexcept {
    updateCartesian();
    return y_;
  }
  inline double r() const noexcept {
    updatePolar();
    return r_;
  }
  inline double phi() const noexcept {
    updatePolar();
    return phi_;
  }

  inline double mag2() const noexcept {
    updatePolar();
    return r_ * r_;
  }
  inline double mag() const noexcept {
    updatePolar();
    return r_;
  }

  // --- Vector algebra ---
  inline vec2d operator+(const vec2d& v) const noexcept {
    // Addition is cheap in Cartesian
    return vec2d(x() + v.x(), y() + v.y());
  }

  inline vec2d operator-(const vec2d& v) const noexcept {
    return vec2d(x() - v.x(), y() - v.y());
  }

  inline vec2d operator*(double s) const noexcept {
    // Scaling is cheap in Polar
    return vec2d::fromPolar(r() * s, phi());
  }

  inline friend vec2d operator*(double s, const vec2d& v) noexcept {
    return v * s;
  }

  inline vec2d& operator+=(const vec2d& v) noexcept {
    setCartesian(x() + v.x(), y() + v.y());
    return *this;
  }

  inline vec2d& operator-=(const vec2d& v) noexcept {
    setCartesian(x() - v.x(), y() - v.y());
    return *this;
  }

  inline vec2d& operator*=(double s) noexcept {
    setPolar(r() * s, phi());
    return *this;
  }

  inline vec2d operator-() const noexcept {
    return vec2d::fromPolar(r(), normalize_angle(phi() + M_PI));
  }

  // --- Dot and cross ---
  inline double dot(const vec2d& v) const noexcept {
    return x() * v.x() + y() * v.y();
  }
  inline double cross(const vec2d& v) const noexcept {
    return x() * v.y() - y() * v.x();  // scalar in 2D
  }
};
