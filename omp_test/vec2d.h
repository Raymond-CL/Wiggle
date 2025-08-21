#pragma once
#include <cmath>

class vec2d {
 public:
  double x{}, y{};  // Cartesian components

  // Constructors
  vec2d() noexcept = default;
  vec2d(double x_, double y_) noexcept : x(x_), y(y_) {}

  // Setters
  inline void setCartesian(double x_, double y_) noexcept {
    x = x_;
    y = y_;
  }

  void setPolar(double mag, double phi) noexcept {
    x = mag * std::cos(phi);
    y = mag * std::sin(phi);
  }

  // Getters
  inline double mag2() const noexcept { return x * x + y * y; }
  double mag() const noexcept { return std::sqrt(mag2()); }
  double phi() const noexcept { return std::atan2(y, x); }

  // Vector algebra
  inline vec2d operator+(const vec2d& v) const noexcept {
    return {x + v.x, y + v.y};
  }
  inline vec2d operator-(const vec2d& v) const noexcept {
    return {x - v.x, y - v.y};
  }
  inline vec2d operator*(double s) const noexcept { return {x * s, y * s}; }
  inline friend vec2d operator*(double s, const vec2d& v) noexcept {
    return v * s;
  }
  inline vec2d& operator+=(const vec2d& v) noexcept {
    x += v.x;
    y += v.y;
    return *this;
  }
  inline vec2d& operator-=(const vec2d& v) noexcept {
    x -= v.x;
    y -= v.y;
    return *this;
  }
  inline vec2d& operator*=(double s) noexcept {
    x *= s;
    y *= s;
    return *this;
  }
  inline vec2d operator-() const noexcept { return {-x, -y}; }

  // Dot and cross products
  inline double dot(const vec2d& v) const noexcept { return x * v.x + y * v.y; }
  inline double cross(const vec2d& v) const noexcept {
    return x * v.y - y * v.x;
  }  // scalar in 2D
};