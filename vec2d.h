#pragma once
#include <cmath>

class vec2d {
 public:
  double x{}, y{};  // Cartesian components

  // Constructors
  vec2d() = default;
  vec2d(double x_, double y_) : x(x_), y(y_) {}

  // Setters
  inline void setCartesian(double x_, double y_) {
    x = x_;
    y = y_;
  }

  inline void setPolar(double mag, double phi) {
    x = mag * std::cos(phi);
    y = mag * std::sin(phi);
  }

  // Getters
  inline double mag2() const { return x * x + y * y; }
  inline double mag() const { return std::hypot(x, y); }
  inline double phi() const {
    return (x == 0.0 && y == 0.0) ? 0.0 : std::atan2(y, x);
  }

  // Vector algebra
  inline vec2d operator+(const vec2d& v) const { return {x + v.x, y + v.y}; }
  inline vec2d operator-(const vec2d& v) const { return {x - v.x, y - v.y}; }
  inline vec2d operator*(double s) const { return {x * s, y * s}; }
  inline friend vec2d operator*(double s, const vec2d& v) { return v * s; }

  // Dot and cross products
  inline double dot(const vec2d& v) const { return x * v.x + y * v.y; }
  inline double cross(const vec2d& v) const {
    return x * v.y - y * v.x;
  }  // scalar in 2D
};