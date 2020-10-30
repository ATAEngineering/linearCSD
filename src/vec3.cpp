// Copyright (C) 2020, ATA Engineering, Inc.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <iostream>
#include "vec3.hpp"

using std::cout;
using std::cerr;
using std::endl;


// -------------------------------------------------------------------------
// function definitions
vec3 &vec3::operator+=(const vec3 &v) {
  std::transform(this->begin(), this->end(), v.begin(), this->begin(),
                 std::plus<double>());
  return *this;
}
vec3 &vec3::operator-=(const vec3 &v) {
  std::transform(this->begin(), this->end(), v.begin(), this->begin(),
                 std::minus<double>());
  return *this;
}
vec3 &vec3::operator*=(const vec3 &v) {
  std::transform(this->begin(), this->end(), v.begin(), this->begin(),
                 std::multiplies<double>());
  return *this;
}
vec3 &vec3::operator/=(const vec3 &v) {
  std::transform(this->begin(), this->end(), v.begin(), this->begin(),
                 std::divides<double>());
  return *this;
}

vec3 operator+(vec3 lhs, const vec3 &rhs) { return lhs += rhs; }
vec3 operator-(vec3 lhs, const vec3 &rhs) { return lhs -= rhs; }
vec3 operator*(vec3 lhs, const vec3 &rhs) { return lhs *= rhs; }
vec3 operator/(vec3 lhs, const vec3 &rhs) { return lhs /= rhs; }

// overloads for scalar addition
vec3 & vec3::operator+=(const double &s) {
  std::for_each(this->begin(), this->end(), [&s](auto &val) { val += s; });
  return *this;
}
vec3 & vec3::operator-=(const double &s) {
  std::for_each(this->begin(), this->end(), [&s](auto &val) { val -= s; });
  return *this;
}
vec3 & vec3::operator*=(const double &s) {
  std::for_each(this->begin(), this->end(), [&s](auto &val) { val *= s; });
  return *this;
}
vec3 & vec3::operator/=(const double &s) {
  std::for_each(this->begin(), this->end(), [&s](auto &val) { val /= s; });
  return *this;
}

vec3 operator+(const double &lhs, vec3 rhs) { return rhs += lhs; }
vec3 operator-(const double &lhs, vec3 rhs) {
  std::for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs - val; });
  return rhs;
}
vec3 operator*(const double &lhs, vec3 rhs) { return rhs *= lhs; }
vec3 operator/(const double &lhs, vec3 rhs) {
  std::for_each(rhs.begin(), rhs.end(), [&lhs](auto &val) { val = lhs / val; });
  return rhs;
}

std::ostream &operator<<(std::ostream &os, const vec3 &v) {
  os << v[0] << " " << v[1] << " " << v[2] << " ";
  return os;
}

double dot(const vec3 &v1, const vec3 &v2) {
  return std::inner_product(v1.begin(), v1.end(), v2.begin(), 0.0);
}

vec3 cross(const vec3 &a, const vec3 &b) { 
  vec3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

double norm(const vec3 &v) { return std::sqrt(dot(v, v)); }
