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

#ifndef VEC3_HEADER_DEF
#define VEC3_HEADER_DEF

#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

class vec3 {
  double data_[3];

 public:
  // constructor
  vec3(const double &x, const double &y, const double &z) : data_{x, y, z} {}
  vec3() : vec3(0, 0, 0) {}

  // move constructor and assignment operator
  vec3(vec3 &&) = default;
  vec3 &operator=(vec3 &&) = default;

  // copy constructor and assignment operator
  vec3(const vec3 &) = default;
  vec3 &operator=(const vec3 &) = default;

  // member functions
  const double &operator[](const int &ii) const { return data_[ii]; }
  double &operator[](const int &ii) { return data_[ii]; }
  const double &X() const { return data_[0]; }
  const double &Y() const { return data_[1]; }
  const double &Z() const { return data_[2]; }

  // provide bein and end so std::begin and std::end can be used
  auto begin() { return std::begin(data_); }
  const auto begin() const { return std::begin(data_); }
  auto end() { return std::end(data_); }
  const auto end() const { return std::end(data_); }

  // overload math operations
  vec3 &operator+=(const vec3 &);
  vec3 &operator-=(const vec3 &);
  vec3 &operator*=(const vec3 &);
  vec3 &operator/=(const vec3 &);

  vec3 &operator+=(const double &);
  vec3 &operator-=(const double &);
  vec3 &operator*=(const double &);
  vec3 &operator/=(const double &);

  vec3 operator+(const double &s) const { 
    auto lhs = *this;
    return lhs += s;
  }
  vec3 operator-(const double &s) const { 
    auto lhs = *this;
    return lhs -= s;
  }
  vec3 operator*(const double &s) const { 
    auto lhs = *this;
    return lhs *= s;
  }
  vec3 operator/(const double &s) const { 
    auto lhs = *this;
    return lhs /= s;
  }

};

// function definitions
vec3 operator+(vec3 lhs, const vec3 &rhs);
vec3 operator-(vec3 lhs, const vec3 &rhs);
vec3 operator*(vec3 lhs, const vec3 &rhs);
vec3 operator/(vec3 lhs, const vec3 &rhs);

vec3 operator+(const double &lhs, vec3 rhs);
vec3 operator-(const double &lhs, vec3 rhs);
vec3 operator*(const double &lhs, vec3 rhs);
vec3 operator/(const double &lhs, vec3 rhs);

std::ostream &operator<<(std::ostream &os, const vec3 &v);

double dot(const vec3 &v1, const vec3 &v2);
vec3 cross(const vec3 &a, const vec3 &b);
double norm(const vec3 &v1);

#endif
