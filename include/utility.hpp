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

#ifndef UTILITY_HEADER_DEF
#define UTILITY_HEADER_DEF

#include <fstream>
#include <string>
#include <vector>
#include "vec3.hpp"

// forward declarations
class wetted_surface;

// function declarations
void GetForceData(const int &inc, const std::string &chemCase,
                  const double &fsf, wetted_surface &fem);
bool IsFaceValid(const wetted_surface &fem, const int &fid);
void WriteMasterFile(const wetted_surface &fem, const std::string &name,
                     const bool &mapNodes);
std::string GetFileName(const std::string &rootName, const int &ii,
                        const bool &mapNodes);
std::vector<vec3> GetFaceForce(const wetted_surface &fem,
                               const std::vector<int> &adjacentFaces);
std::vector<vec3> GetNodeForce(const wetted_surface &fem);
std::vector<vec3> GetNodeMoment(const wetted_surface &fem);
std::vector<int> GetAdjacentFaces(const wetted_surface &fem);
vec3 GetForce(const wetted_surface &fem, const int &id);
vec3 GetMoment(const wetted_surface &fem, const int &id);
void PrintStats(const wetted_surface &fem, const wetted_surface &cfd,
                const std::vector<vec3> &femFace, const vec3 &momentCenter);
void PrintStatsNodes(const wetted_surface &fem, const wetted_surface &cfd,
                     const vec3 &momentCenter);
vec3 SumForces(const wetted_surface &surf);
vec3 SumMoments(const wetted_surface &surf, const vec3 &center);
vec3 SumMomentsFace(const wetted_surface &surf, const std::vector<vec3> &force,
                    const vec3 &center);
vec3 GetNode(const wetted_surface &fem, const int &id);
vec3 GetVecFromString(const std::string &str);
void BroadcastString(std::string &str);
void CreateDirectory(const std::string &dname, const int &rank);

#endif
