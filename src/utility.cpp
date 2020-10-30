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

#include <string>
#include <iostream>
#include <vector>
#include <numeric>
#include <filesystem>
#include <system_error>
#include "wetted_surface_types.h"
#include "utility.hpp"
#include "inputs.hpp"
#include "vec3.hpp"
#include "mpi.h"

using std::cout;
using std::cerr;
using std::endl;


// -------------------------------------------------------------------------
// function definitions
void GetForceData(const int &inc, const std::string &chemCase,
                  const double &fsf, wetted_surface &cfd) {
  constexpr auto width = 7;
  const auto sinc = std::to_string(inc);
  const auto forceFilePath = chemCase + "_" +
                             std::string(width - sinc.length(), '0') + sinc +
                             ".dat";
  // open force file
  cout << "Reading force history file " << forceFilePath << endl;
  std::ifstream inFile(forceFilePath, std::ios::in);
  if (inFile.fail()) {
    cerr << "ERROR: Force file " << forceFilePath << " did not open correctly!"
         << endl;
    exit(EXIT_FAILURE);
  }

  auto id = 0U;
  std::string line;
  while (std::getline(inFile, line)) {
    // split line into tokens
    line = Trim(line);        // trim whitespace and comments
    if (line.length() > 0) {  // only continue if line has data
      auto tokens = Tokenize(line, " ");
      if (tokens.size() != 3) {
        cerr << "ERROR: Expecting three entries per line in force file" << endl;
        cerr << line << endl;
        exit(EXIT_FAILURE);
      }
      cfd.force[id][0] = std::stod(tokens[0]) * fsf;
      cfd.force[id][1] = std::stod(tokens[1]) * fsf;
      cfd.force[id][2] = std::stod(tokens[2]) * fsf;
      id++;
    }
  }
  inFile.close();
  // sanity checks
  if (id != cfd.force.size()) {
    cerr << "ERROR: Number of entries in force file " << forceFilePath << " is "
         << id << " and does not match number of nodes " << cfd.force.size()
         << endl;
    exit(EXIT_FAILURE);
  }
}

void WriteMasterFile(const wetted_surface &fem, const std::string &name,
                     const bool &writeNodes) {
  // Open up mapper master file
  const auto fName = "loadMap/" + name + ".map";
  std::ofstream masterFile(fName, std::ios::out);
  if (masterFile.fail()) {
    cerr << "ERROR: Mapper master file " << fName << " did not open correctly!"
         << endl;
    exit(EXIT_FAILURE);
  }

  auto getNode = [&fem](const int &id) {
    vec3 node(0, 0, 0);
    node[0] = fem.node_coords[id - 1][0];
    node[1] = fem.node_coords[id - 1][1];
    node[2] = fem.node_coords[id - 1][2];
    return node;
  };

  auto areaTri = [&fem, &getNode](const std::vector<int> &conn) {
    auto r1 = getNode(conn[1]) - getNode(conn[0]);
    auto r2 = getNode(conn[2]) - getNode(conn[0]);
    return 0.5 * cross(r1, r2);
  };

  auto areaQuad = [&fem, &getNode](const std::vector<int> &conn) {
    auto r1 = getNode(conn[3]) - getNode(conn[1]);
    auto r2 = getNode(conn[2]) - getNode(conn[0]);
    return 0.5 * cross(r1, r2);
  };

  masterFile << std::scientific;
  masterFile.precision(8);
  if (writeNodes) {
    for (auto ii = 0U; ii < fem.node_ids.size(); ++ii) {
      masterFile << fem.node_ids[ii] << " " << GetFileName(name, ii, writeNodes)
                 << endl;
    }
  } else {  // write data at element centroids
    for (auto ii = 0U; ii < fem.face_conn.size(); ++ii) {
      masterFile << ii + 1 << " " << GetFileName(name, ii, writeNodes) << " ";
      const auto numNodes = fem.face_conn[ii].size();
      if (IsFaceValid(fem, ii)) {
        // calculate centroid
        vec3 centroid(0, 0, 0);
        for (const auto &id : fem.face_conn[ii]) {
          centroid += getNode(id);
        }
        centroid /= fem.face_conn[ii].size();
        masterFile << centroid;
        // calculate element normal
        vec3 area(0, 0, 0);
        if (numNodes == 3 || numNodes == 6) {
          area = areaTri(fem.face_conn[ii]);
        } else {
          area = areaQuad(fem.face_conn[ii]);
        }
        masterFile << area / norm(area) << endl;
      } else {
        cerr << "WARNING: Found FEM element with " << numNodes
             << " nodes - not adding to " << fName << endl;
      }
    }
  }
  masterFile.close();
}

std::string GetFileName(const std::string &rootName, const int &ii,
                        const bool &mapNodes) {
  constexpr auto width = 7;
  auto eid = std::to_string(ii + 1);
  const auto type = mapNodes ? "_node_" : "_elem_";
  auto fileName = rootName + type + std::string(width - eid.length(), '0') +
                  eid + ".csv";
  return fileName;
}

vec3 GetForce(const wetted_surface &fem, const int &id) {
  vec3 force(0, 0, 0);
  force[0] = fem.force[id - 1][0];
  force[1] = fem.force[id - 1][1];
  force[2] = fem.force[id - 1][2];
  return force;
}

vec3 GetMoment(const wetted_surface &fem, const int &id) {
  vec3 moment(0, 0, 0);
  moment[0] = fem.force[id - 1][3];
  moment[1] = fem.force[id - 1][4];
  moment[2] = fem.force[id - 1][5];
  return moment;
}

vec3 GetNode(const wetted_surface &fem, const int &id) {
  vec3 node(0, 0, 0);
  node[0] = fem.node_coords[id - 1][0];
  node[1] = fem.node_coords[id - 1][1];
  node[2] = fem.node_coords[id - 1][2];
  return node;
}

std::vector<vec3> GetNodeMoment(const wetted_surface &fem) {
  std::vector<vec3> nodeMoments;
  nodeMoments.reserve(fem.node_ids.size());
  for (auto id = 1U; id <= fem.node_ids.size(); ++id) {
    nodeMoments.emplace_back(GetMoment(fem, id));
  }
  return nodeMoments;
}

std::vector<vec3> GetNodeForce(const wetted_surface &fem) {
  std::vector<vec3> nodeForces;
  nodeForces.reserve(fem.node_ids.size());
  for (auto id = 1U; id <= fem.node_ids.size(); ++id) {
    nodeForces.emplace_back(GetForce(fem, id));
  }
  return nodeForces;
}

std::vector<vec3> GetFaceForce(const wetted_surface &fem,
                               const std::vector<int> &adjacentFaces) {
  std::vector<vec3> faceForces;
  faceForces.reserve(fem.face_conn.size());
  for (const auto &conn : fem.face_conn) {
    faceForces.emplace_back(0, 0, 0);
    for (const auto &id : conn) {
      faceForces.back() += GetForce(fem, id) / adjacentFaces[id - 1];
    }
  }
  return faceForces;
}

std::vector<int> GetAdjacentFaces(const wetted_surface &fem) {
  std::vector<int> adjacentFaces(fem.node_ids.size(), 0);
  for (const auto &conn : fem.face_conn) {
    for (const auto &id : conn) {
      adjacentFaces[id - 1]++;
    }
  }
  return adjacentFaces;
}

bool IsFaceValid(const wetted_surface &fem, const int &fid) {
  const auto numNodes = fem.face_conn[fid].size();
  const auto isParabolic = fem.isParabolic[fid];
  auto isValid = false;
  if (isParabolic) {
    isValid = numNodes == 6 || numNodes == 8;
  } else {
    isValid = numNodes == 3 || numNodes == 4;
  }
  return isValid;
}

void PrintStats(const wetted_surface &fem, const wetted_surface &cfd,
                const std::vector<vec3> &femFace, const vec3 &momentCenter) {
  const auto cfdForceNode = SumForces(cfd);
  const auto cfdMomentNode = SumMoments(cfd, momentCenter);
  const auto femForceNode = SumForces(fem);
  const auto femMomentNode = SumMoments(fem, momentCenter);
  const auto femForceFace =
      std::accumulate(std::begin(femFace), std::end(femFace), vec3(0, 0, 0));
  const auto femMomentFace = SumMomentsFace(fem, femFace, momentCenter);

  cout << "Sum of forces and moments at CFD nodes: Fx=" << cfdForceNode.X()
       << ", Fy=" << cfdForceNode.Y() << ", Fz=" << cfdForceNode.Z()
       << ", Mx=" << cfdMomentNode.X() << ", My=" << cfdMomentNode.Y()
       << ", Mz=" << cfdMomentNode.Z() << endl;
  cout << "Sum of forces and moments at FEM nodes: Fx=" << femForceNode.X()
       << ", Fy=" << femForceNode.Y() << ", Fz=" << femForceNode.Z()
       << ", Mx=" << femMomentNode.X() << ", My=" << femMomentNode.Y()
       << ", Mz=" << femMomentNode.Z() << endl;
  cout << "Sum of forces and moments at FEM faces: Fx=" << femForceFace.X()
       << ", Fy=" << femForceFace.Y() << ", Fz=" << femForceFace.Z()
       << ", Mx=" << femMomentFace.X() << ", My=" << femMomentFace.Y()
       << ", Mz=" << femMomentFace.Z() << endl;
  constexpr auto eps = 1.0e-30;
  cout << "Error between FEM faces and CFD nodes is: Force="
       << (norm(femForceFace) - norm(cfdForceNode)) /
              (norm(cfdForceNode) + eps) * 100
       << "%, Moment="
       << (norm(femMomentFace) - norm(cfdMomentNode)) /
              (norm(cfdMomentNode) + eps)
       << "%" << endl;
}

void PrintStatsNodes(const wetted_surface &fem, const wetted_surface &cfd,
                     const vec3 &momentCenter) {
  const auto cfdForceNode = SumForces(cfd);
  const auto cfdMomentNode = SumMoments(cfd, momentCenter);
  const auto femForceNode = SumForces(fem);
  const auto femMomentNode = SumMoments(fem, momentCenter);

  cout << "Sum of forces and moments at CFD nodes: Fx=" << cfdForceNode.X()
       << ", Fy=" << cfdForceNode.Y() << ", Fz=" << cfdForceNode.Z()
       << ", Mx=" << cfdMomentNode.X() << ", My=" << cfdMomentNode.Y()
       << ", Mz=" << cfdMomentNode.Z() << endl;
  cout << "Sum of forces and moments at FEM nodes: Fx=" << femForceNode.X()
       << ", Fy=" << femForceNode.Y() << ", Fz=" << femForceNode.Z()
       << ", Mx=" << femMomentNode.X() << ", My=" << femMomentNode.Y()
       << ", Mz=" << femMomentNode.Z() << endl;
  constexpr auto eps = 1.0e-30;
  cout << "Error between FEM nodes and CFD nodes is: Force="
       << (norm(femForceNode) - norm(cfdForceNode)) /
              (norm(cfdForceNode) + eps) * 100
       << "%, Moment="
       << (norm(femMomentNode) - norm(cfdMomentNode)) /
              (norm(cfdMomentNode) + eps)
       << "%" << endl;
}



vec3 SumForces(const wetted_surface &surf) {
  vec3 sum(0, 0, 0);
  for (auto ii = 0U; ii < surf.force.size(); ++ii) {
    sum += vec3(surf.force[ii][0], surf.force[ii][1], surf.force[ii][2]);
  }
  return sum;
}

vec3 SumMoments(const wetted_surface &surf, const vec3 &center) {
  vec3 sum(0, 0, 0);
  for (auto ii = 0U; ii < surf.force.size(); ++ii) {
    const auto force =
        vec3(surf.force[ii][0], surf.force[ii][1], surf.force[ii][2]);
    const auto r = vec3(surf.node_coords[ii][0], surf.node_coords[ii][1],
                        surf.node_coords[ii][2]) -
                   center;
    sum -= cross(r, force);
    const auto moment =
        vec3(surf.force[ii][3], surf.force[ii][4], surf.force[ii][5]);
    sum += moment;
  }
  return sum;
}

vec3 SumMomentsFace(const wetted_surface &surf, const std::vector<vec3> &force,
                    const vec3 &center) {
  vec3 sum(0, 0, 0);
  for (auto ii = 0U; ii < surf.face_conn.size(); ++ii) {
    vec3 centroid(0, 0, 0);
    for (const auto &id : surf.face_conn[ii]) {
      centroid += GetNode(surf, id);
    }
    centroid /= surf.face_conn[ii].size();
    const auto r = centroid - center;
    sum -= cross(r, force[ii]);
  }
  return sum;
}

vec3 GetVecFromString(const std::string &str) {
  const auto strTrim = Trim(str);
  const auto tokens = Tokenize(strTrim, ",");
  return vec3(std::stod(tokens[0]), std::stod(tokens[1]), std::stod(tokens[2]));
}

// function to broadcast a string from ROOT to all processors. This is needed
// because it is not garunteed in the MPI standard that the commmand
// line arguments will be on any processor but ROOT
void BroadcastString(std::string &str) {
  // str -- string to broadcast to all processors

  // get size of string (+1 for c_str end character)
  auto strSize = static_cast<int>(str.size() + 1);
  MPI_Bcast(&strSize, 1, MPI_INT, 0,
            MPI_COMM_WORLD);  // broadcast string size

  // allocate a char buffer of string size
  auto buf = std::make_unique<char[]>(strSize);
  snprintf(buf.get(), strSize, "%s", str.c_str());  // copy string into buffer
  MPI_Bcast(buf.get(), strSize, MPI_CHAR, 0,
            MPI_COMM_WORLD);  // broadcast string as char

  // create new string and assign to old string
  std::string newStr(buf.get(), strSize - 1);  // -1 to not include c_str end character
  str = newStr;
}

// check to see if directory exists and create it if it doesn't
void CreateDirectory(const std::string &dname, const int &rank) {
  if (rank != 0) {
    return;
  }
  std::error_code ec;
  std::filesystem::create_directory(dname, ec);
  if (ec) {
    cerr << "ERROR: Could not create directory: " << dname << endl;
    cerr << "Error code: " << ec.value() << endl;
    cerr << "Error message: " << ec.message() << endl;
    exit(EXIT_FAILURE);
  }
  std::filesystem::permissions(dname, std::filesystem::perms::owner_all,
                               std::filesystem::perm_options::replace);
  std::filesystem::permissions(dname, std::filesystem::perms::group_read,
                               std::filesystem::perm_options::add);
  std::filesystem::permissions(dname, std::filesystem::perms::group_exec,
                               std::filesystem::perm_options::add);
}
