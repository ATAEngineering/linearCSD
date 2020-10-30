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

#ifndef INPUTS_HEADER_DEF
#define INPUTS_HEADER_DEF

#include <string>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include "linearCSD.h"
#include "mappedCoupling.h"
#include "plaintext_output.h"
#include "vec3.hpp"


class inputs {
  std::string cfdFsiFile_ = "undefined";
  std::string femFsiFile_ = "undefined";
  fsi::linearCSD_options csdOpts_;
  fsi::mappedCoupling_options coupOpts_;
  std::map<std::string, std::pair<int, fsi::probeType>> probes_;
  double forceScaleFactor_ = 1.0;
  double lengthScaleFactor_ = 1.0;
  double denseMatrixFillFactor_ = 1.0;
  int startIterCFD_ = 0;
  int endIterCFD_ = 0;
  int incrementCFD_ = 0;
  int numCFDSubIter_ = 1;
  int probeFreq_ = 1;
  int plotFreq_ = 1;
  int plotModulo_ = std::numeric_limits<int>::max();
  int restartFreq_ = 1;
  int restartModulo_ = std::numeric_limits<int>::max();
  int debugLevel_ = pow(2, 14) + 1;
  std::string restartDir_ = "./restart";
  std::string outputDir_ = "./output";
  std::string probeDir_ = "./probes";
  bool printEnergySum_ = true;
  int rank_ = 0;
  std::string readRestart_ = "";

 public:
  // constructor
  inputs(const std::string &inFile, const std::string &resDir, 
         const int &rank);

  // member functions
  const std::string &CfdFsiFile() const { return cfdFsiFile_; }
  const std::string &FemFsiFile() const { return femFsiFile_; }
  std::string CfdRootName() const;
  std::string FemRootName() const;
  bool DoPlot(const int &ii) const { return ii % plotFreq_ == 0; }
  bool DoProbe(const int &ii) const { return ii % probeFreq_ == 0; }
  bool DoRestart(const int &ii) const { 
    return ii % restartFreq_ == 0 && ii != (startIterCFD_ / incrementCFD_); 
  }
  const std::string &RestartDir() const { return restartDir_; }
  bool IsRestart() const { return !readRestart_.empty(); }
  std::string ReadRestartDirPath() const {
    return this->RestartDir() + "/" + readRestart_ + "/";
  }
  const std::string &ProbeDir() const { return probeDir_; }
  std::string RestartModDir(const int &inc) const;
  const std::string &OutputDir() const { return outputDir_; }
  const fsi::mappedCoupling_options &CouplingOptions() const { 
    return coupOpts_; 
  }
  const fsi::linearCSD_options &CsdOptions() const { return csdOpts_; }
  const bool &PrintEnergySum() const { return printEnergySum_; }
  const double &ForceScale() const { return forceScaleFactor_; }
  const double &LengthScale() const { return lengthScaleFactor_; }
  int NumIterations() const { 
    return (endIterCFD_ - startIterCFD_) / incrementCFD_; 
  }
  const int &NumCFDSubIter() const { return numCFDSubIter_; }
  const int &StartIterCFD() const { return startIterCFD_; }
  const int &EndIterCFD() const { return endIterCFD_; }
  const int &IncCFD() const { return incrementCFD_; }
  const int &DebugLevel() const { return debugLevel_; }
  double SimulationTime(const int &ii) const { 
    return csdOpts_.startTime + csdOpts_.timeStep * ii; 
  }
  void EchoInputs() const;
  bool UsingModalDamping() const { 
    return !csdOpts_.modalDampingSchedule.empty(); 
  }
  bool UsingRayleighDamping() const {
    return csdOpts_.rayleighParams[0] > 0.0 || 
           csdOpts_.rayleighParams[1] > 0.0;
  }
  void AddProbes(const std::unique_ptr<fsi::linearCSD> &csd) const;
  void ReadTimeRestartFile(const std::string &fpath);
  void WriteTimeRestartFile(const int &csdIter, const int &cfdIter) const;

  // destructor
  ~inputs(){}
};

// function definitions
std::string Trim(const std::string &s, const std::string & = " \t\r\n");
std::vector<std::string> Tokenize(std::string line, 
                                  const std::string &delimiter,
                                  const unsigned int = 0);
std::string RemoveTrailing(const std::string &str, 
                           const std::string &delimiter);
std::map<std::string, std::string> ReadOptionsList(std::ifstream &inFile,
                                                   std::string &str);
std::set<std::string> ReadList(std::ifstream &inFile,
                               std::string &str);

#endif
