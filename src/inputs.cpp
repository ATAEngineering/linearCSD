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

#include <utility>
#include <map>
#include <set>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include "linearCSD.h"
#include "mappedCoupling.h"
#include "plaintext_output.h"
#include "fsimat_petsc.h"
#include "inputs.hpp"

using std::cout;
using std::cerr;
using std::endl;

// -------------------------------------------------------------------------
// constructor
inputs::inputs(const std::string &inFile, const std::string &resDir, 
               const int &rank) {
  if (rank == 0) {
    cout << "Parsing input file " << inFile << "..." << endl;
  }
  rank_ = rank;
  readRestart_ = resDir;

  if (this->IsRestart()) {
    this->ReadTimeRestartFile(this->ReadRestartDirPath() + "linearFSItime_restart");
  }

  // set any variables w/o inputs here
  coupOpts_.map_data.resize(1);
  // coupOpts_.map_data[0].name = "WS_DISP";
  // coupOpts_.map_data[0].type = fsi::MAP_FEM_DISP;
  coupOpts_.map_data[0].name = "WS_FORCE";
  coupOpts_.map_data[0].type = fsi::MAP_CFD_FORCE;
  //coupOpts_.outputOP2files = true;
  coupOpts_.map_options.splittingAlgorithm = SPLIT_EAR;

  std::ifstream in(inFile, std::ios::in);
  if (in.fail()) {
    cerr << "ERROR: Error in inputs::inputs(). Input file '"
	 << inFile << "' did not open correctly!" << endl;
    exit(EXIT_FAILURE);
  }

  auto openBraces = false;
  auto closeBraces = false;
  auto insideBraces = false;
  std::string line = "";
  while (getline(in, line)) {
    // remove leading and trailing whitespace, ignore comments
    line = Trim(line);

    // determine if we are inside {}
    if (line.find("{") != std::string::npos) {
      if (openBraces) {
	cerr << "ERROR: found second '{' in input file" << endl;
        exit(EXIT_FAILURE);
      }
      openBraces = true;
    }
    if (line.find("}") != std::string::npos) {
      if (closeBraces) {
	cerr << "ERROR: found second '}' in input file" << endl;
        exit(EXIT_FAILURE);
      }
      closeBraces = true;
    }
    insideBraces = openBraces && !closeBraces && 
                   line.find("{") == std::string::npos &&
                   line.find("}") == std::string::npos;

    // only proceed if line has data and inside {}
    if (line.size() > 0 && insideBraces) {
      // tokenize line
      const auto tokens = Tokenize(line, ":", 1);
      const auto key = tokens[0];
      const auto value = tokens[1];
      // see if key matches a valid input variable
      if (key == "fsiFiles") {
        auto options = ReadOptionsList(in, line);
        if (options.find("fsiFile") != options.end()) {
	  femFsiFile_ = options["fsiFile"];
	} else {
	  cerr << "ERROR: fem fsiFile not specified" << endl;
          exit(EXIT_FAILURE);
	}

        if (options.find("massMatrix") != options.end()) {
	  csdOpts_.Mfile = options["massMatrix"];
	} else {
	  cerr << "ERROR: fem massMatrix not specified" << endl;
          exit(EXIT_FAILURE);
	}

        if (options.find("stiffnessMatrix") != options.end()) {
	  csdOpts_.Kfile = options["stiffnessMatrix"];
	} else {
	  cerr << "ERROR: fem stiffnessMatrix not specified" << endl;
          exit(EXIT_FAILURE);
	}

        if (options.find("dampingMatrix") != options.end()) {
	  csdOpts_.Cfile = options["dampingMatrix"];
	}

        if (options.find("dispOTM") != options.end()) {
	  csdOpts_.DispOTMfile = options["dispOTM"];
	} else {
	  cerr << "ERROR: fem dispOTM not specified" << endl;
          exit(EXIT_FAILURE);
	}

        if (options.find("eforceOTM") != options.end()) {
	  csdOpts_.EForceOTMfile = options["eforceOTM"];
	}

        if (options.find("stressOTM") != options.end()) {
	  csdOpts_.StressOTMfile = options["stressOTM"];
	}

        if (options.find("ampFile") != options.end()) {
	  csdOpts_.Ampfile = options["ampFile"];
	}

        if (options.find("spcFile") != options.end()) {
	  csdOpts_.SPCfile = options["spcFile"];
	}

        if (options.find("icFile") != options.end()) {
	  csdOpts_.ICfile = options["icFile"];
	}

      } else if (key == "cfdFsiFile") {
	cfdFsiFile_ = value;
      } else if (key == "fsiModalDamping") {
        auto options = ReadOptionsList(in, line);
        for (const auto &opts : options) {
	  csdOpts_.modalDampingSchedule.push_back(
            std::make_pair(std::stod(opts.first), std::stod(opts.second)));
	}
      } else if (key == "fsiRayleighDamping") {
        auto options = ReadOptionsList(in, line);
        auto specifyMScale = false;
        if (options.find("Mscale") != options.end() ||
            options.find("Kscale") != options.end()) {
	  specifyMScale = true;
	}
        auto calcMScale = false;
        if (options.find("freqA") != options.end() ||
            options.find("freqB") != options.end() ||
            options.find("dampA") != options.end() ||
            options.find("dampB") != options.end()) {
	  calcMScale = true;
	}
        if (specifyMScale && calcMScale) {
	  cerr << "ERROR: in 'fsiRayleighDamping' cannot specify 'Mscale' "
	       << "or 'Kscale' with 'freqA', 'freqB', dampA', and 'dampB'"
	       << endl;
	  exit(EXIT_FAILURE);
	} else if (!specifyMScale && !calcMScale) {
	  cerr << "ERROR: in 'fsiRayleighDamping' must specify 'Mscale', "
               << "'Kscale' or 'freqA', 'dampA', etc" << endl;
	  exit(EXIT_FAILURE);
	}

        if (specifyMScale) {
          if (options.find("Mscale") != options.end()) {
    	    csdOpts_.rayleighParams[0] = std::stod(options["Mscale"]);
	  }
          if (options.find("Kscale") != options.end()) {
    	    csdOpts_.rayleighParams[1] = std::stod(options["Kscale"]);
	  }
	}

        if (calcMScale) {
	  auto freqA = 0.0;
          auto freqB = 0.0;
          auto dampA = 0.0;
          auto dampB = 0.0;
          if (options.find("freqA") != options.end()) {
    	    freqA = std::stod(options["freqA"]);
	  } else {
	    cerr << "ERROR: in 'fsiRayleighDamping', must specify 'freqA'"
  	         << endl;
            exit(EXIT_FAILURE);
	  }
          if (options.find("freqB") != options.end()) {
    	    freqB = std::stod(options["freqB"]);
	  } else {
	    cerr << "ERROR: in 'fsiRayleighDamping', must specify 'freqB'"
  	         << endl;
            exit(EXIT_FAILURE);
	  }
          if (options.find("dampA") != options.end()) {
    	    dampA = std::stod(options["dampA"]);
	  } else {
	    cerr << "ERROR: in 'fsiRayleighDamping', must specify 'dampA'"
  	         << endl;
            exit(EXIT_FAILURE);
	  }
          if (options.find("dampB") != options.end()) {
    	    dampB = std::stod(options["dampB"]);
	  } else {
	    cerr << "ERROR: in 'fsiRayleighDamping', must specify 'dampB'"
  	         << endl;
            exit(EXIT_FAILURE);
	  }
          auto a = freqA * 2.0 * M_PI;
          auto b = freqB * 2.0 * M_PI;
          auto det = 1.0 / (b/a - a/b);
	  csdOpts_.rayleighParams[0] = 2.0 * det * (b * dampA - a * dampB);
          csdOpts_.rayleighParams[1] = 2.0 * det * (-dampA / b + dampB / a);
	}
      } else if (key == "fsiIntegrationScheme") {
        if (value == "newmarkBeta") {
          csdOpts_.integrationType = fsi::INT_NEWMARK;
	} else if (value == "hht") {
	  csdOpts_.integrationType = fsi::INT_HHT;
	} else if (value == "modal") {
	  csdOpts_.integrationType = fsi::INT_MODAL;
	} else if (value == "sol109") {
	  csdOpts_.integrationType = fsi::INT_SOL109;
	} else {
          cerr << "ERROR: fsiIntegrationScheme '" << value
	       << "' is not recognized" << endl;
          exit(EXIT_FAILURE);
	}
      } else if (key == "fsiPrintEnergySum") {
	printEnergySum_ = std::stoi(value);
      } else if (key == "fsiScatterMatricesUsingFile") {
	csdOpts_.scatterMatricesUsingFile = std::stoi(value);
      } else if (key == "fsiComputeResidual") {
	csdOpts_.doResidCalc = std::stoi(value);
      } else if (key == "fsiElementForceProbe") {
        auto options = ReadOptionsList(in, line);
        for (const auto &opt : options) {
	  probes_.insert(std::make_pair(
              opt.first, std::make_pair(std::stoi(opt.second), 
                                        fsi::PROBE_EFORCE)));
	}
      } else if (key == "fsiNodeProbe") {
        auto options = ReadOptionsList(in, line);
        for (const auto &opt : options) {
	  probes_.insert(std::make_pair(
              opt.first, std::make_pair(std::stoi(opt.second), 
                                        fsi::PROBE_UVA)));
	}
      } else if (key == "fsiOutputDir") {
	outputDir_ = value;
      } else if (key == "fsiRestartDir") {
	restartDir_ = value;
      } else if (key == "fsiProbeDir") {
	probeDir_ = value;
      } else if (key == "fsiProbePrecision") {
	csdOpts_.probePrecision = std::stoi(value);
      } else if (key == "fsiOutput") {
        auto lst = ReadList(in, line);
        if (lst.find("disp") != lst.end()) {
	  csdOpts_.outputDisp = true;
	}
        if (lst.find("vel") != lst.end()) {
	  csdOpts_.outputVel = true;
	}
        if (lst.find("accel") != lst.end()) {
	  csdOpts_.outputAccel = true;
	}
        if (lst.find("stress") != lst.end()) {
	  csdOpts_.outputStress = true;
	}
        if (lst.find("strain") != lst.end()) {
	  csdOpts_.outputStrain = true;
	}
        if (lst.find("eforce") != lst.end()) {
	  csdOpts_.outputEForce = true;
	}
        if (lst.find("FSIsurf") != lst.end()) {
	  coupOpts_.outputCFDSurf = true;
          coupOpts_.outputFEMSurf = true;
	}
        if (lst.find("FSIsurf_proj") != lst.end()) {
	  coupOpts_.outputCFDSurf_proj = true;
	}
        if (lst.find("FSIproj_disp") != lst.end()) {
	  coupOpts_.outputProjectionOP2 = true;
	}
        if (lst.find("FSIsurf_map") != lst.end()) {
	  coupOpts_.outputCFDSurf_map = true;
          coupOpts_.outputFEMSurf_map = true;
	}
        if (lst.find("FSIforce") != lst.end()) {
	  coupOpts_.map_data[0].doFEMOutput = true;
	  coupOpts_.map_data[0].doCFDOutput = true;
	}
      } else if (key == "dtmax") {
	csdOpts_.timeStep = std::stod(value);
      } else if (key == "forceScaleFactor") {
	forceScaleFactor_ = std::stod(value);
      } else if (key == "lengthScaleFactor") {
	lengthScaleFactor_ = std::stod(value);
      } else if (key == "probe_freq") {
	probeFreq_ = std::stoi(value);
      } else if (key == "plot_freq") {
	plotFreq_ = std::stoi(value);
      } else if (key == "plot_modulo") {
	plotModulo_ = std::stoi(value);
      } else if (key == "restart_freq") {
	restartFreq_ = std::stoi(value);
      } else if (key == "restart_modulo") {
	restartModulo_ = std::stoi(value);
      } else if (key == "cfdStartIter") {
        if (!this->IsRestart()) {
  	  startIterCFD_ = std::stoi(value);
	}
      } else if (key == "cfdEndIter") {
	endIterCFD_ = std::stoi(value);
      } else if (key == "cfdIncrement") {
	incrementCFD_ = std::stoi(value);
      } else if (key == "fsiDebugLevel") {
	debugLevel_ = std::stoi(value);
      } else if (key == "fsiDenseMatrixFillFactor") {
        denseMatrixFillFactor_ = std::stod(value);
        setDenseMatrixFillFactor(denseMatrixFillFactor_);
      } else if (key == "fsiStaticInit") {
        csdOpts_.doStaticInit = std::stoi(value);
      } else if (key == "fsiOutputMatlabMatrices") {
        csdOpts_.outputMatlabMatrices = std::stoi(value);
      } else if (key == "fsiRestartHDF5") {
        csdOpts_.outputHDF5files = std::stoi(value);
      } else if (key == "fsiIgnoreBoundsCheck") {
        coupOpts_.ignoreBoundsCheck = std::stoi(value);
      } else if (key == "fsiBoundsCheckEpsilon") {
        coupOpts_.boundsCheckEpsilon = std::stod(value);
      } else if (key == "fsiBoundsCheckTolerance") {
        coupOpts_.boundsCheckTolerance = std::stod(value);
      } else if (key == "fsiProjectUnmappedNodes") {
        coupOpts_.map_options.projectNearNodeFlag = std::stoi(value);
      } else if (key == "fsiUsePolyShapeFunctions") {
        coupOpts_.map_options.usePolyShapeFunctions = std::stoi(value);
      } else if (key == "fsiProjectPolyEdges") {
        coupOpts_.map_options.projectPolyEdges = std::stoi(value);
      } else if (key == "fsiPolyEdgeProjectThreshold") {
        coupOpts_.map_options.polyEdgeProjectThreshold = std::stod(value);
      } else if (key == "fsiPolyShapeFunctionExponent") {
        coupOpts_.map_options.IDWShapeFunctionExponent = std::stod(value);
      } else if (key == "fsiMapVisibleFaceDeltaAngle") {
        coupOpts_.map_options.visibleFaceDeltaAngle = std::stod(value);
      } else if (key == "fsiFixInconsistentFEMNormals") {
        coupOpts_.map_options.autoFixInconsistentNormals = std::stoi(value);
      } else if (key == "solverPackage") {
        if (value != "mumps" && value != "pastix" &&
            value != "superlu" && value != "lu") {
	  cerr << "ERROR: solverPackage '" << value 
               << "' is not recognized." << endl;
          cerr << "Choose 'mumps', 'pastix', 'superlu', or 'lu'" << endl;
          exit(EXIT_FAILURE);
	}
	csdOpts_.solverPackage = value;
      } else {
	cerr << "ERROR: input variable '" << key << "' is not recognized"
             << endl;
	exit(EXIT_FAILURE);
      }
    }
  }
  in.close();
}

// -------------------------------------------------------------------------
// function definitions
std::string inputs::RestartModDir(const int &inc) const {
  return this->RestartDir() + "/" + std::to_string(inc % restartModulo_);
}

void inputs::EchoInputs() const {
  if (rank_ == 0) {
    cout << "Echoing inputs..." << endl;
    cout << "cfdFsiFile: " << cfdFsiFile_ << endl;
    cout << "cfdStartIter: " << startIterCFD_ << endl;
    cout << "cfdEndIter: " << endIterCFD_ << endl;
    cout << "cfdIncrement: " << incrementCFD_ << endl;
    cout << "femFsiFile: " << cfdFsiFile_ << endl;
    cout << "massMatrix: " << csdOpts_.Mfile << endl;
    cout << "stiffnessMatrix: " << csdOpts_.Kfile << endl;
    cout << "dampingMatrix: " << csdOpts_.Cfile << endl;
    cout << "dispOTM: " << csdOpts_.DispOTMfile << endl;
    cout << "eforceOTM: " << csdOpts_.EForceOTMfile << endl;
    cout << "stressOTM: " << csdOpts_.StressOTMfile << endl;
    cout << "ampFile: " << csdOpts_.Ampfile << endl;
    cout << "spcFile: " << csdOpts_.SPCfile << endl;
    cout << "icFile: " << csdOpts_.ICfile << endl;
    if (this->UsingRayleighDamping()) {
      cout << "fsiRayleighDamping: <Mscale=" << csdOpts_.rayleighParams[0]
           << ", Kscale=" << csdOpts_.rayleighParams[0] << ">" << endl;
    }
    if (this->UsingModalDamping()) {
      cout << "fsiModalDamping: <";
      csdOpts_.modalDampingSchedule.Print(cout);
      cout << ">" << endl;
    }
    cout << "solverPackage: " << csdOpts_.solverPackage << endl;
    cout << "dtmax: " << csdOpts_.timeStep << endl;
    cout << "plot_freq: " << plotFreq_ << endl;
    cout << "plot_modulo: " << plotModulo_ << endl;
    cout << "fsiOutputDir: " << outputDir_ << endl;
    cout << "restart_freq: " << restartFreq_ << endl;
    cout << "restart_modulo: " << restartModulo_ << endl;
    cout << "fsiDebugLevel: " << debugLevel_ << endl;
    cout << "fsiRestartDir: " << restartDir_ << endl;
    cout << "fsiProbeDir: " << probeDir_ << endl;
    cout << "probe_freq: " << probeFreq_ << endl;
    cout << "forceScaleFactor: " << forceScaleFactor_ << endl;
    cout << "lengthScaleFactor: " << lengthScaleFactor_ << endl;
    cout << "forceScaleFactor: " << forceScaleFactor_ << endl;
    cout << "fsiScatterMatricesUsingFile: " << csdOpts_.scatterMatricesUsingFile << endl;
    cout << "fsiOutputMatlabMatrices: " << csdOpts_.outputMatlabMatrices << endl;
    cout << "fsiDenseMatrixFillFactor: " << denseMatrixFillFactor_ << endl;
    cout << "fsiStaticInit: " << csdOpts_.doStaticInit << endl;
    cout << "fsiRestartHDF5: " << csdOpts_.outputHDF5files << endl;
    cout << "fsiIgnoreBoundsCheck: " << coupOpts_.ignoreBoundsCheck << endl;
    cout << "fsiBoundsCheckEpsilon: " << coupOpts_.boundsCheckEpsilon << endl;
    cout << "fsiBoundsCheckTolerance: " << coupOpts_.boundsCheckTolerance << endl;
    cout << "fsiProjectUnmappedNodes: " << coupOpts_.map_options.projectNearNodeFlag << endl;
    cout << "fsiUsePolyShapeFunctions: " << coupOpts_.map_options.usePolyShapeFunctions << endl;
    cout << "fsiProjectPolyEdges: " << coupOpts_.map_options.projectPolyEdges << endl;
    cout << "fsiPolyEdgeProjectThreshold: " << coupOpts_.map_options.polyEdgeProjectThreshold << endl;
    cout << "fsiPolyShapeFunctionExponent: " << coupOpts_.map_options.IDWShapeFunctionExponent << endl;
    cout << "fsiMapVisibleFaceDeltaAngle: " << coupOpts_.map_options.visibleFaceDeltaAngle << " degrees" << endl;
    cout << "fsiFixInconsistentFEMNormals: " << coupOpts_.map_options.autoFixInconsistentNormals << endl;
    cout << "probes: <" << endl;
    for (const auto &p : probes_) {
      cout << "name=" << p.first << ", id=" << p.second.first 
           << ", type=" << p.second.second << endl;
    }
    cout << endl;
  }
}

void inputs::AddProbes(const std::unique_ptr<fsi::linearCSD> &csd) const {
  for (const auto &p : probes_) {
    csd->addProbeFile(p.first, p.second.first, p.second.second);
  }
}

std::string Trim(const std::string &s, const std::string &whitespace) {
  const std::string comment = "//";
  if (s.empty()) {
    return "";
  } else {
    // find index of first non whitespace character
    const auto sBegin = s.find_first_not_of(whitespace);
    if (sBegin == std::string::npos) {
      return "";  // string is only whitespace
    }

    // find index of last non whitespace character
    const auto sEnd = s.find_last_not_of(whitespace);
    const auto sRange = sEnd - sBegin + 1;  // range to trim to
    auto temp = s.substr(sBegin, sRange);

    // find index of first comment
    auto tempComment = temp.find(comment);
    return temp.substr(0, tempComment);
  }
}

std::vector<std::string> Tokenize(std::string line, 
                                  const std::string &delimiter,
                                  const unsigned int maxSplits) {
  std::vector<std::string> tokens;
  auto reachedMax = false;
  auto pos = line.find(delimiter);
  while (pos != std::string::npos && !reachedMax) {
    auto token = Trim(line.substr(0, pos));
    if (!token.empty()) {tokens.push_back(token);}
    // treat consecutive delimiters as single delimiter
    auto end = line.find_first_not_of(delimiter, pos);
    line.erase(0, end);
    if (maxSplits > 0 && maxSplits == tokens.size()) {
      reachedMax = true;
    }
    pos = line.find(delimiter);
  }
  // add in remainder if not empty
  auto token = Trim(line);
  if (!token.empty()) {tokens.push_back(token);}
  return tokens;
}

// function to remove delimiter if it is last character
std::string RemoveTrailing(const std::string &str, 
                           const std::string &delimiter) {
  auto pos = str.rfind(delimiter);
  return (pos == str.length() - 1) ? str.substr(0, pos) : str;
}

// function to read an options list of the form below
// var = <key=value, key=value,
//        key=value, ...>
std::map<std::string, std::string> ReadOptionsList(std::ifstream &inFile,
                                                   std::string &str) {
  std::map<std::string, std::string> options;
  auto openList = false;
  do {
    const auto start = openList ? 0 : str.find("<");
    const auto listOpened = str.find("<") == std::string::npos ? false : true;
    const auto end = str.find(">");
    openList = (end == std::string::npos) ? true : false;

    // get list info on current line
    std::string list = "";
    if (listOpened && openList) {  // list opened on current line, stays open
      list = str.substr(start + 1, std::string::npos);
    } else if (listOpened && !openList) {  // list opened/closed on current line
      const auto range = end - start -1;
      list = str.substr(start + 1, range);  // +/- 1 to ignore <>
    } else if (!listOpened && openList) {  // list was open and remains open
      list = str.substr(start, std::string::npos);
    } else { // list was open and is now closed
      const auto range = end - start;
      list = str.substr(start, range);
    }
    list = RemoveTrailing(list, ",");

    // tokenize
    const auto tokens = Tokenize(list, ",");
    for (const auto &token : tokens) {
      const auto opt = Tokenize(token, "=");
      if (opt.size() != 2) {
	cerr << "ERROR: Problem reading options list on line '" << str 
             << "'. Key and value should be on same line, separated by '='" 
             << endl;
        exit(EXIT_FAILURE);
      }
      options.insert(std::make_pair(opt[0], opt[1]));
    }

    if (openList) {  // if list is open, get next line
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);
  return options;
}

std::set<std::string> ReadList(std::ifstream &inFile,
                               std::string &str) {
  std::set<std::string> options;
  auto openList = false;
  do {
    const auto start = openList ? 0 : str.find("<");
    const auto listOpened = str.find("<") == std::string::npos ? false : true;
    const auto end = str.find(">");
    openList = (end == std::string::npos) ? true : false;

    // get list info on current line
    std::string list = "";
    if (listOpened && openList) {  // list opened on current line, stays open
      list = str.substr(start + 1, std::string::npos);
    } else if (listOpened && !openList) {  // list opened/closed on current line
      const auto range = end - start -1;
      list = str.substr(start + 1, range);  // +/- 1 to ignore <>
    } else if (!listOpened && openList) {  // list was open and remains open
      list = str.substr(start, std::string::npos);
    } else { // list was open and is now closed
      const auto range = end - start;
      list = str.substr(start, range);
    }
    list = RemoveTrailing(list, ",");

    // tokenize
    const auto tokens = Tokenize(list, ",");
    for (const auto &token : tokens) {
      options.insert(token);
    }

    if (openList) {  // if list is open, get next line
      getline(inFile, str);
      str = Trim(str);
    }
  } while (openList);
  return options;
}

std::string inputs::CfdRootName() const {
  // find index of first comment
  auto idx = cfdFsiFile_.find(".");
  if (idx == std::string::npos) {
    cerr << "ERROR: cannot find root for file name '"
	 << cfdFsiFile_ << "'" << endl;
    exit(EXIT_FAILURE);
  }
  return cfdFsiFile_.substr(0, idx);
}

std::string inputs::FemRootName() const {
  // find index of first comment
  auto idx = femFsiFile_.find(".");
  if (idx == std::string::npos) {
    cerr << "ERROR: cannot find root for file name '"
	 << femFsiFile_ << "'" << endl;
    exit(EXIT_FAILURE);
  }
  return femFsiFile_.substr(0, idx);
}

void inputs::ReadTimeRestartFile(const std::string &fpath) {
  std::ifstream tf(fpath, std::ios::in);
  if (tf.fail()) {
    cerr << "ERROR: Error in inputs::ReadTimeRestartFile()." << endl;
    cerr << "Time restart file '" << fpath << "' did not open correctly!" 
         << endl;
    exit(EXIT_FAILURE);
  }

  tf >> startIterCFD_;
  tf >> csdOpts_.startIter;
  tf >> csdOpts_.startTime;
  tf.close();
}

void inputs::WriteTimeRestartFile(const int &csdIter, const int &cfdIter) const {
  const auto fpath = this->RestartModDir(csdIter) + "/linearFSItime_restart";
  std::ofstream tf(fpath, std::ios::out);
  if (tf.fail()) {
    cerr << "ERROR: Error in inputs::WriteTimeRestartFile()." << endl;
    cerr << "Time restart file '" << fpath << "' did not open correctly!" 
         << endl;
    exit(EXIT_FAILURE);
  }

  tf << cfdIter << " ";
  tf << csdIter << " ";
  tf << std::scientific << std::setprecision(8) << this->SimulationTime(csdIter);
  tf.close();
}
