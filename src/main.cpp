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
#include <memory>
#include <vector>
#include <numeric>
#include <stdexcept>
#include "linearCSD.h"
#include "mappedCoupling.h"
#include "wetted_surface.h"
#include "mpi.h"
#include "petsc.h"

#include "tcout.h"

#include "macros.hpp"
#include "utility.hpp"
#include "inputs.hpp"

using std::cout;
using std::cerr;
using std::endl;


void OutputData(const inputs &inp, const std::unique_ptr<fsi::linearCSD> &csd,
                const std::unique_ptr<fsi::mappedCoupling> &cosim,
                const int &rank, const int &csdIter, const int &cfdIter);

int main(int argc, char** argv) {
  // initialize MPI
  MPI_Init(&argc, &argv);
  PetscInitializeNoArguments();

  // get MPI info
  auto numProcs = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  auto rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto version = 0;
  auto subversion = 0;
  MPI_Get_version(&version, &subversion);

  // check command line arguments
  if (rank == 0) {
    const auto argv1 = argc > 2 ? std::string(argv[1]) : "";
    if (!(argc == 2 || argc == 3) || argv1 == "-h" || argv1 == "--help") {
      cout << "Usage: " << endl;
      cout << "mpirun -np <n> linearFSI <input> <restart>\n" << endl;
      cout << "<n> -- number of processors" << endl;
      cout << "<input> -- input file name" << endl;
      cout << "<restart> -- restart directory name (optional)" << endl;
      PetscFinalize();
      MPI_Finalize();
      return EXIT_FAILURE;
    }
  }

  if (rank == 0) {
    cout << "------------------------------------------------" << endl;
    cout << "               linearFSI stand alone            " << endl;
    cout << "------------------------------------------------" << endl;
    cout << "Version " << MAJORVERSION << "." << MINORVERSION << "."
         << PATCHVERSION << endl;
    cout << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
    cout << "------------------------------------------------" << endl;
    cout << endl << endl;
  }

  try {
    // get inputs
    auto inputFile = std::string(argv[1]);
    BroadcastString(inputFile);
    auto restartDir = argc == 3 ? std::string(argv[2]) : "";
    BroadcastString(restartDir);
    const auto inp = inputs(inputFile, restartDir, rank);
    inp.EchoInputs();

    // create needed directories
    CreateDirectory(inp.OutputDir(), rank);
    CreateDirectory(inp.RestartDir(), rank);
    CreateDirectory(inp.ProbeDir(), rank);

    // initialize solver
    tcout.setDebugLevel(inp.DebugLevel(), false);
    if (rank == 0) {
      cout << "MSG: Initializing CSD solver..." << endl;
    }
    auto csd = std::make_unique<fsi::linearCSD>(MPI_COMM_WORLD, inp.CsdOptions());
    // add probes
    if (rank == 0) {
      cout << "MSG: Adding probes to CSD solver..." << endl;
    }
    inp.AddProbes(csd);

    // read in fem mesh
    wetted_surface femSurf;
    if (rank == 0) {
      cout << "MSG: Reading FEM fsi file: " << inp.FemFsiFile() << endl;
      readFile(inp.FemFsiFile(), femSurf);
      std::vector<double> zero(6, 0.0);
      femSurf.force.resize(femSurf.node_ids.size(), zero);
      std::vector<int> femFaceIds(femSurf.face_conn.size(), 0);
      std::iota(std::begin(femFaceIds), std::end(femFaceIds), 1);
      writeNasSurface(inp.FemRootName() + "_orig.nas", femSurf, femFaceIds);
    }

    // read in cfd mesh
    wetted_surface cfdSurf;
    if (rank == 0) {
      cout << "MSG: Reading CFD fsi file: " << inp.CfdFsiFile() << endl;
      readFile(inp.CfdFsiFile(), cfdSurf);
      // scale cfd by length factor
      for (auto &coords : cfdSurf.node_coords) {
        coords[0] *= inp.LengthScale();
        coords[1] *= inp.LengthScale();
        coords[2] *= inp.LengthScale();
      }
      std::vector<double> zero(3, 0.0);
      cfdSurf.force.resize(cfdSurf.node_ids.size(), zero);
      std::vector<int> cfdFaceIds(cfdSurf.face_conn.size(), 0);
      std::iota(std::begin(cfdFaceIds), std::end(cfdFaceIds), 1);
      writeNasSurface(inp.CfdRootName() + "_orig.nas", cfdSurf, cfdFaceIds);
    }

    // create the cosim data type
    if (rank == 0) {
      cout << "MSG: Creating mapped coupling interface..." << endl;
    }
    auto cosim = std::make_unique<fsi::mappedCoupling>(MPI_COMM_WORLD, MPI_COMM_WORLD, 
                                                       femSurf, cfdSurf, 
                                                       inp.CouplingOptions());
    csd->addCoupledNodeSet_rank0("EXT_WSURF", femSurf.inode_ids, fsi::FORCE_COUP);

    csd->doSetup();

    // load restart data
    if (inp.IsRestart()) {
      if (rank == 0) {
	cout << "MSG: This simulation is a restart, loading restart files..." << endl;
      }
      csd->loadIntegratorRestartFile(inp.ReadRestartDirPath() + "linearFSI_restart");
      csd->loadWorkRestartFile(inp.ReadRestartDirPath() + "linearFSIwork_restart");
      cosim->loadCouplingRestartFile(inp.ReadRestartDirPath() + "linearFSIcoupling_restart");
    }

    csd->doFactor();
    if (rank == 0) {
      cout << "MSG: Output mapped coupling surface files..." << endl;
    }
    cosim->outputSurfaceFiles(inp.OutputDir());

    // solve fem at all times
    for (auto ii = inp.StartIterCFD(); ii <= inp.EndIterCFD(); ii += inp.IncCFD()) {
      auto csdIter = ii / inp.IncCFD();
      if (rank == 0) {
        cout << "MSG: Solving FEM at increment " << csdIter << "..." << endl;
        cout << "MSG: Do plot: " << inp.DoPlot(csdIter) << "; Do restart: " 
             << inp.DoRestart(csdIter) << "; Do probe: " << inp.DoProbe(csdIter) << endl;
      }

      // subiteration loop
      for (auto jj = 0; jj < inp.NumCFDSubIter(); ++jj) {
        // get CFD force data
        if (rank == 0) {
          GetForceData(ii, inp.CfdRootName(), inp.ForceScale(), cfdSurf);
        }
        cosim->setCFDvector_cfdComm(cfdSurf.force, "WS_FORCE");
        csd->setCoupledDOFForce_petscVec("EXT_WSURF", cosim->getFEMvector_petscVec("WS_FORCE"));
        // solve FEM at current step
        csd->doPredict();
  
        // cosim->setFEMvector_petscVec(csd->getCoupledDOFDisp_petscVec("EXT_WSURF"), "WS_DISP");
        // cfdSurf.disp = cosim->getFEMvector_petscVec("WS_DISP");
        // output wetted force sums
        cosim->printMapSum_rank0(std::cout, "WS_FORCE", ii, jj);
        //cosim->outputNewtonFiles(inp.OutputDir(), ii, jj);
      }

      // advance FEM to next time step
      csd->doAdvanceTime();
      // ouput data (restarts, probes, etc)
      OutputData(inp, csd, cosim, rank, csdIter, ii);
    }

    // need to release memory allocated by csd/cosim before MPI_Finalize
    csd.reset(nullptr);
    cosim.reset(nullptr);
    if (rank == 0) {
      cout << "MSG: Program Complete" << endl;
    }
  
  } catch (std::runtime_error &err) {
    cout << "Error running standalone linearFSI: " << err.what() << endl;
    exit(EXIT_FAILURE);
  }

  PetscFinalize();
  MPI_Finalize();
  return EXIT_SUCCESS;
}


void OutputData(const inputs &inp, const std::unique_ptr<fsi::linearCSD> &csd,
                const std::unique_ptr<fsi::mappedCoupling> &cosim,
                const int &rank, const int &csdIter, const int &cfdIter) {
  const auto restartModDir = inp.RestartModDir(csdIter);
  if (inp.DoRestart(csdIter)) {
    CreateDirectory(restartModDir, rank);
  }
  if (inp.PrintEnergySum()) {
    csd->outputEnergySum(inp.OutputDir());
  }
  if (inp.DoPlot(csdIter)) {
    csd->outputPlotFiles(inp.OutputDir());
    if (rank == 0) {
      cout << "MSG: cosim->outputPlotFiles()" << endl;
    }
    cosim->outputPlotFiles(inp.OutputDir(), csdIter, inp.SimulationTime(csdIter));
  }
  if (inp.DoProbe(csdIter)) {
    csd->outputProbeFiles(inp.ProbeDir());
  }

  if (csd->getIntegratorDifferencingType() == 0) {  // central differencing
    auto outputIter = csdIter - 1;
    if (inp.DoPlot(outputIter) && outputIter > 0) {
      csd->outputPlotFiles(inp.OutputDir());
    }
    if (inp.DoRestart(outputIter)) {
      csd->outputWorkRestartFile(restartModDir + "/linearFSIwork_restart");
    }
  } else {
    if (inp.DoPlot(csdIter)) {
      csd->outputPlotFiles(inp.OutputDir());
    }
    if (inp.DoRestart(csdIter)) {
      csd->outputWorkRestartFile(restartModDir + "/linearFSIwork_restart");
    }
  }

  if (inp.DoRestart(csdIter)) {
    inp.WriteTimeRestartFile(csdIter, cfdIter);
    csd->outputIntegratorRestartFile(restartModDir + "/linearFSI_restart");
    cosim->outputCouplingRestartFile(restartModDir + "/linearFSIcoupling_restart");
  }
}
