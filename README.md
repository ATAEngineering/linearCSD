# linearCSD
This project was developed by [ATA Engineering](http://www.ata-e.com) as 
a way to solve for the linear structural response of a finite element
model (FEM) under a long transient loading, typically from a compuational 
fluid dynamics (CFD) simulation. The project uses ATA's linearFSI library 
to provide the mapping and structural solving capabilities. Mass and 
stiffness matrices can be used from commercial structural solvers such as 
NASTRAN. The advantage of using this tool over directly using NASTRAN is 
that long duration loading time histories for large FEMs can approach 
hundreds of GBs in NASTRAN format. This can be problematic since NASTRAN 
reads this data into memory at the start of the simulation. This tool does 
not have the same limitation as the loading is read on a time step by time 
step basis.

# Dependencies
This module depends on the [linearFSI](https://github.com/ATAEngineering/linearFSI) 
library developed by ATA Engineering. The linearFSI module itself depends 
on PETSc and MPI. A C++ compiler supporting the C++17 standard and cmake
3.10 or newer is also required.

# Installation Instructions
First install the dependencies. For very large problems, it can be helpful 
to install PETSc with 64-bit indices. An example of configure options to 
build PETSc for use with linearFSI are shown below. On Cray systems, the 
additional arguments should be used to specify the C, C++, and FORTRAN 
compilers: **--with-cc**, **--with-CC**, and **--with-fc**.

```bash
./configure --prefix=/path/to/install --download-fblaslapack=yes --download-superlu_dist=yes --download-scalapack=yes --with-64-bit-indices --with-debugging=0  --with-pic=1 --with-cxx-dialect=C++11 --with-shared-libraries=1 --download-hdf5=yes --with-mpi-dir=/path/to/mpi --download-hdf5-configure-arguments=--enable-fortran=no
```

The installation follows the standard cmake process. It is recommend to 
build from outside the source directory. The example below shows how to 
build the code. When installing on Cray systems, the additional flag 
**-DCRAY** should be passed to cmake during the configuring process.

```bash
cmake -DCMAKE_INSTALL_PREFIX=/path/to/installation -DCMAKE_BUILD_TYPE=release -DLINEARFSI_LIB_DIR=/path/to/linearFSI/lib -DPETSc_DIR=/path/to/petsc -DMPI_DIR=/path/to/mpi /path/to/source
make
make install
```

# To Do
* Add examples
* Add unit tests
* Improve documentation