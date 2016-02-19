/* This file is for testing hgf development without
   swigging a Python script. This allows for easy use
   of tools like Valgrind, or for easier benchmarking
   and bottleneck finding */
#include <cstdio>
#include <vector>
#include <cstdlib>

#include "hgf.hpp"

int
main( int argc, const char* argv[] )
{

  int nx, ny, nz, nThreads, prec, numSims, simNum, size1, ldi1, ldi2, direction, solver;
  double length, width, height, visc, relax;

  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  int MX = 4;
  int MY = 4;
  int MZ = 4;

  ldi1 = ny;
  ldi2 = nz;

  numSims = 1;
  simNum = 1;
  direction = atoi(argv[4]);
  nThreads = atoi(argv[5]);
  prec = atoi(argv[6]);

  solver = atoi(argv[7]);
  relax = atof(argv[8]);

  length = 1;
  width = 1;
  height = 1;
  visc = 1.0;

  if (nz)
  {
    size1 = nx * ny * nz;
  }
  else
  {
    size1 = nx * ny;
  }
  unsigned long *gridin = new unsigned long [size1];

  for (int cl = 0; cl < size1; cl++)
  {
    gridin[cl] = 0;
  }

  hgfDrive( gridin, size1, ldi1, ldi2, nx, ny, nz, \
                 length, width, height, direction, visc, \
                 nThreads, prec, numSims, simNum, 1e-12, 1e-12, 5000, MX, MY, MZ, solver, relax );

  delete[] gridin;

}