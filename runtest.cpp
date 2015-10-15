/* This file is for testing hgf development without
   swigging a Python script. This allows for easy use
   of tools like Valgrind, or for easier benchmarking
   and bottleneck finding */
#include <cstdio>
#include <vector>

#include "hgf.hpp"
#include "hgfMesh.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

int
main( int argc, const char* argv[] )
{

  int nx, ny, nz, nThreads, prec, numSims, simNum, size1, ldi1, ldi2, direction;
  double length, width, height, visc;

  nx = 30;
  ny = 30;
  nz = 30;

  ldi1 = ny;
  ldi2 = nz;

  prec = 1;
  nThreads = 4;
  numSims = 1;
  simNum = 1;
  direction = 0;

  length = 1;
  width = 1;
  height = 1;
  visc = 1;

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

  hgfStokesDrive( gridin, size1, ldi1, ldi2, nx, ny, nz, \
                 length, width, height, direction, visc, \
                 nThreads, prec, numSims, simNum );

}
