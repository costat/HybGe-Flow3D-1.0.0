#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

#include "hgfM.hpp"
#include "hgfUtils.hpp"

using namespace boost::filesystem;

// main driver
int
main( int argc, const char* argv[] )
{

  int nx, ny, nz, nThreads, prec;
  int ldi1, ldi2, direction, solver, output;
  double length, width, height, visc, relax;

  path ProblemPath( argv[1] );
  path MeshPath;

  // check if mesh already exists
  std::string meshSaved = "Mesh.dat";
  bool isMesh = find_file( ProblemPath, meshSaved, MeshPath );
  if (isMesh) {
    nx = 0; ny = 0; nz = 0;
  }
  else {

  }

}
