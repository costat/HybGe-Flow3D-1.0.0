#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

#include "hgfM.hpp"

using namespace boost::filesystem;

/* function that checks if a file is present in a directory.
   used to determine if the mesh has already been built
   for the problem in this directory */
bool
find_file( const path& ProblemPath, \
           const std::string& fileName, \
           path& MeshPath )
{
  if ( !exists( ProblemPath ) ) return false;
  directory_iterator end_itr;
  for ( directory_iterator itr( ProblemPath );
        itr != end_itr;
        ++itr ) {
    if ( is_directory(itr->status()) ) {
      if ( find_file(itr->path(), fileName, MeshPath ) ) return true;
    }
    else if ( itr->path().leaf() == fileName ) {
      MeshPath = itr->path();
      return true;
    }
  }
  return false;
}

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

  // grab mesh path if yes, set nx = 0 to flag mesh exists

}
