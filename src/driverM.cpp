#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>

#include "hgfM.hpp"
#include "hgfUtils.hpp"
#include "hgfAuxTools.hpp"

namespace bfs = boost::filesystem;

// main driver
int
main( int argc, const char* argv[] )
{
  // struct containing problem parameters
  ProbParam Par;
  // path to problem directory from input
  bfs::path ProblemPath( argv[1] );
  // if it exists, path to Mesh.dat within ProblemPath. if it does not exist: ProblemPath/Mesh.dat.
  bfs::path MeshPath;

  // check if mesh already exists in the problem directory
  std::string meshSaved = "Mesh.dat";
  bool isMesh = find_file( ProblemPath, meshSaved, MeshPath );
  // if the mesh wasn't built, we set the path for the mesh to ProblemPath/Mesh.dat
  if (!isMesh) MeshPath = ProblemPath / "Mesh.dat";

  // whether or not the full mesh was built, we load the data from the geometry file
  bfs::path Geo( ProblemPath / "Geometry.dat" );
  geoFromFile( Par, Geo );
  // load problem parameters
  problemParameters( Par, ProblemPath );

  // send problem to hgf
  //hgfDrive( Par, ProblemPath, Geo );
}
