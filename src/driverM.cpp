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
  std::cout << "\nSolving problem in directory: \t " << ProblemPath << "\n";

  // check if mesh already exists in the problem directory
  std::string meshSaved = "Mesh.dat";
  bfs::path MeshPath;
  Par.isMesh = find_file( ProblemPath, meshSaved, MeshPath );
  // if the mesh wasn't built, we set the path for the mesh to ProblemPath/Mesh.dat
  if (!Par.isMesh) {
    std::cout << "\nMesh file not found. Mesh will be constructed from Geometry file.\n";
    MeshPath = ProblemPath / "Mesh.dat";
  }
  else std::cout << "\nMesh file found. For full domain simulations mesh will be loaded from file.\n";

  // find geometry file, exit if not present
  std::string geom = "Geometry.dat";
  bfs::path Geo;
  bool isGeo = find_file( ProblemPath, geom, Geo);
  if (!isGeo) {
    std::cout << "\nGeometry file not present in problem folder. Exiting\n";
    exit(0);
  }
  // grab geometry data
  geoFromFile( Par, Geo );

  // find parameter file, exit if not present
  std::string param = "Parameters.dat";
  bfs::path Parameters;
  bool isParam = find_file( ProblemPath, param, Parameters );
  if (!isParam) {
    std::cout << "\nParameter file not present in problem folder. Exiting\n";
    exit(0);
  }
  // grab parameter data
  problemParameters( Par, Parameters );
  printParams( Par );

  // send problem to hgf
  hgfDrive( ProblemPath, MeshPath, Par );
}
