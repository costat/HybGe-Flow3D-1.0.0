#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include "hgfUtils.hpp"
#include "hgfAuxTools.hpp"

/* Define a 2d -> 1d array index, uses row major ordering */
#define idx2(i, j, ldi) ((i * ldi) + j)
/* Define a 3d -> 1d array index, uses row major ordering */
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

namespace bfs = boost::filesystem;

/* function that checks if a file is present in a directory.
   used to determine if the mesh has already been built
   for the problem in this directory */
bool
find_file( const bfs::path& ProblemPath, \
           const std::string& fileName, \
           bfs::path& FilePath )
{
  if ( !exists( ProblemPath ) ) return false;
  bfs::directory_iterator end_itr;
  for ( bfs::directory_iterator itr( ProblemPath );
        itr != end_itr;
        ++itr ) {
    if ( bfs::is_directory(itr->status()) ) {
      if ( find_file(itr->path(), fileName, FilePath ) ) return true;
    }
    else if ( itr->path().leaf() == fileName ) {
      FilePath = itr->path();
      return true;
    }
  }
  return false;
}

void
geoFromFile( ProbParam& Par, const bfs::path& Geo )
{
  // read first line to get nx, ny, nz
  std::string line;
  std::string str;
  bfs::ifstream ifs( Geo );
  // grab nx
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ix(line);
  ix >> str >> Par.nx;
  // grab ny
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iy(line);
  iy >> str >> Par.ny;
  // grab nz
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iz(line);
  iz >> str >> Par.nz;

  std::string hold;
  // read remaining lines of geometry into Par.gridin
  if ( Par.nz ) {
    for (int nslices = 0; nslices < Par.nz; nslices++) {
      for (int nrows = 0; nrows < Par.ny; nrows++) {
        if (ifs.good()) {
          std::getline(ifs,line);
        }
        // place line entries into gridin
        std::stringstream stream( line );
        while(1) {
          int n;
          stream >> n;
          if (!stream) break;
          Par.gridin.push_back( n );
        }
      }
    }
  }
  else {
    for (int nrows = 0; nrows < Par.ny; nrows++) {
      // grab row
      if (ifs.good()) {
        std::getline(ifs,line);
      }
      // place line entries into gridin
      std::stringstream stream( line );
      while(1) {
        int n;
        stream >> n;
        if (!stream) break;
        Par.gridin.push_back( n );
      }
    }
  }
}

void
problemParameters( ProbParam& Par, const bfs::path& ProblemPath )
{
  std::string line;
  std::string str;
  bfs::ifstream ifs( ProblemPath / "Paramters.dat" );

}
