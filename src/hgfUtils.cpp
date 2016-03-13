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

  if (!Par.nx || !Par.ny) {
    std::cout << "\nSomething went wrong readying mesh. Make sure formatting is correct.\n";
    exit(0);
  }

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
      if (ifs.good()) {
        std::getline(ifs, line);
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
problemParameters( ProbParam& Par, const bfs::path& ParameterPath )
{
  std::string line;
  std::string str;
  bfs::ifstream ifs( ParameterPath );
  // grab length
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ilength(line);
  ilength >> str >> Par.length;
  // grab width
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iwidth(line);
  iwidth >> str >> Par.width;
  // grab height
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iheight(line);
  iheight >> str >> Par.height;
  // grab viscocity
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ivisc(line);
  ivisc >> str >> Par.visc;
  // grab direction
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream idirec(line);
  idirec >> str >> Par.direction;
  // grab output switch
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iout(line);
  iout >> str >> Par.output;
  // grab nThreads
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream inth(line);
  inth >> str >> Par.nThreads;
  // grab prec
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream iprec(line);
  iprec >> str >> Par.prec;
  // grab tolAbs
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream ita(line);
  ita >> str >> Par.tolAbs;
  // grab tolRel
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream itr(line);
  itr >> str >> Par.tolRel;
  // grab maxIt
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream imaxit(line);
  imaxit >> str >> Par.maxIt;
  // grab nCuts
  if (ifs.good())
  {
    std::getline(ifs, line);
  }
  std::istringstream inc(line);
  inc >> str >> Par.nCuts;
}
void printParams( ProbParam& Par )
{
  std::cout << "\n" << "Length= " << Par.length << "\n";
  std::cout << "Width= " << Par.width << "\n";
  std::cout << "Height= " << Par.height << "\n";
  std::cout << "Viscosity= " << Par.visc << "\n";
  std::cout << "Direction= " << Par.direction << "\n";
  std::cout << "Output= " << Par.output << "\n";
  std::cout << "OMP Threads= " << Par.nThreads << "\n";
  std::cout << "ILU P-Level= " << Par.prec << "\n";
  std::cout << "Absolute Tolernace= " << Par.tolAbs << "\n";
  std::cout << "Relative Tolerance= " << Par.tolRel << "\n";
  std::cout << "Maximum Iterations= " << Par.maxIt << "\n";
  std::cout << "Cuts for Subdomains= " << Par.nCuts << "\n\n";
}
