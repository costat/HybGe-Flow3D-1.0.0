#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <iterator>

#include <paralution.hpp>

#include "hgf.hpp"
#include "hgfMesh.hpp"
#include "hgfStokes.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

/* hgfStokesDrive solves the Stokes equations on a cartesian grid
   with an immersed boundary. The model is solved in 2d or 3d, either
   once for upscaling a constant conductivity, or 2/3 times for
   upscaling a conductivity tensor. */
void
hgfDrive( unsigned long *gridin, int size1, int ldi1, int ldi2, \
          int nx, int ny, int nz, \
          double length, double width, double height, int direction, \
          double visc, int nThreads, int prec, int numSims, int simNum, \
          double tolAbs, double tolRel, int maxIt )
{

  if (simNum == 1) init_paralution();

  switch ( direction )
  {
    case 3 : // Solve all 3 flow directions for upscaled tensor
    {
      std::string outNameX;
      std::string outNameY;
      std::string outNameZ;
      outNameX = "flowrunX.dat";
      outNameY = "flowrunY.dat";
      outNameZ = "flowrunZ.dat";
      double mesh_duration, stokes_duration, total_duration;
      double start, stokes_start;

      start = omp_get_wtime();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );

      // Mesh Timer
      mesh_duration = ( omp_get_wtime() - start );
      stokes_start = omp_get_wtime();

      switch ( Mesh.DIM )
      {
        case 3 :
        {
          // call Stokes' solver for each axis flow problem
          std::vector< double > solX;
          std::vector< double > solY;
          std::vector< double > solZ;
          solX.resize( Mesh.dofTotal );
          solY.resize( Mesh.dofTotal );
          solZ.resize( Mesh.dofTotal );
          StokesSolve( Mesh, visc, 0, solX, tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolve( Mesh, visc, 1, solY, tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolve( Mesh, visc, 2, solZ, tolAbs, tolRel, maxIt, nThreads, prec );

          // stokes timer
          stokes_duration = ( omp_get_wtime() - stokes_start );

          // Write solution to file
          writeSolutionTP ( Mesh, solX, outNameX );
          writeSolutionTP ( Mesh, solY, outNameY );
          writeSolutionTP ( Mesh, solZ, outNameZ );
          computeKTensorL ( Mesh, solX, solY, solZ );

          std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
          std::cout << "Stokes problems solved in " << stokes_duration << "seconds\n";
          std::cout  << stokes_duration << "seconds\n";
          // Total timers
          total_duration = ( omp_get_wtime() - start );
          std::cout << "Total time: " << total_duration << "seconds\n";
          break;
        }
        case 2 :
        {
          // call Stokes' solver for each axis flow problem
          std::vector< double > solX;
          std::vector< double > solY;
          std::vector< double > solZ;
          solX.resize( Mesh.dofTotal );
          solY.resize( Mesh.dofTotal );
          solZ.resize( 1 );
          StokesSolve( Mesh, visc, 0, solX, tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolve( Mesh, visc, 1, solY, tolAbs, tolRel, maxIt, nThreads, prec );

          // stokes timer
          stokes_duration = ( omp_get_wtime() - stokes_start );

          // Write solution to file
          writeSolutionTP ( Mesh, solX, outNameX );
          writeSolutionTP ( Mesh, solY, outNameY );
          computeKTensorL ( Mesh, solX, solY, solZ );

          std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
          std::cout << "Stokes problems solved in " << stokes_duration;
          std::cout  << "seconds\n";
          // Total timers
          total_duration = ( omp_get_wtime() - start );
          std::cout << "Total time: " << total_duration << "seconds\n";
          break;
        }
      }
      break;
    }
    default : // Solve a single flow direction, upscale constant conductivity
    {
      std::string outName;
      outName = "flowrun.dat";
      double mesh_duration, stokes_duration, total_duration;
      double start, stokes_start;

      start = omp_get_wtime();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );
      // Mesh Timer
      mesh_duration = ( omp_get_wtime() - start );
      stokes_start = omp_get_wtime();

      // call Stokes' solver
      std::vector< double > sol;
      sol.resize( Mesh.dofTotal );
      StokesSolve( Mesh, visc, direction, sol, tolAbs, tolRel, maxIt, nThreads, prec );

      // Linear solve timer
      stokes_duration = ( omp_get_wtime() - stokes_start );

      // Write solution to file
      writeSolutionTP ( Mesh, sol, outName );
      computeKConstantDrive ( Mesh, sol, direction );

      std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
      std::cout << "Stokes problems solved in " << stokes_duration << "seconds\n";
      // Total timers
      total_duration = ( omp_get_wtime() - start );
      std::cout << "Total time: " << total_duration << "seconds\n";
      break;
    }
  }

  if (simNum == numSims) stop_paralution();
}

