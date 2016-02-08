#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <math.h>
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
          double tolAbs, double tolRel, int maxIt, \
          int MX, int MY, int MZ )
{

  int type;

  if (simNum == 1) init_paralution();

  if (direction < 3) type = 0;
  else if (direction == 3) type = 1;
  else if (direction == 4) type = 2;
  else
  {
    std::cout << "\nProblem type not properly defined. Solving x-flow problem.\n\n";
    type = 0;
    direction = 0;
  }

  switch ( type )
  {
    case 0 : // Solve a single flow direction, upscale constant conductivity
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
      StokesSolveDirect( Mesh, visc, direction, sol, tolAbs, tolRel, maxIt, nThreads, prec );

      // Linear solve timer
      stokes_duration = ( omp_get_wtime() - stokes_start );

      // Write solution to file
      writeSolutionTP ( Mesh, sol, outName );
      double K;
      computeKConstantDrive ( Mesh, sol, K, direction, 1 );

      std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
      std::cout << "Stokes problems solved in " << stokes_duration << "seconds\n";
      // Total timers
      total_duration = ( omp_get_wtime() - start );
      std::cout << "Total time: " << total_duration << "seconds\n";
      break;
    }
    case 1 : // Solve all 2/3 flow directions for upscaled tensor
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

      // call Stokes' solver for each axis flow problem
      std::vector< double > solX;
      std::vector< double > solY;
      std::vector< double > solZ;
      solX.resize( Mesh.dofTotal );
      solY.resize( Mesh.dofTotal );
      StokesSolveDirect( Mesh, visc, 0, solX, tolAbs, tolRel, maxIt, nThreads, prec );
      StokesSolveDirect( Mesh, visc, 1, solY, tolAbs, tolRel, maxIt, nThreads, prec );
      if ( Mesh.DIM == 3) {
        solZ.resize( Mesh.dofTotal );
        StokesSolveDirect( Mesh, visc, 2, solZ, tolAbs, tolRel, maxIt, nThreads, prec );
      }
      else solZ.resize( 1 );

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
    { // compute conductivities on subdomains, then solve constructed pore-network problem

      // determine array slices for subdomains
      std::vector< std::vector< unsigned long > > slices;
      std::vector< double > lengths, widths, heights;
      std::vector< int > nxs, nys, nzs;
      MeshSubdivide( gridin, ldi1, ldi2, nx, ny, nz, \
                     length, width, height, MX, MY, MZ, \
                     slices, lengths, widths, heights, \
                     nxs, nys, nzs );
      for (int sd = 0; sd < slices.size(); sd++)
      {
        std::cout << "\nMesh # " << sd << "; size = " << slices[sd].size() << ".\n";
        for (int jj = 0; jj < slices[sd].size(); jj++)
        {
          std::cout << slices[sd][jj] << "\n";
        }
      }

      // build meshes
      std::vector< FluidMesh > Meshes;
      Meshes.resize( slices.size() );
      for (int sd = 0; sd < slices.size(); sd++) {
        Meshes[sd].BuildUniformMesh( slices[sd].data(), nys[sd], nzs[sd], nxs[sd], nys[sd], nzs[sd], lengths[sd], widths[sd], heights[sd] );
      }

      // solve porescale problems
      std::vector< std::vector< double > > Solutions;
      Solutions.resize( Meshes[0].DIM*slices.size() );
      if ( nz ) {
        for (int sd = 0; sd < slices.size(); sd++) {
          Solutions[ idx2( sd, 0, 3 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 1, 3 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 2, 3 ) ].resize( Meshes[sd].dofTotal );
          StokesSolveDirect( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 3 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolveDirect( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 3 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolveDirect( Meshes[sd], visc, 2, Solutions[ idx2( sd, 2, 3 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
        }
      }
      else {
        for (int sd = 0; sd < slices.size(); sd++) {
          Solutions[ idx2( sd, 0, 2 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 1, 2 ) ].resize( Meshes[sd].dofTotal );
          StokesSolveDirect( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 2 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
          StokesSolveDirect( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 2 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
        }
      }
      // post-process porescale solutions
      std::vector< double > Ks;
      Ks.resize( Meshes[0].DIM*slices.size() );
      if ( nz ) {
        for (int sd = 0; sd < slices.size(); sd++)
        {
          computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 0, 3 ) ], Ks[ idx2( sd, 0, 3 ) ], 0, 0 );
          computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 1, 3 ) ], Ks[ idx2( sd, 1, 3 ) ], 1, 0 );
          computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 2, 3 ) ], Ks[ idx2( sd, 2, 3 ) ], 2, 0 );
          std::cout << "\nDomain " << sd << ": Mesh size: " << Meshes[sd].DOF[0] << " cells \t Kxx = " << Ks[idx2(sd,0,3)] << "\t Kyy = " << Ks[idx2(sd,1,3)];
          std::cout << "\t Kzz = " << Ks[idx2(sd,2,3)] << "\n";
        }
      }
      else {
        for (int sd = 0; sd < slices.size(); sd++) {
          computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 0, 2 ) ], Ks[ idx2( sd, 0, 2 ) ], 0, 0 );
          computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 1, 2 ) ], Ks[ idx2( sd, 1, 2 ) ], 1, 0 );
          std::cout << "\nDomain " << sd << ": Mesh size: " << Meshes[sd].DOF[0] << " cells \t Kxx = " << Ks[idx2(sd,0,2)] << "\t Kyy = " << Ks[idx2(sd,1,2)] << "\n";
        }
      }

      // solve pore-network model with porescale computed Ks

    }
  }

  if (simNum == numSims) stop_paralution();
}

