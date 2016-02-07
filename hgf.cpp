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
      StokesSolve( Mesh, visc, direction, sol, tolAbs, tolRel, maxIt, nThreads, prec );

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
    case 1 : // Solve all 3 flow directions for upscaled tensor
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
    case 2 :
    { // compute conductivities on subdomains, then solve constructed pore-network problem
      switch ( nz )
      {
        case 0 :
        {
          int nSubDomains = 0;
          int xStart, nxRemainder, yStart, nyRemainder, nyG, nxG, xCount, yCount;
          // determine array slices for subdomains
          std::vector< std::vector< unsigned long > > slices;
          std::vector< double > lengths, widths;
          std::vector< int > nxs, nys;
          slices.resize( MX * MY );
          lengths.resize( MX * MY );
          widths.resize( MX * MY );
          nxs.resize( MX * MY );
          nys.resize( MX * MY );
          yCount = 0;
          yStart = 0;
          for (int yy = 0; yy < MY; yy++) {
            xStart = 0;
            nxRemainder = nx;
            xCount = 0;
            nyG = (int)(round((double)(nyRemainder)/(MY-yCount)));
            for (int xx = 0; xx < MX; xx++) {
              nSubDomains++;
              nxG = (int)(round((double)(nxRemainder)/(MX-xCount)));
              for (int cx = 0; cx < nxG; cx++) {
                for (int cy = 0; cy < nyG; cy++) {
                  slices[ nSubDomains-1 ].push_back( gridin[ idx2( (cx+xStart), (cy+yStart), ldi1 ) ] );
                }
              }
              lengths[ nSubDomains-1 ] = length * ((double)nxG / nx);
              widths[ nSubDomains-1 ] = width * ((double)nyG / ny);
              nxs[ nSubDomains-1 ] = nxG;
              nys[ nSubDomains-1 ] = nyG;
              xStart = xStart + nxG;
              nxRemainder = nx - xStart;
              xCount++;
            }
            yStart = yStart + nyG;
            nyRemainder = ny - yStart;
            yCount++;
          }

          // build meshes
          std::vector< FluidMesh > Meshes;
          Meshes.resize(nSubDomains);
          for (int sd = 0; sd < nSubDomains; sd++) {
            Meshes[sd].BuildUniformMesh( slices[sd].data(), nys[sd], 0, nxs[sd], nys[sd], 0, lengths[sd], widths[sd], 0 );
          }

          // solve porescale problems
          std::vector< std::vector< double > > Solutions;
          Solutions.resize( 2*nSubDomains );
          for (int sd = 0; sd < nSubDomains; sd++) {
            Solutions[ idx2( sd, 0, 2 ) ].resize( Meshes[sd].dofTotal );
            Solutions[ idx2( sd, 1, 2 ) ].resize( Meshes[sd].dofTotal );
            StokesSolve( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 2 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
            StokesSolve( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 2 ) ], tolAbs, tolRel, maxIt, nThreads, prec );
          }

          // post-process porescale solutions
          std::vector< double > Ks;
          Ks.resize( 2*nSubDomains );
          for (int sd = 0; sd < nSubDomains; sd++)
          {
            computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 0, 2 ) ], Ks[ idx2( sd, 0, 2 ) ], 0, 0 );
            computeKConstantDrive( Meshes[sd], Solutions[ idx2( sd, 1, 2 ) ], Ks[ idx2( sd, 1, 2 ) ], 1, 0 );
            std::cout << "\nDomain " << sd << ": Mesh size: " << Meshes[sd].DOF[0] << " cells \t Kxx = " << Ks[idx2(sd,0,2)] << "\t Kyy = " << Ks[idx2(sd,1,2)] << "\n";
          }

          // solve pore-network model with porescale computed Ks


          break; // break for 2d dim switch
        }
        default :
        {
          // determine array slices for subdomains

          // build meshes

          // solve porescale problems

          // solve pore-network model with porescale computed Ks

          break; // break for 3d dim switch
        }
      }
      break; // break for type 2 problem solve
    }
  }

  if (simNum == numSims) stop_paralution();
}

