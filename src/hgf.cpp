#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>
#include <math.h>
#include <iterator>
#include <string>
#include <sstream>

#include <paralution.hpp>

#include "hgf.hpp"
#include "hgfMesh.hpp"
#include "hgfIB.hpp"
#include "hgfStokes.hpp"
#include "hgfPoreNetwork.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

/* current g++ has a bug claiming to_string() is not a memeber of std.
   when this is fixed this "patch" can be deleted. called as patch::to_string()*/
namespace patch
{
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm;
    stm << n;
    return stm.str();
  }
}

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
          int MX, int MY, int MZ, int solver, double relax, int output )
{
  std::string outExt;
  if ( output == 0 ) outExt = ".dat";
  else outExt = ".vtk";
  int type;

  if (simNum == 1) init_paralution();

  if (direction < 3) type = 0;
  else if (direction == 3) type = 1;
  else if (direction == 4) type = 2;
  else if (direction == 5) type = 3;
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
      outName = "./output/flowrun" + outExt;
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
      if (solver == 0) StokesSolveDirect( Mesh, visc, direction, sol, tolAbs, tolRel, maxIt, nThreads, prec );
      else StokesSolveRich( Mesh, visc, direction, sol, tolAbs, tolRel, maxIt, nThreads, prec, relax );

      // Linear solve timer
      stokes_duration = ( omp_get_wtime() - stokes_start );

      // Write solution to file
      if (output) writeSolutionPV ( Mesh, sol, outName );
      else writeSolutionTP( Mesh, sol, outName );
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
      outNameX = "./output/flowrunX" + outExt;
      outNameY = "./output/flowrunY" + outExt;
      outNameZ = "./output/flowrunZ" + outExt;
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
      if (solver == 0) {
        StokesSolveDirect( Mesh, visc, 0, solX, tolAbs, tolRel, maxIt, nThreads, prec );
        StokesSolveDirect( Mesh, visc, 1, solY, tolAbs, tolRel, maxIt, nThreads, prec );
        if ( Mesh.DIM == 3) {
          solZ.resize( Mesh.dofTotal );
          StokesSolveDirect( Mesh, visc, 2, solZ, tolAbs, tolRel, maxIt, nThreads, prec );
        }
        else solZ.resize( 1 );
      }
      else {
        StokesSolveRich( Mesh, visc, 0, solX, tolAbs, tolRel, maxIt, nThreads, prec, relax );
        StokesSolveRich( Mesh, visc, 1, solY, tolAbs, tolRel, maxIt, nThreads, prec, relax );
        if ( Mesh.DIM == 3) {
          solZ.resize( Mesh.dofTotal );
          StokesSolveRich( Mesh, visc, 2, solZ, tolAbs, tolRel, maxIt, nThreads, prec, relax );
        }
        else solZ.resize( 1 );
      }

      // stokes timer
      stokes_duration = ( omp_get_wtime() - stokes_start );

      // Write solution to file
      if (output) {
        writeSolutionPV ( Mesh, solX, outNameX );
        writeSolutionPV ( Mesh, solY, outNameY );
        writeSolutionPV ( Mesh, solZ, outNameZ );
      }
      else {
        writeSolutionTP ( Mesh, solX, outNameX );
        writeSolutionTP ( Mesh, solY, outNameY );
        writeSolutionTP ( Mesh, solZ, outNameZ );
      }
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

      double mesh_duration, stokes_duration, pn_duration, total_duration;
      double start, stokes_start, pn_start;

      start = omp_get_wtime();
      // determine array slices for subdomains
      std::vector< std::vector< unsigned long > > slices;
      std::vector< double > lengths, widths, heights;
      std::vector< int > nxs, nys, nzs;
      MeshSubdivide( gridin, ldi1, ldi2, nx, ny, nz, \
                     length, width, height, MX, MY, MZ, \
                     slices, lengths, widths, heights, \
                     nxs, nys, nzs );

      // build meshes
      std::vector< FluidMesh > Meshes;
      Meshes.resize( slices.size() );
      for (unsigned long sd = 0; sd < slices.size(); sd++) {
        Meshes[sd].BuildUniformMesh( slices[sd].data(), nys[sd], nzs[sd], nxs[sd], nys[sd], nzs[sd], \
                                     lengths[sd], widths[sd], heights[sd] );
      }

      mesh_duration = ( omp_get_wtime() - start );
      stokes_start = omp_get_wtime();

      // solve porescale problems
      std::vector< std::vector< double > > Solutions;
      Solutions.resize( Meshes[0].DIM*slices.size() );

      if ( nz ) {
        for (unsigned long sd = 0; sd < slices.size(); sd++) {
          Solutions[ idx2( sd, 0, 3 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 1, 3 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 2, 3 ) ].resize( Meshes[sd].dofTotal );
          if (solver == 0) {
            std::cout << "\nSolving x-flow problem on submesh " << sd << "\n\n";
            StokesSolveDirect( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 3 ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
            std::cout << "\nSolving y-flow problem on submesh " << sd << "\n\n";
            StokesSolveDirect( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 3 ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
            std::cout << "\nSolving z-flow problem on submesh " << sd << "\n\n";
            StokesSolveDirect( Meshes[sd], visc, 2, Solutions[ idx2( sd, 2, 3 ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
          }
          else {
            std::cout << "\nSolving x-flow problem on submesh " << sd << "\n\n";
            StokesSolveRich( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 3 ) ], \
                             tolAbs, tolRel, maxIt, nThreads, prec, relax );
            std::cout << "\nSolving y-flow problem on submesh " << sd << "\n\n";
            StokesSolveRich( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 3 ) ], \
                             tolAbs, tolRel, maxIt, nThreads, prec, relax );
            std::cout << "\nSolving z-flow problem on submesh " << sd << "\n\n";
            StokesSolveRich( Meshes[sd], visc, 2, Solutions[ idx2( sd, 2, 3 ) ], \
                             tolAbs, tolRel, maxIt, nThreads, prec, relax );
          }
        }
      }
      else {
        for (unsigned long sd = 0; sd < slices.size(); sd++) {
          Solutions[ idx2( sd, 0, 2 ) ].resize( Meshes[sd].dofTotal );
          Solutions[ idx2( sd, 1, 2 ) ].resize( Meshes[sd].dofTotal );
          if (solver == 0) {
            std::cout << "\nSolving x-flow problem on submesh " << sd << "\n\n";
            StokesSolveDirect( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 2 ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
            std::cout << "\nSolving y-flow problem on submesh " << sd << "\n\n";
            StokesSolveDirect( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 2 ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
          }
          else {
            std::cout << "\nSolving x-flow problem on submesh " << sd << "\n\n";
            StokesSolveRich( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, 2 ) ], \
                             tolAbs, tolRel, maxIt, nThreads, prec, relax );
            std::cout << "\nSolving y-flow problem on submesh " << sd << "\n\n";
            StokesSolveRich( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, 2 ) ], \
                             tolAbs, tolRel, maxIt, nThreads, prec, relax );
          }
        }
      }

      stokes_duration = ( omp_get_wtime() - stokes_start );
      pn_start = omp_get_wtime();

      // post-process porescale solutions
      std::vector< double > Ks;
      Ks.resize( Meshes[0].DIM*slices.size() );
      if ( nz ) {
        for (unsigned long sd = 0; sd < slices.size(); sd++) {
          computeKConstantDrive( Meshes[sd], \
            Solutions[ idx2( sd, 0, 3 ) ], Ks[ idx2( sd, 0, 3 ) ], 0, 0 );
          computeKConstantDrive( Meshes[sd], \
            Solutions[ idx2( sd, 1, 3 ) ], Ks[ idx2( sd, 1, 3 ) ], 1, 0 );
          computeKConstantDrive( Meshes[sd], \
            Solutions[ idx2( sd, 2, 3 ) ], Ks[ idx2( sd, 2, 3 ) ], 2, 0 );
        }
      }
      else {
        for (unsigned long sd = 0; sd < slices.size(); sd++) {
          computeKConstantDrive( Meshes[sd], \
            Solutions[ idx2( sd, 0, 2 ) ], Ks[ idx2( sd, 0, 2 ) ], 0, 0 );
          computeKConstantDrive( Meshes[sd], \
            Solutions[ idx2( sd, 1, 2 ) ], Ks[ idx2( sd, 1, 2 ) ], 1, 0 );
        }
      }
      if (Meshes[0].DIM == 2) MZ = 0;

      // solve pore-network model with porescale computed Ks
      PoreNetwork pn;
      pn.UniformPN( length, width, height, MX, MY, MZ );
      std::vector<double> PNSolution;
      PNSolution.resize( pn.nPores );
      double KPN;
      std::cout << "\nSolving the Pore-Network problem\n\n";
      PoreNetworkSolveDirect( pn, Ks, PNSolution, 0 );

      pn_duration = ( omp_get_wtime() - pn_start );
      total_duration = ( omp_get_wtime() - start );

      // Write PN
      std::string outName;
      outName = "./output/PNSolution.dat";
      writePoreNetworkSolutionTP ( pn, PNSolution, outName );

      // compute and save K
      computeKPoreNetwork( pn, PNSolution, Ks, KPN, 0, 1 );

      std::cout << "\nPorescale meshes constructed in " << mesh_duration << "seconds\n";
      std::cout << "Stokes problems solved in " << stokes_duration << "seconds\n";
      std::cout << "Pore-Network problem solved in " << pn_duration << "seconds\n";
      // Total timers
      total_duration = ( omp_get_wtime() - start );
      std::cout << "Total time: " << total_duration << "seconds\n";
      break;
    }
    case 3 :
    { // compute empirical distributions on subdomains, sample for constructed pore-network problem
      std::string outNameX;
      std::string outNameY;
      std::string outNameZ;
      int NVF = 2;
      int nSims = 2;
      double mesh_duration, stokes_duration, pn_duration, total_duration;
      double start, stokes_start, pn_start;
      double objectVolumeFrac = 0.05;

      start = omp_get_wtime();
      // determine array slices for subdomains
      std::vector< std::vector< unsigned long > > slices;
      std::vector< double > lengths, widths, heights;
      std::vector< int > nxs, nys, nzs;
      MeshSubdivide( gridin, ldi1, ldi2, nx, ny, nz, \
                     length, width, height, MX, MY, MZ, \
                     slices, lengths, widths, heights, \
                     nxs, nys, nzs );

      // build meshes
      std::vector< FluidMesh > Meshes;
      Meshes.resize( slices.size() );
      for (unsigned long sd = 0; sd < slices.size(); sd++) {
        Meshes[sd].BuildUniformMesh( slices[sd].data(), nys[sd], nzs[sd], nxs[sd], nys[sd], nzs[sd], \
                                     lengths[sd], widths[sd], heights[sd] );
      }

      mesh_duration = ( omp_get_wtime() - start );
      stokes_start = omp_get_wtime();

      // solve porescale problems
      std::vector< std::vector< double > > Solutions;
      Solutions.resize( Meshes[0].DIM*slices.size() );

      // loop for subdomains
      for (int sd = 0; sd < Meshes.size(); sd++) {
        Solutions[ idx2( sd, 0, Meshes[sd].DIM ) ].resize( Meshes[sd].dofTotal );
        Solutions[ idx2( sd, 1, Meshes[sd].DIM ) ].resize( Meshes[sd].dofTotal );
        if (Meshes[sd].DIM == 3) Solutions[ idx2( sd, 2, Meshes[sd].DIM ) ].resize( Meshes[sd].dofTotal );
        // loop for volume franction
        for (int vf = 0; vf < NVF; vf++) {
          // loop for samples
          for (int spl = 0; spl < nSims; spl++) {

            // console info
            std::cout << "\n\n";
            std::cout << "---------------------------------------------------------------------\n";
            std::cout << "Solving sample " << spl << " for volume fraction ";
            std::cout << (vf+1)*objectVolumeFrac << " on subdomain " << sd << "\n.";
            std::cout << "---------------------------------------------------------------------\n";

            // build and set immersed boundary
            BuildImmersedBoundary( Meshes[sd], objectVolumeFrac, (vf+1) );
            // solve problems
            StokesSolveDirect( Meshes[sd], visc, 0, Solutions[ idx2( sd, 0, Meshes[sd].DIM ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
            StokesSolveDirect( Meshes[sd], visc, 1, Solutions[ idx2( sd, 1, Meshes[sd].DIM ) ], \
                               tolAbs, tolRel, maxIt, nThreads, prec );
            if (Meshes[sd].DIM == 3) {
              StokesSolveDirect( Meshes[sd], visc, 2, Solutions[ idx2( sd, 2, Meshes[sd].DIM ) ], \
                                 tolAbs, tolRel, maxIt, nThreads, prec );
            }
            // store solutions
            outNameX = "./output/domain" + patch::to_string(sd) + "_" + patch::to_string((vf+1)) + "obs_sample" + patch::to_string(spl) + "_X" + outExt;
            outNameY = "./output/domain" + patch::to_string(sd) + "_" + patch::to_string((vf+1)) + "obs_sample" + patch::to_string(spl) + "_Y" + outExt;
            if (Meshes[sd].DIM == 3) {
              outNameZ = "./output/domain" + patch::to_string(sd) + "_" + patch::to_string((vf+1)) + "obs_sample" + patch::to_string(spl) + "_Z" + outExt;
            }
            if (output) {
              writeSolutionPV ( Meshes[sd], Solutions[ idx2( sd, 0, Meshes[sd].DIM ) ], outNameX );
              writeSolutionPV ( Meshes[sd], Solutions[ idx2( sd, 1, Meshes[sd].DIM ) ], outNameY );
              if (Meshes[sd].DIM == 3) {
                writeSolutionPV ( Meshes[sd], Solutions[ idx2( sd, 2, Meshes[sd].DIM ) ], outNameZ );
              }
            }
            else {
              writeSolutionTP ( Meshes[sd], Solutions[ idx2( sd, 0, Meshes[sd].DIM ) ], outNameX );
              writeSolutionTP ( Meshes[sd], Solutions[ idx2( sd, 1, Meshes[sd].DIM ) ], outNameY );
              if (Meshes[sd].DIM == 3) {
                writeSolutionTP ( Meshes[sd], Solutions[ idx2( sd, 2, Meshes[sd].DIM ) ], outNameZ );
              }
            }
            // compute K

            // save K in appropriate vector

          }
        }
      }

      break; // break type 3
    }
  }
  if (simNum == numSims) stop_paralution();
}
