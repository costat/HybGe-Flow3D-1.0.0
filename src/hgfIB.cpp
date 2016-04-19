#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

/* immersedBoundary adds a penalty to the diagonal entry in rows associated
   to immersed boundary voxels. */
void
immersedBoundary ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                   std::vector<int>& matJs, std::vector<double>& matVals )
{

  double pen = 10000;

  switch ( Mesh.DIM )
  {
    case 3 :
    {
      /* Immersed boundary voxels are tracked from the original grid, which
         is given by pressure cells. Loop through pressure nodes
         to check IB status */
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        // Check if pressure cell is an IB voxel
        if (Mesh.ImmersedBoundary[cl] == 2) {
          /* Add entries in the COO matrix storage
             vectors. Note these will be repeated indices. Paralution handles
             this by summing repeated entries. If a different solver is used that
             does not allow repeated entries, care should be taken to
             overcome this issue. */
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matJs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matIs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matJs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
        }
      }
      break;
    }
    case 2 :
    {
      /* Immersed boundary voxels are tracked from the original grid, which
         is given by pressure cells. Loop through pressure nodes
         to check IB status */
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        // Check if pressure cell is an IB voxel
        if (Mesh.ImmersedBoundary[cl] == 2) {
          /* Add entries in the COO matrix storage
             vectors. Note these will be repeated indices. Paralution handles
             this by summing repeated entries. If a different solver is used that
             does not allow repeated entries, care should be taken to
             overcome this issue. */
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] \
              + Mesh.DOF[1]);
          matVals.push_back(pen);
        }
      }
      break;
    }
  }
}
void
immersedBoundarySingleComponent ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                                  std::vector<int>& matJs, std::vector<double>& matVals, \
                                  int component )
{

  double pen = 10000;

  switch ( component )
  {
    case 0 :
    {
       for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        // Check if pressure cell is an IB voxel
        if (Mesh.ImmersedBoundary[cl] == 2) {
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellUNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
        }
      }
      break;
    }
    case 1 :
    {
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        // Check if pressure cell is an IB voxel
        if (Mesh.ImmersedBoundary[cl] == 2) {
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
        }
      }
      break;
    }
    case 2 :
    {
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        // Check if pressure cell is an IB voxel
        if (Mesh.ImmersedBoundary[cl] == 2) {
          matIs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
          matIs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matJs.push_back( Mesh.PressureCellWNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] );
          matVals.push_back(pen);
        }
      }
      break;
    }
  }
}
int
BuildImmersedBoundary( FluidMesh& Mesh, double vf, int nObs )
{
  double dvfrac = nObs * vf + 1;
  int nIBCells = std::max((int)round( dvfrac * Mesh.ImmersedBoundary.size() ));

  // set IB to zeros before randomly filling blocks
  std::fill( Mesh.ImmersedBoundary.begin(), Mesh.ImmersedBoundary.end(), 0 );
  int nObsPlaced = 0;
  int nFails = 0;
  int failMax = 1000;
  int seedCell;
  // generate a random location and test if that spot is free for an IB cell
  do
  {
    seedCell = rand() % (int)Mesh.ImmersedBoundary.size();
    if (!Mesh.ImmersedBoundary[ seedCell ]) {
      Mesh.ImmersedBoundary[ seedCell ] = 2;
      nFails = 0;
      nObsPlaced++;
    }
    else {
      nFails++;
    }
  } while( nObsPlaced < nIBCells && nFails < failMax )

  if (nFails > failMax) return 1;
  else return 0;
}
