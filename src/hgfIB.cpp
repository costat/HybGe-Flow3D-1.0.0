#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <math.h>

// hgf includes
#ifndef CUDA_BUILD
# define CUDA_BUILD 0
#endif

#if CUDA_BUILD
#include "hgfMeshCu.cuh"
#else
#include "hgfMesh.hpp"
#endif

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
void
BuildImmersedBoundary( FluidMesh& Mesh, double vf, int nObs )
{
  int radSize = std::max((int)(round( pow( (vf * Mesh.ImmersedBoundary.size()), (double)(1.0/Mesh.DIM) ) )), 1);
  // set IB to zeros before randomly filling blocks
  std::fill( Mesh.ImmersedBoundary.begin(), Mesh.ImmersedBoundary.end(), 0 );
  std::vector< unsigned long > cellNums;
  cellNums.resize( pow(radSize, Mesh.DIM) );
  double blockNorm;
  int nObsPlaced = 0;
  int seedCell;
  // generate a random location and test if a block of appropriate size 'fits' from that anchor
  if (Mesh.DIM == 2) goto newseed2;
  else goto newseed3;
  // 2d seed generation
  newseed2 :
  {
    seedCell = rand() % (int)Mesh.ImmersedBoundary.size();
    for (int ii = 0; ii < radSize; ii++) {
      if (ii == 0) {
        cellNums[ idx2( ii, 0, radSize) ] = seedCell;
      }
      else {
        if (Mesh.PFaceConnectivity[ idx2( cellNums[ idx2( (ii-1), 0, radSize ) ], 1, 4 ) ]) {
          cellNums[ idx2( ii, 0, radSize) ] = Mesh.PFaceConnectivity[ idx2( cellNums[ idx2( (ii-1), 0, radSize ) ], 1, 4 ) ]-1;
        }
        else {
          goto newseed2;
        }
      }
      for (int jj = 1; jj < radSize; jj++) {
        if (Mesh.PFaceConnectivity[ idx2( cellNums[ idx2( ii, (jj-1), radSize ) ], 2, 4 ) ]) {
          cellNums[ idx2( ii, jj, radSize ) ] = Mesh.PFaceConnectivity[ idx2( cellNums[ idx2( ii, (jj-1), radSize ) ], 2, 4 ) ]-1;
        }
        else {
          goto newseed2;
        }
      }
    }
    goto blockcheck;
  }
  // 3d seed generation
  newseed3 :
  {
    seedCell = rand() % (int)Mesh.ImmersedBoundary.size();

    // first we set the anchor entries (ii, jj, 0)
    for (int ii = 0; ii < radSize; ii++) {
      if ( ii == 0 ) {
        cellNums[ idx3( ii, 0, 0, radSize, radSize ) ] = seedCell;
      }
      else {
        if (Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( (ii-1), 0, 0, radSize, radSize ) ], 1, 6 ) ]) {
          cellNums[ idx3( ii, 0, 0, radSize, radSize ) ] = \
            Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( (ii-1), 0, 0, radSize, radSize ) ], 1, 6 ) ]-1;
        }
        else {
          goto newseed3;
        }
      }
      for (int jj = 1; jj < radSize; jj++) {
        if (Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( ii, (jj-1), 0, radSize, radSize ) ], 4, 6 ) ]) {
          cellNums[ idx3(ii, jj, 0, radSize, radSize ) ] = \
            Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( ii, (jj-1), 0, radSize, radSize ) ], 4, 6 ) ]-1;
        }
        else {
          goto newseed3;
        }
      }
    }
    // now we set the entries (*,*,kk)
    for (int ii = 0; ii < radSize; ii++) {
      for (int jj = 0; jj < radSize; jj++) {
        for (int kk = 1; kk < radSize; kk++) {
          if (Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( ii, jj, (kk-1), radSize, radSize ) ], 2, 6 ) ]) {
            cellNums[ idx3( ii, jj, kk, radSize, radSize ) ] = \
              Mesh.PFaceConnectivity[ idx2( cellNums[ idx3( ii, jj, (kk-1), radSize, radSize ) ], 2, 6 ) ]-1;
          }
          else {
            goto newseed3;
          }
        }
      }
    }
    goto blockcheck;
  }
  blockcheck :
  {
    blockNorm = 0;
    // define blockNorm by summing current IB entries. if empty, set cells to IB;
    for (int ii = 0; ii < cellNums.size(); ii++) {
      blockNorm = blockNorm + Mesh.ImmersedBoundary[ cellNums[ ii ] ];
    }
    // if blockNorm is 0 then the location is void space. we then set the IB. otherwise return to newseed and start over
    if (!blockNorm) {
      for (int ii = 0; ii < cellNums.size(); ii++) {
        Mesh.ImmersedBoundary[ cellNums[ ii ] ] = 2;
      }
      nObsPlaced++;
      if (nObsPlaced < nObs) {
        if (Mesh.DIM == 2) goto newseed2;
        else goto newseed3;
      }
    }
    else {
      if (Mesh.DIM == 2) goto newseed2;
      else goto newseed3;
    }
  }
}
