#include <iostream>
#include <cstdlib>

#include <paralution.hpp>

#include "hgfMesh.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

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
  int radSize = max((int)(round(sqrt(vf * Mesh.ImmersedBoundary.size()))), 1);
  std::vector< unsigned long > cellNums;
  cellNums.resize( radSize * radSize );
  double blockNorm = 1;
  int seedCell;
  while (blockNorm) {
    seedCell = rand() % (int)Mesh.ImmersedBoundary.size();
    for (int ii = 0; ii < radSize; ii++) {
      if (ii == 0) {
        cellNums[ idx2( ii, 0, radSize) ] = seedCell;
      }
      else {
        if (Mesh.FaceConnectivity[ idx2( cellNums[ idx2( (ii-1), 0, radSize ) ], 1, 4 ) ]) {

        }
        else {
          // break out to while level, new seed
        }
      }
      for (int jj = 1; jj < radSize; jj++) {
        if (Mesh.FaceConnectivity[ idx2( cellNums[ idx2( ii, (jj-1), radSize ) ], 2, 4 ) ]) {

        }
        else {
           // break out to while loop level, new seed
        }
      }
    }
    // define blockNorm by summing current IB entries. if empty, set cells to IB;
    for (int ii = 0; ii < cellNums.size(); ii++) {
      blockNorm = blockNorm + Mesh.ImmersedBoundary[ cellNums[ ii ] ];
    }
    if (!blockNorm) {
      for (int ii = 0; ii < cellNum.size(); ii++) {
        Mesh.ImmersedBoundary[ cellNums[ ii ] ] = 1;
      }
    }
  }
}
