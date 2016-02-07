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
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 0, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matIs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
          matJs.push_back( Mesh.PressureCellVNeighbor[ \
              idx2( cl, 1, Mesh.PressureCellVelocityNeighborLDI ) ] + Mesh.DOF[1] + Mesh.DOF[2]);
        }
      }
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
    }
  }
}
