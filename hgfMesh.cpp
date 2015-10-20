#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <omp.h>
#include "hgfMesh.hpp"
#include "hgf.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index, uses row major ordering */
#define idx2(i, j, ldi) ((i * ldi) + j)
/* Define a 3d -> 1d array index, uses row major ordering */
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

// Function to construct the mesh from voxel array input
void FluidMesh::BuildUniformMesh( unsigned long *gridin, int ldi1, int ldi2, \
                                  int nx, int ny, int nz, \
                                  double length, double width, double height )
{
  int checkVert;
  int maxPNodes;
  int numPCells = 0;
  xLim[0] = 0;
  xLim[1] = length;
  yLim[0] = 0;
  yLim[1] = width;
  zLim[0] = 0;
  zLim[1] = height;

  switch ( nz )
  {
    case 0 : // 2D problem
    {
      // Constants determined by dimension alone
      DIM = 2;
      NodesLDI = 2;
      CellCentersLDI = 2;
      CellWidthsLDI = 2;
      FaceConnectivityLDI = 4;
      PressureCellVelocityNeighborLDI = 2;
      VelocityCellPressureNeighborLDI = 2;
      NX = nx;
      NY = ny;
      NZ = nz;

      // Max pressure nodes
      std::vector<double> nodeHold;
      nodeHold.resize(2);
      for (int yi = 0; yi < ny; yi++)
      {
        for (int xi = 0; xi < nx; xi++)
        {
          if (gridin[ idx2(xi, yi, ldi1) ] != 1)
          {
            numPCells++;
          }
        }
      }
      maxPNodes = numPCells * 4;
      int nNodes = 0;
      double dx = length / nx;
      double dy = width / ny;
      double dz = 0;
      int countCell = -1;
      mv.resize( (numPCells * 4) );
      double cellVert [ 8 ];

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((nx * ny));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((maxPNodes * 2));

      for (int yi = 0; yi < ny; yi++)
      {
        for (int xi = 0; xi < nx; xi++)
        {
          FullGrid.push_back(gridin[ idx2(xi, yi, ldi1) ]);
          if (gridin[ idx2( xi, yi, ldi1 ) ] != 1)
          {
            countCell++;
            ImmersedBoundary.push_back(gridin[ idx2(xi, yi, ldi1) ]);
            cellVert[0] = (xi + 1) * dx - dx;
            cellVert[1] = (yi + 1) * dy - dy;
            cellVert[2] = (xi + 1) * dx;
            cellVert[3] = cellVert[1];
            cellVert[4] = cellVert[2];
            cellVert[5] = (yi + 1) * dy;
            cellVert[6] = cellVert[0];
            cellVert[7] = cellVert[5];

            for (int pcount = 0; pcount < 4; pcount++)
            {
              nodeHold[0] = cellVert[ idx2( pcount, 0, 2 ) ];
              nodeHold[1] = cellVert[ idx2( pcount, 1, 2 ) ];
              if ( !countCell ) // countCell = 0 -> first cell so no possible
                                // node duplicates
              {
                checkVert = -1;
              }
              else
              {
                checkVert = isNear( nodeHold, Nodes, dx, dy, dz, nNodes, DIM );
              }
              if (checkVert == -1) // node is not a duplicate
              {
                nNodes++;
                Nodes.push_back(nodeHold[0]);
                Nodes.push_back(nodeHold[1]);
                mv[ idx2( countCell, pcount, 4 ) ] = nNodes-1;
              }
              else // node is a duplicate
              {
                mv[ idx2( countCell, pcount, 4 ) ] = checkVert;
              }
            }
          }
        }
      }
      Nodes.resize(Nodes.size());

      // Next we compute cell centers for pressures
      PCellCenters.resize((numPCells * 2));
      for (int cl = 0; cl < numPCells; cl++)
      {
        PCellCenters[ idx2( cl, 0, 2 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 0, 4 ) ], 0, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 1, 4) ], 0, NodesLDI ) ]);
        PCellCenters[ idx2( cl, 1, 2 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 1, 4 ) ], 1, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 2, 4) ], 1, NodesLDI ) ]);
      }

      DOF.resize(3);

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, independent
         of other velocity components */
      #pragma omp parallel sections
      {
        { // Final P computations
          DOF[0] = numPCells;
          PFaceConnectivity.resize((numPCells * 4));
          innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells );
          PCellWidths.resize((DOF[0] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[0]; cl++) {
            PCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            PCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
          }
        }
        #pragma omp section
        { // U grid
          int checkVertU;
          std::vector<double> nodeHoldU;
          nodeHoldU.resize(2);
          int countUCells = 0;
          double uStep = 0.5*dx;
          int maxUCells = nx * ny + ny;
          UCellCenters.reserve((maxUCells * 2));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++)
          {
            cellVert[0] = PCellCenters[ idx2( cl, 0, 2 ) ] - uStep;
            cellVert[1] = PCellCenters[ idx2( cl, 1, 2 ) ];
            cellVert[2] = PCellCenters[ idx2( cl, 0, 2 ) ] + uStep;
            cellVert[3] = PCellCenters[ idx2( cl, 1, 2 ) ];
            for (int pcount = 0; pcount < 2; pcount++)
            {
              // U Component
              nodeHoldU[0] = cellVert[ idx2( pcount, 0, 2 ) ];
              nodeHoldU[1] = cellVert[ idx2( pcount, 1, 2 ) ];
              if ( !countUCells ) // First U node
              {
                checkVertU = -1;
              }
              else
              {
                checkVertU = isNear( nodeHoldU, UCellCenters, \
                                     dx, dy, dz, countUCells, DIM );
              }
              if (checkVertU == -1) // cell center location is not a duplicate
              {
                countUCells++;
                UCellCenters.push_back(nodeHoldU[0]);
                UCellCenters.push_back(nodeHoldU[1]);
                if (pcount == 0) {
                  UCellPressureNeighbor[ idx2( (countUCells-1), 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countUCells;
                }
                else if (pcount == 1) {
                  UCellPressureNeighbor[ idx2( (countUCells-1), 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countUCells;
                }
              }
              if (checkVertU != -1)
              {
                if (pcount == 0) {
                  UCellPressureNeighbor[ idx2( checkVertU, 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertU+1;
                }
                else if (pcount== 1) {
                  UCellPressureNeighbor[ idx2( checkVertU, 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertU+1;
                }
              }
            }
          }
          UCellCenters.resize(UCellCenters.size());
          DOF[1] = countUCells;
          UCellPressureNeighbor.resize((DOF[1] * 2));
          UFaceConnectivity.resize((countUCells * 4));
          innerFaceConnectivity( UFaceConnectivity, UCellCenters, \
                                 dx, dy, dz, countUCells );
          UInteriorCells.reserve(DOF[1]);
          UBoundaryCells.reserve(DOF[1]);
          int nbrsu;
          for (int cl = 0; cl < DOF[1]; cl++) {
            nbrsu = 0;
            for (int position = 0; position < 4; position++) {
              if (UFaceConnectivity[ idx2( cl, position, FaceConnectivityLDI ) ] != 0) {
                nbrsu++;
              }
            }
            if (nbrsu == 4) {
              UInteriorCells.push_back(cl);
            }
            else {
              UBoundaryCells.push_back(cl);
            }
          }
          UInteriorCells.resize(UInteriorCells.size());
          UBoundaryCells.resize(UBoundaryCells.size());
          UCellWidths.resize((DOF[1] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[1]; cl++) {
            UCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            UCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
          }
        }
        #pragma omp section
        { // V grid
          int checkVertV;
          std::vector<double> nodeHoldV;
          nodeHoldV.resize(2);
          int countVCells = 0;
          double vStep = 0.5*dy;
          int maxVCells = nx * ny + nx;
          VCellCenters.reserve((maxVCells * 2));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int cl = 0; cl < numPCells; cl++)
          {
            cellVert[4] = PCellCenters[ idx2( cl, 0, 2 ) ];
            cellVert[5] = PCellCenters[ idx2( cl, 1, 2 ) ] - vStep;
            cellVert[6] = PCellCenters[ idx2( cl, 0, 2 ) ];
            cellVert[7] = PCellCenters[ idx2( cl, 1, 2 ) ] + vStep;
            for (int pcount = 0; pcount < 2; pcount++)
            {
              // V Component
              nodeHoldV[0] = cellVert[ idx2( (pcount + 2), 0, 2 ) ];
              nodeHoldV[1] = cellVert[ idx2( (pcount + 2), 1, 2 ) ];
              if ( !countVCells ) // First V node
              {
                checkVertV = -1;
              }
              else
              {
                checkVertV = isNear( nodeHoldV, VCellCenters, \
                                    dx, dy, dz, countVCells, DIM );
              }
              if (checkVertV == -1)
              {
                countVCells++;
                VCellCenters.push_back(nodeHoldV[0]);
                VCellCenters.push_back(nodeHoldV[1]);
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
                else if (pcount == 1) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
              }
              if (checkVertV != -1)
              {
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( checkVertV, 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
                else if (pcount== 1) {
                  VCellPressureNeighbor[ idx2( checkVertV, 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
              }
            }
          }
          UCellCenters.resize(UCellCenters.size());
          DOF[2] = countVCells;
          VCellPressureNeighbor.resize((DOF[2] * 2));
          VFaceConnectivity.resize((countVCells * 4));
          innerFaceConnectivity( VFaceConnectivity, VCellCenters, \
                                 dx, dy, dz, countVCells );
          VInteriorCells.reserve(DOF[2]);
          VBoundaryCells.reserve(DOF[2]);
          int nbrsv;
          for (int cl = 0; cl < DOF[2]; cl++) {
            nbrsv = 0;
            for (int position = 0; position < 4; position++) {
              if (VFaceConnectivity[ idx2( cl, position, FaceConnectivityLDI ) ] != 0) {
                nbrsv++;
              }
            }
            if (nbrsv == 4) {
              VInteriorCells.push_back(cl);
            }
            else {
              VBoundaryCells.push_back(cl);
            }
          }
          VInteriorCells.resize(VInteriorCells.size());
          VBoundaryCells.resize(VBoundaryCells.size());
          VCellWidths.resize((DOF[2] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[2]; cl++) {
            VCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            VCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
          }
        }
      }
      break;
    }
    default : // 3D problem
    {
      // Constants determined by dimension alone
      DIM = 3;
      NodesLDI = 3;
      CellCentersLDI = 3;
      CellWidthsLDI = 3;
      FaceConnectivityLDI = 6;
      PressureCellVelocityNeighborLDI = 2;
      VelocityCellPressureNeighborLDI = 2;

      // Max pressure nodes
      std::vector<double> nodeHold;
      nodeHold.resize(3);
      for (int yi = 0; yi < ny; yi++)
      {
        for (int xi = 0; xi < nx; xi++)
        {
          for (int zi = 0; zi < nz; zi++)
          {
            if (gridin[ idx3(xi, yi, zi, ldi1, ldi2) ] != 1)
            {
              numPCells++;
            }
          }
        }
      }
      maxPNodes = numPCells * 8;

      int nNodes = 0;
      double dx = length / nx;
      double dy = width / ny;
      double dz = height / nz;
      int countCell = -1;
      mv.resize( (numPCells * 8) );
      double cellVert [ 24 ];

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((nx * ny * nz));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((maxPNodes * 3));
      for (int zi = 0; zi < nz; zi++)
      {
        for (int yi = 0; yi < ny; yi++)
        {
          for (int xi = 0; xi < nx; xi++)
          {
            FullGrid.push_back(gridin[ idx3(xi, yi, zi, ldi1, ldi2) ]);
            if (gridin[ idx3( xi, yi, zi, ldi1, ldi2 ) ] != 1)
            {
              countCell++;
              ImmersedBoundary.push_back(gridin[ idx3(xi, yi, zi, ldi1, ldi2) ]);

              cellVert[0] = (xi + 1) * dx - dx;
              cellVert[1] = (yi + 1) * dy - dy;
              cellVert[2] = (zi + 1) * dz - dz;

              cellVert[3] = (xi + 1) * dx;
              cellVert[4] = cellVert[1];
              cellVert[5] = cellVert[2];

              cellVert[6] = cellVert[3];
              cellVert[7] = (yi + 1) * dy;
              cellVert[8] = cellVert[2];

              cellVert[9] = cellVert[0];
              cellVert[10] = cellVert[7];
              cellVert[11] = cellVert[2];

              cellVert[12] = cellVert[0];
              cellVert[13] = cellVert[7];
              cellVert[14] = (zi + 1) * dz;

              cellVert[15] = cellVert[3];
              cellVert[16] = cellVert[7];
              cellVert[17] = cellVert[14];

              cellVert[18] = cellVert[3];
              cellVert[19] = cellVert[1];
              cellVert[20] = cellVert[14];

              cellVert[21] = cellVert[0];
              cellVert[22] = cellVert[1];
              cellVert[23] = cellVert[14];

              for (int pcount = 0; pcount < 8; pcount++)
              {
                nodeHold[0] = cellVert[ idx2( pcount, 0, 3 ) ];
                nodeHold[1] = cellVert[ idx2( pcount, 1, 3 ) ];
                nodeHold[2] = cellVert[ idx2( pcount, 2, 3 ) ];
                if ( !countCell ) // countCell = 0 -> first cell so no possible
                                  // node duplicates
                {
                  checkVert = -1;
                }
                else
                {
                  checkVert = isNear( nodeHold, Nodes, dx, dy, dz, nNodes, DIM );
                }
                if (checkVert == -1) // node is not a duplicate
                {
                  nNodes++;
                  Nodes.push_back(nodeHold[0]);
                  Nodes.push_back(nodeHold[1]);
                  Nodes.push_back(nodeHold[2]);
                  mv[ idx2( countCell, pcount, 8 ) ] = nNodes-1;
                }
                else // node is a duplicate
                {
                  mv[ idx2( countCell, pcount, 8 ) ] = checkVert;
                }
              }
            }
          }
        }
      }
      Nodes.resize(Nodes.size());

      // Next we compute cell centers for pressures
      PCellCenters.resize((numPCells * 3));
      for (int cl = 0; cl < numPCells; cl++)
      {
        PCellCenters[ idx2( cl, 0, 3 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 0, 8 ) ], 0, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 1, 8 ) ], 0, NodesLDI ) ]);
        PCellCenters[ idx2( cl, 1, 3 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 1, 8 ) ], 1, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 2, 8 ) ], 1, NodesLDI ) ]);
        PCellCenters[ idx2( cl, 2, 3 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 0, 8 ) ], 2, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 7, 8 ) ], 2, NodesLDI ) ]);
      }

      DOF.resize(4);

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, indendent of other
         velocity components.*/
      #pragma omp parallel sections
      {
        { // Final P computations
          DOF[0] = numPCells;
          PFaceConnectivity.resize((numPCells * 6));
          innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells );
          PCellWidths.resize((DOF[0] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[0]; cl++) {
            PCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            PCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
            PCellWidths[ idx2( cl, 2, CellWidthsLDI ) ] = dz;
          }
        } // End P section
        #pragma omp section
        { // U section
          int checkVertU;
          std::vector<double> nodeHoldU;
          nodeHoldU.resize(3);
          int countUCells = 0;
          double uStep = 0.5*dx;
          int maxUCells = nx * ny * nz + ny * nz;
          UCellCenters.reserve((maxUCells * 3));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++)
          {
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[0] = PCellCenters[ idx2( cl, 0, 3 ) ] - uStep;
            cellVert[1] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[2] = PCellCenters[ idx2( cl, 2, 3 ) ];

            cellVert[3] = PCellCenters[ idx2( cl, 0, 3 ) ] + uStep;
            cellVert[4] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[5] = PCellCenters[ idx2( cl, 2, 3 ) ];

            for (int pcount = 0; pcount < 2; pcount++)
            {
              // U Component
              nodeHoldU[0] = cellVert[ idx2( pcount, 0, 3 ) ];
              nodeHoldU[1] = cellVert[ idx2( pcount, 1, 3 ) ];
              nodeHoldU[2] = cellVert[ idx2( pcount, 2, 3 ) ];
              if ( !countUCells ) // First U node
              {
                checkVertU = -1;
              }
              else
              {
                checkVertU = isNear( nodeHoldU, UCellCenters, \
                                    dx, dy, dz, countUCells, DIM );
              }
              if (checkVertU == -1) // cell center location is not a duplicate
              {
                countUCells++;
                UCellCenters.push_back(nodeHoldU[0]);
                UCellCenters.push_back(nodeHoldU[1]);
                UCellCenters.push_back(nodeHoldU[2]);
                if (pcount == 0) {
                  UCellPressureNeighbor[ idx2( (countUCells-1), 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countUCells;
                }
                else if (pcount == 1) {
                  UCellPressureNeighbor[ idx2( (countUCells-1), 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countUCells;
                }
              }
              if (checkVertU != -1)
              {
                if (pcount == 0) {
                  UCellPressureNeighbor[ idx2( checkVertU, 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertU+1;
                }
                else if (pcount== 1) {
                  UCellPressureNeighbor[ idx2( checkVertU, 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellUNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertU+1;
                }
              }
            }
          }
          UCellCenters.resize(UCellCenters.size());
          DOF[1] = countUCells;
          UCellPressureNeighbor.resize((DOF[1] * 2));
          UFaceConnectivity.resize((countUCells * 6));
          innerFaceConnectivity( UFaceConnectivity, UCellCenters, \
                                 dx, dy, dz, countUCells );
          UInteriorCells.reserve(DOF[1]);
          UBoundaryCells.reserve(DOF[1]);
          int nbrsu;
          for (int cl = 0; cl < DOF[1]; cl++) {
            nbrsu = 0;
            for (int position = 0; position < 6; position++) {
              if (UFaceConnectivity[ idx2( cl, position, FaceConnectivityLDI ) ] != 0) {
                nbrsu++;
              }
            }
            if (nbrsu == 6) {
              UInteriorCells.push_back(cl);
            }
            else {
              UBoundaryCells.push_back(cl);
            }
          }
          UInteriorCells.resize(UInteriorCells.size());
          UBoundaryCells.resize(UBoundaryCells.size());
          UCellWidths.resize((DOF[1] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[1]; cl++) {
            UCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            UCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
            UCellWidths[ idx2( cl, 2, CellWidthsLDI ) ] = dz;
          }
        } // end U section
        #pragma omp section
        { // V section
          int checkVertV;
          std::vector<double> nodeHoldV;
          nodeHoldV.resize(3);
          int countVCells = 0;
          double vStep = 0.5*dy;
          int maxVCells = nx * ny * nz + nx * nz;
          VCellCenters.reserve((maxVCells * 3));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int cl = 0; cl < numPCells; cl++)
          {
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[6] = PCellCenters[ idx2( cl, 0, 3 ) ];
            cellVert[7] = PCellCenters[ idx2( cl, 1, 3 ) ] - vStep;
            cellVert[8] = PCellCenters[ idx2( cl, 2, 3 ) ];

            cellVert[9] = PCellCenters[ idx2( cl, 0, 3 ) ];
            cellVert[10] = PCellCenters[ idx2( cl, 1, 3 ) ] + vStep;
            cellVert[11] = PCellCenters[ idx2( cl, 2, 3 ) ];

            for (int pcount = 0; pcount < 2; pcount++)
            {
              // V Component
              nodeHoldV[0] = cellVert[ idx2( (pcount + 2), 0, 3 ) ];
              nodeHoldV[1] = cellVert[ idx2( (pcount + 2), 1, 3 ) ];
              nodeHoldV[2] = cellVert[ idx2( (pcount + 2), 2, 3 ) ];
              if ( !countVCells ) // First V node
              {
                checkVertV = -1;
              }
              else
              {
                checkVertV = isNear( nodeHoldV, VCellCenters, \
                                    dx, dy, dz, countVCells, DIM );
              }
              if (checkVertV == -1)
              {
                countVCells++;
                VCellCenters.push_back(nodeHoldV[0]);
                VCellCenters.push_back(nodeHoldV[1]);
                VCellCenters.push_back(nodeHoldV[2]);
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
                else if (pcount == 1) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
              }
              if (checkVertV != -1)
              {
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( checkVertV, 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
                else if (pcount== 1) {
                  VCellPressureNeighbor[ idx2( checkVertV, 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellVNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
              }
            }
          }
          UCellCenters.resize(UCellCenters.size());
          DOF[2] = countVCells;
          VCellPressureNeighbor.resize((DOF[2] * 2));
          VFaceConnectivity.resize((countVCells * 6));
          innerFaceConnectivity( VFaceConnectivity, VCellCenters, \
                                 dx, dy, dz, countVCells );
          VInteriorCells.reserve(DOF[2]);
          VBoundaryCells.reserve(DOF[2]);
          int nbrsv;
          for (int cl = 0; cl < DOF[2]; cl++) {
            nbrsv = 0;
            for (int position = 0; position < 6; position++) {
              if (VFaceConnectivity[ idx2( cl, position, FaceConnectivityLDI ) ] != 0) {
                nbrsv++;
              }
            }
            if (nbrsv == 6) {
              VInteriorCells.push_back(cl);
            }
            else {
              VBoundaryCells.push_back(cl);
            }
          }
          VInteriorCells.resize(VInteriorCells.size());
          VBoundaryCells.resize(VBoundaryCells.size());
          VCellWidths.resize((DOF[2] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[2]; cl++) {
            VCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            VCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
            VCellWidths[ idx2( cl, 2, CellWidthsLDI ) ] = dz;
          }
        } // end V Section
        #pragma omp section
        { // W Section
          std::vector<double> nodeHoldW;
          nodeHoldW.resize(3);
          int countWCells = 0;
          double wStep = 0.5*dz;
          int maxWCells = nx * ny * nz + nx * ny;
          WCellCenters.reserve((maxWCells * 3));
          PressureCellWNeighbor.resize((numPCells * 2));
          WCellPressureNeighbor.resize((maxWCells * 2));
          for (int cl = 0; cl < numPCells; cl++)
          {
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[12] = PCellCenters[ idx2( cl, 0, 3 ) ];
            cellVert[13] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[14] = PCellCenters[ idx2( cl, 2, 3 ) ] - wStep;

            cellVert[15] = PCellCenters[ idx2( cl, 0, 3 ) ];
            cellVert[16] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[17] = PCellCenters[ idx2( cl, 2, 3 ) ] + wStep;

            for (int pcount = 0; pcount < 2; pcount++)
            {
              int checkVertW;
              // W Component
              nodeHoldW[0] = cellVert[ idx2( (pcount + 4), 0, 3 ) ];
              nodeHoldW[1] = cellVert[ idx2( (pcount + 4), 1, 3 ) ];
              nodeHoldW[2] = cellVert[ idx2( (pcount + 4), 2, 3 ) ];
              if ( !countWCells ) // First W node
              {
                checkVertW = -1;
              }
              else
              {
                checkVertW = isNear( nodeHoldW, WCellCenters, \
                                    dx, dy, dz, countWCells, DIM );
              }
              if (checkVertW == -1)
              {
                countWCells++;
                WCellCenters.push_back(nodeHoldW[0]);
                WCellCenters.push_back(nodeHoldW[1]);
                WCellCenters.push_back(nodeHoldW[2]);
                if (pcount == 0) {
                  WCellPressureNeighbor[ idx2( (countWCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellWNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countWCells;
                }
                else if (pcount == 1) {
                  WCellPressureNeighbor[ idx2( (countWCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = cl+1;
                  PressureCellWNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countWCells;
                }
              }
              if (checkVertW != -1)
              {
                if (pcount == 0) {
                  WCellPressureNeighbor[ idx2( checkVertW, 1, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellWNeighbor[ idx2( cl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertW+1;
                }
                else if (pcount== 1) {
                  WCellPressureNeighbor[ idx2( checkVertW, 0, \
                               VelocityCellPressureNeighborLDI) ] = cl+1;
                  PressureCellWNeighbor[ idx2( cl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertW+1;
                }
              }
            }
          }
          WCellCenters.resize(WCellCenters.size());
          DOF[3] = countWCells;
          WCellPressureNeighbor.resize((DOF[3] * 2));
          WFaceConnectivity.resize((countWCells * 6));
          innerFaceConnectivity( WFaceConnectivity, WCellCenters, \
                                 dx, dy, dz, countWCells );
          WInteriorCells.reserve(DOF[3]);
          WBoundaryCells.reserve(DOF[3]);
          int nbrsw;
          for (int cl = 0; cl < DOF[3]; cl++) {
            nbrsw = 0;
            for (int position = 0; position < 6; position++) {
              if (WFaceConnectivity[ idx2( cl, position, FaceConnectivityLDI ) ] != 0) {
                nbrsw++;
              }
            }
            if (nbrsw == 6) {
              WInteriorCells.push_back(cl);
            }
            else {
              WBoundaryCells.push_back(cl);
            }
          }
          WInteriorCells.resize(WInteriorCells.size());
          WBoundaryCells.resize(WBoundaryCells.size());
          WCellWidths.resize((DOF[3] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[3]; cl++) {
            WCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            WCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
            WCellWidths[ idx2( cl, 2, CellWidthsLDI ) ] = dz;
          }
        }
      }
      break;
    }
  } // End of dimension switch
}

int FluidMesh::isNear( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                       double dx, double dy, double dz, int nNodes, int DIM )
{
  int cl = nNodes;
  int prox = -1;

  double xr = 0;
  double yr = 0;
  double zr = 0;
  double epsx = 0.2 * dx;
  double epsy = 0.2 * dy;
  double epsz = 0.2 * dz;

  switch ( DIM )
  {
    case 2 :
    {
      do
      {
        cl--;
        xr = fabs( Vector1[0] - Vector2[ idx2( cl, 0, 2 ) ]);
        yr = fabs( Vector1[1] - Vector2[ idx2( cl, 1, 2 ) ]);
        if (xr < epsx)
        {
          if (yr < epsy)
          {
            prox = cl;
          }
        }
      } while (prox == -1 && cl > 0);
      break;
    }
    default :
    {
      do
      {
        cl--;
        xr = fabs( Vector1[0] - Vector2[ idx2( cl, 0, 3 ) ]);
        yr = fabs( Vector1[1] - Vector2[ idx2( cl, 1, 3 ) ]);
        zr = fabs( Vector1[2] - Vector2[ idx2( cl, 2, 3 ) ]);
        if (xr < epsx)
        {
          if (yr < epsy)
          {
            if (zr < epsz)
            {
              prox = cl;
            }
          }
        }
      } while (prox == -1 && cl > 0);
      break;
    }
  }
  return prox;
}
// Function to compute cell face connectivity information
void FluidMesh::innerFaceConnectivity( \
                std::vector<unsigned long>& ComponentFaceConnectivity, \
                std::vector<double> ComponentCellCenters, \
                double dx, double dy, double dz, int nCells )
{
  int incr, numNeighbors, nl;
  double xr, yr, zr, epsx, epsy, epsz, xtol, ytol, ztol;

  epsx = 0.2 * dx;
  epsy = 0.2 * dy;
  epsz = 0.2 * dz;
  xtol = 1.2 * dx;
  ytol = 1.2 * dy;
  ztol = 1.2 * dz;

  switch ( DIM )
  {
    case 2 :
    {
      double cccl0, cccl1;
      int int0, int1, int2, int3;
      for (int cl = 0; cl < nCells; cl++) {
        numNeighbors = 0;
        nl = 0;
        cccl0 = ComponentCellCenters[ idx2( cl, 0, 2 ) ];
        cccl1 = ComponentCellCenters[ idx2( cl, 1, 2 ) ];
        int0 = 0;
        int1 = 0;
        int2 = 0;
        int3 = 0;
        do
        {
          incr = nl + 1;
          xr = ComponentCellCenters[ idx2(nl, 0, 2) ] \
               - cccl0;
          yr = ComponentCellCenters[ idx2(nl, 1, 2) ] \
               - cccl1;
          if (fabs(xr) < epsx) {
            if (fabs(yr) < ytol) {
              if (yr < 0)
              {
                int0 = incr;
                numNeighbors++;
              }
              else if (yr > 0)
              {
                int2 = incr;
                numNeighbors++;
              }
            }
          }
          else if (fabs(yr) < epsy) {
            if (fabs(xr) < xtol) {
              if (xr < 0)
              {
                int3 = incr;
                numNeighbors++;
              }
              else if (xr > 0)
              {
                int1 = incr;
                numNeighbors++;
              }
            }
          }
          nl++;
        } while (nl < nCells && numNeighbors < 4);
        ComponentFaceConnectivity[ idx2(cl, 0, 4) ] = int0;
        ComponentFaceConnectivity[ idx2(cl, 1, 4) ] = int1;
        ComponentFaceConnectivity[ idx2(cl, 2, 4) ] = int2;
        ComponentFaceConnectivity[ idx2(cl, 3, 4) ] = int3;
      }
      break;
    }
    default :
    {
      double cccl0, cccl1, cccl2;
      int int0, int1, int2, int3, int4, int5;
      for (int cl = 0; cl < nCells; cl++) {
        numNeighbors = 0;
        cccl0 = ComponentCellCenters[ idx2( cl, 0, 3 ) ];
        cccl1 = ComponentCellCenters[ idx2( cl, 1, 3 ) ];
        cccl2 = ComponentCellCenters[ idx2( cl, 2, 3 ) ];
        int0 = 0;
        int1 = 0;
        int2 = 0;
        int3 = 0;
        int4 = 0;
        int5 = 0;
        nl = 0;
        do
        {
          incr = nl + 1;
          xr = ComponentCellCenters[ idx2(nl, 0, 3) ] \
               - cccl0;
          yr = ComponentCellCenters[ idx2(nl, 1, 3) ] \
               - cccl1;
          zr = ComponentCellCenters[ idx2(nl, 2, 3) ] \
               - cccl2;

          if (fabs(xr) < epsx) {
            if (fabs(yr) < epsy) {
              if (fabs(zr) < ztol) {
                if (zr < 0)
                {
                  int0 = incr;
                  numNeighbors++;
                }
                else if (zr > 0)
                {
                  int2 = incr;
                  numNeighbors++;
                }
              }
            }
            if (fabs(zr) < epsz) {
              if (fabs(yr) < ytol) {
                if (yr < 0)
                {
                  int5 = incr;
                  numNeighbors++;
                }
                else if (yr > 0)
                {
                  int4 = incr;
                  numNeighbors++;
                }
              }
            }
          }
          else if (fabs(yr) < epsy) {
            if (fabs(zr) < epsz) {
              if (fabs(xr) < xtol) {
                if (xr < 0)
                {
                  int3 = incr;
                  numNeighbors++;
                }
                else if (xr > 0)
                {
                  int1 = incr;
                  numNeighbors++;
                }
              }
            }
          }
          nl++;
        } while (nl < nCells && numNeighbors < 6);
        ComponentFaceConnectivity[ idx2(cl, 0, 6) ] = int0;
        ComponentFaceConnectivity[ idx2(cl, 1, 6) ] = int1;
        ComponentFaceConnectivity[ idx2(cl, 2, 6) ] = int2;
        ComponentFaceConnectivity[ idx2(cl, 3, 6) ] = int3;
        ComponentFaceConnectivity[ idx2(cl, 4, 6) ] = int4;
        ComponentFaceConnectivity[ idx2(cl, 5, 6) ] = int5;
       }
      break;
    }
  }
}
// Compute total DOF
int FluidMesh::TotalDOF( void )
{
  int outVal = 0;
  switch ( DIM )
  {
    case 2 :
      outVal =  DOF[0] + DOF[1] + DOF[2];
      break;
    case 3 :
      outVal = DOF[0] + DOF[1] + DOF[2] + DOF[3];
      break;
  }
  return outVal;
}
// Compute DOF for velocities
int FluidMesh::VelocityDOF( void )
{
  int outVal = 0;
  switch ( DIM )
  {
    case 2 :
      outVal = DOF[1] + DOF[2];
      break;
    case 3 :
      outVal = DOF[1] + DOF[2] + DOF[3];
      break;
  }
  return outVal;
}
// Compute maximum possible nonzero entries in linear system
int FluidMesh::MaxNonZero( void )
{
  int outVal = 0;
  switch ( DIM )
  {
    case 2 :
      outVal = 4 * DOF[0] + 8 * DOF[1] + 8 * DOF[2];
      break;
    case 3 :
      outVal = 6 * DOF[0] + 10 * DOF[1] + 10 * DOF[2] + 10 * DOF[3];
      break;
  }
  return outVal;
}
