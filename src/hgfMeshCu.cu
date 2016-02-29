#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <cuda_runtime.h>
#include <fstream>
#include <boost/filesystem.hpp>

#include "hgfMeshCu.hpp"
#include "hgf.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"
#include "hgfAuxTools.hpp"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

/* Define a 2d -> 1d array index, uses row major ordering */
#define idx2(i, j, ldi) ((i * ldi) + j)
/* Define a 3d -> 1d array index, uses row major ordering */
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

__global__ void ifcKernel2D ( unsigned long *d_CFC, const double *d_CCC, \
                              double epsx, double epsy, \
                              double xtol, double ytol, int nCells )
{
  int cl = blockIdx.x * blockDim.x + threadIdx.x;
  int incr;
  double xr, yr;
  int numNeighbors = 0;
  int nl = 0;
  double cccl0 = d_CCC[ idx2( cl, 0, 2 ) ];
  double cccl1 = d_CCC[ idx2( cl, 1, 2 ) ];
  int int0 = 0;
  int int1 = 0;
  int int2 = 0;
  int int3 = 0;
  do
  {
    incr = nl + 1;
    xr = d_CCC[ idx2(nl, 0, 2) ] \
         - cccl0;
    yr = d_CCC[ idx2(nl, 1, 2) ] \
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
  d_CFC[ idx2(cl, 0, 4) ] = int0;
  d_CFC[ idx2(cl, 1, 4) ] = int1;
  d_CFC[ idx2(cl, 2, 4) ] = int2;
  d_CFC[ idx2(cl, 3, 4) ] = int3;
}
__global__ void ifcKernel3D( unsigned long *d_CFC, const double *d_CCC, \
                             double epsx, double epsy, double epsz, \
                             double xtol, double ytol, double ztol, \
                             int nCells )
{
  int cl = blockIdx.x * blockDim.x + threadIdx.x;
  int incr;
  double xr, yr, zr;
  int numNeighbors = 0;
  double cccl0 = d_CCC[ idx2( cl, 0, 3 ) ];
  double cccl1 = d_CCC[ idx2( cl, 1, 3 ) ];
  double cccl2 = d_CCC[ idx2( cl, 2, 3 ) ];
  int int0 = 0;
  int int1 = 0;
  int int2 = 0;
  int int3 = 0;
  int int4 = 0;
  int int5 = 0;
  int nl = 0;
  do
  {
    incr = nl + 1;
    xr = d_CCC[ idx2(nl, 0, 3) ] \
         - cccl0;
    yr = d_CCC[ idx2(nl, 1, 3) ] \
         - cccl1;
    zr = d_CCC[ idx2(nl, 2, 3) ] \
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
  d_CFC[ idx2(cl, 0, 6) ] = int0;
  d_CFC[ idx2(cl, 1, 6) ] = int1;
  d_CFC[ idx2(cl, 2, 6) ] = int2;
  d_CFC[ idx2(cl, 3, 6) ] = int3;
  d_CFC[ idx2(cl, 4, 6) ] = int4;
  d_CFC[ idx2(cl, 5, 6) ] = int5;
}
void
MeshSubdivide( unsigned long *gridin, int ldi1, int ldi2, \
               int nx, int ny, int nz, \
               double length, double width, double height, \
               int MX, int MY, int MZ, \
               std::vector< std::vector<unsigned long> >& slices, \
               std::vector<double>& lengths, \
               std::vector<double>& widths, \
               std::vector<double>& heights, \
               std::vector<int>& nxs, \
               std::vector<int>& nys, \
               std::vector<int>& nzs )
{
  switch ( nz )
  {
    case 0 :
    {
      int nSubDomains = 0;
      int xStart, nxRemainder, yStart, nyRemainder, nyG, nxG, xCount, yCount;
      slices.resize( MX * MY );
      lengths.resize( MX * MY );
      widths.resize( MX * MY );
      heights.resize( MX * MY );
      nxs.resize( MX * MY );
      nys.resize( MX * MY );
      nzs.resize( MX * MY );
      yCount = 0;
      yStart = 0;
      nyRemainder = ny;
      for (int yy = 0; yy < MY; yy++) {
        xStart = 0;
        nxRemainder = nx;
        xCount = 0;
        nyG = (int)(round(((double)nyRemainder)/(MY-yCount)));
        for (int xx = 0; xx < MX; xx++) {
          nSubDomains++;
          nxG = (int)(round(((double)nxRemainder)/(MX-xCount)));
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
      break;
    }
    default :
    {
      int nSubDomains = 0;
      int xStart, nxRemainder, yStart, nyRemainder, zStart, nzRemainder, nxG, nyG, nzG, xCount, yCount, zCount;
      slices.resize( MX * MY * MZ );
      lengths.resize( MX * MY * MZ );
      widths.resize( MY * MY * MZ );
      heights.resize( MX * MY * MZ );
      nxs.resize( MX * MY * MZ );
      nys.resize( MX * MY * MZ );
      nzs.resize( MX * MY * MZ );
      zCount = 0;
      zStart = 0;
      nzRemainder = nz;
      for (int zz = 0; zz < MZ; zz++) {
        yStart = 0;
        nyRemainder = ny;
        yCount = 0;
        nzG = (int)(round(((double)nzRemainder)/(MZ-zCount)));
        for (int yy = 0; yy < MY; yy++) {
          xStart = 0;
          nxRemainder = nx;
          xCount = 0;
          nyG = (int)(round(((double)nyRemainder)/(MY-yCount)));
          for (int xx = 0; xx < MX; xx++) {
            nSubDomains++;
            nxG = (int)(round(((double)nxRemainder)/(MX-xCount)));
            for (int cx = 0; cx < nxG; cx++) {
              for (int cy = 0; cy < nyG; cy++) {
                for (int cz = 0; cz < nzG; cz++) {
                  slices[ nSubDomains-1 ].push_back( gridin[ idx3( (cx+xStart), (cy+yStart), (cz+zStart), ldi1, ldi2 ) ] );
                }
              }
            }
            lengths[ nSubDomains-1 ] = length * ((double)nxG / nx);
            widths[ nSubDomains-1 ] = width * ((double)nyG / ny);
            heights[ nSubDomains-1 ] = height * ((double)nzG / nz);
            nxs[ nSubDomains-1 ] = nxG;
            nys[ nSubDomains-1 ] = nyG;
            nzs[ nSubDomains-1 ] = nzG;
            xStart = xStart + nxG;
            nxRemainder = nx - xStart;
            xCount++;
          }
          yStart = yStart + nyG;
          nyRemainder = ny - yStart;
          yCount++;
        }
        zStart = zStart + nzG;
        nzRemainder = nz - zStart;
        zCount++;
      }
      break;
    }
  }
}
// Function to compute cell face connectivity information
void innerFaceConnectivity( \
       std::vector<unsigned long>& ComponentFaceConnectivity, \
       const std::vector<double>& ComponentCellCenters, \
       double dx, double dy, double dz, int nCells, int DIM )
{

  double epsx = 0.2 * dx;
  double epsy = 0.2 * dy;
  double epsz = 0.2 * dz;
  double xtol = 1.2 * dx;
  double ytol = 1.2 * dy;
  double ztol = 1.2 * dz;

  // initialize device memory for faceconnectivity and cellcenters
  unsigned long *d_CFC = NULL;
  double *d_CCC = NULL;

  gpuErrchk( cudaMalloc( (void **)&d_CFC, ComponentFaceConnectivity.size() * sizeof(unsigned long) ) );

  gpuErrchk( cudaMalloc( (void **)&d_CCC, ComponentCellCenters.size() * sizeof(double) ) );

  // copy cell centers data from host to device
  gpuErrchk( cudaMemcpy( d_CCC, ComponentCellCenters.data(), \
    ComponentCellCenters.size() * sizeof(double), cudaMemcpyHostToDevice ) );

  int blockSize;
  int minGridSize;
  int gridSize;

  // compute
  if (DIM == 2)
  {
    gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, ifcKernel2D, 0, nCells ) );
    gridSize = (nCells + blockSize - 1) / blockSize;
    ifcKernel2D<<< gridSize, blockSize >>>( d_CFC, d_CCC, epsx, epsy, xtol, ytol, nCells );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
  }
  else if (DIM == 3)
  {
    gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, ifcKernel3D, 0, nCells ) );
    gridSize = (nCells + blockSize - 1) / blockSize;
    ifcKernel3D<<< gridSize, blockSize >>>( d_CFC, d_CCC, epsx, epsy, epsz, xtol, ytol, ztol, nCells );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
  }

  // copy results back to host
  gpuErrchk( cudaMemcpy( ComponentFaceConnectivity.data(), d_CFC, \
    ComponentFaceConnectivity.size() * sizeof(unsigned long), \
    cudaMemcpyDeviceToHost ) );

  // free device memory
  gpuErrchk( cudaFree( d_CFC ) );
  gpuErrchk( cudaFree( d_CCC ) );

}
// Function to construct the mesh from voxel array input
void FluidMesh::BuildUniformMesh( const ProbParam& Par )
{
  int checkVert;
  int maxPNodes;
  int numPCells = 0;
  int numVoid = 0;
  xLim[0] = 0;
  xLim[1] = Par.length;
  yLim[0] = 0;
  yLim[1] = Par.width;
  zLim[0] = 0;
  zLim[1] = Par.height;

  switch ( Par.nz ) {
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
      NX = Par.nx;
      NY = Par.ny;
      NZ = Par.nz;

      // Max pressure nodes
      std::vector<double> nodeHold;
      nodeHold.resize(2);
      for (int yi = 0; yi < Par.ny; yi++) {
        for (int xi = 0; xi < Par.nx; xi++) {
          if (Par.gridin[ idx2(yi, xi, Par.nx) ] != 1) {
            numPCells++;
            if (Par.gridin[ idx2(yi, xi, Par.nx) ] == 0) numVoid++;
          }
        }
      }
      maxPNodes = numPCells * 4;
      int nNodes = 0;
      double dx = Par.length / Par.nx;
      double dy = Par.width / Par.ny;
      double dz = 0;
      int countCell = -1;
      mv.resize( (numPCells * 4) );
      double cellVert [ 8 ];
      porosity = numVoid/(double)(Par.nx * Par.ny);

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((Par.nx * Par.ny));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((maxPNodes * 2));

      for (int yi = 0; yi < Par.ny; yi++) {
        for (int xi = 0; xi < Par.nx; xi++) {
          FullGrid.push_back(Par.gridin[ idx2(yi, xi, Par.nx) ]);
          if (Par.gridin[ idx2( yi, xi, Par.nx ) ] != 1) {
            countCell++;
            ImmersedBoundary.push_back(Par.gridin[ idx2(yi, xi, Par.nx) ]);
            cellVert[0] = (xi + 1) * dx - dx;
            cellVert[1] = (yi + 1) * dy - dy;
            cellVert[2] = (xi + 1) * dx;
            cellVert[3] = cellVert[1];
            cellVert[4] = cellVert[2];
            cellVert[5] = (yi + 1) * dy;
            cellVert[6] = cellVert[0];
            cellVert[7] = cellVert[5];

            for (int pcount = 0; pcount < 4; pcount++) {
              nodeHold[0] = cellVert[ idx2( pcount, 0, 2 ) ];
              nodeHold[1] = cellVert[ idx2( pcount, 1, 2 ) ];
              if ( !countCell ) { // countCell = 0 -> first cell so no possible
                                // node duplicates
                checkVert = -1;
              }
              else {
                checkVert = isNear2d( nodeHold, Nodes, dx, dy, dz, nNodes );
              }
              if (checkVert == -1) { // node is not a duplicate
                nNodes++;
                Nodes.push_back(nodeHold[0]);
                Nodes.push_back(nodeHold[1]);
                mv[ idx2( countCell, pcount, 4 ) ] = nNodes-1;
              }
              else { // node is a duplicate
                mv[ idx2( countCell, pcount, 4 ) ] = checkVert;
              }
            }
          }
        }
      }
      Nodes.resize(Nodes.size());

      // Next we compute cell centers for pressures
      PCellCenters.resize((numPCells * 2));
      for (int cl = 0; cl < numPCells; cl++) {
        PCellCenters[ idx2( cl, 0, 2 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 0, 4 ) ], 0, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 1, 4) ], 0, NodesLDI ) ]);
        PCellCenters[ idx2( cl, 1, 2 ) ] = 0.5 \
          * (Nodes[ idx2( mv[ idx2( cl, 1, 4 ) ], 1, NodesLDI ) ] \
             + Nodes[ idx2( mv[ idx2( cl, 2, 4) ], 1, NodesLDI ) ]);
      }

      DOF.resize(3);
      sortPV();

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, independent
         of other velocity components */
      #pragma omp parallel sections
      {
        { // Final P computations
          DOF[0] = numPCells;
          PFaceConnectivity.resize((numPCells * 4));
          innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells, DIM );
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
          int maxUCells = Par.nx * Par.ny + Par.ny;
          UCellCenters.reserve((maxUCells * 2));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            cellVert[0] = PCellCenters[ idx2( cl, 0, 2 ) ] - uStep;
            cellVert[1] = PCellCenters[ idx2( cl, 1, 2 ) ];
            cellVert[2] = PCellCenters[ idx2( cl, 0, 2 ) ] + uStep;
            cellVert[3] = PCellCenters[ idx2( cl, 1, 2 ) ];
            for (int pcount = 0; pcount < 2; pcount++) {
              // U Component
              nodeHoldU[0] = cellVert[ idx2( pcount, 0, 2 ) ];
              nodeHoldU[1] = cellVert[ idx2( pcount, 1, 2 ) ];
              if ( !countUCells ) {// First U node
                checkVertU = -1;
              }
              else {
                checkVertU = isNear2d( nodeHoldU, UCellCenters, \
                                     dx, dy, dz, countUCells );
              }
              if (checkVertU == -1) { // cell center location is not a duplicate
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
              if (checkVertU != -1) {
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
                                 dx, dy, dz, countUCells, DIM );
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
          int maxVCells = Par.nx * Par.ny + Par.nx;
          int yl;
          VCellCenters.reserve((maxVCells * 2));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int xl = 0; xl < numPCells; xl++) {
            yl = PresListByY[ xl ];
            cellVert[4] = PCellCenters[ idx2( yl, 0, 2 ) ];
            cellVert[5] = PCellCenters[ idx2( yl, 1, 2 ) ] - vStep;
            cellVert[6] = PCellCenters[ idx2( yl, 0, 2 ) ];
            cellVert[7] = PCellCenters[ idx2( yl, 1, 2 ) ] + vStep;
            for (int pcount = 0; pcount < 2; pcount++) {
              // V Component
              nodeHoldV[0] = cellVert[ idx2( (pcount + 2), 0, 2 ) ];
              nodeHoldV[1] = cellVert[ idx2( (pcount + 2), 1, 2 ) ];
              if ( !countVCells ) { // First V node
                checkVertV = -1;
              }
              else {
                checkVertV = isNear2d( nodeHoldV, VCellCenters, \
                                    dx, dy, dz, countVCells );
              }
              if (checkVertV == -1) {
                countVCells++;
                VCellCenters.push_back(nodeHoldV[0]);
                VCellCenters.push_back(nodeHoldV[1]);
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
                else if (pcount == 1) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
              }
              if (checkVertV != -1) {
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( checkVertV, 1, \
                               VelocityCellPressureNeighborLDI) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
                else if (pcount== 1) {
                  VCellPressureNeighbor[ idx2( checkVertV, 0, \
                               VelocityCellPressureNeighborLDI) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
              }
            }
          }
          VCellCenters.resize(VCellCenters.size());
          DOF[2] = countVCells;
          VCellPressureNeighbor.resize((DOF[2] * 2));
          VFaceConnectivity.resize((countVCells * 4));
          innerFaceConnectivity( VFaceConnectivity, VCellCenters, \
                                 dx, dy, dz, countVCells, DIM );
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
      for (int zi = 0; zi < Par.nz; zi++) {
        for (int yi = 0; yi < Par.ny; yi++) {
          for (int xi = 0; xi < Par.nx; xi++) {
            if (Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ] != 1) {
              numPCells++;
              if (Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx ) ] == 0) numVoid++;
            }
          }
        }
      }
      maxPNodes = numPCells * 8;

      int nNodes = 0;
      double dx = Par.length / Par.nx;
      double dy = Par.width / Par.ny;
      double dz = Par.height / Par.nz;
      int countCell = -1;
      mv.resize( (numPCells * 8) );
      double cellVert [ 24 ];
      porosity = numVoid /(double)(Par.nx * Par.ny * Par.nz);

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((Par.nx * Par.ny * Par.nz));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((maxPNodes * 3));
      for (int zi = 0; zi < Par.nz; zi++) {
        for (int yi = 0; yi < Par.ny; yi++) {
          for (int xi = 0; xi < Par.nx; xi++) {
            FullGrid.push_back(Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ]);
            if (Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] != 1) {
              countCell++;
              ImmersedBoundary.push_back(Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ]);

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

              for (int pcount = 0; pcount < 8; pcount++) {
                nodeHold[0] = cellVert[ idx2( pcount, 0, 3 ) ];
                nodeHold[1] = cellVert[ idx2( pcount, 1, 3 ) ];
                nodeHold[2] = cellVert[ idx2( pcount, 2, 3 ) ];
                if ( !countCell ) { // countCell = 0 -> first cell so no possible
                                  // node duplicates
                  checkVert = -1;
                }
                else {
                  checkVert = isNear3d( nodeHold, Nodes, dx, dy, dz, nNodes );
                }
                if (checkVert == -1) { // node is not a duplicate
                  nNodes++;
                  Nodes.push_back(nodeHold[0]);
                  Nodes.push_back(nodeHold[1]);
                  Nodes.push_back(nodeHold[2]);
                  mv[ idx2( countCell, pcount, 8 ) ] = nNodes-1;
                }
                else { // node is a duplicate
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
      for (int cl = 0; cl < numPCells; cl++) {
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
      sortPV();

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, indendent of other
         velocity components.*/
      #pragma omp parallel sections
      {
        { // Final P computations
          DOF[0] = numPCells;
          PFaceConnectivity.resize((numPCells * 6));
          innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells, DIM );
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
          int maxUCells = Par.nx * Par.ny * Par.nz + Par.ny * Par.nz;
          UCellCenters.reserve((maxUCells * 3));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[0] = PCellCenters[ idx2( cl, 0, 3 ) ] - uStep;
            cellVert[1] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[2] = PCellCenters[ idx2( cl, 2, 3 ) ];

            cellVert[3] = PCellCenters[ idx2( cl, 0, 3 ) ] + uStep;
            cellVert[4] = PCellCenters[ idx2( cl, 1, 3 ) ];
            cellVert[5] = PCellCenters[ idx2( cl, 2, 3 ) ];

            for (int pcount = 0; pcount < 2; pcount++) {
              // U Component
              nodeHoldU[0] = cellVert[ idx2( pcount, 0, 3 ) ];
              nodeHoldU[1] = cellVert[ idx2( pcount, 1, 3 ) ];
              nodeHoldU[2] = cellVert[ idx2( pcount, 2, 3 ) ];
              if ( !countUCells ) { // First U node
                checkVertU = -1;
              }
              else {
                checkVertU = isNear3d( nodeHoldU, UCellCenters, \
                                    dx, dy, dz, countUCells );
              }
              if (checkVertU == -1) { // cell center location is not a duplicate
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
              if (checkVertU != -1) {
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
                                 dx, dy, dz, countUCells, DIM );
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
          int maxVCells = Par.nx * Par.ny * Par.nz + Par.nx * Par.nz;
          int yl;
          VCellCenters.reserve((maxVCells * 3));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            yl = PresListByY[ cl ];
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[6] = PCellCenters[ idx2( yl, 0, 3 ) ];
            cellVert[7] = PCellCenters[ idx2( yl, 1, 3 ) ] - vStep;
            cellVert[8] = PCellCenters[ idx2( yl, 2, 3 ) ];

            cellVert[9] = PCellCenters[ idx2( yl, 0, 3 ) ];
            cellVert[10] = PCellCenters[ idx2( yl, 1, 3 ) ] + vStep;
            cellVert[11] = PCellCenters[ idx2( yl, 2, 3 ) ];

            for (int pcount = 0; pcount < 2; pcount++) {
              // V Component
              nodeHoldV[0] = cellVert[ idx2( (pcount + 2), 0, 3 ) ];
              nodeHoldV[1] = cellVert[ idx2( (pcount + 2), 1, 3 ) ];
              nodeHoldV[2] = cellVert[ idx2( (pcount + 2), 2, 3 ) ];
              if ( !countVCells ) { // First V node
                checkVertV = -1;
              }
              else {
                checkVertV = isNear3d( nodeHoldV, VCellCenters, \
                                    dx, dy, dz, countVCells );
              }
              if (checkVertV == -1) {
                countVCells++;
                VCellCenters.push_back(nodeHoldV[0]);
                VCellCenters.push_back(nodeHoldV[1]);
                VCellCenters.push_back(nodeHoldV[2]);
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
                else if (pcount == 1) {
                  VCellPressureNeighbor[ idx2( (countVCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countVCells;
                }
              }
              if (checkVertV != -1) {
                if (pcount == 0) {
                  VCellPressureNeighbor[ idx2( checkVertV, 1, \
                               VelocityCellPressureNeighborLDI) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
                else if (pcount== 1) {
                  VCellPressureNeighbor[ idx2( checkVertV, 0, \
                               VelocityCellPressureNeighborLDI) ] = yl+1;
                  PressureCellVNeighbor[ idx2( yl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertV+1;
                }
              }
            }
          }
          VCellCenters.resize(VCellCenters.size());
          DOF[2] = countVCells;
          VCellPressureNeighbor.resize((DOF[2] * 2));
          VFaceConnectivity.resize((countVCells * 6));
          innerFaceConnectivity( VFaceConnectivity, VCellCenters, \
                                 dx, dy, dz, countVCells, DIM );
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
          int maxWCells = Par.nx * Par.ny * Par.nz + Par.nx * Par.ny;
          int zl;
          WCellCenters.reserve((maxWCells * 3));
          PressureCellWNeighbor.resize((numPCells * 2));
          WCellPressureNeighbor.resize((maxWCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            zl = PresListByZ[ cl ];
            // Each 3 value block is a 'row' of the 2d celLVert array
            cellVert[12] = PCellCenters[ idx2( zl, 0, 3 ) ];
            cellVert[13] = PCellCenters[ idx2( zl, 1, 3 ) ];
            cellVert[14] = PCellCenters[ idx2( zl, 2, 3 ) ] - wStep;

            cellVert[15] = PCellCenters[ idx2( zl, 0, 3 ) ];
            cellVert[16] = PCellCenters[ idx2( zl, 1, 3 ) ];
            cellVert[17] = PCellCenters[ idx2( zl, 2, 3 ) ] + wStep;

            for (int pcount = 0; pcount < 2; pcount++) {
              int checkVertW;
              // W Component
              nodeHoldW[0] = cellVert[ idx2( (pcount + 4), 0, 3 ) ];
              nodeHoldW[1] = cellVert[ idx2( (pcount + 4), 1, 3 ) ];
              nodeHoldW[2] = cellVert[ idx2( (pcount + 4), 2, 3 ) ];
              if ( !countWCells ) { // First W node
                checkVertW = -1;
              }
              else {
                checkVertW = isNear3d( nodeHoldW, WCellCenters, \
                                    dx, dy, dz, countWCells );
              }
              if (checkVertW == -1) {
                countWCells++;
                WCellCenters.push_back(nodeHoldW[0]);
                WCellCenters.push_back(nodeHoldW[1]);
                WCellCenters.push_back(nodeHoldW[2]);
                if (pcount == 0) {
                  WCellPressureNeighbor[ idx2( (countWCells-1), 1, \
                               VelocityCellPressureNeighborLDI ) ] = zl+1;
                  PressureCellWNeighbor[ idx2( zl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = countWCells;
                }
                else if (pcount == 1) {
                  WCellPressureNeighbor[ idx2( (countWCells-1), 0, \
                               VelocityCellPressureNeighborLDI ) ] = zl+1;
                  PressureCellWNeighbor[ idx2( zl, 1, \
                               PressureCellVelocityNeighborLDI ) ] = countWCells;
                }
              }
              if (checkVertW != -1) {
                if (pcount == 0) {
                  WCellPressureNeighbor[ idx2( checkVertW, 1, \
                               VelocityCellPressureNeighborLDI) ] = zl+1;
                  PressureCellWNeighbor[ idx2( zl, 0, \
                               PressureCellVelocityNeighborLDI ) ] = checkVertW+1;
                }
                else if (pcount== 1) {
                  WCellPressureNeighbor[ idx2( checkVertW, 0, \
                               VelocityCellPressureNeighborLDI) ] = zl+1;
                  PressureCellWNeighbor[ idx2( zl, 1, \
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
                                 dx, dy, dz, countWCells, DIM );
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
  TotalDOF();
  MaxNonZero();
}
// Function to find duplicate node in 2d
int FluidMesh::isNear2d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                         double dx, double dy, double dz, int nNodes )
{
  int cl = nNodes;
  int prox = -1;

  double xr = 0;
  double yr = 0;
  double epsx = 0.2 * dx;
  double epsy = 0.2 * dy;

  do
  {
    cl--;
    xr = fabs( Vector1[0] - Vector2[ idx2( cl, 0, 2 ) ]);
    yr = fabs( Vector1[1] - Vector2[ idx2( cl, 1, 2 ) ]);
    if (xr < epsx) {
      if (yr < epsy) {
        prox = cl;
      }
    }
  } while (prox == -1 && cl > 0);
  return prox;
}
// Function to find duplicate node in 3d
int FluidMesh::isNear3d( std::vector<double>& Vector1, std::vector<double>& Vector2, \
                         double dx, double dy, double dz, int nNodes )
{
  int cl = nNodes;
  int prox = -1;

  double xr = 0;
  double yr = 0;
  double zr = 0;
  double epsx = 0.2 * dx;
  double epsy = 0.2 * dy;
  double epsz = 0.2 * dz;

  do
  {
    cl--;
    xr = fabs( Vector1[0] - Vector2[ idx2( cl, 0, 3 ) ]);
    yr = fabs( Vector1[1] - Vector2[ idx2( cl, 1, 3 ) ]);
    zr = fabs( Vector1[2] - Vector2[ idx2( cl, 2, 3 ) ]);
    if (xr < epsx) {
      if (yr < epsy) {
        if (zr < epsz) {
          prox = cl;
        }
      }
    }
  } while (prox == -1 && cl > 0);
  return prox;
}
// Compute total DOF
void FluidMesh::TotalDOF( void )
{
  switch ( DIM ) {
    case 2 :
      dofTotal =  DOF[0] + DOF[1] + DOF[2];
      break;
    case 3 :
      dofTotal = DOF[0] + DOF[1] + DOF[2] + DOF[3];
      break;
  }
}
// Compute DOF for velocities
int FluidMesh::VelocityDOF( void )
{
  int outVal = 0;
  switch ( DIM ) {
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
void FluidMesh::MaxNonZero( void )
{
  switch ( DIM ) {
    case 2 :
      maxNNZ = 4 * DOF[0] + 8 * DOF[1] + 8 * DOF[2];
      break;
    case 3 :
      maxNNZ = 6 * DOF[0] + 10 * DOF[1] + 10 * DOF[2] + 10 * DOF[3];
      break;
  }
}
// create sorted pressure index for y grid
void FluidMesh::sortPV( void )
{
  switch ( DIM ) {
    case 2 :
    {
      std::vector< sortStruc2 > pYtrans2( (PCellCenters.size()/2) );
      PresListByY.reserve( PCellCenters.size()/2 );
      for (unsigned long cl = 0; cl < (PCellCenters.size()/2); cl++) {
        pYtrans2[cl].xx = PCellCenters[ idx2( cl, 0, 2 ) ];
        pYtrans2[cl].yy = PCellCenters[ idx2( cl, 1, 2 ) ];
        pYtrans2[cl].ind = cl;
      }
      std::sort(pYtrans2.begin(), pYtrans2.end(), byXbyY());
      for (unsigned long cl = 0; cl < (PCellCenters.size()/2); cl++) {
        PresListByY.push_back( pYtrans2[cl].ind );
      }
      break;
    }
    case 3 :
    {
      std::vector< sortStruc3 > pYtrans3( (PCellCenters.size()/3) );
      std::vector< sortStruc3 > pZtrans3( (PCellCenters.size()/3) );
      PresListByY.reserve( PCellCenters.size()/3 );
      PresListByZ.reserve( PCellCenters.size()/3 );
      for (unsigned long cl = 0; cl < (PCellCenters.size()/3); cl++) {
        pYtrans3[cl].xx = PCellCenters[ idx2( cl, 0, 3 ) ];
        pYtrans3[cl].yy = PCellCenters[ idx2( cl, 1, 3 ) ];
        pYtrans3[cl].zz = PCellCenters[ idx2( cl, 2, 3 ) ];
        pYtrans3[cl].ind = cl;
        pZtrans3[cl].xx = PCellCenters[ idx2( cl, 0, 3 ) ];
        pZtrans3[cl].yy = PCellCenters[ idx2( cl, 1, 3 ) ];
        pZtrans3[cl].zz = PCellCenters[ idx2( cl, 2, 3 ) ];
        pZtrans3[cl].ind = cl;
      }
      std::sort(pYtrans3.begin(), pYtrans3.end(), byZbyXbyY());
      std::sort(pZtrans3.begin(), pZtrans3.end(), byYbyXbyZ());
      for (unsigned long cl = 0; cl < (PCellCenters.size()/3); cl++) {
        PresListByY.push_back( pYtrans3[cl].ind );
        PresListByZ.push_back( pZtrans3[cl].ind );
      }
      break;
    }
  }
}
// create the pore-network from porescale meshes
void PoreNetwork::UniformPN( double length, double width, double height, int nx, int ny, int nz )
{
  if (nz) {
    DIM = 3;
    nPores = nx * ny * nz;
    dx = length/nx;
    dy = width/ny;
    dz = height/nz;
  }
  else {
    DIM = 2;
    nPores = nx * ny;
    dx = length/nx;
    dy = width/ny;
  }
  psLength = length;
  psWidth = width;
  psHeight = height;
  PoresXYZ.resize( nPores * DIM );
  Throats.resize( nPores * DIM * 2 );
  // set pore locations
  if (DIM == 2) {
    for (int porey = 0; porey < ny; porey++) {
      for (int porex = 0; porex < nx; porex++) {
        PoresXYZ[ idx2( idx2( porey, porex, ny ), 0, 2 ) ] = 0.5 * dx + dx * ( porex );
        PoresXYZ[ idx2( idx2( porey, porex, ny ), 1, 2 ) ] = 0.5 * dy + dy * ( porey );
      }
    }
  }
  else {
    for (int porez = 0; porez < nz; porez++) {
      for (int porey = 0; porey < ny; porey++) {
        for (int porex = 0; porex < nx; porex++) {
          PoresXYZ[ idx2( idx3( porez, porey, porex, ny, nz ), 0, 3 ) ] = 0.5 * dx + dx * ( porex );
          PoresXYZ[ idx2( idx3( porez, porey, porex, ny, nz ), 1, 3 ) ] = 0.5 * dy + dy * ( porey );
          PoresXYZ[ idx2( idx3( porez, porey, porex, ny, nz ), 2, 3 ) ] = 0.5 * dz + dz * ( porez );
        }
      }
    }
  }
  int nConnections;
  innerFaceConnectivity( Throats, PoresXYZ, dx, dy, dz, nPores, DIM );
  for (int pore = 0; pore < nPores; pore++) {
    nConnections = 0;
    for (int side = 0; side < (2*DIM); side++) {
      if (Throats[ idx2( pore, side, 2*DIM ) ]) nConnections++;
    }
    if (nConnections == 2*DIM) {
      InteriorPores.push_back( pore );
    }
    else {
      BoundaryPores.push_back( pore );
    }
  }
}
void SaveFluidMesh( const FluidMesh& Mesh, const std::string& outName )
{
  {
    std::ofstream ofs(outName.c_str());
    boost::archive::text_oarchive oa(ofs);
    oa << Mesh;
  }
}
void LoadFluidMesh( FluidMesh& Mesh, const std::string& inName )
{
  // load vectors
  {
    std::ifstream ifs(inName.c_str());
    boost::archive::text_iarchive ia(ifs);
    ia >> Mesh;
  }
}
