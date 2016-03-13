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
  if (cl < nCells) {
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
}
__global__ void ifcKernel3D( unsigned long *d_CFC, const double *d_CCC, \
                             double epsx, double epsy, double epsz, \
                             double xtol, double ytol, double ztol, \
                             int nCells )
{
  int cl = blockIdx.x * blockDim.x + threadIdx.x;
  if (cl < nCells) {
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
}
void
MeshSubdivide( const ProbParam& Par, \
               std::vector< ProbParam >& SubPar  )
{
  switch ( Par.nz )
  {
    case 0 :
    {
      int nSubDomains = 0;
      int xStart, nxRemainder, yStart, nyRemainder, nyG, nxG, xCount, yCount;
      yCount = 0;
      yStart = 0;
      nyRemainder = Par.ny;
      for (int yy = 0; yy < Par.nCuts; yy++) {
        xStart = 0;
        nxRemainder = Par.nx;
        xCount = 0;
        nyG = (int)(round(((double)nyRemainder)/(Par.nCuts-yCount)));
        for (int xx = 0; xx < Par.nCuts; xx++) {
          nSubDomains++;
          nxG = (int)(round(((double)nxRemainder)/(Par.nCuts-xCount)));
          for (int cy = 0; cy < nyG; cy++) {
            for (int cx = 0; cx < nxG; cx++) {
              SubPar[ nSubDomains-1 ].gridin.push_back( Par.gridin[ idx2( (cy+yStart), (cx+xStart), Par.nx ) ] );
            }
          }
          SubPar[ nSubDomains-1 ].length = Par.length * ((double)nxG / Par.nx);
          SubPar[ nSubDomains-1 ].width = Par.width * ((double)nyG / Par.ny);
          SubPar[ nSubDomains-1 ].nx = nxG;
          SubPar[ nSubDomains-1 ].ny = nyG;
          xStart = xStart + nxG;
          nxRemainder = Par.nx - xStart;
          xCount++;
        }
        yStart = yStart + nyG;
        nyRemainder = Par.ny - yStart;
        yCount++;
      }
      break;
    }
    default :
    {
      int nSubDomains = 0;
      int xStart, nxRemainder, yStart, nyRemainder, zStart, nzRemainder, nxG, nyG, nzG, xCount, yCount, zCount;
      zCount = 0;
      zStart = 0;
      nzRemainder = Par.nz;
      for (int zz = 0; zz < Par.nCuts; zz++) {
        yStart = 0;
        nyRemainder = Par.ny;
        yCount = 0;
        nzG = (int)(round(((double)nzRemainder)/(Par.nCuts-zCount)));
        for (int yy = 0; yy < Par.nCuts; yy++) {
          xStart = 0;
          nxRemainder = Par.nx;
          xCount = 0;
          nyG = (int)(round(((double)nyRemainder)/(Par.nCuts-yCount)));
          for (int xx = 0; xx < Par.nCuts; xx++) {
            nSubDomains++;
            nxG = (int)(round(((double)nxRemainder)/(Par.nCuts-xCount)));
            for (int cz = 0; cz < nzG; cz++) {
              for (int cy = 0; cy < nyG; cy++) {
                for (int cx = 0; cx < nxG; cx++) {
                  SubPar[ nSubDomains-1 ].gridin.push_back( Par.gridin[ idx3( (cz+zStart), (cy+yStart), (cx+xStart), Par.ny, Par.nx ) ] );
                }
              }
            }
            SubPar[ nSubDomains-1 ].length = Par.length * ((double)nxG / Par.nx);
            SubPar[ nSubDomains-1 ].width = Par.width * ((double)nyG / Par.ny);
            SubPar[ nSubDomains-1 ].height = Par.height * ((double)nzG / Par.nz);
            SubPar[ nSubDomains-1 ].nx = nxG;
            SubPar[ nSubDomains-1 ].ny = nyG;
            SubPar[ nSubDomains-1 ].nz = nzG;
            xStart = xStart + nxG;
            nxRemainder = Par.nx - xStart;
            xCount++;
          }
          yStart = yStart + nyG;
          nyRemainder = Par.ny - yStart;
          yCount++;
        }
        zStart = zStart + nzG;
        nzRemainder = Par.nz - zStart;
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
    gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, (void*)ifcKernel2D, 0, nCells ) );
    gridSize = (nCells + blockSize - 1) / blockSize;
    ifcKernel2D<<< gridSize, blockSize >>>( d_CFC, d_CCC, epsx, epsy, xtol, ytol, nCells );
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
  }
  else if (DIM == 3)
  {
    gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize, (void*)ifcKernel3D, 0, nCells ) );
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
void FluidMesh::BuildUniformMesh( ProbParam& Par )
{
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

      // # pressure nodes
      for (int yi = 0; yi < Par.ny; yi++) {
        for (int xi = 0; xi < Par.nx; xi++) {
          if (Par.gridin[ idx2(yi, xi, Par.nx) ] != 1) {
            numPCells++;
            if (Par.gridin[ idx2(yi, xi, Par.nx) ] == 0) numVoid++;
          }
        }
      }
      int nNodes = numPCells * 4;
      double dx = Par.length / Par.nx;
      double dy = Par.width / Par.ny;
      double dz = 0;
      mv.resize( (numPCells * 4) );
      porosity = numVoid/(double)(Par.nx * Par.ny);

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((Par.nx * Par.ny));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((nNodes * 2));
      int countCell = -1;
      for (int yi = 0; yi < Par.ny; yi++) {
        for (int xi = 0; xi < Par.nx; xi++) {
          FullGrid.push_back(Par.gridin[ idx2(yi, xi, Par.nx) ]);
          if (Par.gridin[ idx2( yi, xi, Par.nx ) ] != 1) {
            countCell++;
            ImmersedBoundary.push_back(Par.gridin[ idx2(yi, xi, Par.nx) ]);

            Nodes.push_back( (xi + 1) * dx - dx );
            Nodes.push_back( (yi + 1) * dy - dy );

            Nodes.push_back( (xi + 1) * dx );
            Nodes.push_back( (yi + 1) * dy - dy );

            Nodes.push_back( (xi + 1) * dx );
            Nodes.push_back( (yi + 1) * dy );

            Nodes.push_back( (xi + 1) * dx - dx );
            Nodes.push_back( (yi + 1) * dy );

            mv[ idx2( countCell, 0, 4 ) ] = countCell*4;
            mv[ idx2( countCell, 1, 4 ) ] = countCell*4 + 1;
            mv[ idx2( countCell, 2, 4 ) ] = countCell*4 + 2;
            mv[ idx2( countCell, 3, 4 ) ] = countCell*4 + 3;
          }
        }
      }

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

      DOF[0] = numPCells;
      PFaceConnectivity.resize((numPCells * 4));
      innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells, DIM );

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, independent
         of other velocity components */
      #pragma omp parallel sections
      {
        { // Final P computations
          PCellWidths.resize((DOF[0] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[0]; cl++) {
            PCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            PCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
          }
        }
        #pragma omp section
        { // U grid
          int pl, ul;
          int countUCells = 0;
          double uStep = 0.5*dx;
          int maxUCells = Par.nx * Par.ny + Par.ny;
          UCellCenters.reserve((maxUCells * 2));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            if (PFaceConnectivity[ idx2( cl, 3, 4 ) ]) {
              countUCells++;
              pl = PFaceConnectivity[ idx2( cl, 3, 4 ) ]-1;
              ul = PressureCellUNeighbor[ idx2( pl, 1, 2 ) ]-1;
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 2 ) ] + uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 2 ) ]) );
              PressureCellUNeighbor[ idx2( cl, 0, 2 ) ] = ul+1;
              PressureCellUNeighbor[ idx2( cl, 1, 2 ) ] = countUCells;
              UCellPressureNeighbor[ idx2( ul, 1, 2 ) ] = cl+1;
              UCellPressureNeighbor[ idx2( (countUCells-1), 0, 2 ) ] = cl+1;
            }
            else
            {
              countUCells += 2;
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 2 ) ] - uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 2 ) ]) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 2 ) ] + uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 2 ) ]) );
              PressureCellUNeighbor[ idx2( cl, 0, 2 ) ] = countUCells-1;
              PressureCellUNeighbor[ idx2( cl, 1, 2 ) ] = countUCells;
              UCellPressureNeighbor[ idx2( (countUCells-2), 1, 2 ) ] = cl+1;
              UCellPressureNeighbor[ idx2( (countUCells-1), 0, 2 ) ] = cl+1;
            }
          }
          UCellCenters.shrink_to_fit();
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
          int pd, vd;
          int countVCells = 0;
          double vStep = 0.5*dy;
          int maxVCells = Par.nx * Par.ny + Par.nx;
          int yl;
          VCellCenters.reserve((maxVCells * 2));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int xl = 0; xl < numPCells; xl++) {
            yl = PresListByY[ xl ];
            if (PFaceConnectivity[ idx2( yl, 0, 4 ) ]) {
              countVCells++;
              pd = PFaceConnectivity[ idx2( yl, 0, 4 ) ]-1;
              vd = PressureCellVNeighbor[ idx2( pd, 1, 2 ) ]-1;
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 2 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 2 ) ] + vStep) );
              PressureCellVNeighbor[ idx2( yl, 0, 2 ) ] = vd+1;
              PressureCellVNeighbor[ idx2( yl, 1, 2 ) ] = countVCells;
              VCellPressureNeighbor[ idx2( vd, 1, 2 ) ] = yl+1;
              VCellPressureNeighbor[ idx2( (countVCells-1), 0, 2 ) ] = yl+1;
            }
            else {
              countVCells += 2;
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 2 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 2 ) ] - vStep) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 2 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 2 ) ] + vStep) );
              PressureCellVNeighbor[ idx2( yl, 0, 2 ) ] = countVCells - 1;
              PressureCellVNeighbor[ idx2( yl, 1, 2 ) ] = countVCells;
              VCellPressureNeighbor[ idx2( (countVCells-2), 1, 2 ) ] = yl+1;
              VCellPressureNeighbor[ idx2( (countVCells-1), 0, 2 ) ] = yl+1;
            }
          }
          VCellCenters.shrink_to_fit();
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
      NX = Par.nx;
      NY = Par.ny;
      NZ = Par.nz;

      // call to trim dead pores and incompatible cells
      Sanity( Par );

      // # P cells
      for (int zi = 0; zi < Par.nz; zi++) {
        for (int yi = 0; yi < Par.ny; yi++) {
          for (int xi = 0; xi < Par.nx; xi++) {
            if (Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] != 1) {
              numPCells++;
              if (Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] == 0) {
                numVoid++;
              }
            }
          }
        }
      }

      int nNodes = numPCells * 8;
      double dx = Par.length / Par.nx;
      double dy = Par.width / Par.ny;
      double dz = Par.height / Par.nz;
      mv.resize( (numPCells * 8) );
      porosity = numVoid /(double)(Par.nx * Par.ny * Par.nz);

      // First we buil FullGrid, ImmersedBoundary, and Nodes.
      FullGrid.reserve((Par.nx * Par.ny * Par.nz));
      ImmersedBoundary.reserve(numPCells);
      Nodes.reserve((nNodes * 3));
      int countCell = -1;
      for (int zi = 0; zi < Par.nz; zi++) {
        for (int yi = 0; yi < Par.ny; yi++) {
          for (int xi = 0; xi < Par.nx; xi++) {
            FullGrid.push_back(Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ]);
            if (Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] != 1) {
              countCell++;
              ImmersedBoundary.push_back(Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ]);

              Nodes.push_back( (xi + 1) * dx - dx );
              Nodes.push_back( (yi + 1) * dy - dy );
              Nodes.push_back( (zi + 1) * dz - dz );

              Nodes.push_back( (xi + 1) * dx );
              Nodes.push_back( (yi + 1) * dy - dy );
              Nodes.push_back( (zi + 1) * dz - dz );

              Nodes.push_back( (xi + 1) * dx );
              Nodes.push_back( (yi + 1) * dy );
              Nodes.push_back( (zi + 1) * dz - dz );

              Nodes.push_back( (xi + 1) * dx - dx );
              Nodes.push_back( (yi + 1) * dy );
              Nodes.push_back( (zi + 1) * dz - dz );

              Nodes.push_back( (xi + 1) * dx - dx);
              Nodes.push_back( (yi + 1) * dy );
              Nodes.push_back( (zi + 1) * dz );

              Nodes.push_back( (xi + 1) * dx );
              Nodes.push_back( (yi + 1) * dy );
              Nodes.push_back( (zi + 1) * dz );

              Nodes.push_back( (xi + 1) * dx );
              Nodes.push_back( (yi + 1) * dy - dy );
              Nodes.push_back( (zi + 1) * dz );

              Nodes.push_back( (xi + 1) * dx - dx );
              Nodes.push_back( (yi + 1) * dy - dy );
              Nodes.push_back( (zi + 1) * dz );

              mv[ idx2( countCell, 0, 8 ) ] = countCell*8;
              mv[ idx2( countCell, 1, 8 ) ] = countCell*8+1;
              mv[ idx2( countCell, 2, 8 ) ] = countCell*8+2;
              mv[ idx2( countCell, 3, 8 ) ] = countCell*8+3;
              mv[ idx2( countCell, 4, 8 ) ] = countCell*8+4;
              mv[ idx2( countCell, 5, 8 ) ] = countCell*8+5;
              mv[ idx2( countCell, 6, 8 ) ] = countCell*8+6;
              mv[ idx2( countCell, 7, 8 ) ] = countCell*8+7;
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
      DOF[0] = numPCells;
      PFaceConnectivity.resize((numPCells * 6));
      innerFaceConnectivity( PFaceConnectivity, PCellCenters, dx, dy, dz, numPCells, DIM );

      /* We finish mesh construction concurrently, since staggered grids
         for each component are constructed from the P grid, indendent of other
         velocity components.*/
      #pragma omp parallel sections
      {
        { // Final P computations
          PCellWidths.resize((DOF[0] * CellWidthsLDI));
          for (int cl = 0; cl < DOF[0]; cl++) {
            PCellWidths[ idx2( cl, 0, CellWidthsLDI ) ] = dx;
            PCellWidths[ idx2( cl, 1, CellWidthsLDI ) ] = dy;
            PCellWidths[ idx2( cl, 2, CellWidthsLDI ) ] = dz;
          }
        } // End P section
        #pragma omp section
        { // U section
          int pl, ul;
          int countUCells = 0;
          double uStep = 0.5*dx;
          int maxUCells = Par.nx * Par.ny * Par.nz + Par.ny * Par.nz;
          UCellCenters.reserve((maxUCells * 3));
          PressureCellUNeighbor.resize((numPCells * 2));
          UCellPressureNeighbor.resize((maxUCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            if (PFaceConnectivity[ idx2( cl, 3, 6 ) ]) {
              countUCells++;
              pl = PFaceConnectivity[ idx2( cl, 3, 6 ) ]-1;
              ul = PressureCellUNeighbor[ idx2( pl, 1, 2 ) ]-1;
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 3 ) ] + uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 3 ) ]) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 2, 3 ) ]) );
              PressureCellUNeighbor[ idx2( cl, 0, 2 ) ] = ul+1;
              PressureCellUNeighbor[ idx2( cl, 1, 2 ) ] = countUCells;
              UCellPressureNeighbor[ idx2( ul, 1, 2 ) ] = cl+1;
              UCellPressureNeighbor[ idx2( (countUCells-1), 0, 2 ) ] = cl+1;
            }
            else
            {
              countUCells += 2;
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 3 ) ] - uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 3 ) ]) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 2, 3 ) ]) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 0, 3 ) ] + uStep) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 1, 3 ) ]) );
              UCellCenters.push_back( (PCellCenters[ idx2( cl, 2, 3 ) ]) );
              PressureCellUNeighbor[ idx2( cl, 0, 2 ) ] = countUCells-1;
              PressureCellUNeighbor[ idx2( cl, 1, 2 ) ] = countUCells;
              UCellPressureNeighbor[ idx2( (countUCells-2), 1, 2 ) ] = cl+1;
              UCellPressureNeighbor[ idx2( (countUCells-1), 0, 2 ) ] = cl+1;
            }
          }
          UCellCenters.shrink_to_fit();
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
          int pb, vb;
          int countVCells = 0;
          double vStep = 0.5*dy;
          int maxVCells = Par.nx * Par.ny * Par.nz + Par.nx * Par.nz;
          int yl;
          VCellCenters.reserve((maxVCells * 3));
          PressureCellVNeighbor.resize((numPCells * 2));
          VCellPressureNeighbor.resize((maxVCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            yl = PresListByY[ cl ];
            if (PFaceConnectivity[ idx2( yl, 5, 6 ) ]) {
              countVCells++;
              pb = PFaceConnectivity[ idx2( yl, 5, 6 ) ]-1;
              vb = PressureCellVNeighbor[ idx2( pb, 1, 2 ) ]-1;
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 3 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 3 ) ] + vStep) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 2, 3 ) ]) );
              PressureCellVNeighbor[ idx2( yl, 0, 2 ) ] = vb+1;
              PressureCellVNeighbor[ idx2( yl, 1, 2 ) ] = countVCells;
              VCellPressureNeighbor[ idx2( vb, 1, 2 ) ] = yl+1;
              VCellPressureNeighbor[ idx2( (countVCells-1), 0, 2 ) ] = yl+1;
            }
            else {
              countVCells += 2;
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 3 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 3 ) ] - vStep) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 2, 3 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 0, 3 ) ]) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 1, 3 ) ] + vStep) );
              VCellCenters.push_back( (PCellCenters[ idx2( yl, 2, 3 ) ]) );
              PressureCellVNeighbor[ idx2( yl, 0, 2 ) ] = countVCells - 1;
              PressureCellVNeighbor[ idx2( yl, 1, 2 ) ] = countVCells;
              VCellPressureNeighbor[ idx2( (countVCells-2), 1, 2 ) ] = yl+1;
              VCellPressureNeighbor[ idx2( (countVCells-1), 0, 2 ) ] = yl+1;
            }
          }
          VCellCenters.shrink_to_fit();
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
          int pd, wd;
          int countWCells = 0;
          double wStep = 0.5*dz;
          int maxWCells = Par.nx * Par.ny * Par.nz + Par.nx * Par.ny;
          int zl;
          WCellCenters.reserve((maxWCells * 3));
          PressureCellWNeighbor.resize((numPCells * 2));
          WCellPressureNeighbor.resize((maxWCells * 2));
          for (int cl = 0; cl < numPCells; cl++) {
            zl = PresListByZ[ cl ];
            if (PFaceConnectivity[ idx2( zl, 0, 6 ) ]) {
              countWCells++;
              pd = PFaceConnectivity[ idx2( zl, 0, 6 ) ]-1;
              wd = PressureCellWNeighbor[ idx2( pd, 1, 2 ) ]-1;
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 0, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 1, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 2, 3 ) ] + wStep) );
              PressureCellWNeighbor[ idx2( zl, 0, 2 ) ] = wd+1;
              PressureCellWNeighbor[ idx2( zl, 1, 2 ) ] = countWCells;
              WCellPressureNeighbor[ idx2( wd, 1, 2 ) ] = zl+1;
              WCellPressureNeighbor[ idx2( (countWCells-1), 0, 2 ) ] = zl+1;
            }
            else {
              countWCells += 2;
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 0, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 1, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 2, 3 ) ] - wStep) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 0, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 1, 3 ) ]) );
              WCellCenters.push_back( (PCellCenters[ idx2( zl, 2, 3 ) ] + wStep) );
              PressureCellWNeighbor[ idx2( zl, 0, 2 ) ] = countWCells - 1;
              PressureCellWNeighbor[ idx2( zl, 1, 2 ) ] = countWCells;
              WCellPressureNeighbor[ idx2( (countWCells-2), 1, 2 ) ] = zl+1;
              WCellPressureNeighbor[ idx2( (countWCells-1), 0, 2 ) ] = zl+1;
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
// function removes cells that are boundaries in opposite directions
int FluidMesh::Sanity( ProbParam& Par )
{
  int totalChanged = 0;
  int nChanged;

  if (DIM == 3) goto sanityCheck3;
  else goto sanityCheck2;

  sanityCheck3 :
  {
    nChanged = 0;
    // sanity
    for (int zi = 0; zi < Par.nz; zi++) {
      for (int yi = 0; yi < Par.ny; yi++) {
        for (int xi = 0; xi < Par.nx; xi++) {
          if (Par.gridin[ idx3(zi, yi, xi, Par.ny, Par.nx) ] != 1) {
            // xi sanity
            if (xi==0) {
              if (Par.gridin[ idx3( zi, yi, (xi+1), Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else if (xi==Par.nx-1) {
              if (Par.gridin[ idx3( zi, yi, (xi-1), Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else {
              if (Par.gridin[ idx3( zi, yi, (xi+1), Par.ny, Par.nx ) ] == 1 \
                  && Par.gridin[ idx3( zi, yi, (xi-1), Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            // yi sanity
            if (yi==0) {
              if (Par.gridin[ idx3( zi, (yi+1), xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else if (yi==Par.ny-1) {
              if (Par.gridin[ idx3( zi, (yi-1), xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else {
              if (Par.gridin[ idx3( zi, (yi+1), xi, Par.ny, Par.nx ) ] == 1 \
                  && Par.gridin[ idx3( zi, (yi-1), xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            // zi sanity
            if (zi==0) {
              if (Par.gridin[ idx3( (zi+1), yi, xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else if (zi==Par.nz-1) {
              if (Par.gridin[ idx3( (zi-1), yi, xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
            else {
              if (Par.gridin[ idx3( (zi+1), yi, xi, Par.ny, Par.nx ) ] == 1 \
                  && Par.gridin[ idx3( (zi-1), yi, xi, Par.ny, Par.nx ) ] == 1) {
                Par.gridin[ idx3( zi, yi, xi, Par.ny, Par.nx ) ] = 1;
                nChanged++;
              }
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck3;
    else goto cleanup;
  }

  sanityCheck2 :
  {
    nChanged = 0;
    // sanity
    for (int yi = 0; yi < Par.ny; yi++) {
      for (int xi = 0; xi < Par.nx; xi++) {
        if (Par.gridin[ idx2( yi, xi, Par.nx ) ] != 1) {
          // xi sanity
          if (xi==0) {
            if (Par.gridin[ idx2( yi, (xi+1), Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
          else if (xi==Par.nx-1) {
            if (Par.gridin[ idx2( yi, (xi-1), Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
          else {
            if (Par.gridin[ idx2( yi, (xi+1), Par.nx ) ] == 1 \
                && Par.gridin[ idx2( yi, (xi-1), Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
          // yi sanity
          if (yi==0) {
            if (Par.gridin[ idx2( (yi+1), xi, Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
          else if (yi==Par.ny-1) {
            if (Par.gridin[ idx2( (yi-1), xi, Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
          else {
            if (Par.gridin[ idx2( (yi+1), xi, Par.nx ) ] == 1 \
                && Par.gridin[ idx2( (yi-1), xi, Par.nx ) ] == 1) {
              Par.gridin[ idx2( yi, xi, Par.nx ) ] = 1;
              nChanged++;
            }
          }
        }
      }
    }
    totalChanged += nChanged;
    if (nChanged != 0) goto sanityCheck2;
    else goto cleanup;
  }

  cleanup :
    std::cout << "\nWarning, input geometry was incompatible.\n";
    std::cout << totalChanged << " cells, representing ";
    if (DIM == 3) std::cout << (double)100*totalChanged/(Par.nx*Par.ny*Par.nz);
    else std::cout << (double)100*totalChanged/(Par.nx*Par.ny);
    std::cout << "% of the input geometry, with boundaries on opposite faces \nwere found and removed from void space.\n\n";
    return totalChanged;
}
// create the pore-network from porescale meshes
void PoreNetwork::UniformPN( const ProbParam& Par )
{
  if (Par.nz) {
    DIM = 3;
    nPores = Par.nCuts * Par.nCuts * Par.nCuts;
    dx = Par.length/Par.nCuts;
    dy = Par.width/Par.nCuts;
    dz = Par.height/Par.nCuts;
  }
  else {
    DIM = 2;
    nPores = Par.nCuts * Par.nCuts;
    dx = Par.length/Par.nCuts;
    dy = Par.width/Par.nCuts;
  }
  psLength = Par.length;
  psWidth = Par.width;
  psHeight = Par.height;
  PoresXYZ.resize( nPores * DIM );
  Throats.resize( nPores * DIM * 2 );
  // set pore locations
  if (DIM == 2) {
    for (int porey = 0; porey < Par.nCuts; porey++) {
      for (int porex = 0; porex < Par.nCuts; porex++) {
        PoresXYZ[ idx2( idx2( porey, porex, Par.nCuts ), 0, 2 ) ] = 0.5 * dx + dx * ( porex );
        PoresXYZ[ idx2( idx2( porey, porex, Par.nCuts ), 1, 2 ) ] = 0.5 * dy + dy * ( porey );
      }
    }
  }
  else {
    for (int porez = 0; porez < Par.nCuts; porez++) {
      for (int porey = 0; porey < Par.nCuts; porey++) {
        for (int porex = 0; porex < Par.nCuts; porex++) {
          PoresXYZ[ idx2( idx3( porez, porey, porex, Par.nCuts, Par.nCuts ), 0, 3 ) ] = 0.5 * dx + dx * ( porex );
          PoresXYZ[ idx2( idx3( porez, porey, porex, Par.nCuts, Par.nCuts ), 1, 3 ) ] = 0.5 * dy + dy * ( porey );
          PoresXYZ[ idx2( idx3( porez, porey, porex, Par.nCuts, Par.nCuts ), 2, 3 ) ] = 0.5 * dz + dz * ( porez );
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
