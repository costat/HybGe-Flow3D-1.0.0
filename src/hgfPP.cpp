#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfPP.hpp"

#define idx2(i, j, ldi) ((i * ldi) + j)

/* computeKTensorL computes an upscaled conductivity tensor using
   flow solutions in each principal axis direction */
void
computeKTensorL ( const FluidMesh& Mesh, \
                  const std::vector<double>& xSolution, \
                  const std::vector<double>& ySolution, \
                  const std::vector<double>& zSolution, \
                  const std::vector<std::string>& KoutNames )
{
  switch ( Mesh.DIM )
  {
    case 3 :
    {
      double GVals[27], Vels[9];
      int GIs[27], GJs[27];

      // Assign Indices
      for (int i = 0; i < 3; i++) {
        // First block
        GIs[i] = 0;
        GJs[i] = i;
        GIs[3+i] = 1;
        GJs[3+i] = i;
        GIs[6+i] = 2;
        GJs[6+i] = i;
        // Second block
        GIs[9+i] = 3;
        GJs[9+i] = 3 + i;
        GIs[12+i] = 4;
        GJs[12+i] = 3 + i;
        GIs[15+i] = 5;
        GJs[15+i] = 3 + i;
        // Third block
        GIs[18+i] = 6;
        GJs[18+i] = 6 + i;
        GIs[21+i] = 7;
        GJs[21+i] = 6 + i;
        GIs[24+i] = 8;
        GJs[24+i] = 6 + i;
      }

      // Compute averages
      computeAveragesX ( Mesh, xSolution, Vels[0], GVals[0], 1, KoutNames[0] );
      computeAveragesY ( Mesh, xSolution, Vels[1], GVals[1], 0, KoutNames[0] );
      computeAveragesZ ( Mesh, xSolution, Vels[2], GVals[2], 0, KoutNames[0] );
      computeAveragesX ( Mesh, ySolution, Vels[3], GVals[3], 0, KoutNames[1] );
      computeAveragesY ( Mesh, ySolution, Vels[4], GVals[4], 1, KoutNames[1] );
      computeAveragesZ ( Mesh, ySolution, Vels[5], GVals[5], 0, KoutNames[1] );
      computeAveragesX ( Mesh, zSolution, Vels[6], GVals[6], 0, KoutNames[2] );
      computeAveragesY ( Mesh, zSolution, Vels[7], GVals[7], 0, KoutNames[2] );
      computeAveragesZ ( Mesh, zSolution, Vels[8], GVals[8], 1, KoutNames[2] );

      for (int i = 0; i < 9; i++) {
        GVals[i+9] = GVals[i];
        GVals[i+18] = GVals[i];
      }

      // solve linear system for K tensor

      break;
    }
    case 2 :
    {
      double GVals[8], Vels[4];
      int GIs[8], GJs[8];

      // Assign Indices
      for (int i = 0; i < 2; i++) {
        // First block
        GIs[i] = 0;
        GJs[i] = i;
        GIs[2+i] = 1;
        GJs[2+i] = i;
        // Second block
        GIs[4+i] = 2;
        GJs[4+i] = 2 + i;
        GIs[6+i] = 3;
        GJs[6+i] = 2 + i;
      }

      // Compute averages
      computeAveragesX ( Mesh, xSolution, Vels[0], GVals[0], 1, KoutNames[0] );
      computeAveragesY ( Mesh, xSolution, Vels[1], GVals[1], 0, KoutNames[0] );
      computeAveragesX ( Mesh, ySolution, Vels[2], GVals[2], 0, KoutNames[1] );
      computeAveragesY ( Mesh, ySolution, Vels[3], GVals[3], 1, KoutNames[1] );

      for (int i = 0; i < 4; i++) {
        GVals[i+4] = GVals[i];
      }

      // solve linear system for K tensor

      break;
    }
  }
}

/* computeAveragesX computes appropriate pressure and velocity averages
   for upscaling conductivity in the X direction.
   If requested a constant upscaled conductivity
   is computed and written to an output file */
void
computeAveragesX ( const FluidMesh& Mesh, \
                   const std::vector< double >& Solution, \
                   double& V, double& G, int print, \
                   const std::string& KoutName )
{
  double xmin, xmax, xmid, midRangex, ymin, ymax, zmin, zmax, K, pNode1, pNode2;
  double P1 = 0;
  double P2 = 0;
  double r = 0.09;
  double pressure;
  int vInd = 0;
  int p1Ind = 0;
  int p2Ind = 0;

  G = 0;
  V = 0;

  switch ( Mesh.DIM )
  {
    case 3 :
    {
      int vDOF = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      // Geometric properties for volume averaging
      xmin = Mesh.xLim[0];
      xmax = Mesh.xLim[1];
      ymin = Mesh.yLim[0];
      ymax = Mesh.yLim[1];
      zmin = Mesh.zLim[0];
      zmax = Mesh.zLim[1];
      double L = xmax - xmin;
      double W = ymax - ymin;
      double H = zmax - zmin;

      // Adjust min and max to avoid potential recirculation near boundaries.
      // Then compute mid and midRange.
      xmin = xmin + r*L;
      xmax = xmax - r*L;
      xmid = 0.5 * (xmin + xmax);
      zmin = zmin + r*H;
      zmax = zmax - r*H;
      ymin = ymin + r*W;
      ymax = ymax - r*W;
      midRangex = 0.5 * (xmax + xmid) - 0.5 * (xmid + xmin);

      // Compute averages
      for (int cl = 0; cl < Mesh.DOF[1]; cl++)
      {
        pNode1 = Mesh.UCellPressureNeighbor[ \
                   idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ];
        pNode2 = Mesh.UCellPressureNeighbor[ \
                   idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ];
        if (pNode1 && pNode2)
        {
          pNode1--;
          pNode2--;
          if (!Mesh.ImmersedBoundary[ pNode1 ] && !Mesh.ImmersedBoundary[ pNode2 ])
          {
            if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmin)
            {
              if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmax) // Inside x limits
              {
                if (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymin)
                {
                  if (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymax) // Inside y limits
                  {
                    if (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] > zmin)
                    {
                      if (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] < zmax) // Inside z limits
                      {
                        pressure = 0.5 * (Solution[ pNode1 + vDOF ] \
                                   + Solution[ pNode2 + vDOF ] );
                        if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmid)
                        {
                          P1 = P1 + pressure;
                          p1Ind++;
                          V = V + Solution[ cl ];
                          vInd++;
                        }
                        else if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmid)
                        {
                          P2 = P2 + pressure;
                          p2Ind++;
                          V = V + Solution[ cl ];
                          vInd++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangex;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantX;
        KConstantX.open( KoutName.c_str() );
        KConstantX << K*Mesh.porosity;
        KConstantX.close();
      }
      break;
    }
    case 2 :
    {
      int vDOF = Mesh.DOF[1] + Mesh.DOF[2];
      // Geometric properties for volume averaging
      xmin = Mesh.xLim[0];
      xmax = Mesh.xLim[1];
      ymin = Mesh.yLim[0];
      ymax = Mesh.yLim[1];
      double L = xmax - xmin;
      double W = ymax - ymin;

      // Adjust min and max to avoid potential recirculation near boundaries.
      // Then compute mid and midRange.
      xmin = xmin + r*L;
      xmax = xmax - r*L;
      xmid = 0.5 * (xmin + xmax);
      ymin = ymin + r*W;
      ymax = ymax - r*W;
      midRangex = 0.5 * (xmax + xmid) - 0.5 * (xmid + xmin);
      // Compute averages
      for (int cl = 0; cl < Mesh.DOF[1]; cl++)
      {
        pNode1 = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ];
        pNode2 = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ];
        if (pNode1 && pNode2)
        {
          pNode1--;
          pNode2--;
          if (!Mesh.ImmersedBoundary[ pNode1 ] && !Mesh.ImmersedBoundary[ pNode2 ])
          {
            if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmin)
            {
              if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmax) // Inside x limits
              {
                if (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymin)
                {
                  if (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymax) // Inside y limits
                  {
                    pressure = 0.5 * (Solution[ pNode1 + vDOF ] + Solution[ pNode2 + vDOF ] );
                    if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmid)
                    {
                      P1 = P1 + pressure;
                      p1Ind++;
                      V = V + Solution[ cl ];
                      vInd++;
                    }
                    else if (Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmid)
                    {
                      P2 = P2 + pressure;
                      p2Ind++;
                      V = V + Solution[ cl ];
                      vInd++;
                    }
                  }
                }
              }
            }
          }
        }
      }
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangex;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantX;
        KConstantX.open( KoutName.c_str() );
        KConstantX << K*Mesh.porosity;
        KConstantX.close();
      }
      break;
    }
  }

}

/* computeAveragesY computes appropriate pressure and velocity averages
   for upscaling conductivity in the Y direction.
   If requested a constant upscaled conductivity
   is computed and written to an output file */
void
computeAveragesY ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print, \
                   const std::string& KoutName )
{
  double xmin, xmax, midRangey, ymin, ymax, ymid, zmin, zmax, K, pNode1, pNode2;
  double P1 = 0;
  double P2 = 0;
  double r = 0.09;
  double pressure;
  int vInd = 0;
  int p1Ind = 0;
  int p2Ind = 0;

  V = 0;
  G = 0;

  switch ( Mesh.DIM )
  {
    case 3 :
    {
      int vDOF = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      // Geometric properties for volume averaging
      xmin = Mesh.xLim[0];
      xmax = Mesh.xLim[1];
      ymin = Mesh.yLim[0];
      ymax = Mesh.yLim[1];
      zmin = Mesh.zLim[0];
      zmax = Mesh.zLim[1];
      double L = xmax - xmin;
      double W = ymax - ymin;
      double H = zmax - zmin;

      // Adjust min and max to avoid potential recirculation near boundaries.
      // Then compute mid and midRange.
      xmin = xmin + r*L;
      xmax = xmax - r*L;
      zmin = zmin + r*H;
      zmax = zmax - r*H;
      ymin = ymin + r*W;
      ymax = ymax - r*W;
      ymid = 0.5 * (ymin + ymax);
      midRangey = 0.5 * (ymax + ymid) - 0.5 * (ymid + ymin);

      // Compute averages
      for (int cl = 0; cl < Mesh.DOF[2]; cl++)
      {
        pNode1 = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ];
        pNode2 = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ];
        if (pNode1 && pNode2)
        {
          pNode1--;
          pNode2--;
          if (!Mesh.ImmersedBoundary[ pNode1 ] && !Mesh.ImmersedBoundary[ pNode2 ])
          {
            if (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmin)
            {
              if (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmax) // Inside x limits
              {
                if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymin)
                {
                  if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymax) // Inside y limits
                  {
                    if (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] > zmin)
                    {
                      if (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] < zmax) // Inside z limits
                      {
                        pressure = 0.5 * (Solution[ pNode1 + vDOF ] + Solution[ pNode2 + vDOF ] );
                        if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymid)
                          {
                          P1 = P1 + pressure;
                          p1Ind++;
                          V = V + Solution[ cl + Mesh.DOF[1] ];
                          vInd++;
                        }
                        else if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymid)
                        {
                          P2 = P2 + pressure;
                          p2Ind++;
                          V = V + Solution[ cl + Mesh.DOF[1] ];
                          vInd++;
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangey;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantY;
        KConstantY.open( KoutName.c_str() );
        KConstantY << K*Mesh.porosity;
        KConstantY.close();
      }
      break;
    }
    case 2 :
    {
      int vDOF = Mesh.DOF[1] + Mesh.DOF[2];
      // Geometric properties for volume averaging
      xmin = Mesh.xLim[0];
      xmax = Mesh.xLim[1];
      ymin = Mesh.yLim[0];
      ymax = Mesh.yLim[1];
      double L = xmax - xmin;
      double W = ymax - ymin;

      // Adjust min and max to avoid potential recirculation near boundaries.
      // Then compute mid and midRange.
      xmin = xmin + r*L;
      xmax = xmax - r*L;
      ymin = ymin + r*W;
      ymax = ymax - r*W;
      ymid = 0.5 * (ymin + ymax);
      midRangey = 0.5 * (ymax + ymid) - 0.5 * (ymid + ymin);

      // Compute averages
      for (int cl = 0; cl < Mesh.DOF[2]; cl++)
      {
        pNode1 = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ];
        pNode2 = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ];
        if (pNode1 && pNode2)
        {
          pNode1--;
          pNode2--;
          if (!Mesh.ImmersedBoundary[ pNode1 ] && !Mesh.ImmersedBoundary[ pNode2 ])
          {
            if (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmin)
            {
              if (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmax) // Inside x limits
              {
                if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymin)
                {
                  if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymax) // Inside y limits
                  {
                    pressure = 0.5 * (Solution[ pNode1 + vDOF ] + Solution[ pNode2 + vDOF ] );
                    if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymid)
                    {
                      P1 = P1 + pressure;
                      p1Ind++;
                      V = V + Solution[ cl + Mesh.DOF[1] ];
                      vInd++;
                    }
                    else if (Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymid)
                    {
                      P2 = P2 + pressure;
                      p2Ind++;
                      V = V + Solution[ cl + Mesh.DOF[1] ];
                      vInd++;
                    }
                  }
                }
              }
            }
          }
        }
      }
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangey;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantY;
        KConstantY.open( KoutName.c_str() );
        KConstantY << K*Mesh.porosity;
        KConstantY.close();
      }
      break;
    }
  }
}

/* computeAveragesZ computes appropriate pressure and velocity averages
   for upscaling conductivity in the Z direction.
   If requested a constant upscaled conductivity
   is computed and written to an output file */
void
computeAveragesZ ( const FluidMesh& Mesh, \
                   const std::vector<double>& Solution, \
                   double& V, double& G, int print, \
                   const std::string& KoutName )
{
  double xmin, xmax, midRangez, ymin, ymax, zmid, zmin, zmax, K, pNode1, pNode2;
  double P1 = 0;
  double P2 = 0;
  double r = 0.09;
  double pressure;
  int vInd = 0;
  int p1Ind = 0;
  int p2Ind = 0;

  V = 0;
  G = 0;

  int vDOF = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
  // Geometric properties for volume averaging
  xmin = Mesh.xLim[0];
  xmax = Mesh.xLim[1];
  ymin = Mesh.yLim[0];
  ymax = Mesh.yLim[1];
  zmin = Mesh.zLim[0];
  zmax = Mesh.zLim[1];
  double L = xmax - xmin;
  double W = ymax - ymin;
  double H = zmax - zmin;

  // Adjust min and max to avoid potential recirculation near boundaries.
  // Then compute mid and midRange.
  xmin = xmin + r*L;
  xmax = xmax - r*L;
  zmin = zmin + r*H;
  zmax = zmax - r*H;
  ymin = ymin + r*W;
  ymax = ymax - r*W;
  zmid = 0.5 * (zmin + zmax);
  midRangez = 0.5 * (zmax + zmid) - 0.5 * (zmid + zmin);

  // Compute averages
  for (int cl = 0; cl < Mesh.DOF[3]; cl++)
  {
    pNode1 = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ];
    pNode2 = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ];
    if (pNode1 && pNode2)
    {
      pNode1--;
      pNode2--;
      if (!Mesh.ImmersedBoundary[ pNode1 ] && !Mesh.ImmersedBoundary[ pNode2 ])
      {
        if (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] > xmin)
        {
          if (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] < xmax) // Inside x limits
          {
            if (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] > ymin)
            {
              if (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] < ymax) // Inside y limits
              {
                if (Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] > zmin)
                {
                  if (Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] < zmax) // Inside z limits
                  {
                    pressure = 0.5 * (Solution[ pNode1 + vDOF ] + Solution[ pNode2 + vDOF ] );
                    if (Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] < zmid)
                    {
                      P1 = P1 + pressure;
                      p1Ind++;
                      V = V + Solution[ cl + Mesh.DOF[1] + Mesh.DOF[2] ];
                      vInd++;
                    }
                    else if (Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] > zmid)
                    {
                      P2 = P2 + pressure;
                      p2Ind++;
                      V = V + Solution[ cl + Mesh.DOF[1] + Mesh.DOF[2] ];
                      vInd++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  V = V / vInd;
  P1 = P1 / p1Ind;
  P2 = P2 / p2Ind;
  G = (P1 - P2) / midRangez;

  if (print)
  {
    K = V / G;
    // Write out K
    std::ofstream KConstantZ;
    KConstantZ.open( KoutName.c_str() );
    KConstantZ << K*Mesh.porosity;
    KConstantZ.close();
  }
}

/* computeKConstantDrive drives the computation and write
   of a constant upscaled conductivity */
void
computeKConstantDrive ( const FluidMesh & Mesh, \
                        const std::vector<double>& Solution, \
                        double& K, \
                        int direction, \
                        int print, \
                        const std::string& KoutName )
{
  double V, G;
  switch ( direction )
  {
    case 0 :
      computeAveragesX ( Mesh, Solution, V, G, print, KoutName );
      break;
    case 1 :
      computeAveragesY ( Mesh, Solution, V, G, print, KoutName );
      break;
    case 2 :
      computeAveragesZ ( Mesh, Solution, V, G, print, KoutName );
      break;
  }
  K = (V/G) * Mesh.porosity;
}

/* writeSolutionTP writes the solution to an output file appropriate for tecplot
   visualization */
void
writeSolutionTP ( const FluidMesh& Mesh, const std::vector<double>& sol, \
                  std::string& outName )
{
  double uval, vval, wval;
  int nNodes;
  int nEls;
  int horizCount = 0;
  std::ofstream flowrun;
  flowrun.open ( outName.c_str() );
  flowrun << "Title = Stokes Solution\n";
  switch ( Mesh.DIM )
  {
    case 2 :
    {
      nNodes = Mesh.Nodes.size() / 2;
      nEls = Mesh.DOF[0];
      flowrun << "VARIABLES = X, Y, P, U, V, IB\n";
      flowrun << "ZONE N=" << nNodes << ", E=" << nEls << ", DATAPACKING=BLOCK, ZONETYPE=FEQuadrilateral\n";
      flowrun << "VARLOCATION=( [3-6]=CELLCENTERED )\n";
      for (int row = 0; row < nNodes; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.Nodes[ idx2( row, 0, 2 ) ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.Nodes[ idx2( row, 0, 2 ) ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nNodes; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.Nodes[ idx2( row, 1, 2 ) ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.Nodes[ idx2( row, 1, 2 ) ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << sol[ row  + Mesh.DOF[1] + Mesh.DOF[2] ] << "\t";
        }
        else
        {
          flowrun << "\n" << sol[ row + Mesh.DOF[1] + Mesh.DOF[2] ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        uval = 0.5 * (sol[ Mesh.PressureCellUNeighbor[ idx2( row, 0, 2 ) ] - 1 ] \
                    + sol[ Mesh.PressureCellUNeighbor[ idx2( row, 1, 2 ) ] - 1 ] );
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << uval << "\t";
        }
        else
        {
          flowrun << "\n" << uval << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        vval = 0.5 * (sol[ Mesh.PressureCellVNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] ] \
                    + sol[ Mesh.PressureCellVNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] ] );
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << vval << "\t";
        }
        else
        {
          flowrun << "\n" << vval << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.ImmersedBoundary[ row ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.ImmersedBoundary[ row ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 200)
        {
          flowrun << Mesh.mv[ idx2( row, 0, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 1, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 2, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 3, 4 ) ] + 1 << "\t";
          horizCount = horizCount + 4;
        }
        else
        {
          flowrun << "\n";
          flowrun << Mesh.mv[ idx2( row, 0, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 1, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 2, 4 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 3, 4 ) ] + 1 << "\t";
          horizCount = 4;
        }
      }
      break;
    }
    case 3 :
    {
      nNodes = Mesh.Nodes.size() / 3;
      nEls = Mesh.DOF[0];
      flowrun << "VARIABLES = X, Y, Z, P, U, V, W, IB\n";
      flowrun << "ZONE N=" << nNodes << ", E=" << nEls << ", DATAPACKING=BLOCK, ZONETYPE=FEBrick\n";
      flowrun << "VARLOCATION=( [4-8]=CELLCENTERED)\n";
      for (int row = 0; row < nNodes; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.Nodes[ idx2( row, 0, 3 ) ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.Nodes[ idx2( row, 0, 3 ) ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nNodes; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.Nodes[ idx2( row, 1, 3 ) ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.Nodes[ idx2( row, 1, 3 ) ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nNodes; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.Nodes[ idx2( row, 2, 3 ) ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.Nodes[ idx2( row, 2, 3 ) ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << sol[ row  + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3] ] << "\t";
        }
        else
        {
          flowrun << "\n" << sol[ row + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3]  ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        uval = 0.5 * (sol[ Mesh.PressureCellUNeighbor[ idx2( row, 0, 2 ) ] - 1 ] \
                    + sol[ Mesh.PressureCellUNeighbor[ idx2( row, 1, 2 ) ] - 1 ] );
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << uval << "\t";
        }
        else
        {
          flowrun << "\n" << uval << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        vval = 0.5 * (sol[ Mesh.PressureCellVNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] ] \
                    + sol[ Mesh.PressureCellVNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] ] );
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << vval << "\t";
        }
        else
        {
          flowrun << "\n" << vval << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        wval = 0.5 * (sol[ Mesh.PressureCellWNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] + Mesh.DOF[2] ] \
                    + sol[ Mesh.PressureCellWNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] + Mesh.DOF[2] ] );
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << wval << "\t";
        }
        else
        {
          flowrun << "\n" << wval << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";

      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 1000)
        {
          flowrun << Mesh.ImmersedBoundary[ row ] << "\t";
        }
        else
        {
          flowrun << "\n" << Mesh.ImmersedBoundary[ row ] << "\t";
          horizCount = 1;
        }
      }
      horizCount = 0;
      flowrun << "\n";
      for (int row = 0; row < nEls; row++)
      {
        horizCount++;
        if (horizCount < 200)
        {
          flowrun << Mesh.mv[ idx2( row, 0, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 1, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 2, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 3, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 7, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 6, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 5, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 4, 8 ) ] + 1 << "\t";
          horizCount = horizCount + 8;
        }
        else
        {
          flowrun << "\n";
          flowrun << Mesh.mv[ idx2( row, 0, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 1, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 2, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 3, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 7, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 6, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 5, 8 ) ] + 1 << "\t";
          flowrun << Mesh.mv[ idx2( row, 4, 8 ) ] + 1 << "\t";
          horizCount = 8;
        }
      }
      break;
    }
  }
  flowrun.close();
}
/* writeSolutionPV writes the solution to an output file appropriate for paraview
   visualization */
void
writeSolutionPV ( const FluidMesh& Mesh, const std::vector<double>& sol, \
                  std::string& outName )
{
  double uval, vval, wval;
  int nNodes;
  int nEls;
  std::ofstream flowrun;
  flowrun.open( outName.c_str() );

  switch ( Mesh.DIM )
  {
    case 2 :
    {
      nNodes = Mesh.Nodes.size() / 2;
      nEls = Mesh.DOF[0];
      flowrun << "# vtk DataFile Version 3.0\n";
      flowrun << "vtk output\n";
      flowrun << "ASCII\n\n";
      flowrun << "DATASET UNSTRUCTURED_GRID\n";
      flowrun << "POINTS " << nNodes << " double\n";
      for (int row = 0; row < nNodes; row++) {
        flowrun << Mesh.Nodes[ idx2( row, 0, 2 ) ] << "\t";
        flowrun << Mesh.Nodes[ idx2( row, 1, 2 ) ] << "\t";
        flowrun << 0.0 << "\n";
      }
      flowrun << "\n";
      flowrun << "CELLS " << nEls << " " << 5*nEls << "\n";
      for (int row = 0; row < nEls; row++) {
        flowrun << 4 << "\t";
        flowrun << Mesh.mv[ idx2( row, 0, 4 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 1, 4 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 2, 4 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 3, 4 ) ] << "\t";
      }
      flowrun << "\n";
      flowrun << "CELL_TYPES " << nEls << "\n";
      for (int row = 0; row < nEls; row++) {
        flowrun << 9 << "\n";
      }
      flowrun << "\n";
      flowrun << "CELL_DATA " << nEls << "\n";
      flowrun << "SCALARS pressure double\n";
      flowrun << "LOOKUP_TABLE default\n";
      for (int row = 0; row < nEls; row++) {
        flowrun << sol[ Mesh.DOF[1] + Mesh.DOF[2] + row ] << "\n";
      }
      flowrun << "\n";
      flowrun << "VECTORS velocity double\n";
      for (int row = 0; row < nEls; row++) {
        uval = 0.5 * (sol[ Mesh.PressureCellUNeighbor[ idx2( row, 0, 2 ) ] - 1 ] \
                   + sol[ Mesh.PressureCellUNeighbor[ idx2( row, 1, 2 ) ] - 1 ] );
        vval = 0.5 * (sol[ Mesh.PressureCellVNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] ] \
                   + sol[ Mesh.PressureCellVNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] ] );
        flowrun << uval << "\t" << vval << "\t" << 0 << "\n";
      }
      flowrun << "\n";
      flowrun << "SCALARS immersedboundary int\n";
      flowrun << "LOOKUP_TABLE default\n";
      for (int row = 0; row < nEls; row++)
      {
        flowrun << Mesh.ImmersedBoundary[ row ] << "\n";
      }
      break;
    }
    case 3 :
    {
      nNodes = Mesh.Nodes.size() / 3;
      nEls = Mesh.DOF[0];
      flowrun << "# vtk DataFile Version 3.0\n";
      flowrun << "vtk output\n";
      flowrun << "ASCII\n\n";
      flowrun << "DATASET UNSTRUCTURED_GRID\n";
      flowrun << "POINTS " << nNodes << " double\n";
      for (int row = 0; row < nNodes; row++)
      {
        flowrun << Mesh.Nodes[ idx2( row, 0, 3 ) ] << "\t";
        flowrun << Mesh.Nodes[ idx2( row, 1, 3 ) ] << "\t";
        flowrun << Mesh.Nodes[ idx2( row, 2, 3 ) ] << "\n";
      }
      flowrun << "\n";
      flowrun << "CELLS " << nEls << " " << 9*nEls << "\n";
      for (int row = 0; row < nEls; row++)
      {
        flowrun << 8 << "\t";
        flowrun << Mesh.mv[ idx2( row, 0, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 1, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 2, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 3, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 7, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 6, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 5, 8 ) ] << "\t";
        flowrun << Mesh.mv[ idx2( row, 4, 8 ) ] << "\n";
      }
      flowrun << "\n";
      flowrun << "CELL_TYPES " << nEls << "\n";
      for (int row = 0; row < nEls; row++)
      {
        flowrun << 12 << "\n";
      }
      flowrun << "\n";
      flowrun << "CELL_DATA " << nEls << "\n";
      flowrun << "SCALARS pressure double\n";
      flowrun << "LOOKUP_TABLE default\n";
      for (int row = 0; row < nEls; row++)
      {
        flowrun << sol[ Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3] + row ] << "\n";
      }
      flowrun << "\n";
      flowrun << "VECTORS velocity double\n";
      for (int row = 0; row < nEls; row++)
      {
        uval = 0.5 * (sol[ Mesh.PressureCellUNeighbor[ idx2( row, 0, 2 ) ] - 1 ] \
                   + sol[ Mesh.PressureCellUNeighbor[ idx2( row, 1, 2 ) ] - 1 ] );
        vval = 0.5 * (sol[ Mesh.PressureCellVNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] ] \
                   + sol[ Mesh.PressureCellVNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] ] );
        wval = 0.5 * (sol[ Mesh.PressureCellWNeighbor[ idx2( row, 0, 2 ) ] - 1 + Mesh.DOF[1] + Mesh.DOF[2] ] \
                   + sol[ Mesh.PressureCellWNeighbor[ idx2( row, 1, 2 ) ] - 1 + Mesh.DOF[1] + Mesh.DOF[2] ] );
        flowrun << uval << "\t" << vval << "\t" << wval << "\n";
      }
      flowrun << "\n";
      flowrun << "SCALARS immersedboundary int\n";
      flowrun << "LOOKUP_TABLE default\n";
      for (int row = 0; row < nEls; row++)
      {
        flowrun << Mesh.ImmersedBoundary[ row ] << "\n";
      }
      break;
    }
  }
  flowrun.close();
}
void
computeKPoreNetwork( const PoreNetwork& pn, const std::vector<double>& Solution, \
                     const std::vector<double>& Ks, double& K, int direction, int print, \
                     const std::string& KoutName )
{
  int dirIn;
  double A, L;
  if (direction == 0) {
    dirIn = 3;
    if (pn.DIM == 2) {
      L = pn.psLength;
      A = (pn.psWidth-pn.dy) * pn.dx;
    }
    else if (pn.DIM == 3) {
      L = pn.psLength;
      A = (pn.psWidth-pn.dy) * (pn.psHeight-pn.dz);
    }
  }
  else if (direction == 1) {
    if (pn.DIM == 2) {
      dirIn = 0;
      L = pn.psWidth;
      A = (pn.psLength-pn.dx) * pn.dy;
    }
    else if (pn.DIM == 3) {
      dirIn = 5;
      L = pn.psWidth;
      A = (pn.psLength-pn.dx) * (pn.psHeight-pn.dz);
    }
  }
  else {
    dirIn = 0;
    L = pn.psHeight;
    A = (pn.psWidth-pn.dy) * (pn.psLength-pn.dx);
  }

  int pore;
  std::vector<unsigned long> lPores;
  lPores.reserve( pn.BoundaryPores.size() );
  for (int ii = 0; ii < pn.BoundaryPores.size(); ii++) {
    pore = pn.BoundaryPores[ ii ];
    if (!pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]) lPores.push_back( pore );
  }
  double KHOLD = 0;
  double TK;
  for (int ii = 0; ii < lPores.size(); ii++) {
    pore = lPores[ ii ];
    for (int dir = 0; dir < (pn.DIM*2); dir++) {
      if ( dir == 1 || dir == 3 ) {
        TK = 0.5 * Ks[ idx2( pore, 0, 2 ) ];
      }
      else if ( dir == 0 || dir == 2 ) {
        if (pn.DIM == 2) {
          TK = 0.5 * Ks[ idx2( pore, 1, (pn.DIM) ) ];
        }
        else {
          TK = 0.5 * Ks[ idx2( pore, 2, (pn.DIM) ) ];
        }
      }
      else {
        TK = 0.5 * Ks[ idx2( pore, 1, (pn.DIM) ) ];
      }
      KHOLD = KHOLD + TK * ( 1 - Solution[ pore ] );
    }
  }
  K = (KHOLD * L) / A;

  if (print) {
    // Write out K
    std::ofstream KPNout;
    KPNout.open(KoutName.c_str());
    KPNout << K;
    KPNout.close();
  }
}
void
writePoreNetworkSolutionTP ( const PoreNetwork& pn, const std::vector<double>& sol, \
                             std::string& outName )
{

  std::vector<unsigned long> ThroatConns;
  ThroatConns.reserve( pn.nPores*4 );
  for (int ii = 0; ii < pn.nPores; ii++) {
    for (int jj = 0; jj < (pn.DIM*2); jj++) {
      if ((int)(pn.Throats[ idx2( ii, jj, (pn.DIM*2) ) ]-1) > ii) {
        ThroatConns.push_back(ii);
        ThroatConns.push_back(pn.Throats[ idx2( ii, jj, (pn.DIM*2) ) ]-1);
      }
    }
  }
  int nThroats = ThroatConns.size()/2;

  std::ofstream flowrun;
  flowrun.open ( outName.c_str() );
  flowrun << "Title = Pore-Network Solution\n";
  flowrun << "Variables = X, Y, Z, P\n";
  flowrun << "ZONE N=" << pn.nPores << ", E=" << nThroats << ", F=FEPOINT\n";

  switch (pn.DIM)
  {
    case 2 :
    {
      for (int ii = 0; ii < pn.nPores; ii++) {
        flowrun << pn.PoresXYZ[ idx2( ii, 0, 2 ) ] << "\t" << pn.PoresXYZ[ idx2( ii, 1, 2 ) ] << "\t" << 0.0 << "\t" << sol[ ii ] << "\n";
      }
      break;
    }
    case 3 :
    {
      for (int ii = 0; ii < pn.nPores; ii++) {
        flowrun << pn.PoresXYZ[ idx2( ii, 0, 3 ) ] << "\t" << pn.PoresXYZ[ idx2( ii, 1, 3 ) ] << "\t";
        flowrun << pn.PoresXYZ[ idx2( ii, 2, 3 ) ] << "\t" << sol[ ii ] << "\n";
      }
      break;
    }
  }
  for (int ii = 0; ii < nThroats; ii++) {
    flowrun << ThroatConns[ idx2( ii, 0, 2 ) ]+1 << "\t" << ThroatConns[ idx2( ii, 0, 2 ) ]+1 << "\t" << ThroatConns[ idx2( ii, 0, 2 ) ]+1;
    flowrun << "\t" << ThroatConns[ idx2( ii, 1, 2 ) ]+1 << "\n";
  }

  flowrun.close();
}

void
writeKCollection( const std::vector<double>& Kvec, int dim, std::string& outName )
{
  std::ofstream kv;
  kv.open( outName.c_str() );
  for (int ii = 0; ii < (Kvec.size() / dim); ii++) {
    kv << Kvec[ idx2( ii, 0, dim ) ];
    kv << "\t" << Kvec[ idx2( ii, 1, dim ) ];
    if (dim == 3) kv << "\t" << Kvec[ idx2( ii, 2, dim ) ];
    kv << "\n";
  }
  kv.close();
}
