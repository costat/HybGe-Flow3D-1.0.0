#include <iostream>
#include <cstdlib>
#include <vector>
#include <paralution.hpp>

#include "hgfMesh.hpp"
#include "hgf.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

void
computeKTensorL ( const FluidMesh& Mesh, \
                  const paralution::LocalVector<double>& xSolution, \
                  const paralution::LocalVector<double>& ySolution, \
                  const paralution::LocalVector<double>& zSolution )
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
      computeAveragesX ( Mesh, xSolution, Vels[0], GVals[0], 1 );
      computeAveragesY ( Mesh, xSolution, Vels[1], GVals[1], 0 );
      computeAveragesZ ( Mesh, xSolution, Vels[2], GVals[2], 0 );
      computeAveragesX ( Mesh, ySolution, Vels[3], GVals[3], 0 );
      computeAveragesY ( Mesh, ySolution, Vels[4], GVals[4], 1 );
      computeAveragesZ ( Mesh, ySolution, Vels[5], GVals[5], 0 );
      computeAveragesX ( Mesh, zSolution, Vels[6], GVals[6], 0 );
      computeAveragesY ( Mesh, zSolution, Vels[7], GVals[7], 0 );
      computeAveragesZ ( Mesh, zSolution, Vels[8], GVals[8], 1 );

      for (int i = 0; i < 9; i++) {
        GVals[i+9] = GVals[i];
        GVals[i+18] = GVals[i];
      }

      // Initialize paralution arrays
      LocalVector<double> K;
      LocalVector<double> MacroVels;
      LocalMatrix<double> MacroPresGrads;

      // Conductivity tensor
      K.Allocate("conductivity tensor", 9);
      K.Zeros();

      // Macro velocities
      MacroVels.Allocate("macroscopic velocities", 9);
      MacroVels.Zeros();
      for (int i = 0; i < 9; i++) {
        MacroVels[i] = Vels[i];
      }

      // Pressure gradients matrix
      MacroPresGrads.Assemble(GIs, GJs, GVals, 27, "pressure gradients", 9, 9);

      // Solver setup
      GMRES<LocalMatrix<double>, LocalVector<double>, double > ls;
      ls.SetOperator(MacroPresGrads);
      ls.Verbose(2);
      ls.Build();

      ls.Solve(MacroVels, &K);

      // Write out solution to file
      std::ofstream KTensor;
      KTensor.open ("KTensor.dat");
      KTensor << K[0] << "\t" << K[1] << "\t" << K[2] << "\n";
      KTensor << K[3] << "\t" << K[4] << "\t" << K[5] << "\n";
      KTensor << K[6] << "\t" << K[7] << "\t" << K[8];
      KTensor.close();
      MacroVels.Clear();
      K.Clear();
      MacroPresGrads.Clear();
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
      computeAveragesX ( Mesh, xSolution, Vels[0], GVals[0], 1 );
      computeAveragesY ( Mesh, xSolution, Vels[1], GVals[1], 0 );
      computeAveragesX ( Mesh, ySolution, Vels[2], GVals[2], 0 );
      computeAveragesY ( Mesh, ySolution, Vels[3], GVals[3], 1 );

      for (int i = 0; i < 4; i++) {
        GVals[i+4] = GVals[i];
      }

      // Initialize paralution arrays
      LocalVector<double> K;
      LocalVector<double> MacroVels;
      LocalMatrix<double> MacroPresGrads;

      // Conductivity tensor
      K.Allocate("conductivity tensor", 4);
      K.Zeros();

      // Macro velocities
      MacroVels.Allocate("macroscopic velocities", 4);
      MacroVels.Zeros();
      for (int i = 0; i < 4; i++) {
        MacroVels[i] = Vels[i];
      }

      // Pressure gradients matrix
      MacroPresGrads.Assemble(GIs, GJs, GVals, 8, "pressure gradients", 4, 4);

      // Solver setup
      GMRES<LocalMatrix<double>, LocalVector<double>, double > ls;
      ls.SetOperator(MacroPresGrads);
      ls.Verbose(2);
      ls.Build();

      ls.Solve(MacroVels, &K);

      // Write out solution to file
      std::ofstream KTensor;
      KTensor.open ("KTensor.dat");
      KTensor << K[0] << "\t" << K[1] << "\n";
      KTensor << K[2] << "\t" << K[3] << "\n";
      KTensor.close();
      MacroVels.Clear();
      K.Clear();
      MacroPresGrads.Clear();
      break;
    }
  }
}

void
computeAveragesX ( const FluidMesh& Mesh, \
                   const paralution::LocalVector<double>& Solution, \
                   double& V, double& G, int print )
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
        pNode1 = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
        pNode2 = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
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
        KConstantX.open("KConstantX.dat");
        KConstantX << K;
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
        pNode1 = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
        pNode2 = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
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
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangex;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantX;
        KConstantX.open("KConstantX.dat");
        KConstantX << K;
        KConstantX.close();
      }
      break;
    }
  }

}

void
computeAveragesY ( const FluidMesh& Mesh, \
                   const paralution::LocalVector<double>& Solution, \
                   double& V, double& G, int print )
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
        pNode1 = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
        pNode2 = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
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
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangey;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantY;
        KConstantY.open("KConstantY.dat");
        KConstantY << K;
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
        pNode1 = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
        pNode2 = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
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
      V = V / vInd;
      P1 = P1 / p1Ind;
      P2 = P2 / p2Ind;
      G = (P1 - P2) / midRangey;

      if (print)
      {
        K = V / G;
        // Write out K
        std::ofstream KConstantY;
        KConstantY.open("KConstantY.dat");
        KConstantY << K;
        KConstantY.close();
      }
      break;
    }
  }
}

void
computeAveragesZ ( const FluidMesh& Mesh, \
                    const paralution::LocalVector<double>& Solution, \
                    double& V, double& G, int print )
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
    pNode1 = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
    pNode2 = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] + 1;
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
  V = V / vInd;
  P1 = P1 / p1Ind;
  P2 = P2 / p2Ind;
  G = (P1 - P2) / midRangez;

  if (print)
  {
    K = V / G;
    // Write out K
    std::ofstream KConstantZ;
    KConstantZ.open("KConstantZ.dat");
    KConstantZ << K;
    KConstantZ.close();
  }
}

void
computeKConstantDrive ( const FluidMesh & Mesh, \
                        const paralution::LocalVector<double>& Solution,
                        int direction )
{
  double V, G;
  switch ( direction )
  {
    case 0 :
      computeAveragesX ( Mesh, Solution, V, G, 1 );
      break;
    case 1 :
      computeAveragesY ( Mesh, Solution, V, G, 1 );
      break;
    case 2 :
      computeAveragesZ ( Mesh, Solution, V, G, 1 );
      break;
  }
}


void
writeSolutionL ( const FluidMesh& Mesh, const paralution::LocalVector<double>& sol, \
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
                    + sol[ Mesh.PressureCellUNeighbor[ idx2( row, 1, 2 ) ] - 1] );
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
