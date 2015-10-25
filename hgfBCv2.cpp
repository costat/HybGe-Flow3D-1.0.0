#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>

#include "hgf.hpp"
#include "hgfMesh.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index,
   uses row marjor indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

void
FlowComponent3d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, const std::vector<unsigned long>& componentBoundary, \
                  const std::vector<double> componentCellWidths, \
                  const std::vector<double> CellWidthsLR, \
                  const std::vector<double> CellWidthsUD, \
                  const std::vector<unsigned long> componentConnectivity, \
                  double visc, int direction, int component )
{

  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  int cl, cl2, fc, colId[6], nbrfaces[6], dirIn, dirOut, dirLeft, dirRight, dirUp, dirDown, colShift;
  double val[7], dx[6], dy[6], dz[6];

  switch (component)
  {
    case 0 :
    {
      dirIn = 3;
      dirOut = 1;
      dirLeft = 4;
      dirRight = 5;
      dirDown = 0;
      dirUp = 2;
      break;
    }
    case 1 :
    {
      dirIn = 5;
      dirOut = 4;
      dirLeft = 3;
      dirRight = 1;
      dirDown = 0;
      dirUp = 2;
      break;
    }
    case 2 :
    {
      dirIn = 0;
      dirOut = 2;
      dirLeft = 3;
      dirRight = 1;
      dirDown = 4;
      dirUp = 5;
      break;
    }
  }

  if (component == 0) colShift = 0;
  else if (component == 1) colShift = Mesh.DOF[1];
  else if (component == 2) colShift = Mesh.DOF[1] + Mesh.DOF[2];

  for (unsigned long arrayIndex = 0; arrayIndex < componentBoundary.size(); arrayIndex++) {
    cl = componentBoundary[ arrayIndex ];
    cl2 = cl + colShift;
    for (fc = 0; fc < 6; fc++)
    {
      nbrfaces[fc] = componentConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
      if (nbrfaces[fc]) {
        dx[fc] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 0, Mesh.CellWidthsLDI ) ];
        dy[fc] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 1, Mesh.CellWidthsLDI ) ];
        dz[fc] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 2, Mesh.CellWidthsLDI ) ];
      }
      else {
        dx[fc] = 0;
        dy[fc] = 0;
        dz[fc] = 0;
      }
    }
    if (!nbrfaces[dirIn]) // !dirIn section
    {
      if (component == direction) // inflow, nonzero force
      {

      }
    }
    else if (!nbrfaces[dirOut]) // !dirOut section
    {
      if (component == direction) // outflow
      {

      }
      else // noslip
      {

      }
    }
    else if (!nbrfaces[dirDown]) // !dirDown section
    {
      if (!nbrfaces[dirUp])
      {
        if (!nbrfaces[dirLeft])
        {
          if (!nbrfaces[dirRight]) // !dirDown, !dirUp, !dirLeft, !dirRight
          {

          }
          else // !dirDown, !dirUp, !dirLeft
          {

          }
        }
        else if (!nbrfaces[dirRight]) // !dirDown, !dirUp, !dirRight
        {

        }
        else // !dirDown, !dirUp
        {

        }
      }
      else if (!nbrfaces[dirLeft])
      {
        if (!nbrfaces[dirRight]) // !dirDown, !dirLeft, !dirRight
        {

        }
        else // !dirDown, !dirLeft
        {

        }
      }
      else if (!nbrfaces[dirRight]) // !dirDown, !dirRight
      {

      }
      else // !dirDown
      {

      }
    } // end !dirDown primary section
    else if (!nbrfaces[dirUp]) // !dirUp section
    {
      if (!nbrfaces[dirLeft])
      {
        if (!nbrfaces[dirRight]) // !dirUp, !dirLeft, !dirLeft
        {

        }
        else // !dirUp, !dirLeft
        {

        }
      }
      else if (!nbrfaces[dirRight]) // !dirUp, !dirRight
      {

      }
      else // !dirUp
      {

      }
    } // end !dirUp primary section
    else if (!nbrfaces[dirLeft]) // !dirLeft section
    {
      if (!nbrfaces[dirRight]) // !dirLeft, !dirRight
      {

      }
      else // !dirLeft
      {

      }
    } // end !dirLeft primary section
    else if (!nbrfaces[dirRight]) // !dirRight section
    {

    } // end !dirRight primary section
  }
}

void
AxisFlowDrive ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                std::vector<int>& matJs, std::vector<double>& matVals, \
                std::vector<double>& force, double visc, int direction )
{

  switch (Mesh.DIM)
  {
    case 2 :
    {

      break;
    }
    default :
    {
      const std::vector<unsigned long>& componentBoundary = Mesh.UBoundaryCells;
      const std::vector<double>& componentCellWidths = Mesh.UCellWidths;
      const std::vector<unsigned long>& componentConnectivity = Mesh.UFaceConnectivity;
      const std::vector<double>& cellWidthsLR = Mesh.VCellWidths;
      const std::vector<double>& cellWidthsUD = Mesh.WCellWidths;
      FlowComponent3D( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary, componentCellWidths, \
                       CellWidthsLR, CellWidthsUD, componentConnectivity, \
                       visc, direction, component );
      break;
    }
}

