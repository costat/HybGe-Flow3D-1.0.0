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
                  double visc, int direction, int component )
{

  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  int cl, cl2, fc, colId[5], nbrfaces[5], dirIn, dirOut, dirLeft, dirRight, dirUp, dirDown;
  double val[5], dx[5], dy[5];

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
    }
    case 1 :
    {
      dirIn = 5;
      dirOut = 4;
      dirLeft = 3;
      dirRight = 1;
      dirDown = 0;
      dirUp = 2;
    }
    case 2 :
    {
      dirIn = 0;
      dirOut = 2;
      dirLeft = 3;
      dirRight = 1;
      dirDown = 4;
      dirUp = 5;
    }
  }

  for (unsigned long arrayIndex = 0; arrayIndex < componentBoundary.size(); arrayIndex++) {
    cl = componentBoundary[ arrayIndex ];
    for (fc = 0; fc < 6; fc++)
    {
      nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
      if (nbrfaces[fc]) {
        dx[fc] = Mesh.UCellWidths[ idx2( (nbrfaces[fc] - 1), 0, Mesh.CellWidthsLDI ) ];
        dy[fc] = Mesh.UCellWidths[ idx2( (nbrfaces[fc] - 1), 1, Mesh.CellWidthsLDI ) ];
        dz[fc] = Mesh.UCellWidths[ idx2( (nbrfaces[fc] - 1), 2, Mesh.CellWidthsLDI ) ];
      }
      else {
        dx[fc] = 0;
        dy[fc] = 0;
        dz[fc] = 0;
      }
    }
    if (!nbrfaces[dirIn]) // !dirIn section
    {

    }
    else if (!nbrfaces[dirOut]) // !dirOut section
    {

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
AxisFlowDrive3d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, double visc, int direction, int component )
{

}

