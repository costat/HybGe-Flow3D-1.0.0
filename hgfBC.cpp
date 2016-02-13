#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>
#include <math.h>

#include "hgfMesh.hpp"
#include "hgfBC.hpp"

/* Define a 2d -> 1d array index,
   uses row marjor indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

void
FlowComponent3d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, \
                  const std::vector<unsigned long>& componentBoundary, \
                  const std::vector<double>& componentCellCenters, \
                  const std::vector<double>& componentCellWidths, \
                  const std::vector<double>& CellWidthsLR, \
                  const std::vector<double>& CellWidthsUD, \
                  const std::vector<unsigned long>& componentConnectivity, \
                  const std::vector<unsigned long>& PressureNeighbor, \
                  double visc, int direction, int component, int shift, int pressure )
{
  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  double zmin = Mesh.zLim[0];
  double zmax = Mesh.zLim[1];
  int entries = 0;
  int cl, cl2, fc, colId[8], nbrfaces[6], LR, UD;
  int dirIn, dirOut, dirLeft, dirRight, dirUp, dirDown, colShift;
  double val[8], dxyz[18], componentMin, componentMax, minLR, maxLR, minUD, maxUD;

  if (component == 0) // sets various constants according to which velocity component is being considered
  {
    dirIn = 3;
    dirOut = 1;
    dirLeft = 4;
    dirRight = 5;
    dirDown = 0;
    dirUp = 2;
    LR = 1;
    UD = 2;
    componentMin = xmin;
    componentMax = xmax;
    minLR = ymin;
    maxLR = ymax;
    minUD = zmin;
    maxUD = zmax;
    colShift = 0;
  }
  else if (component == 1)
  {
    dirIn = 5;
    dirOut = 4;
    dirLeft = 3;
    dirRight = 1;
    dirDown = 0;
    dirUp = 2;
    LR = 0;
    UD = 2;
    componentMin = ymin;
    componentMax = ymax;
    minLR = xmin;
    maxLR = xmax;
    minUD = zmin;
    maxUD = zmax;
    if (shift) colShift = Mesh.DOF[1];
    else colShift = 0;
  }
  else
  {
    dirIn = 0;
    dirOut = 2;
    dirLeft = 3;
    dirRight = 1;
    dirDown = 4;
    dirUp = 5;
    LR = 0;
    UD = 1;
    componentMin = zmin;
    componentMax = zmax;
    minLR = xmin;
    maxLR = xmax;
    minUD = ymin;
    maxUD = ymax;
    if (shift) colShift = Mesh.DOF[1] + Mesh.DOF[2];
    else colShift = 0;
  }

  for (unsigned long arrayIndex = 0; arrayIndex < componentBoundary.size(); arrayIndex++)
  {
    cl = componentBoundary[ arrayIndex ];
    cl2 = cl + colShift;
    for (fc = 0; fc < 6; fc++)
    {
      nbrfaces[fc] = componentConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
      if (nbrfaces[fc])
      {
        dxyz[ idx2( fc, 0, 3 ) ] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 0, Mesh.CellWidthsLDI ) ];
        dxyz[ idx2( fc, 1, 3 ) ] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 1, Mesh.CellWidthsLDI ) ];
        dxyz[ idx2( fc, 2, 3 ) ] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 2, Mesh.CellWidthsLDI ) ];
      }
      else {
        dxyz[ idx2( fc, 0, 3 ) ] = 0;
        dxyz[ idx2( fc, 1, 3 ) ] = 0;
        dxyz[ idx2( fc, 2, 3 ) ] = 0;
      }
    }
    if (!nbrfaces[dirIn]) // !dirIn section, no slip regardless of flow direction
    {
      if (component == direction && (componentCellCenters[ idx2( cl, direction, Mesh.CellCentersLDI ) ] \
                                    - 0.5 * componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI ) ]) \
                                    < componentMin) // inflow, nonzero force
      {
        force[cl2] = maxin \
          * (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] - minLR) \
          * (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] - maxLR) \
          * (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] - minUD) \
          * (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] - maxUD);
      }
      colId[0] = cl2;
      val[0] = 1;
      entries = 1;
    }
    else if (!nbrfaces[dirOut]) // !dirOut section
    {
      if (component == direction && (componentCellCenters[ idx2( cl, direction, Mesh.CellCentersLDI ) ] \
                                    + 0.5 * componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI) ]) \
                                    > componentMax) // outflow
      {
        // Set values
        val[0] = -1. / componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI ) ];
        colId[0] = componentConnectivity[ idx2( cl, dirIn, Mesh.FaceConnectivityLDI ) ] - 1 + colShift;
        colId[1] = cl2;
        entries = 2;
      }
      else // noslip
      {
        colId[0] = cl2;
        val[0] = 1;
        entries = 1;
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
            val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                           / (0.5 * (dxyz[ idx2( dirOut, component, 3 )] \
                                  + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
            val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                           / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                                  + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
            colId[0] = nbrfaces[dirOut] - 1 + colShift;
            colId[1] = nbrfaces[dirIn] - 1 + colShift;
            colId[2] = cl2;
            entries = 3;
             // check for outflow possibility in LR direction
            if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                    + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
            {
              val[2] = -(val[0] + val[1]) \
                       + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
            // check for outflow possibility in UD direction
            else if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                         + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
            {
                val[2] = -(val[0] + val[1]) \
                       + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
            else // full no slip
            {
                val[2] = -(val[0] + val[1]) \
                       + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
          }
          else // !dirDown, !dirUp, !dirLeft
          {
            val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                           / (0.5 * (dxyz[ idx2( dirOut, component, 3 )] \
                                  + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
            val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                           / (0.5 * (dxyz[ idx2( dirIn, component, 3 )] \
                                  + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
            val[2] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                          / (0.5 * (dxyz[ idx2( dirRight, LR, 3 )] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
            colId[0] = nbrfaces[dirOut] - 1 + colShift;
            colId[1] = nbrfaces[dirIn] - 1 + colShift;
            colId[2] = nbrfaces[dirRight] - 1 + colShift;
            colId[3] = cl2;
            // check for outflow possibility in UD direction
            if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                    + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]) > maxUD)
            {
              val[3] = -(val[0] + val[1] + val[2]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
            else if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                    + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
            {
              val[3] = -(val[0] + val[1] + val[2]) \
                       + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                  * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                  / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
            else
            {
              val[3] = -(val[0] + val[1] + val[2]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
            }
            entries = 4;
          }
        }
        else if (!nbrfaces[dirRight]) // !dirDown, !dirUp, !dirRight
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                        / (0.5 * (dxyz[ idx2( dirLeft, LR, 3 )] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirLeft] - 1 + colShift;
          colId[3] = cl2;
          // check for outflow possibility in LR direction
          if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                     + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          // check for outflow possibility in UD direction
          else if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                       + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
          {
              val[3] = -(val[0] + val[1] + val[2]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          else // full no slip
          {
              val[3] = -(val[0] + val[1] + val[2]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          entries = 4;
        }
        else // !dirDown, !dirUp
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                        / (0.5 * (dxyz[ idx2( dirLeft, LR, 3 )] \
                               + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
          val[3] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirRight, LR, 3 ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirLeft] - 1 + colShift;
          colId[3] = nbrfaces[dirRight] - 1 + colShift;
          colId[4] = cl2;
          // check for outflow possibility in UD direction
          if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                 + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          else
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          entries = 5;
        }
      }
      else if (!nbrfaces[dirLeft])
      {
        if (!nbrfaces[dirRight]) // !dirDown, !dirLeft, !dirRight
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 )] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                        / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                               + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirUp] - 1 + colShift;
          colId[3] = cl2;
          // check for outflow possibility in LR direction
          if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
          }
          else
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                   + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
          }
          entries = 4;
        }
        else // !dirDown, !dirLeft
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
          val[3] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirRight, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
          if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirUp] - 1 + colShift;
          colId[3] = nbrfaces[dirRight] - 1 + colShift;
          colId[4] = cl2;
          entries = 5;
        }
      }
      else if (!nbrfaces[dirRight]) // !dirDown, !dirRight
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirLeft, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        // check for outflow possibility in LR direction
        if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellCentersLDI ) ]) > maxLR)
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        }
        else
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirUp] - 1 + colShift;
        colId[3] = nbrfaces[dirLeft] - 1 + colShift;
        colId[4] = cl2;
        entries = 5;
      }
      else // !dirDown
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirRight, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirLeft, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                          * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                          / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirUp] - 1 + colShift;
        colId[3] = nbrfaces[dirRight] - 1 + colShift;
        colId[4] = nbrfaces[dirLeft] - 1 + colShift;
        colId[5] = cl2;
        entries = 6;
      }
      // Pressure component to BC
      if (pressure)
      {
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        entries = entries + 2;
      }
    } // end !dirDown primary section
    else if (!nbrfaces[dirUp]) // !dirUp section
    {
      if (!nbrfaces[dirLeft])
      {
        if (!nbrfaces[dirRight]) // !dirUp, !dirLeft, !dirRight
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirDown] - 1 + colShift;
          colId[3] = cl2;
          // check for outflow possibility in LR direction
          if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI)] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI)] \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI)];
          }
          // check for outflow possibility in UD direction
          else if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                       + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                     + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
          }
          else
          {
            val[3] = -(val[0] + val[1] + val[2]) \
                     + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          entries = 4;
        }
        else // !dirUp, !dirLeft
        {
          val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
          val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
          val[3] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                         / (0.5 * (dxyz[ idx2( dirRight, component, 3 ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
          colId[0] = nbrfaces[dirOut] - 1 + colShift;
          colId[1] = nbrfaces[dirIn] - 1 + colShift;
          colId[2] = nbrfaces[dirDown] - 1 + colShift;
          colId[3] = nbrfaces[dirRight] - 1 + colShift;
          colId[4] = cl2;
          // Check for outflow possibility in UD direction
          if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                       + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
          }
          // check for outflow possibility in LR direction
          else if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                  + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI)] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI)];
          }
          else
          {
            val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]
                     + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                                * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                                / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
          }
          entries = 5;
        }
      }
      else if (!nbrfaces[dirRight]) // !dirUp, !dirRight
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirLeft, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirDown] - 1 + colShift;
        colId[3] = nbrfaces[dirLeft] - 1 + colShift;
        colId[4] = cl2;
        // Check for outflow possibility in UD direction
        if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                     + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        else
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        }
        entries = 5;
      }
      else // !dirUp
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirLeft, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirRight, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirDown] - 1 + colShift;
        colId[3] = nbrfaces[dirLeft] - 1 + colShift;
        colId[4] = nbrfaces[dirRight] - 1 + colShift;
        colId[5] = cl2;
        // Check for outflow possibility in UD direction
        if (direction == UD && (componentCellCenters[ idx2( cl, UD, Mesh.CellCentersLDI ) ] \
                                     + componentCellWidths[ idx2( cl, UD, Mesh.CellCentersLDI ) ]) > maxUD)
        {
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
        }
        else
        {
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        }
        entries = 6;
      }
      if (pressure)
      {
        // Pressure component to BC
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        entries = entries + 2;
      }
    } // end !dirUp primary section
    else if (!nbrfaces[dirLeft]) // !dirLeft section
    {
      if (!nbrfaces[dirRight]) // !dirLeft, !dirRight
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirDown] - 1 + colShift;
        colId[3] = nbrfaces[dirUp] - 1 + colShift;
        colId[4] = cl2;
        // Check for outflow possibility in LR direction
        if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        else
        {
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                   + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        entries = 5;
      }
      else // !dirLeft
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                       * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirRight, component, 3 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        // Check for outflow possibility in LR direction
        if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
        {
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
        }
        else
        {
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirDown] - 1 + colShift;
        colId[3] = nbrfaces[dirUp] - 1 + colShift;
        colId[4] = nbrfaces[dirRight] - 1 + colShift;
        colId[5] = cl2;
        entries = 6;
      }
      if (pressure)
      {
        // Pressure component to BC
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        entries = entries + 2;
      }
    } // end !dirLeft primary section
    else if (!nbrfaces[dirRight]) // !dirRight section
    {
      val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirOut, component, 3 ) ]  \
                            + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirIn, component, 3 ) ] \
                            + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                     * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirDown, component, 3 ) ] \
                            + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirUp, component, 3 ) ] \
                            + componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                     * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirLeft, component, 3 ) ] \
                            + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
      colId[0] = nbrfaces[dirOut] - 1 + colShift;
      colId[1] = nbrfaces[dirIn] - 1 + colShift;
      colId[2] = nbrfaces[dirDown] - 1 + colShift;
      colId[3] = nbrfaces[dirUp] - 1 + colShift;
      colId[4] = nbrfaces[dirLeft] - 1 + colShift;
      colId[5] = cl2;
      // Check for outflow possibility in LR direction
      if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
      {
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
      }
      else
      {
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
      }
      entries = 6;
      if (pressure)
      {
        // Pressure component to BC
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                           * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
        entries = entries + 2;
      }
    } // end !dirRight primary section
    // Enter values
    for (int col = 0; col < entries; col++)
    {
      matIs.push_back(cl2);
      matJs.push_back(colId[col]);
      matVals.push_back(val[col]);
    }
  }
}

void
FlowComponent2d ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                  std::vector<int>& matJs, std::vector<double>& matVals, \
                  std::vector<double>& force, \
                  const std::vector<unsigned long>& componentBoundary, \
                  const std::vector<double>& componentCellCenters, \
                  const std::vector<double>& componentCellWidths, \
                  const std::vector<double>& CellWidthsLR, \
                  const std::vector<unsigned long>& componentConnectivity, \
                  const std::vector<unsigned long>& PressureNeighbor, \
                  double visc, int direction, int component, int shift, int pressure
                )
{
  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  int entries = 0;
  int cl, cl2, fc, colId[6], nbrfaces[4], LR;
  int dirIn, dirOut, dirLeft, dirRight, colShift;
  double val[8], dxyz[8], componentMin, componentMax, minLR, maxLR;

  if (component == 0) // sets various constants according to which velocity component is being considered
  {
    dirIn = 3;
    dirOut = 1;
    dirLeft = 2;
    dirRight = 0;
    LR = 1;
    componentMin = xmin;
    componentMax = xmax;
    minLR = ymin;
    maxLR = ymax;
    colShift = 0;
  }
  else
  {
    dirIn = 0;
    dirOut = 2;
    dirLeft = 3;
    dirRight = 1;
    LR = 0;
    componentMin = ymin;
    componentMax = ymax;
    minLR = xmin;
    maxLR = xmax;
    if (shift) colShift = Mesh.DOF[1];
    else colShift = 0;
  }

  for (unsigned long arrayIndex = 0; arrayIndex < componentBoundary.size(); arrayIndex++)
  { // loop over component boundary cells
    cl = componentBoundary[ arrayIndex ];
    cl2 = cl + colShift;
    for (fc = 0; fc < 4; fc++)
    {
      nbrfaces[fc] = componentConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
      if (nbrfaces[fc])
      {
        dxyz[ idx2( fc, 0, 2 ) ] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 0, Mesh.CellWidthsLDI ) ];
        dxyz[ idx2( fc, 1, 2 ) ] = componentCellWidths[ idx2( (nbrfaces[fc] - 1), 1, Mesh.CellWidthsLDI ) ];
      }
      else
      {
        dxyz[ idx2( fc, 0, 2 ) ] = 0;
        dxyz[ idx2( fc, 1, 2 ) ] = 0;
      }
    }
    if (!nbrfaces[dirIn]) // !dirIn section, no slip regardless of flow directions
    {
      if (component == direction && (componentCellCenters[ idx2( cl, direction, Mesh.CellCentersLDI ) ] \
                                    - 0.6 * componentCellWidths[ idx2( cl, direction, Mesh.CellCentersLDI ) ]) \
                                    < componentMin) // inflow, nonzero force
      {
        force[cl2] = -maxin \
          * (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] - minLR) \
          * (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] - maxLR);
      }
      colId[0] = cl2;
      val[0] = 1;
      entries = 1;
    }
    else if (!nbrfaces[dirOut]) // !dirOut section, no slip regardless of flow direction
    {
      if (component == direction && (componentCellCenters[ idx2( cl, direction, Mesh.CellCentersLDI ) ] \
                                    + 0.5 * componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI ) ]) \
                                    > componentMax) // outflow
      {
        val[0] = -1. / componentCellWidths[ idx2( cl, direction, Mesh.CellWidthsLDI ) ];
        val[1] = -val[0];
        colId[0] = componentConnectivity[ idx2( cl, dirIn, Mesh.FaceConnectivityLDI ) ] - 1 + colShift;
        colId[1] = cl2;
        entries = 2;
      }
      else // no slip
      {
        colId[0] = cl2;
        val[0] = 1;
        entries = 1;
      }
    }
    else if (!nbrfaces[dirLeft]) // !dirLeft section
    {
      if (!nbrfaces[dirRight]) // !dirLeft, !dirRight
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 2 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 2 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = cl2;
        // check for outflow possibility in LR direction
        if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                                + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
        {
          val[2] = -(val[0] + val[1]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        else
        {
          val[2] = -(val[0] + val[1]) \
                   + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        entries = 3;
      }
      else // !dirLeft
      {
        val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirOut, component, 2 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                       / (0.5 * (dxyz[ idx2( dirIn, component, 2 ) ] \
                              + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI )] \
                       / (0.5 * (dxyz[ idx2( dirRight, component, 2 ) ] \
                              + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
        // check for outflow possibility in LR direction
        if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                               + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
        {
          val[3] = -(val[0] + val[1] + val[2]);
        }
        else
        {
          val[3] = -(val[0] + val[1] + val[2]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        colId[0] = nbrfaces[dirOut] - 1 + colShift;
        colId[1] = nbrfaces[dirIn] - 1 + colShift;
        colId[2] = nbrfaces[dirRight] - 1 + colShift;
        colId[3] = cl2;
        entries = 4;
      }
      if (pressure)
      {
        // Pressure component
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        entries = entries + 2;
      }
    }
    else if (!nbrfaces[dirRight]) // !dirRight section
    {
      val[0] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirOut, component, 2 ) ] \
                            + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirIn, component, 2 ) ] \
                            + componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dxyz[ idx2( dirLeft, component, 2 ) ] \
                            + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]));
      colId[0] = nbrfaces[dirOut] - 1 + colShift;
      colId[1] = nbrfaces[dirIn] - 1 + colShift;
      colId[2] = nbrfaces[dirLeft] - 1 + colShift;
      colId[3] = cl2;
      entries = 4;
      // check for outflow possibility in LR direction
      if (direction == LR && (componentCellCenters[ idx2( cl, LR, Mesh.CellCentersLDI ) ] \
                             + componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ]) > maxLR)
      {
        val[3] = -(val[0] + val[1] + val[2]);
      }
      else
      {
        val[3] = -(val[0] + val[1] + val[2]) \
                 + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                            / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
      }
      if (pressure)
      {
        // Pressure component
        val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        val[entries+1] = -val[entries];
        colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                         - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        entries = entries + 2;
      }
    }
    // Enter values
    for (int col = 0; col < entries; col++)
    {
      matIs.push_back(cl2);
      matJs.push_back(colId[col]);
      matVals.push_back(val[col]);
    }
  } // end loop over component boundary cells
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
      // Set U boundary conditions
      FlowComponent2d( Mesh, matIs, matJs, matVals, \
                       force, Mesh.UBoundaryCells, Mesh.UCellCenters, Mesh.UCellWidths, \
                       Mesh.VCellWidths, Mesh.UFaceConnectivity, Mesh.UCellPressureNeighbor, \
                       visc, direction, 0, 1, 1 );
      // Set V boundary conditions
      FlowComponent2d( Mesh, matIs, matJs, matVals, \
                       force, Mesh.VBoundaryCells, Mesh.VCellCenters, Mesh.VCellWidths, \
                       Mesh.UCellWidths, Mesh.VFaceConnectivity, Mesh.VCellPressureNeighbor, \
                       visc, direction, 1, 1, 1 );
      break;
    }
    default :
    {
      // Set U boundary conditions
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, Mesh.UBoundaryCells, Mesh.UCellCenters, Mesh.UCellWidths, \
                       Mesh.VCellWidths, Mesh.WCellWidths, Mesh.UFaceConnectivity, Mesh.UCellPressureNeighbor, \
                       visc, direction, 0, 1, 1 );
      // Set V boundary conditions
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, Mesh.VBoundaryCells, Mesh.VCellCenters, Mesh.VCellWidths, \
                       Mesh.UCellWidths, Mesh.WCellWidths, Mesh.VFaceConnectivity, Mesh.VCellPressureNeighbor, \
                       visc, direction, 1, 1, 1 );
      // Set W boundary conditions
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, Mesh.WBoundaryCells, Mesh.WCellCenters, Mesh.WCellWidths, \
                       Mesh.UCellWidths, Mesh.VCellWidths, Mesh.WFaceConnectivity, Mesh.WCellPressureNeighbor, \
                       visc, direction, 2, 1, 1 );
      break;
    }
  }
}

void
AxisFlowSingleComponent ( const FluidMesh& Mesh, std::vector<int>& matIs, \
                          std::vector<int>& matJs, std::vector<double>& matVals, \
                          std::vector<double>& force, double visc, int direction, int component )
{

  switch (Mesh.DIM)
  {
    case 2 :
    {
      switch ( component )
      {
        case 0 :
        {
          FlowComponent2d( Mesh, matIs, matJs, matVals, \
                           force, Mesh.UBoundaryCells, Mesh.UCellCenters, Mesh.UCellWidths, \
                           Mesh.VCellWidths, Mesh.UFaceConnectivity, Mesh.UCellPressureNeighbor, \
                           visc, direction, 0, 0, 0 );

          break; // u component break;
        }
        case 1 :
        {
          FlowComponent2d( Mesh, matIs, matJs, matVals, \
                           force, Mesh.VBoundaryCells, Mesh.VCellCenters, Mesh.VCellWidths, \
                           Mesh.UCellWidths, Mesh.VFaceConnectivity, Mesh.VCellPressureNeighbor, \
                           visc, direction, 1, 0, 0 );
          break; // v component break;
        }
      }
      break; // dim 2 break;
    }
    case 3 :
    {
      switch ( component )
      {
        case 0 :
        {
          // Set U boundary conditions
          FlowComponent3d( Mesh, matIs, matJs, matVals, \
                           force, Mesh.UBoundaryCells, Mesh.UCellCenters, Mesh.UCellWidths, \
                           Mesh.VCellWidths, Mesh.WCellWidths, Mesh.UFaceConnectivity, Mesh.UCellPressureNeighbor, \
                           visc, direction, 0, 0, 0 );
          break; // u component break;
        }
        case 1 :
        {
          // Set V boundary conditions
          FlowComponent3d( Mesh, matIs, matJs, matVals, \
                           force, Mesh.VBoundaryCells, Mesh.VCellCenters, Mesh.VCellWidths, \
                           Mesh.UCellWidths, Mesh.WCellWidths, Mesh.VFaceConnectivity, Mesh.VCellPressureNeighbor, \
                           visc, direction, 1, 0, 0 );
          break; // v component break;
        }
        case 2 :
        {
          // Set U boundary conditions
          FlowComponent3d( Mesh, matIs, matJs, matVals, \
                           force, Mesh.WBoundaryCells, Mesh.WCellCenters, Mesh.WCellWidths, \
                           Mesh.UCellWidths, Mesh.VCellWidths, Mesh.WFaceConnectivity, Mesh.WCellPressureNeighbor, \
                           visc, direction, 2, 0, 0 );
          break; // w component break;
        }
      }
      break; // dim 3 break;
    }
  }

}
int
FindPeriodicPair( const PoreNetwork& pn, int pore, int dir, int side )
{
  int potpp;
  double dt;
  if (dir == 0) dt = pn.dx;
  else if (dir == 1) dt = pn.dy;
  else dt = pn.dz;
  for (int potp = 0; potp < (pn.BoundaryPores.size()); potp++)
  {
    potpp = pn.BoundaryPores[ potp ];
    if (!pn.Throats[ idx2( potpp, side, pn.DIM*2 ) ]) {
      if (fabs(pn.PoresXYZ[ idx2( pore, dir, pn.DIM ) ] - pn.PoresXYZ[ idx2( potpp, dir, pn.DIM ) ]) < 0.2*dt) {
        return potpp;
      }
    }
  }
}
void
PoreNetworkBoundary( const PoreNetwork& pn, std::vector<int>& matIs, \
                     std::vector<int>& matJs, std::vector<double>& matVals, \
                     std::vector<double>& force, \
                     const std::vector<double>& Ks, int direction )
{
  double *val = new double[ 2*pn.DIM + 1 ];
  int *colId = new int[ 2*pn.DIM + 1 ];
  switch (pn.DIM)
  {
    case 2 :
    {
      int dirIn, dirOut, dirLeft, dirRight, dirLR;
      if (direction == 0)
      {
        dirIn = 3;
        dirOut = 1;
        dirRight = 0;
        dirLeft = 2;
        dirLR = 1;
      }
      else if (direction == 1)
      {
        dirIn = 0;
        dirOut = 2;
        dirLeft = 3;
        dirRight = 1;
        dirLR = 0;
      }

      int pore, ppore, entries;
      for (int pi = 0; pi < pn.BoundaryPores.size(); pi++)
      {
        pore = pn.BoundaryPores[ pi ];
        if (!pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
            // this is partly a periodic pore, need to find pair node
            ppore = FindPeriodicPair( pn, pore, direction, dirRight );
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[2] = ppore;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
          else if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
            // this is partly a periodic pore, need to find pair node
            ppore = FindPeriodicPair( pn, pore, direction, dirLeft );
            colId[0] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[2] = ppore;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
          else {
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[2] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
        }
        else if (!pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
            // this is partly a periodic pore, need to find pair node
            ppore = FindPeriodicPair( pn, pore, direction, dirRight );
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]-1;
            colId[2] = ppore;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
          else if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
            // this is partly a periodic pore, need to find pair node
            ppore = FindPeriodicPair( pn, pore, direction, dirLeft );
            colId[0] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]-1;
            colId[2] = ppore;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
          else {
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]-1;
            colId[2] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[3] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
            val[3] = -(val[0] + val[1] + val[2]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 4;
          }
        }
        else if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
          // this is a periodic pore, need to find pair node
          ppore = FindPeriodicPair( pn, pore, direction, dirLeft );
          colId[0] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
          colId[1] = pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]-1;
          colId[2] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
          colId[3] = ppore;
          colId[4] = pore;
          val[0] = -0.5 * ( Ks[ idx2( colId[0], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
          val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
          val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
          val[3] = -0.5 * ( Ks[ idx2( colId[3], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
          val[4] = -(val[0] + val[1] + val[2] + val[3]);
          entries = 5;
        }
        else if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
          // this is a periodic pore, need to find pair node
          ppore = FindPeriodicPair( pn, pore, direction, dirRight );
          colId[0] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
          colId[1] = pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]-1;
          colId[2] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
          colId[3] = ppore;
          colId[4] = pore;
          val[0] = -0.5 * ( Ks[ idx2( colId[0], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
          val[1] = -0.5 * ( Ks[ idx2( colId[1], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
          val[2] = -0.5 * ( Ks[ idx2( colId[2], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
          val[3] = -0.5 * ( Ks[ idx2( colId[3], dirLR, pn.DIM ) ] + Ks[ idx2( pore, dirLR, pn.DIM ) ] );
          val[4] = -(val[0] + val[1] + val[2] + val[3]);
          entries = 5;
        }
        for (int dir = 0; dir < entries; dir++) {
          matIs.push_back( pore );
          matJs.push_back( colId[dir] );
          matVals.push_back( val[dir] );
        }
      }
      break;
    }
    case 3 :
    {
      int dirIn, dirOut, dirLeft, dirRight, dirDown, dirUp, LR, UD;
      if (direction == 0) {
        dirIn = 3;
        dirOut = 1;
        dirLeft = 4;
        dirRight = 5;
        dirDown = 0;
        dirUp = 2;
        LR = 1;
        UD = 2;
      }
      else if (direction == 1) {
        dirIn = 5;
        dirOut = 4;
        dirLeft = 3;
        dirRight = 1;
        dirDown = 0;
        dirUp = 2;
        LR = 0;
        UD = 2;
      }
      else {
        dirIn = 0;
        dirOut = 2;
        dirLeft = 3;
        dirRight = 1;
        dirDown = 4;
        dirUp = 5;
        LR = 0;
        UD = 1;
      }
      int pore, ppore1, ppore2, entries;
      for (int pi = 0; pi < pn.BoundaryPores.size(); pi++) {
        pore = pn.BoundaryPores[ pi ];
        if (!pn.Throats[ idx2( pore, dirIn, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
            if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
              // !in, !right, !down
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );
              colId[0] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[3] = ppore1;
              colId[4] = ppore2;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
            else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
              // !in, !right, !up
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );
              colId[0] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[3] = ppore1;
              colId[4] = ppore2;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
            else {
              // !in !right
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              colId[0] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
              colId[3] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[4] = ppore1;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
          }
          else if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
            if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
              // !in, !left, !down
              ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );
              colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[3] = ppore1;
              colId[4] = ppore2;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
            else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
              // !in, !left, !up
              ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );
              colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[3] = ppore1;
              colId[4] = ppore2;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
            else {
              // !in !left
              ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );
              colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
              colId[1] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
              colId[2] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
              colId[3] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
              colId[4] = ppore1;
              colId[5] = pore;
              val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[1] = -0.5 * ( Ks[ idx2( colId[1], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
              val[2] = -0.5 * ( Ks[ idx2( colId[2], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] ):
              val[3] = -0.5 * ( Ks[ idx2( colId[3], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
              val[4] = -0.5 * ( Ks[ idx2( colId[4], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
              val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                       0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
              entries = 6;
            }
          }
          else if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
            // !in !down
            ppore1 = FindPeriodicPair( pn, pore, direction, dirUp );
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[2] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
            colId[3] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[4] = ppore1;
            colId[5] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] ):
            val[3] = -0.5 * ( Ks[ idx2( colId[3], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
            val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 6;
          }
          else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
            // !in !up
            ppore1 = FindPeriodicPair( pn, pore, direction, dirDown );
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[2] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
            colId[3] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[4] = ppore1;
            colId[5] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] ):
            val[3] = -0.5 * ( Ks[ idx2( colId[3], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
            val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 6;
          }
          else {
            // !in
            colId[0] = pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]-1;
            colId[1] = pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]-1;
            colId[2] = pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]-1;
            colId[3] = pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]-1;
            colId[4] = pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]-1;
            colId[5] = pore;
            val[0] = -0.5 * ( Ks[ idx2( colId[0], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[1] = -0.5 * ( Ks[ idx2( colId[1], LR, pn.DIM ) ] + Ks[ idx2( pore, LR, pn.DIM ) ] );
            val[2] = -0.5 * ( Ks[ idx2( colId[2], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] ):
            val[3] = -0.5 * ( Ks[ idx2( colId[3], direction, pn.DIM ) ] + Ks[ idx2( pore, direction, pn.DIM ) ] );
            val[4] = -0.5 * ( Ks[ idx2( colId[4], UD, pn.DIM ) ] + Ks[ idx2( pore, UD, pn.DIM ) ] );
            val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) + \
                     0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            force[pore] = force[pore] + 0.5 * Ks[ idx2( pore, direction, pn.DIM ) ];
            entries = 6;
          }
        } // end dirIn primary section
        else if (!pn.Throats[ idx2( pore, dirOut, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
            if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
              // !out, !right, !down
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );

            }
            else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
              // !out, !right, !up
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );

            }
            else {
              // !out !right
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );

            }
          }
          else if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
            if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
              // !out, !right, !down
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );

            }
            else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
              // !out, !right, !up
              ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
              ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );

            }
            else {
              // !out !left
              pore1 = FindPeriodicPair( pn, pore, direction, dirRight );

            }
          }
          else if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
            // !out !down
            ppore1 = FindPeriodicPair( pn, pore, direction, dirUp );

          }
          else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
            // !out !up
            ppore1 = FindPeriodicPair( pn, pore, direction, dirDown );

          }
        } // end dirOut primary section
        else if (!pn.Throats[ idx2( pore, dirRight, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
            // !right !down
            ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
            ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );

          }
          else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
            // !right !up
            ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );
            ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );

          }
          else {
            // !right
            ppore1 = FindPeriodicPair( pn, pore, direction, dirLeft );

          }
        } // end dirRight primary section
        else if (!pn.Throats[ idx2( pore, dirLeft, (pn.DIM*2) ) ]) {
          if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
            // !left !down
            ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );
            ppore2 = FindPeriodicPair( pn, pore, direction, dirUp );

          }
          else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
            // !left !up
            ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );
            ppore2 = FindPeriodicPair( pn, pore, direction, dirDown );

          }
          else {
            // !left
            ppore1 = FindPeriodicPair( pn, pore, direction, dirRight );

          }
        } // end dirLeft primary section
        else if (!pn.Throats[ idx2( pore, dirDown, (pn.DIM*2) ) ]) {
          ppore1 = FindPeriodicPair( pn, pore, direction, dirUp );

        }
        else if (!pn.Throats[ idx2( pore, dirUp, (pn.DIM*2) ) ]) {
          ppore1 = FindPeriodicPair( pn, pore, direction, dirDown );

        }
        for (int dir = 0; dir < entries; dir++) {
          matIs.push_back( pore );
          matJs.push_back( colId[dir] );
          matVals.push_back( val[dir] );
        }
      } // end loop bounary pores

      break;
    }
  }
  delete[] val;
  delete[] colId;
}
