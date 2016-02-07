#include <iostream>
#include <vector>
#include <cstdlib>
#include <string.h>

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
                  double visc, int direction, int component )
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
    colShift = Mesh.DOF[1];
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
    colShift = Mesh.DOF[1] + Mesh.DOF[2];
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
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      entries = entries + 2;
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
      // Pressure component to BC
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      entries = entries + 2;
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
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                   + 2 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
                              * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ] \
                              / componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
        }
        else
        {
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                   + 4 * visc * componentCellWidths[ idx2( cl, component, Mesh.CellWidthsLDI ) ] \
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
      // Pressure component to BC
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      entries = entries + 2;
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
      // Pressure component to BC
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ] \
                         * componentCellWidths[ idx2( cl, UD, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      entries = entries + 2;
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
                  double visc, int direction, int component
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
    colShift = Mesh.DOF[1];
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
      // Pressure component
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      entries = entries + 2;
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
      // Pressure component
      val[entries] = -1. * componentCellWidths[ idx2( cl, LR, Mesh.CellWidthsLDI ) ];
      val[entries+1] = -val[entries];
      colId[entries] = PressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[entries+1] = PressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                       - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      entries = entries + 2;
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
      const std::vector<unsigned long>& componentBoundary = Mesh.UBoundaryCells;
      const std::vector<double>& componentCellWidths = Mesh.UCellWidths;
      const std::vector<double>& componentCellCenters = Mesh.UCellCenters;
      const std::vector<unsigned long>& componentConnectivity = Mesh.UFaceConnectivity;
      const std::vector<double>& cellWidthsLR = Mesh.VCellWidths;
      const std::vector<unsigned long>& PressureNeighbor = Mesh.UCellPressureNeighbor;
      FlowComponent2d( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary, componentCellCenters, componentCellWidths, \
                       cellWidthsLR, componentConnectivity, PressureNeighbor, \
                       visc, direction, 0 );
      // Set V boundary conditions
      const std::vector<unsigned long>& componentBoundary1 = Mesh.VBoundaryCells;
      const std::vector<double>& componentCellWidths1 = Mesh.VCellWidths;
      const std::vector<double>& componentCellCenters1 = Mesh.VCellCenters;
      const std::vector<unsigned long>& componentConnectivity1 = Mesh.VFaceConnectivity;
      const std::vector<double>& cellWidthsLR1 = Mesh.UCellWidths;
      const std::vector<unsigned long>& PressureNeighbor1 = Mesh.VCellPressureNeighbor;
      FlowComponent2d( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary1, componentCellCenters1, componentCellWidths1, \
                       cellWidthsLR1, componentConnectivity1, PressureNeighbor1, \
                       visc, direction, 1 );
      break;
    }
    default :
    {
      // Set U boundary conditions
      const std::vector<unsigned long>& componentBoundary = Mesh.UBoundaryCells;
      const std::vector<double>& componentCellWidths = Mesh.UCellWidths;
      const std::vector<double>& componentCellCenters = Mesh.UCellCenters;
      const std::vector<unsigned long>& componentConnectivity = Mesh.UFaceConnectivity;
      const std::vector<double>& cellWidthsLR = Mesh.VCellWidths;
      const std::vector<double>& cellWidthsUD = Mesh.WCellWidths;
      const std::vector<unsigned long>& PressureNeighbor = Mesh.UCellPressureNeighbor;
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary, componentCellCenters, componentCellWidths, \
                       cellWidthsLR, cellWidthsUD, componentConnectivity, PressureNeighbor, \
                       visc, direction, 0 );
      // Set V boundary conditions
      const std::vector<unsigned long>& componentBoundary1 = Mesh.VBoundaryCells;
      const std::vector<double>& componentCellWidths1 = Mesh.VCellWidths;
      const std::vector<double>& componentCellCenters1 = Mesh.VCellCenters;
      const std::vector<unsigned long>& componentConnectivity1 = Mesh.VFaceConnectivity;
      const std::vector<double>& cellWidthsLR1 = Mesh.UCellWidths;
      const std::vector<double>& cellWidthsUD1 = Mesh.WCellWidths;
      const std::vector<unsigned long>& PressureNeighbor1 = Mesh.VCellPressureNeighbor;
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary1, componentCellCenters1, componentCellWidths1, \
                       cellWidthsLR1, cellWidthsUD1, componentConnectivity1, PressureNeighbor1, \
                       visc, direction, 1 );
      // Set W boundary conditions
      const std::vector<unsigned long>& componentBoundary2 = Mesh.WBoundaryCells;
      const std::vector<double>& componentCellWidths2 = Mesh.WCellWidths;
      const std::vector<double>& componentCellCenters2 = Mesh.WCellCenters;
      const std::vector<unsigned long>& componentConnectivity2 = Mesh.WFaceConnectivity;
      const std::vector<double>& cellWidthsLR2 = Mesh.UCellWidths;
      const std::vector<double>& cellWidthsUD2 = Mesh.VCellWidths;
      const std::vector<unsigned long>& PressureNeighbor2 = Mesh.WCellPressureNeighbor;
      FlowComponent3d( Mesh, matIs, matJs, matVals, \
                       force, componentBoundary2, componentCellCenters2, componentCellWidths2, \
                       cellWidthsLR2, cellWidthsUD2, componentConnectivity2, PressureNeighbor2, \
                       visc, direction, 2 );
      break;
    }
  }
}
