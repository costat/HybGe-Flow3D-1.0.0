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

/* Define a 2d . 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

/* ucartflow implements boundary conditions in rows of the linear system
   associated to u components of boundary voxels. It implements
   a flow problem in a principal axis direction. */
void
ucartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction )
{
  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  double zmin = Mesh.zLim[0];
  double zmax = Mesh.zLim[1];
  int cl, fc, colId[6], nbrfaces[6];
  double val[7], dx[6], dy[6], dz[6];

  /* u boundary conditions */
  for (unsigned long arrayIndex = 0; arrayIndex < Mesh.UBoundaryCells.size(); arrayIndex++) {
    cl = Mesh.UBoundaryCells[arrayIndex];
    if (!Mesh.UFaceConnectivity[ idx2( cl, 1, Mesh.FaceConnectivityLDI ) ]) {
      // Right boundary
      if (direction == 0 && Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] == xmax) { // Outflow condition for x principal flow direction
        val[0] = -1. / Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        colId[0] = Mesh.UFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI ) ] - 1;
        colId[1] = cl;
        // Insert values
        matIs.push_back(cl);
        matJs.push_back(colId[0]);
        matVals.push_back(val[0]);
        matIs.push_back(cl);
        matJs.push_back(colId[1]);
        matVals.push_back(val[1]);
      }
      else { // No slip for y and z principal flow directions
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(1);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI ) ]) {
      // Left boundary
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(1);
      if (direction == 0 && Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] == xmin) {
        force[cl] = maxin \
          * (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymin) \
          * (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymax) \
          * (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] - zmin) \
          * (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] - zmax);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ]) {
      // bottom boundary, no slip
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.UCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[4]) {
        /* bottom back corner */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                    * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                    / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                    * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                    / (0.5 * (dz[2] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                    * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                    / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                    * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                    / (0.5 * (dy[5] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 1 && (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] \
                               + Mesh.UCellWidths[ cl ]) > ymax) { // In yflow problem, there is outflow on back
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // For xflow and zflow we just have no slip
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[1] - 1;
        colId[1] = nbrfaces[2] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[5] - 1;
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[4]);
      }
      else if (!nbrfaces[5]) {
        /* bottom front corner */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                       * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                       * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[2] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                       * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                       * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[4] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -(val[0] + val[1] + val[2] + val[3]) \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];

        // Columns
        colId[0] = nbrfaces[1] - 1;
        colId[1] = nbrfaces[2] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[4] - 1;
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[4]);
      }
      else {
        /* bottom boundary, not a corner */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]
                  * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[5] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];

        // Columns
        colId[0] = nbrfaces[1] - 1;
        colId[1] = nbrfaces[2] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[4] - 1;
        colId[4] = nbrfaces[5] - 1;
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[5]);
      }
      // Pressure component to BC
      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 2, Mesh.FaceConnectivityLDI ) ]) {
      /* top boundary, no slip or outflow */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.UCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[4]) {
        /* top back corner */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[5] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0) { // full no slip
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        else if (direction == 1 && (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] \
                                    + Mesh.UCellWidths[ cl ]) > ymax) { // outflow in y direction
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        else if (direction == 2 && (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                                    + Mesh.UCellWidths[ cl ]) > zmax) { // outflow in z direction
           val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1;
        colId[1] = nbrfaces[1] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[5] - 1;
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[4]);
      }
      else if (!nbrfaces[5]) {
        /* top front corner, no slip or outflow */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[4] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 2 && (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                               + Mesh.UCellWidths[ cl ]) > zmax) { // Outflow at top
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // Full no slip
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1;
        colId[1] = nbrfaces[1] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[4] - 1;
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[4]);

      }
      else {
        /* top boundary no corner */
        val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[4] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[5] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 2 && (Mesh.UCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                               + Mesh.UCellWidths[ cl ]) > zmax) { // outflow condition
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
        }
        else { // No slip
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                 + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1;
        colId[1] = nbrfaces[1] - 1;
        colId[2] = nbrfaces[3] - 1;
        colId[3] = nbrfaces[4] - 1;
        colId[4] = nbrfaces[5] - 1;
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(val[5]);

      }
      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      // Insert pressure values
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }

    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 5, Mesh.FaceConnectivityLDI ) ]) {
      /* front boundary, no slip */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.UCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dz[0] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dz[2] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dy[4] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
             + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                        * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                        / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];

      // Columns
      colId[0] = nbrfaces[0] - 1;
      colId[1] = nbrfaces[1] - 1;
      colId[2] = nbrfaces[2] - 1;
      colId[3] = nbrfaces[3] - 1;
      colId[4] = nbrfaces[4] - 1;
      for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(val[5]);

      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];

      // Insert pressure values
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }

    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 4, Mesh.FaceConnectivityLDI ) ]) {
      /* back boundary, no slip or outflow */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.UCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dz[0] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dz[2] + Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
               / (0.5 * (dy[5] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      if (direction == 1 && (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] \
                             + Mesh.UCellWidths[ cl ]) > ymax) { // outflow for yflow problem
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
      }
      else { // no slip for all others
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      }
      // Columns
      colId[0] = nbrfaces[0] - 1;
      colId[1] = nbrfaces[1] - 1;
      colId[2] = nbrfaces[2] - 1;
      colId[3] = nbrfaces[3] - 1;
      colId[4] = nbrfaces[5] - 1;
      for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(val[5]);

      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.UCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];

      // Insert pressure values
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
  }
}

/* vcartflow implements boundary conditions in rows of the linear system
   associated to v components of boundary voxels. It implements
   a flow problem in a principal axis direction. */
void
vcartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction )
{
  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  double zmin = Mesh.zLim[0];
  double zmax = Mesh.zLim[1];
  int cl, cl2, fc, colId[6], nbrfaces[6];
  double val[7], dx[6], dy[6], dz[6];

  /* v boundary conditions */
  for (unsigned long arrayIndex = 0; arrayIndex < Mesh.VBoundaryCells.size(); arrayIndex++) {
    cl = Mesh.VBoundaryCells[arrayIndex];
    cl2 = cl + Mesh.DOF[1];
    if (!Mesh.VFaceConnectivity[ idx2( cl, 4, Mesh.FaceConnectivityLDI ) ]) {
      // back boundary, outflow or noslip
      if (direction == 1 && Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] == ymax) { // Outflow boundary
        val[0] = -1. / Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        colId[0] = Mesh.VFaceConnectivity[ idx2( cl, 5, Mesh.FaceConnectivityLDI ) ] \
                   - 1 + Mesh.DOF[1];
        colId[1] = cl2;
        // Insert values
        matIs.push_back(cl2);
        matJs.push_back(colId[0]);
        matVals.push_back(val[0]);
        matIs.push_back(cl2);
        matJs.push_back(colId[1]);
        matVals.push_back(val[1]);
      }
      else { // no slip
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(1);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2( cl, 5, Mesh.FaceConnectivityLDI ) ]) {
      // front boundary, inflow or no slip
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(1);
      if (direction == 1 && Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] == ymin) { // inflow for yflow problem
        force[cl2] = maxin * (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmin) \
                         * (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmax) \
                         * (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] - zmin) \
                         * (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] - zmax);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2( cl, 1, Mesh.FaceConnectivityLDI ) ]) {
      // right boundary, outflow or no slip
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.VCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[2]) {
        // top right corner
        val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.VCellWidths[ cl ]) > xmax) { // outflow to the right, no slip on top
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        else if (direction == 2 && (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                                    + Mesh.VCellWidths[ cl ]) > zmax) { // no slip to the right, outflow on top
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[3] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else if (!nbrfaces[0]) {
        // bottom right corner, no slip or outflow
        val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] ));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.VCellWidths[ cl ]) > xmax) { // outflow to right, no slip on bottom
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[2] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[3] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else {
        // right boundary, not a corner, no slip or outflow
        val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.VCellWidths[ cl ]) > xmax) { // outflow
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
        }
        else { // no slip
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[5]);
      }
      // Pressure component
      val[0] = -1. * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI ) ]) {
      /* left boundary */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.VCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[2]) {
        // top left corner, no slips or outflows
        val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 2 && (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                              + Mesh.VCellWidths[ cl ]) > zmax) { // outflow on top, no slip on left
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else if (!nbrfaces[0]) {
        // bottom left corner, no slip
        val[0] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -(val[0] + val[1] + val[2] + val[3]) \
               + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        // Columns
        colId[0] = nbrfaces[1] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else {
        /* left boundary no corner */
        val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
        colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1];
        colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1];
        colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1];
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[5]);
      }
      // Pressure components
      val[0] = -1. * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ]) {
      /* bottom boundary */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.VCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[2] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
             + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                        * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                        / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      // Columns
      colId[0] = nbrfaces[1] - 1 + Mesh.DOF[1];
      colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1];
      colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1];
      colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1];
      colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1];
      for (fc = 0; fc < 5; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[5]);
      // Pressure component
      val[0] = -1. * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2( cl, 2, Mesh.FaceConnectivityLDI ) ]) {
      // top boundary, no slip or outflow
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.VCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[0] + Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI) ] \
                / (0.5 * (dy[4] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[5] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      if (direction == 2 && (Mesh.VCellCenters[ idx2( cl, 2, Mesh.CellCentersLDI ) ] \
                              + Mesh.VCellWidths[ cl ]) > zmax) { // outflow for z flow problem
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
      }
      else { // no slip for x and y flows
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      }
      // Columns
      colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
      colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1];
      colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1];
      colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1];
      colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1];
      for (fc = 0; fc < 5; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[5]);
      // Pressure component
      val[0] = -1. * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.VCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
  }
}

/* wcartflow implements boundary conditions in rows of the linear system
   associated to u components of boundary voxels. It implements
   a flow problem in a principal axis direction. */
void
wcartflow ( const FluidMesh& Mesh, std::vector<int>& matIs, \
            std::vector<int>& matJs, std::vector<double>& matVals, \
            std::vector<double>& force, \
            double visc, int direction )
{
  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  double zmin = Mesh.zLim[0];
  double zmax = Mesh.zLim[1];
  int cl, cl2, fc, colId[6], nbrfaces[6];
  double val[7], dx[6], dy[6], dz[6];

  /* w boundary conditions */
  for (unsigned long arrayIndex = 0; arrayIndex < Mesh.WBoundaryCells.size(); arrayIndex++) {
    cl = Mesh.WBoundaryCells[arrayIndex];
    cl2 = cl + Mesh.DOF[1] + Mesh.DOF[2];
    if (!Mesh.WFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ]) {
      // bottom boundary, inflow or no slip
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(1);
      if (direction == 2 && Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellWidthsLDI ) ] == zmin) { // inflow force for z flow problem
        force[cl2] = maxin * (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmin) \
                               * (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmax) \
                               * (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymin) \
                               * (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymax);
      }
    }
    else if (!Mesh.WFaceConnectivity[ idx2( cl, 2, Mesh.FaceConnectivityLDI ) ]) {
      // top boundary, outflow or no slip
      if (direction == 2 && Mesh.WCellCenters[ idx2( cl, 2, Mesh.CellWidthsLDI ) ] == zmax) { // outflow
        val[0] = -1. / Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ];
        colId[0] = Mesh.WFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ] \
                   - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = cl2;
        // Insert values
        matIs.push_back(cl2);
        matJs.push_back(colId[0]);
        matVals.push_back(val[0]);
        matIs.push_back(cl2);
        matJs.push_back(colId[1]);
        matVals.push_back(val[1]);
      }
      else { // no slip
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(1);
      }
    }
    else if (!Mesh.WFaceConnectivity[ idx2( cl, 1, Mesh.FaceConnectivityLDI ) ]) {
      // right boundary
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.WFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.WCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.WCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.WCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[4]) {
        // back right corner, outflows or no slips
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.WCellWidths[ cl ]) > xmax) { // outflow to right, no slip to the back
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        else if (direction == 1 && (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellWidthsLDI ) ] \
                                    + Mesh.WCellWidths[ cl ]) > ymax) { // outflow to the back, no slip on right
           val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                  + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                             * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                             / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);

      }
      else if (!nbrfaces[5]) {
        // front right corner, noslip or outflows
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[3] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.WCellWidths[ cl ]) > xmax) { // outflow to the right
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else {
        // right boundary, not a corner, outflow or no slip
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[3] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[4] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[5] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 0 && (Mesh.WCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.WCellWidths[ cl ]) > xmax) { // outflow
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
        }
        else { // no slip
          val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[5]);
      }
      val[0] = -1. * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.WFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI ) ]) {
      /* left boundary */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.WFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.WCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.WCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.WCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      if (!nbrfaces[4]) {
        // back left corner, no slip or outflows
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[1] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        if (direction == 1 && (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellWidthsLDI ) ] \
                                    + Mesh.WCellWidths[ cl ]) > ymax) { // outflow to back
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        }
        else { // no slip in both directions
          val[4] = -(val[0] + val[1] + val[2] + val[3]) \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                            * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                            / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        }
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[5] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else if (!nbrfaces[5]) {
        /* front left corner */
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dx[1] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -(val[0] + val[1] + val[2] + val[3]) \
               + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
               + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 4; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[4]);
      }
      else {
        /* left boundary no corner */
        val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[1] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                 * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                 / (0.5 * (dx[1] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
        val[2] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
        val[3] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[4] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[4] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                  * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                  / (0.5 * (dy[5] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        // Columns
        colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[3] = nbrfaces[4] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
        for (fc = 0; fc < 5; fc++) {
          matIs.push_back(cl2);
          matJs.push_back(colId[fc]);
          matVals.push_back(val[fc]);
        }
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(val[5]);
      }
      val[0] = -1. * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.WFaceConnectivity[ idx2( cl, 5, Mesh.FaceConnectivityLDI ) ]) {
      /* front boundary */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.WFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.WCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.WCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.WCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[1] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[3] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[4] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
             + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                        * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                        / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      // Columns
      colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[3] = nbrfaces[3] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[4] = nbrfaces[4] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 5; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[5]);
      // Pressure components
      val[0] = -1. * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.WFaceConnectivity[ idx2( cl, 4, Mesh.FaceConnectivityLDI ) ]) {
      /* back boundary, no slip or outflow */
      for (fc = 0; fc < 6; fc++) {
        nbrfaces[fc] = Mesh.WFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.WCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.WCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dz[fc] = Mesh.WCellWidths[ idx2( 2, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
          dz[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[0] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[1] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dz[2] + Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -visc * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dx[3] + Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[4] = -visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                / (0.5 * (dy[5] + Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      if (direction == 1 && (Mesh.WCellCenters[ idx2( cl, 1, Mesh.CellWidthsLDI ) ] \
                                    + Mesh.WCellWidths[ cl ]) > ymax) { // outflow
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]);
      }
      else { // no slip
        val[5] = -(val[0] + val[1] + val[2] + val[3] + val[4]) \
               + 2 * visc * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          * Mesh.WCellWidths[ idx2( 2, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      }
      // Columns
      colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[3] = nbrfaces[3] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[4] = nbrfaces[5] - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 5; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[5]);
      // Pressure components
      val[0] = -1. * Mesh.WCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                   * Mesh.WCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.WCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      colId[1] = Mesh.WCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
  }
}

/* ucartflow2D implements boundary conditions in rows of the linear system
   associated to u components of boundary voxels for the 2D problem.
   It implements a flow problem in a principal axis direction. */
void
ucartflow2D ( const FluidMesh& Mesh, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals, \
              std::vector<double>& force, \
              double visc, int direction )
{

  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  int cl, fc, colId[4], nbrfaces[4];
  double val[5], dx[4], dy[4];

 // u boundary conditions
  for (unsigned long arrayIndex = 0; arrayIndex < Mesh.UBoundaryCells.size(); arrayIndex++) {
    cl = Mesh.UBoundaryCells[arrayIndex];
    if (!Mesh.UFaceConnectivity[ idx2( cl, 1, Mesh.FaceConnectivityLDI) ]) {
      // right boundary
      if (direction == 0 && Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] == xmax) { // xflow problem, this is an outflow boundary
        val[0] = -1. / Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
        colId[0] = Mesh.UFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI ) ] - 1;
        colId[1] = cl;
        // Insert values
        matIs.push_back(cl);
        matJs.push_back(colId[0]);
        matVals.push_back(val[0]);
        matIs.push_back(cl);
        matJs.push_back(colId[1]);
        matVals.push_back(val[1]);
      }
      else { // yflow problem, this is a no slip boundary
        matIs.push_back(cl);
        matJs.push_back(cl);
        matVals.push_back(1);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 3, Mesh.FaceConnectivityLDI) ]) {
      // left boundary, Dirichlet boundary
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(1);
      if (direction == 0 && Mesh.UCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] == xmin) { // xflow problem, force is nonzero here
        force[cl] = -maxin * (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymin) \
                                   * (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] - ymax);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ]) {
      // bottom boundary, 0 Dirichlet value for u in either case
      for (fc = 0; fc < 4; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[2] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -(val[0] + val[1] + val[2]) \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      colId[0] = nbrfaces[1] - 1;
      colId[1] = nbrfaces[2] - 1;
      colId[2] = nbrfaces[3] - 1;
      for (fc = 0; fc < 3; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(val[3]);
      // Pressure component to BC
      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.UFaceConnectivity[ idx2( cl, 2, Mesh.FaceConnectivityLDI ) ]) {
      // top boundary
      for (fc = 0; fc < 4; fc++) {
        nbrfaces[fc] = Mesh.UFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.UCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.UCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[0] + Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[1] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[3] + Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      if (direction == 1 && (Mesh.UCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] \
                               + Mesh.UCellWidths[ cl ]) > ymax ) { // yflow problem, this is an outflow
        val[3] = -(val[0] + val[1] + val[2]);
      }
      else { // xflow problem, this is a no slip boundary
        val[3] = -(val[0] + val[1] + val[2]) \
               + 2 * visc * Mesh.UCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                          / Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      }
      colId[0] = nbrfaces[0] - 1;
      colId[1] = nbrfaces[1] - 1;
      colId[2] = nbrfaces[3] - 1;
      for (fc = 0; fc < 3; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl);
      matJs.push_back(cl);
      matVals.push_back(val[3]);
      // Pressure component to BC
      val[0] = -1. * Mesh.UCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.UCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = Mesh.UCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
  }
}

/* vcartflow2D implements boundary conditions in rows of the linear system
   associated to v components of boundary voxels for the 2D problem.
   It implements a flow problem in a principal axis direction. */
void
vcartflow2D ( const FluidMesh& Mesh, std::vector<int>& matIs, \
              std::vector<int>& matJs, std::vector<double>& matVals, \
              std::vector<double>& force, \
              double visc, int direction )
{

  double maxin = 1.0;
  double xmin = Mesh.xLim[0];
  double xmax = Mesh.xLim[1];
  double ymin = Mesh.yLim[0];
  double ymax = Mesh.yLim[1];
  int cl, cl2, fc, colId[4], nbrfaces[4];
  double val[4], dx[4], dy[4];

  for (unsigned long arrayIndex = 0; arrayIndex < Mesh.VBoundaryCells.size(); arrayIndex++) {
    cl = Mesh.VBoundaryCells[arrayIndex];
    cl2 = cl + Mesh.DOF[1];
    if (!Mesh.VFaceConnectivity[ idx2(cl, 0, Mesh.FaceConnectivityLDI) ]) {
      // bottom boundary, Dirichlet in either case
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(1);
      if (direction == 1 && Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] == ymin) { // yflow problem, force is nonzero here
        force[cl2] = -maxin * (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmin) \
                                 * (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] - xmax);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2(cl, 2, Mesh.FaceConnectivityLDI) ]) {
      // top boundary
      if (direction == 1 && Mesh.VCellCenters[ idx2( cl, 1, Mesh.CellCentersLDI ) ] == ymax) { // outflow boundary for yflow problem
        val[0] = -1. / Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        val[1] = 1. / Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ];
        colId[0] = Mesh.VFaceConnectivity[ idx2( cl, 0, Mesh.FaceConnectivityLDI ) ] - 1 + Mesh.DOF[1];
        colId[1] = cl2;
        // Insert values
        matIs.push_back(cl2);
        matJs.push_back(colId[0]);
        matVals.push_back(val[0]);
        matIs.push_back(cl2);
        matJs.push_back(colId[1]);
        matVals.push_back(val[1]);
      }
       else { // no slip
        matIs.push_back(cl2);
        matJs.push_back(cl2);
        matVals.push_back(1);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2(cl, 1, Mesh.FaceConnectivityLDI) ]) {
      // right boundary
      for (fc = 0; fc < 4; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[0] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[2] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[3] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      if (direction == 0 && (Mesh.VCellCenters[ idx2( cl, 0, Mesh.CellCentersLDI ) ] \
                               + Mesh.VCellWidths[ cl ]) > xmax) {
        val[3] = -(val[0] + val[1] + val[2]);
      }
      else {
        val[3] = -(val[0] + val[1] + val[2]) \
                  + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                             / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      }
      colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
      colId[1] = nbrfaces[2] - 1 + Mesh.DOF[1];
      colId[2] = nbrfaces[3] - 1 + Mesh.DOF[1];
      for (fc = 0; fc < 3; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[3]);
      // Pressure component
      val[0] = -1. * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
    else if (!Mesh.VFaceConnectivity[ idx2(cl, 3, Mesh.FaceConnectivityLDI) ]) {
      // left boundary
      for (fc = 0; fc < 4; fc++) {
        nbrfaces[fc] = Mesh.VFaceConnectivity[ idx2( cl, fc, Mesh.FaceConnectivityLDI ) ];
        if (nbrfaces[fc]) {
          dx[fc] = Mesh.VCellWidths[ idx2( 0, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
          dy[fc] = Mesh.VCellWidths[ idx2( 1, (nbrfaces[fc] - 1), Mesh.CellWidthsLDI ) ];
        }
        else {
          dx[fc] = 0;
          dy[fc] = 0;
        }
      }
      val[0] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[0] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[1] = -visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dx[1] + Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ]));
      val[2] = -visc * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ] \
                     / (0.5 * (dy[2] + Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ]));
      val[3] = -(val[0] + val[1] + val[2]) \
                + 2 * visc * Mesh.VCellWidths[ idx2( 1, cl, Mesh.CellWidthsLDI ) ] \
                           / Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      colId[0] = nbrfaces[0] - 1 + Mesh.DOF[1];
      colId[1] = nbrfaces[1] - 1 + Mesh.DOF[1];
      colId[2] = nbrfaces[2] - 1 + Mesh.DOF[1];
      for (fc = 0; fc < 3; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
      matIs.push_back(cl2);
      matJs.push_back(cl2);
      matVals.push_back(val[3]);
      // Pressure component
      val[0] = -1. * Mesh.VCellWidths[ idx2( 0, cl, Mesh.CellWidthsLDI ) ];
      val[1] = -val[0];
      colId[0] = Mesh.VCellPressureNeighbor[ idx2( cl, 0, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      colId[1] = Mesh.VCellPressureNeighbor[ idx2( cl, 1, Mesh.VelocityCellPressureNeighborLDI ) ] \
                 - 1 + Mesh.DOF[1] + Mesh.DOF[2];
      for (fc = 0; fc < 2; fc++) {
        matIs.push_back(cl2);
        matJs.push_back(colId[fc]);
        matVals.push_back(val[fc]);
      }
    }
  }
}

void
BoundaryConditions ( const FluidMesh& Mesh, double visc, \
                     std::vector<int>& matIs, std::vector<int>& matJs, \
                     std::vector<double>& matVals, \
                     std::vector<double>& force, int direction )
{

  switch ( Mesh.DIM )
  {
    case 3 :
      ucartflow( Mesh, matIs, matJs, matVals, force, visc, direction );
      vcartflow( Mesh, matIs, matJs, matVals, force, visc, direction );
      wcartflow( Mesh, matIs, matJs, matVals, force, visc, direction );
      break;
    case 2 :
      ucartflow2D( Mesh, matIs, matJs, matVals, force, visc, direction );
      vcartflow2D( Mesh, matIs, matJs, matVals, force, visc, direction );
      break;
  }

}

