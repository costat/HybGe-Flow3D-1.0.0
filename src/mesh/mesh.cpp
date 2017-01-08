/* mesh */
// system includes
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cstring>
#include <algorithm>

// nlem includes
#include "../../include/hgflow.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// 1d->3d index
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))

//-------------- quadrilateral elements section ----------------//

void
hgf::mesh::build( parameters& par)
{
  if (par.dimension == 2) { // 2d
    build_from_voxel_quad(par);
  }
  else { // 3d
    build_from_voxel_hex(par);
  }

#ifdef _MESH_VISUAL_DEBUG
  printVTK(par);
#endif
}

void
hgf::mesh::build_from_voxel_quad( parameters& par)
{
  double dx = (double)par.length / par.nx;
  double dy = (double)par.width / par.ny;
  int nCells = 0;
  int nNodes = 0;
  int nEdges = 0;
  std::vector< int > x_index;
  std::vector< int > y_index;
  x_index.resize(par.voxel_geometry.size());
  y_index.resize(par.voxel_geometry.size());
  std::vector< int > cell_numbers(par.voxel_geometry.data(), \
    par.voxel_geometry.data() + par.voxel_geometry.size());

#pragma omp parallel
  {
    // determine the number of elements, edges, nodes
#pragma omp for reduction(+:nCells,nNodes,nEdges)
    for (int yi = 0; yi < par.ny; yi++) {
      for (int xi = 0; xi < par.nx; xi++) {
        if (par.voxel_geometry[idx2(yi, xi, par.nx)] != 1) {
          x_index[idx2(yi, xi, par.nx)] = xi;
          y_index[idx2(yi, xi, par.nx)] = yi;
          nCells++;
          if (yi == 0) { // first y layer
            if (xi == 0) { // first cell
              nNodes += 4;
              nEdges += 4;
            }
            else if (par.voxel_geometry[idx2(yi, (xi-1), par.nx)] == 1) {
              nNodes += 4;
              nEdges += 4;
            }
            else {
              nNodes += 2;
              nEdges += 3;
            }
          }
          else {
            if (xi == 0) { // first cell
              if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) {
                if (par.voxel_geometry[idx2((yi - 1), (xi + 1), par.nx)] == 1) {
                  // all four nodes and edges are new
                  nNodes += 4;
                  nEdges += 4;
                }
                else {
                  // 3 of the nodes are new, 4 of the edges
                  nNodes += 3;
                  nEdges += 4;
                }
              }
              else {
                // there is a cell below, 2 new nodes, 3 new edges
                nNodes += 2;
                nEdges += 3;
              }
            }
            else if (xi < par.nx - 1) { // not the last cell in the y-level
              if (par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) { // there's no cell to the left
                if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) { // there's no cell below
                  if (par.voxel_geometry[idx2((yi - 1), (xi - 1), par.nx)] == 1) { // there's no cell bott left
                    if (par.voxel_geometry[idx2((yi - 1), (xi + 1), par.nx)] == 1) { // there's no cell bottom right
                      // all four nodes and edges are new
                      nNodes += 4;
                      nEdges += 4;
                    }
                    else {
                      // 3 of the nodes are new, 4 of the edges
                      nNodes += 3;
                      nEdges += 4;
                    }
                  }
                  else { // there's a cell bott left
                    if (par.voxel_geometry[idx2((yi - 1), (xi + 1), par.nx)] == 1) { // no cell bott right
                      nNodes += 3;
                      nEdges += 4;
                    }
                    else {
                      nNodes += 2;
                      nEdges += 4;
                    }
                  }
                }
                else { // there's a cell below
                  nNodes += 2;
                  nEdges += 3;
                }
              }
              else { // there is a cell to the left
                if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) { // there's no cell below
                  if (par.voxel_geometry[idx2((yi - 1), (xi + 1), par.nx)] == 1) { // there's no cell bott right
                    nNodes += 2;
                    nEdges += 3;
                  }
                  else {
                    nNodes += 1;
                    nEdges += 3;
                  }
                }
                else { // there is a cell below
                  nNodes += 1;
                  nEdges += 2;
                }
              }
            }
            else { // last cell in the y-level
              if (par.voxel_geometry[idx2(yi, (xi - 1), par.nx)] == 1) { // no cell to the left
                if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) { // no cell below
                  if (par.voxel_geometry[idx2((yi-1),(xi-1),par.nx)] == 1) { // no cell bott left
                    nNodes += 4;
                    nEdges += 4;
                  }
                  else {
                    nNodes += 3;
                    nEdges += 4;
                  }
                }
                else {
                  nNodes += 2;
                  nEdges += 3;
                }
              }
              else { // there is a cell to the left
                if (par.voxel_geometry[idx2((yi - 1), xi, par.nx)] == 1) { // no cell below
                  nNodes += 2;
                  nEdges += 3;
                }
                else {
                  nNodes += 1;
                  nEdges += 2;
                }
              }
            }
          }
        }
        else {
          x_index[idx2(yi, xi, par.nx)] = -1;
          y_index[idx2(yi, xi, par.nx)] = -1;
        }
      }
    }

#pragma omp single
    { // various tasks that need to be performed serially
      gtlNode.resize(nNodes * 4);
      els.resize(nCells);
      std::cout << "\nBuilding " << nCells << " quadrilateral element mesh from voxel data.\n\n";
      int cell_num = -1;
      for (int cell = 0; cell < cell_numbers.size(); cell++) {
        if (cell_numbers[cell] != 1) {
          cell_num++;
          cell_numbers[cell] = cell_num;
        }
        else {
          cell_numbers[cell] = -1;
        }
      }
      int ccount = 0;
      x_index.erase(std::remove(x_index.begin(), x_index.end(), -1), x_index.end());
      y_index.erase(std::remove(y_index.begin(), y_index.end(), -1), y_index.end());
    }

    // x direction loop on elements
#pragma omp for schedule(static,1) nowait
    for (int cell = 0; cell < nCells; cell++) {
      // grab ii, jj for coordinates
      int ii = x_index[cell];
      int jj = y_index[cell];

      els[cell].dx = dx;
      els[cell].dy = dy;

      //----- node information -----//
      // Note: gnums are temporary here, and are unique'd later
      // local node 0, bottom left on element
      els[cell].vtx[0].coords[0] = ii*dx;
      els[cell].vtx[0].coords[1] = jj*dy;

      // local node 1, bottom right on element
      els[cell].vtx[1].coords[0] = (ii + 1)*dx;
      els[cell].vtx[1].coords[1] = jj*dy;

      // local node 2, top right on element
      els[cell].vtx[2].coords[0] = (ii + 1)*dx;
      els[cell].vtx[2].coords[1] = (jj + 1)*dy;

      // local node 3, top left on element
      els[cell].vtx[3].coords[0] = ii*dx;
      els[cell].vtx[3].coords[1] = (jj + 1)*dy;

      //------ edge information ------//
      // local edge 0, bottom on element
      if (jj == 0 || par.voxel_geometry[idx2((jj - 1), ii, par.nx)] == 1) {
        // this is a boundary edge
        els[cell].edg[0].bctype = 1;
        els[cell].edg[0].neighbor = -1;
      }
      else {
        els[cell].edg[0].bctype = 0;
        els[cell].edg[0].neighbor = cell_numbers[idx2((jj - 1), ii, par.nx)];
      }

      // local edge 1, right on element
      if (ii == par.nx-1 || par.voxel_geometry[idx2(jj, (ii + 1), par.nx)] == 1) {
        // this is a boundary edge
        els[cell].edg[1].bctype = 1;
        els[cell].edg[1].neighbor = -1;
      }
      else {
        els[cell].edg[1].bctype = 0;
        els[cell].edg[1].neighbor = cell_numbers[idx2(jj, (ii + 1), par.nx)];
      }

      // local edge 2, top on element
      if (jj == par.ny-1 || par.voxel_geometry[idx2((jj + 1), ii, par.nx)] == 1) {
        // this is a boundary edge
        els[cell].edg[2].bctype = 1;
        els[cell].edg[2].neighbor = -1;
      }
      else {
        els[cell].edg[2].bctype = 0;
        els[cell].edg[2].neighbor = cell_numbers[idx2((jj + 1), ii, par.nx)];
      }

// local edge 3, left on element
if (ii == 0 || par.voxel_geometry[idx2(jj, (ii - 1), par.nx)] == 1) {
  // this is a boundary edge
  els[cell].edg[3].bctype = 1;
  els[cell].edg[3].neighbor = -1;
}
else {
  els[cell].edg[3].bctype = 0;
  els[cell].edg[3].neighbor = cell_numbers[idx2(jj, (ii - 1), par.nx)];
}

    }

#pragma omp sections
    {
#pragma omp section
      {
        int edge_number = -1;
        int neighbor_number;
        int cell_number;
        // unique the edge numbers
        // first pass for x components
        for (int xi = 0; xi < par.nx; xi++) {
          for (int yi = 0; yi < par.ny; yi++) {
            cell_number = cell_numbers[idx2(yi, xi, par.nx)];
            if (cell_number != -1) {
              if (els[cell_number].edg[0].neighbor == -1) {
                edge_number++;
                els[cell_number].edg[0].gnum = edge_number;
              }
              else {
                neighbor_number = els[cell_number].edg[0].neighbor;
                els[cell_number].edg[0].gnum = els[neighbor_number].edg[2].gnum;
              }
              edge_number++;
              els[cell_number].edg[2].gnum = edge_number;
            }
          }
        }
        // second pass for y components
        for (int yi = 0; yi < par.ny; yi++) {
          for (int xi = 0; xi < par.nx; xi++) {
            cell_number = cell_numbers[idx2(yi, xi, par.nx)];
            if (cell_number != -1) {
              if (els[cell_number].edg[3].neighbor == -1) {
                edge_number++;
                els[cell_number].edg[3].gnum = edge_number;
              }
              else {
                neighbor_number = els[cell_number].edg[3].neighbor;
                els[cell_number].edg[3].gnum = els[neighbor_number].edg[1].gnum;
              }
              edge_number++;
              els[cell_number].edg[1].gnum = edge_number;
            }
          }
        }
      }
#pragma omp section
      {
        // unique the node numbers
        int node_number = -1;
        int neighbor_number;
        int cell_number;
        for (int yi = 0; yi < par.ny; yi++) {
          for (int xi = 0; xi < par.nx; xi++) {
            if (par.voxel_geometry[idx2(yi, xi, par.nx)] != 1) { // live cell
              if (xi == 0) {
                // 2 new nodes to add in this cell and appropriate neighbors
                cell_number = cell_numbers[idx2(yi, xi, par.nx)];
                node_number += 2;
                els[cell_number].vtx[0].gnum = node_number - 1;
                els[cell_number].vtx[1].gnum = node_number;
                // bottom neighbors (if applicable)
                if (yi) {
                  if (cell_numbers[idx2((yi - 1), xi, par.nx)] != -1) {
                    // there is a cell below
                    neighbor_number = cell_numbers[idx2((yi - 1), xi, par.nx)];
                    els[neighbor_number].vtx[3].gnum = node_number - 1;
                    els[neighbor_number].vtx[2].gnum = node_number;
                  }
                  if (xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                    // there is a cell bott right
                    neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                    els[neighbor_number].vtx[3].gnum = node_number;
                  }
                }
                // right neighbor (if applicable)
                if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                  neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                  els[neighbor_number].vtx[0].gnum = node_number;
                }
              }
              else { // live, xi != 0
                // condition under which we have 2 nodes to add
                if (!yi && cell_numbers[idx2(yi, (xi - 1), par.nx)] == -1) { // no cells below or left
                  cell_number = cell_numbers[idx2(yi, xi, par.nx)];
                  node_number += 2;
                  els[cell_number].vtx[0].gnum = node_number - 1;
                  els[cell_number].vtx[1].gnum = node_number;
                  // cell to right if applicable
                  if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                    els[neighbor_number].vtx[0].gnum = node_number;
                  }
                }
                else if (cell_numbers[idx2(yi, (xi - 1), par.nx)] == -1 && yi \
                    && cell_numbers[idx2((yi - 1), (xi - 1), par.nx)] == -1 \
                    && cell_numbers[idx2((yi - 1), xi, par.nx)] == -1) {

                  cell_number = cell_numbers[idx2(yi, xi, par.nx)];
                  node_number += 2;
                  els[cell_number].vtx[0].gnum = node_number - 1;
                  els[cell_number].vtx[1].gnum = node_number;
                  // bott right if applicable
                  if (xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                    els[neighbor_number].vtx[3].gnum = node_number;
                  }
                  // cell to right if applicable
                  if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                    els[neighbor_number].vtx[0].gnum = node_number;
                  }
                }
                // else we have 1 node to add
                else {
                  node_number++;
                  cell_number = cell_numbers[idx2(yi, xi, par.nx)];
                  els[cell_number].vtx[1].gnum = node_number;
                  if (yi) {
                    // bottom if applicable
                    if (cell_numbers[idx2((yi - 1), xi, par.nx)] != -1) {
                      neighbor_number = cell_numbers[idx2((yi - 1), xi, par.nx)];
                      els[neighbor_number].vtx[2].gnum = node_number;
                    }
                    // bottom right if applicable
                    if (xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                      neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                      els[neighbor_number].vtx[3].gnum = node_number;
                    }
                  }
                  // right if applicable
                  if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                    els[neighbor_number].vtx[0].gnum = node_number;
                  }
                }
              }
            }
            else {
              if (xi == 0) {
                // dead first cell
                if (yi && cell_numbers[idx2((yi - 1), xi, par.nx)] != -1) {
                  // there's a cell below, need 2 new nodes
                  node_number += 2;
                  neighbor_number = cell_numbers[idx2((yi - 1), xi, par.nx)];
                  els[neighbor_number].vtx[3].gnum = node_number - 1;
                  els[neighbor_number].vtx[2].gnum = node_number;
                  // cell bott right
                  if (xi < par.nx && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                    els[neighbor_number].vtx[3].gnum = node_number;
                  }
                  // cell right
                  if (xi < par.nx && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                    els[neighbor_number].vtx[0].gnum = node_number;
                  }
                }
                else if (yi && xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                  // there's a cell bott right, need 1 new node
                  node_number++;
                  neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                  els[neighbor_number].vtx[3].gnum = node_number;
                  // cell right
                  if (cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                    neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                    els[neighbor_number].vtx[0].gnum = node_number;
                  }
                }
                else if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                  // there's a cell to the right, need 1 new node
                  node_number++;
                  neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                  els[neighbor_number].vtx[0].gnum = node_number;
                }
              }
              else {
                // dead, not first cell
                if (yi) {
                  if (cell_numbers[idx2((yi - 1), xi, par.nx)] != -1) {
                    // cell below, 1 new node
                    node_number++;
                    neighbor_number = cell_numbers[idx2((yi - 1), xi, par.nx)];
                    els[neighbor_number].vtx[2].gnum = node_number;
                    if (xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                      // cell bott right
                      neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                      els[neighbor_number].vtx[3].gnum = node_number;
                    }
                    if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                      // cell right
                      neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                      els[neighbor_number].vtx[0].gnum = node_number;
                    }
                  }
                  else if (xi < par.nx - 1 && cell_numbers[idx2((yi - 1), (xi + 1), par.nx)] != -1) {
                    // cell bott right, 1 new node
                    node_number++;
                    neighbor_number = cell_numbers[idx2((yi - 1), (xi + 1), par.nx)];
                    els[neighbor_number].vtx[3].gnum = node_number;
                    if (cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                      // cell right
                      neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                      els[neighbor_number].vtx[0].gnum = node_number;
                    }
                  }
                }
                else if (xi < par.nx - 1 && cell_numbers[idx2(yi, (xi + 1), par.nx)] != -1) {
                  // cell to right, +1 node
                  node_number++;
                  neighbor_number = cell_numbers[idx2(yi, (xi + 1), par.nx)];
                  els[neighbor_number].vtx[0].gnum = node_number;
                }
              }
            }
          }
        }
        // last sweep for top y-level roofs
        for (int xi = 0; xi < par.nx; xi++) {
          if (cell_numbers[idx2((par.ny - 1), xi, par.nx)] != -1) { // alive top cell
            if (xi == 0 || cell_numbers[idx2((par.ny - 1), (xi - 1), par.nx)] == -1) {
              // 2 nodes to add
              node_number += 2;
              cell_number = cell_numbers[idx2((par.ny - 1), xi, par.nx)];
              els[cell_number].vtx[3].gnum = node_number - 1;
              els[cell_number].vtx[2].gnum = node_number;
              // neighbor right
              if (xi < par.nx - 1 && cell_numbers[idx2((par.ny - 1), (xi + 1), par.nx)] != -1) {
                neighbor_number = cell_numbers[idx2((par.ny - 1), (xi + 1), par.nx)];
                els[neighbor_number].vtx[3].gnum = node_number;
              }
            }
            else {
              node_number++;
              cell_number = cell_numbers[idx2((par.ny - 1), xi, par.nx)];
              els[cell_number].vtx[2].gnum = node_number;
              if (xi < par.nx - 1 && cell_numbers[idx2((par.ny - 1), (xi + 1), par.nx)] != -1) {
                neighbor_number = cell_numbers[idx2((par.ny - 1), (xi + 1), par.nx)];
                els[neighbor_number].vtx[3].gnum = node_number;
              }
            }
          }
        }
      }
    }
#pragma omp barrier
    // need unique'd nodes before next section

#pragma omp for schedule(static,1) nowait
      for (int cell = 0; cell < nCells; cell++) {
      gtlNode[idx2(els[cell].vtx[0].gnum, 3, 4)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[1].gnum, 2, 4)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[2].gnum, 0, 4)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[3].gnum, 1, 4)] = cell + 1;
    }

  }

}

void
hgf::mesh::build_from_voxel_hex(parameters& par)
{

  double dx = (double)par.length / par.nx;
  double dy = (double)par.width / par.ny;
  double dz = (double)par.height / par.nz;
  int nCells = 0;
  int nNodes = 0;
  int nFaces = 0;
  std::vector< int > x_index;
  std::vector< int > y_index;
  std::vector< int > z_index;
  x_index.resize(par.voxel_geometry.size());
  y_index.resize(par.voxel_geometry.size());
  z_index.resize(par.voxel_geometry.size());
  std::vector< int > cell_numbers(par.voxel_geometry.data(), \
    par.voxel_geometry.data() + par.voxel_geometry.size());

#pragma omp parallel
  {
    // determine the number of elements, faces, nodes -- no edges currently, add later if an application demands
#pragma omp for reduction(+:nCells,nNodes,nFaces)
    for (int zi = 0; zi < par.nz; zi++) {
      for (int yi = 0; yi < par.ny; yi++) {
        for (int xi = 0; xi < par.nx; xi++) {
          if (par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] != 1) {
            x_index[idx3(zi, yi, xi, par.ny, par.nx)] = xi;
            y_index[idx3(zi, yi, xi, par.ny, par.nx)] = yi;
            z_index[idx3(zi, yi, xi, par.ny, par.nx)] = zi;
            nCells++;
            int tmp_nFaces = 3;
            int tmp_nNodes = 1;
            // nFaces
            if (xi == 0 || par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) tmp_nFaces++;
            if (yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) tmp_nFaces++;
            if (zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) tmp_nFaces++;
            // nNodes
            // node 1
            if ((xi == 0 || par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) && \
              (xi == 0 || yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), (xi - 1), par.ny, par.nx)] == 1) && \
              (yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) && \
              (yi == 0 || zi == 0 || par.voxel_geometry[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)] == 1) && \
              (xi == 0 || yi == 0 || zi == 0 || par.voxel_geometry[idx3((zi - 1), (yi - 1), (xi - 1), par.ny, par.nx)] == 1) && \
              (xi == 0 || zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, (xi - 1), par.ny, par.nx)] == 1) && \
              (zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1)) \
              tmp_nNodes++;
            // node 2
            if ((zi == 0 || xi == par.nx - 1 || par.voxel_geometry[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)] == 1) && \
              (zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) && \
              (zi == 0 || yi == 0 || par.voxel_geometry[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)] == 1) && \
              (zi == 0 || yi == 0 || xi == par.nx - 1 || par.voxel_geometry[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)] == 1) && \
              (yi == 0 || xi == par.nx - 1 || par.voxel_geometry[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)] == 1) && \
              (yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1)) \
              tmp_nNodes++;
            // node 3
            if ((zi == 0 || xi == par.nx - 1 || par.voxel_geometry[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)] == 1) && \
              (zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1)) \
              tmp_nNodes++;
            // node 4
            if ((zi == 0 || par.voxel_geometry[idx3((zi - 1), yi, xi, par.ny, par.nx)] == 1) && \
              (zi == 0 || xi == 0 || par.voxel_geometry[idx3((zi - 1), yi, (xi - 1), par.ny, par.nx)] == 1) && \
              (xi == 0 || par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1)) \
              tmp_nNodes++;
            // node 5
            if (xi == 0 || par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) tmp_nNodes++;
            // node 6
            if ((xi == 0 || par.voxel_geometry[idx3(zi, yi, (xi - 1), par.ny, par.nx)] == 1) && \
              (xi == 0 || yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), (xi - 1), par.ny, par.nx)] == 1) && \
              (yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1)) \
              tmp_nNodes++;
            // node 7
            if ((yi == 0 || par.voxel_geometry[idx3(zi, (yi - 1), xi, par.ny, par.nx)] == 1) && \
              (yi == 0 || xi == par.nx - 1 || par.voxel_geometry[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)] == 1)) \
              tmp_nNodes++;

            nFaces += tmp_nFaces;
            nNodes += tmp_nNodes;
          }
          else {
            x_index[idx3(zi, yi, xi, par.ny, par.nx)] = -1;
            y_index[idx3(zi, yi, xi, par.ny, par.nx)] = -1;
            z_index[idx3(zi, yi, xi, par.ny, par.nx)] = -1;
          }
        }
      }
    }

#pragma omp barrier

#pragma omp single
    { // various tasks that need to be performed serially
      gtlNode.resize(nNodes * 8);
      els.resize(nCells);
      std::cout << "\nBuilding " << nCells << " hexahedral cell mesh from voxel data.\n\n";
      int cell_num = -1;
      for (int cell = 0; cell < cell_numbers.size(); cell++) {
        if (cell_numbers[cell] != 1) {
          cell_num++;
          cell_numbers[cell] = cell_num;
        }
        else {
          cell_numbers[cell] = -1;
        }
      }
      int ccount = 0;
      x_index.erase(std::remove(x_index.begin(), x_index.end(), -1), x_index.end());
      y_index.erase(std::remove(y_index.begin(), y_index.end(), -1), y_index.end());
      z_index.erase(std::remove(z_index.begin(), z_index.end(), -1), z_index.end());
    } // end of omp single section

#pragma omp barrier

#pragma omp for schedule(static,1) nowait
    for (int cell = 0; cell < nCells; cell++) {
      // grab ii, jj, kk for coordinates
      int ii = x_index[cell];
      int jj = y_index[cell];
      int kk = z_index[cell];

      els[cell].dx = dx;
      els[cell].dy = dy;
      els[cell].dz = dz;

      //----- node information -----//
      // local node 0
      els[cell].vtx[0].coords[0] = ii*dx;
      els[cell].vtx[0].coords[1] = jj*dy;
      els[cell].vtx[0].coords[2] = kk*dz;

      // local node 1
      els[cell].vtx[1].coords[0] = (ii + 1)*dx;
      els[cell].vtx[1].coords[1] = jj*dy;
      els[cell].vtx[1].coords[2] = kk*dz;

      // local node 2
      els[cell].vtx[2].coords[0] = (ii + 1)*dx;
      els[cell].vtx[2].coords[1] = (jj + 1)*dy;
      els[cell].vtx[2].coords[2] = kk*dz;

      // local node 3
      els[cell].vtx[3].coords[0] = ii*dx;
      els[cell].vtx[3].coords[1] = (jj + 1)*dy;
      els[cell].vtx[3].coords[2] = kk*dz;

      // local node 4
      els[cell].vtx[4].coords[0] = ii*dx;
      els[cell].vtx[4].coords[1] = (jj + 1)*dy;
      els[cell].vtx[4].coords[2] = (kk + 1)*dz;

      // local node 5
      els[cell].vtx[5].coords[0] = (ii + 1)*dx;
      els[cell].vtx[5].coords[1] = (jj + 1)*dy;
      els[cell].vtx[5].coords[2] = (kk + 1)*dz;

      // local node 6
      els[cell].vtx[6].coords[0] = (ii + 1)*dx;
      els[cell].vtx[6].coords[1] = jj*dy;
      els[cell].vtx[6].coords[2] = (kk + 1)*dz;

      // local node 7
      els[cell].vtx[7].coords[0] = ii*dx;
      els[cell].vtx[7].coords[1] = jj*dy;
      els[cell].vtx[7].coords[2] = (kk + 1)*dz;

      //------ face information ------//
      // local face 0, negative y direction
      if (jj == 0 || par.voxel_geometry[idx3(kk, (jj - 1), ii, par.ny, par.nx)] == 1) {
        els[cell].fac[0].bctype = 1;
        els[cell].fac[0].neighbor = -1;
      }
      else {
        els[cell].fac[0].bctype = 0;
        els[cell].fac[0].neighbor = cell_numbers[idx3(kk, (jj - 1), ii, par.ny, par.nx)];
      }
      // local face 1, positive x direction
      if (ii == par.nx - 1 || par.voxel_geometry[idx3(kk, jj, (ii + 1), par.ny, par.nx)] == 1) {
        els[cell].fac[1].bctype = 1;
        els[cell].fac[1].neighbor = -1;
      }
      else {
        els[cell].fac[1].bctype = 0;
        els[cell].fac[1].neighbor = cell_numbers[idx3(kk, jj, (ii + 1), par.ny, par.nx)];
      }
      // local face 2, positive y direction
      if (jj == par.ny - 1 || par.voxel_geometry[idx3(kk, (jj + 1), ii, par.ny, par.nx)] == 1) {
        els[cell].fac[2].bctype = 1;
        els[cell].fac[2].neighbor = -1;
      }
      else {
        els[cell].fac[2].bctype = 0;
        els[cell].fac[2].neighbor = cell_numbers[idx3(kk, (jj + 1), ii, par.ny, par.nx)];
      }
      // local face 3, negative x direction
      if (ii == 0 || par.voxel_geometry[idx3(kk, jj, (ii - 1), par.ny, par.nx)] == 1) {
        els[cell].fac[3].bctype = 1;
        els[cell].fac[3].neighbor = -1;
      }
      else {
        els[cell].fac[3].bctype = 0;
        els[cell].fac[3].neighbor = cell_numbers[idx3(kk, jj, (ii - 1), par.ny, par.nx)];
      }
      // local face 4, negative z direction
      if (kk == 0 || par.voxel_geometry[idx3((kk - 1), jj, ii, par.ny, par.nx)] == 1) {
        els[cell].fac[4].bctype = 1;
        els[cell].fac[4].neighbor = -1;
      }
      else {
        els[cell].fac[4].bctype = 0;
        els[cell].fac[4].neighbor = cell_numbers[idx3((kk - 1), jj, ii, par.ny, par.nx)];
      }
      // local face 5, positive z direction
      if (kk == par.nz - 1 || par.voxel_geometry[idx3((kk + 1), jj, ii, par.ny, par.nx)] == 1) {
        els[cell].fac[5].bctype = 1;
        els[cell].fac[5].neighbor = -1;
      }
      else {
        els[cell].fac[5].bctype = 0;
        els[cell].fac[5].neighbor = cell_numbers[idx3((kk + 1), jj, ii, par.ny, par.nx)];
      }
    }

#pragma omp sections
    {
#pragma omp section
      { // unique the faces
        int face_number = -1;
        int neighbor_number;
        int cell_number;
        // unique the face numbers
        // first pass for x faces
        for (int zi = 0; zi < par.nz; zi++) {
          for (int yi = 0; yi < par.ny; yi++) {
            for (int xi = 0; xi < par.nx; xi++) {
              cell_number = cell_numbers[idx3(zi, yi, xi, par.ny, par.nx)];
              if (cell_number != -1) {
                if (els[cell_number].fac[3].neighbor == -1) {
                  face_number++;
                  els[cell_number].fac[3].gnum = face_number;
                }
                else {
                  neighbor_number = els[cell_number].fac[3].neighbor;
                  els[cell_number].fac[3].gnum = els[neighbor_number].fac[1].gnum;
                }
                face_number++;
                els[cell_number].fac[1].gnum = face_number;
              }
            }
          }
        }
        // second pass for y faces
        for (int zi = 0; zi < par.nz; zi++) {
          for (int xi = 0; xi < par.nx; xi++) {
            for (int yi = 0; yi < par.ny; yi++) {
              cell_number = cell_numbers[idx3(zi, yi, xi, par.ny, par.nx)];
              if (cell_number != -1) {
                if (els[cell_number].fac[0].neighbor == -1) {
                  face_number++;
                  els[cell_number].fac[0].gnum = face_number;
                }
                else {
                  neighbor_number = els[cell_number].fac[0].neighbor;
                  els[cell_number].fac[0].gnum = els[neighbor_number].fac[2].gnum;
                }
                face_number++;
                els[cell_number].fac[2].gnum = face_number;
              }
            }
          }
        }
        // third pass for z faces
        for (int yi = 0; yi < par.ny; yi++) {
          for (int xi = 0; xi < par.nx; xi++) {
            for (int zi = 0; zi < par.nz; zi++) {
              cell_number = cell_numbers[idx3(zi, yi, xi, par.ny, par.nx)];
              if (cell_number != -1) {
                if (els[cell_number].fac[4].neighbor == -1) {
                  face_number++;
                  els[cell_number].fac[4].gnum = face_number;
                }
                else {
                  neighbor_number = els[cell_number].fac[4].neighbor;
                  els[cell_number].fac[4].gnum = els[neighbor_number].fac[5].gnum;
                }
                face_number++;
                els[cell_number].fac[5].gnum = face_number;
              }
            }
          }
        }
      }
#pragma omp section
      { // unique the nodes
        int node_number = -1;
        int neighbor_number;
        int cell_number;
        for (int zi = 0; zi < par.nz; zi++) {
          for (int yi = 0; yi < par.ny; yi++) {
            for (int xi = 0; xi < par.nx; xi++) {
              if (par.voxel_geometry[idx3(zi, yi, xi, par.ny, par.nx)] != 1) { // live cell
                if (xi == 0) {
                  // 2 new nodes to add in this cell and appropriate neighbors
                  cell_number = cell_numbers[idx3(zi, yi, xi, par.ny, par.nx)];
                  node_number += 2;
                  els[cell_number].vtx[0].gnum = node_number - 1;
                  els[cell_number].vtx[1].gnum = node_number;
                  goto two_new_nodes;
                } // xi == 0
                else { // live cell, xi != 0
                  // 1 node to add
                  cell_number = cell_numbers[idx3(zi, yi, xi, par.ny, par.nx)];
                  node_number++;
                  els[cell_number].vtx[1].gnum = node_number;
                  goto one_new_node;
                }
              }
              else { // dead cell
                // does this dead cell contribute any nodes? check for neighbors and jump if any of them exist
                if (xi < par.nx - 1 && cell_numbers[idx3(zi, yi, (xi + 1), par.ny, par.nx)] != -1) goto dead_active;
                else if (xi < par.nx - 1 && zi && cell_numbers[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)] != -1) goto dead_active;
                else if (zi && cell_numbers[idx3((zi - 1), yi, xi, par.ny, par.nx)] != -1) goto dead_active;
                else if (yi && xi < par.nx - 1 && cell_numbers[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)] != -1) goto dead_active;
                else if (yi && zi && xi < par.nx - 1 && cell_numbers[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) goto dead_active;
                else if (yi && zi && cell_numbers[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)] != -1) goto dead_active;
                else if (yi && cell_numbers[idx3(zi, (yi - 1), xi, par.ny, par.nx)] != -1) goto dead_active;
                else continue;
              }
            dead_active:
              {
                // 2 new nodes cases
                if (xi == 0) {
                  node_number += 2;
                  goto two_new_nodes;
                }
                // 1 new node cases
                else {
                  node_number++;
                  goto one_new_node;
                }
              }
            two_new_nodes:
              { // assumes node count was already incremented by 2 and adjusts neighboring cells
                if (yi && cell_numbers[idx3(zi, (yi - 1), xi, par.ny, par.nx)] != -1) {
                  // there is a cell backwards in y
                  neighbor_number = cell_numbers[idx3(zi, (yi - 1), xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[3].gnum = node_number - 1;
                  els[neighbor_number].vtx[2].gnum = node_number;
                }
                if (yi && xi < par.nx - 1 && cell_numbers[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                  // there is a cell back in y, forward in x
                  neighbor_number = cell_numbers[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[3].gnum = node_number;
                }
                if (xi < par.nx - 1 && cell_numbers[idx3(zi, yi, (xi + 1), par.ny, par.nx)] != -1) {
                  neighbor_number = cell_numbers[idx3(zi, yi, (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[0].gnum = node_number;
                }
                if (zi && cell_numbers[idx3((zi - 1), yi, xi, par.ny, par.nx)] != -1) {
                  // there's a cell back in z
                  neighbor_number = cell_numbers[idx3((zi - 1), yi, xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[7].gnum = node_number - 1;
                  els[neighbor_number].vtx[6].gnum = node_number;
                }
                if (zi && yi && cell_numbers[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)] != -1) {
                  // there's a cell back in z, back in y
                  neighbor_number = cell_numbers[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[4].gnum = node_number - 1;
                  els[neighbor_number].vtx[5].gnum = node_number;
                }
                if (zi && xi < par.nx - 1 && cell_numbers[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)] != -1) {
                  // there's a cell back in z, forward in x
                  neighbor_number = cell_numbers[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[7].gnum = node_number;
                }
                if (zi && yi && xi < par.nx - 1 && cell_numbers[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                  // there's a cell back in z, back in y, forward in x
                  neighbor_number = cell_numbers[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[4].gnum = node_number;
                }
                continue;
              }
            one_new_node:
              { // assumes node count was already incremented and adjusts neighboring cells
                // neighbors
                if (xi < par.nx - 1 && cell_numbers[idx3(zi, yi, (xi + 1), par.ny, par.nx)] != -1) {
                  // cell forward in x
                  neighbor_number = cell_numbers[idx3(zi, yi, (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[0].gnum = node_number;
                }
                if (xi < par.nx - 1 && zi && cell_numbers[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)] != -1) {
                  // cell back in z, forward in x
                  neighbor_number = cell_numbers[idx3((zi - 1), yi, (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[7].gnum = node_number;
                }
                if (zi && cell_numbers[idx3((zi - 1), yi, xi, par.ny, par.nx)] != -1) {
                  // cell back in z
                  neighbor_number = cell_numbers[idx3((zi - 1), yi, xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[6].gnum = node_number;
                }
                if (yi && xi < par.nx - 1 && cell_numbers[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                  // cell back in y, forward in x
                  neighbor_number = cell_numbers[idx3(zi, (yi - 1), (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[3].gnum = node_number;
                }
                if (yi && zi && xi < par.nx - 1 && cell_numbers[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                  // cell back in z, back in y, forward in x
                  neighbor_number = cell_numbers[idx3((zi - 1), (yi - 1), (xi + 1), par.ny, par.nx)];
                  els[neighbor_number].vtx[4].gnum = node_number;
                }
                if (yi && zi && cell_numbers[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)] != -1) {
                  // cell back in z, back in y
                  neighbor_number = cell_numbers[idx3((zi - 1), (yi - 1), xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[5].gnum = node_number;
                }
                if (yi && cell_numbers[idx3(zi, (yi - 1), xi, par.ny, par.nx)] != -1) {
                  // cell back in y
                  neighbor_number = cell_numbers[idx3(zi, (yi - 1), xi, par.ny, par.nx)];
                  els[neighbor_number].vtx[2].gnum = node_number;
                }
              }
            }
          }
          // sweep for y-roof at this z-level
          for (int xi = 0; xi < par.nx; xi++) {
            if (par.voxel_geometry[idx3(zi, (par.ny - 1), xi, par.ny, par.nx)] != 1) { // live top cell
              if (xi == 0) {
                node_number += 2;
                cell_number = cell_numbers[idx3(zi, (par.ny - 1), xi, par.ny, par.nx)];
                els[cell_number].vtx[3].gnum = node_number - 1;
                els[cell_number].vtx[2].gnum = node_number;
                goto two_new_nodes_y_sweep;
              }
              else {
                node_number++;
                cell_number = cell_numbers[idx3(zi, (par.ny - 1), xi, par.ny, par.nx)];
                els[cell_number].vtx[2].gnum = node_number;
                goto one_new_node_y_sweep;
              }
            }
            else { // dead top cell
              // does this dead cell have anything to contribute?
              if (xi < par.nx - 1 && cell_numbers[idx3(zi, (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) goto dead_active_y_sweep;
              else if (xi < par.nx - 1 && zi && cell_numbers[idx3((zi - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) goto dead_active_y_sweep;
              else if (zi && cell_numbers[idx3((zi - 1), (par.ny - 1), xi, par.ny, par.nx)] != -1) goto dead_active_y_sweep;
              else continue;
            }
          dead_active_y_sweep:
            {
              if (xi == 0) {
                node_number += 2;
                goto two_new_nodes_y_sweep;
              }
              else {
                node_number++;
                goto one_new_node_y_sweep;
              }
            }
          two_new_nodes_y_sweep:
            {
              if (xi < par.nx - 1 && cell_numbers[idx3(zi, (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3(zi, (par.ny - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[3].gnum = node_number;
              }
              if (xi < par.nx - 1 && zi && cell_numbers[idx3((zi - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((zi - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number;
              }
              if (zi && cell_numbers[idx3((zi - 1), (par.ny - 1), xi, par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((zi - 1), (par.ny - 1), xi, par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number - 1;
                els[neighbor_number].vtx[5].gnum = node_number;
              }
              continue;
            }
          one_new_node_y_sweep:
            {
              if (xi < par.nx - 1 && cell_numbers[idx3(zi, (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3(zi, (par.ny - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[3].gnum = node_number;
              }
              if (xi < par.nx - 1 && zi && cell_numbers[idx3((zi - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((zi - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number;
              }
              if (zi && cell_numbers[idx3((zi - 1), (par.ny - 1), xi, par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((zi - 1), (par.ny - 1), xi, par.ny, par.nx)];
                els[neighbor_number].vtx[5].gnum = node_number;
              }
            }
          }
        }
        // last sweep for z ceiling
        for (int yi = 0; yi < par.ny; yi++) {
          for (int xi = 0; xi < par.nx; xi++) {
            if (par.voxel_geometry[idx3((par.nz - 1), yi, xi, par.ny, par.nx)] != 1) { // live top cell
              if (xi == 0) {
                node_number += 2;
                cell_number = cell_numbers[idx3((par.nz - 1), yi, xi, par.ny, par.nx)];
                els[cell_number].vtx[7].gnum = node_number - 1;
                els[cell_number].vtx[6].gnum = node_number;
                goto two_new_nodes_z_sweep;
              }
              else {
                node_number++;
                cell_number = cell_numbers[idx3((par.nz - 1), yi, xi, par.ny, par.nx)];
                els[cell_number].vtx[6].gnum = node_number;
                goto one_new_node_z_sweep;
              }
            }
            else { // dead top cell
              // does this dead cell contribute a node?
              if (xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), yi, (xi + 1), par.ny, par.nx)] != -1) goto dead_active_z_sweep;
              else if (yi && xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) goto dead_active_z_sweep;
              else if (yi && cell_numbers[idx3((par.nz - 1), (yi - 1), xi, par.ny, par.nx)] != -1) goto dead_active_z_sweep;
              else continue;
            }
          dead_active_z_sweep:
            {
              if (xi == 0) {
                node_number += 2;
                goto two_new_nodes_z_sweep;
              }
              else {
                node_number++;
                goto one_new_node_z_sweep;
              }
            }
          two_new_nodes_z_sweep:
            {
              if (xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), yi, (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), yi, (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[7].gnum = node_number;
              }
              if (yi && xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), (yi - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number;
              }
              if (yi && cell_numbers[idx3((par.nz - 1), (yi - 1), xi, par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), (yi - 1), xi, par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number - 1;
                els[neighbor_number].vtx[5].gnum = node_number;
              }
              continue;
            }
          one_new_node_z_sweep:
            {
              if (xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), yi, (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), yi, (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[7].gnum = node_number;
              }
              if (yi && xi < par.nx - 1 && cell_numbers[idx3((par.nz - 1), (yi - 1), (xi + 1), par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), (yi - 1), (xi + 1), par.ny, par.nx)];
                els[neighbor_number].vtx[4].gnum = node_number;
              }
              if (yi && cell_numbers[idx3((par.nz - 1), (yi - 1), xi, par.ny, par.nx)] != -1) {
                neighbor_number = cell_numbers[idx3((par.nz - 1), (yi - 1), xi, par.ny, par.nx)];
                els[neighbor_number].vtx[5].gnum = node_number;
              }
            }
          }
        }
        // sweep for y-roof at z-roof
        for (int xi = 0; xi < par.nx; xi++) {
          if (xi == 0 && par.voxel_geometry[idx3((par.nz - 1), (par.ny - 1), xi, par.ny, par.nx)] != 1) {
            node_number += 2;
            cell_number = cell_numbers[idx3((par.nz - 1), (par.ny - 1), xi, par.ny, par.nx)];
            els[cell_number].vtx[4].gnum = node_number - 1;
            els[cell_number].vtx[5].gnum = node_number;
            if (xi < par.nx - 1 && par.voxel_geometry[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != 1) {
              neighbor_number = cell_numbers[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)];
              els[neighbor_number].vtx[4].gnum = node_number;
            }
          }
          else if (par.voxel_geometry[idx3((par.nz - 1), (par.ny - 1), xi, par.ny, par.nx)] != 1) {
            node_number++;
            cell_number = cell_numbers[idx3((par.nz - 1), (par.ny - 1), xi, par.ny, par.nx)];
            els[cell_number].vtx[5].gnum = node_number;
            if (xi < par.nx - 1 && par.voxel_geometry[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != 1) {
              neighbor_number = cell_numbers[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)];
              els[neighbor_number].vtx[4].gnum = node_number;
            }
          }
          else if (xi < par.nx - 1 && par.voxel_geometry[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)] != 1) {
            node_number++;
            neighbor_number = cell_numbers[idx3((par.nz - 1), (par.ny - 1), (xi + 1), par.ny, par.nx)];
            els[neighbor_number].vtx[4].gnum = node_number;
          }
        }
      }
    }

#pragma omp barrier

#pragma omp for schedule(static,1) nowait
    for (int cell = 0; cell < nCells; cell++) {
      // 5 4 7 6 1 0 3 2
      gtlNode[idx2(els[cell].vtx[0].gnum, 5, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[1].gnum, 4, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[2].gnum, 7, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[3].gnum, 6, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[4].gnum, 1, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[5].gnum, 0, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[6].gnum, 3, 8)] = cell + 1;
      gtlNode[idx2(els[cell].vtx[7].gnum, 2, 8)] = cell + 1;
    }

  } // end of omp parallel region
#ifdef _MESH_PRINT_DEBUG
  for (int ii = 0; ii < els.size(); ii++) {
    std::cout << "el " << ii << ":\n";
    for (int jj = 0; jj < 8; jj++) {
      std::cout << "global node " << els[ii].vtx[jj].gnum << ":\t" << els[ii].vtx[jj].coords[0] << "\t" << els[ii].vtx[jj].coords[1] << "\t" << els[ii].vtx[jj].coords[2] << "\n";
    }
    std::cout << "\n";
  }

#endif
}

void
hgf::mesh::printVTK(const parameters& par)
{

  int nEls = (int)els.size();
  if (par.dimension == 3) {
    int nNodes = (int)gtlNode.size() / 8;
    std::vector<double> nodes(nNodes * 3);
    for (int ii = 0; ii < nNodes; ii++) {
      if (gtlNode[idx2(ii, 0, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[2];
      }
      else if (gtlNode[idx2(ii, 1, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[2];
      }
      else if (gtlNode[idx2(ii, 2, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[2];
      }
      else if (gtlNode[idx2(ii, 3, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[2];
      }
      else if (gtlNode[idx2(ii, 4, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[2];
      }
      else if (gtlNode[idx2(ii, 5, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[2];
      }
      else if (gtlNode[idx2(ii, 6, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[2];
      }
      else if (gtlNode[idx2(ii, 7, 8)]) {
        nodes[idx2(ii, 0, 3)] = els[gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 3)] = els[gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[1];
        nodes[idx2(ii, 2, 3)] = els[gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[2];
      }
    }

    std::ofstream meshvis;
    std::string mesh_file = par.problem_path.string() + "/Mesh.vtk";
    meshvis.open(mesh_file);

    meshvis << "# vtk DataFile Version 3.0\n";
    meshvis << "vtk output\n";
    meshvis << "ASCII\n\n";
    meshvis << "DATASET UNSTRUCTURED_GRID\n";
    meshvis << "POINTS " << nNodes << " double\n";
    for (int ii = 0; ii < nNodes; ii++) {
      meshvis << nodes[idx2(ii, 0, 3)] << "\t";
      meshvis << nodes[idx2(ii, 1, 3)] << "\t";
      meshvis << nodes[idx2(ii, 2, 3)] << "\n";
    }

    meshvis << "\n";
    meshvis << "CELLS " << nEls << " " << 9 * nEls << "\n";
    for (int ii = 0; ii < nEls; ii++) {
      meshvis << 8 << "\t";
      meshvis << els[ii].vtx[0].gnum << "\t";
      meshvis << els[ii].vtx[1].gnum << "\t";
      meshvis << els[ii].vtx[2].gnum << "\t";
      meshvis << els[ii].vtx[3].gnum << "\t";
      meshvis << els[ii].vtx[7].gnum << "\t";
      meshvis << els[ii].vtx[6].gnum << "\t";
      meshvis << els[ii].vtx[5].gnum << "\t";
      meshvis << els[ii].vtx[4].gnum << "\n";
    }
    meshvis << "\n";
    meshvis << "CELL_TYPES " << nEls << "\n";
    for (int ii = 0; ii < nEls; ii++) {
      meshvis << 12 << "\n";
    }

  }
  else {
    int nNodes = (int)gtlNode.size() / 4;
    // build an exclusive nodes vector
    std::vector<double> nodes(nNodes * 2);
    for (int ii = 0; ii < nNodes; ii++) {
      if (gtlNode[idx2(ii, 0, 4)]) {
        nodes[idx2(ii, 0, 2)] = els[gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 2)] = els[gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[1];
      }
      else if (gtlNode[idx2(ii, 1, 4)]) {
        nodes[idx2(ii, 0, 2)] = els[gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 2)] = els[gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[1];
      }
      else if (gtlNode[idx2(ii, 2, 4)]) {
        nodes[idx2(ii, 0, 2)] = els[gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 2)] = els[gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[1];
      }
      else if (gtlNode[idx2(ii, 3, 4)]) {
        nodes[idx2(ii, 0, 2)] = els[gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 2)] = els[gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[1];
      }
    }

    std::ofstream meshvis;
    std::string mesh_file = par.problem_path.string() + "/Mesh.vtk";
    meshvis.open(mesh_file);

    meshvis << "# vtk DataFile Version 3.0\n";
    meshvis << "vtk output\n";
    meshvis << "ASCII\n\n";
    meshvis << "DATASET UNSTRUCTURED_GRID\n";
    meshvis << "POINTS " << nNodes << " double\n";
    for (int ii = 0; ii < nNodes; ii++) {
      meshvis << nodes[idx2(ii, 0, 2)] << "\t";
      meshvis << nodes[idx2(ii, 1, 2)] << "\t";
      meshvis << 0.0 << "\n";
    }
    meshvis << "\n";
    meshvis << "CELLS " << nEls << " " << 5 * nEls << "\n";
    for (int ii = 0; ii < nEls; ii++) {
      meshvis << 4 << "\t";
      meshvis << els[ii].vtx[0].gnum << "\t";
      meshvis << els[ii].vtx[1].gnum << "\t";
      meshvis << els[ii].vtx[2].gnum << "\t";
      meshvis << els[ii].vtx[3].gnum << "\t";
    }
    meshvis << "\n";
    meshvis << "CELL_TYPES " << nEls << "\n";
    for (int ii = 0; ii < nEls; ii++) {
      meshvis << 9 << "\n";
    }
  }
}
