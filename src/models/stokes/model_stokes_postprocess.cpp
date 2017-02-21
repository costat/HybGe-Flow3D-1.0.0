/* stokes postprocess source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// NOTE: flux inserts currently assume a cartesian grid, e.g. flux's can be computed using single component distances. This should be extended.

/** \brief hgf::models::stokes::solution_build builds the solution vector using boundary information and the Stokes interior solution, solution_int.
 *
 */
void
hgf::models::stokes::solution_build(void)
{
  // Get sizes. In 2d, nW and velocity_w.size should both == 0]
  int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int nW = std::accumulate(interior_w.begin(), interior_w.end(), 0);
  int nVel = nU + nV + nW;
  solution.resize(velocity_u.size() + velocity_v.size() + velocity_w.size() + pressure.size());

  if (nW) {
#pragma omp parallel
    {
#pragma omp for schedule(static,1) nowait // u loop
      for (int ii = 0; ii < velocity_u.size(); ii++) {
        double dx;
        // am i an interior node?
        if (interior_u[ii]) solution[ii] = solution_int[interior_u_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii].type == 1) solution[ii] = boundary[ii].value;
        // neumann bc node
        else {
          // check for a node to the right, if yes, calculate using prescribed flux and value to the right
          if (velocity_u[ii].neighbors[1] != -1) {
            dx = velocity_u[velocity_u[ii].neighbors[1]].coords[0] - velocity_u[ii].coords[0];
            solution[ii] = solution_int[interior_u_nums[velocity_u[ii].neighbors[1]]] - boundary[ii].value * dx;
          }
          // else left, calculate using prescribed flux and u value to left
          else {
            dx = velocity_u[ii].coords[0] - velocity_u[velocity_u[ii].neighbors[3]].coords[0];
            solution[ii] = solution_int[interior_u_nums[velocity_u[ii].neighbors[3]]] + boundary[ii].value * dx;
          }
        }
      }
#pragma omp for schedule(static,1) nowait // v loop
      for (int ii = 0; ii < velocity_v.size(); ii++) {
        double dy;
        // am i an interior node?
        if (interior_v[ii]) solution[velocity_u.size() + ii] = solution_int[nU + interior_v_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii + velocity_u.size()].type == 1) solution[ii + velocity_u.size()] = boundary[ii + velocity_u.size()].value;
        // neumann bc node
        else {
          // check for a node in the y+ direction, if yes, calculate value using prescribed flux and y+ v value
          if (velocity_v[ii].neighbors[2] != -1) {
            dy = velocity_v[velocity_v[ii].neighbors[2]].coords[1] - velocity_v[ii].coords[1];
            solution[ii + velocity_u.size()] = solution_int[nU + interior_v_nums[velocity_v[ii].neighbors[2]]] - boundary[ii + velocity_u.size()].value * dy;
          }
          // else y-, calculate using prescribed flux and v value to y- direction
          else {
            dy = velocity_v[ii].coords[1] - velocity_v[velocity_v[ii].neighbors[0]].coords[1];
            solution[ii + velocity_u.size()] = solution_int[nU + interior_v_nums[velocity_v[ii].neighbors[0]]] + boundary[ii + velocity_u.size()].value * dy;
          }
        }
      }
#pragma omp for schedule(static,1) nowait // w loop
      for (int ii = 0; ii < velocity_w.size(); ii++) {
        double dz;
        // am i an interior node?
        if (interior_w[ii]) solution[velocity_u.size() + velocity_v.size() + ii] = solution_int[nU + nV + interior_w_nums[ii]];
        // am i a dirichlet boundary node?
        else if (boundary[ii + velocity_u.size() + velocity_v.size()].type == 1) solution[ii + velocity_u.size() + velocity_v.size()] = boundary[ii + velocity_u.size() + velocity_v.size()].value;
        // neumann bc node
        else {
          // check for a node in the z+ direction, if yes, calculate value using prescribed flux and z+ w value
          if (velocity_w[ii].neighbors[5] != -1) {
            dz = velocity_w[velocity_w[ii].neighbors[5]].coords[2] - velocity_w[ii].coords[2];
            solution[ii + velocity_u.size() + velocity_v.size()] = solution_int[nU + nV + interior_w_nums[velocity_w[ii].neighbors[5]]] - boundary[ii + velocity_u.size() + velocity_v.size()].value * dz;
          }
          // else z-, calculate using prescribed flux and w value in z- direction
          else {
            dz = velocity_w[ii].coords[2] - velocity_w[velocity_w[ii].neighbors[4]].coords[2];
            solution[ii + velocity_u.size() + velocity_v.size()] = solution_int[nU + nV + interior_w_nums[velocity_w[ii].neighbors[4]]] + boundary[ii + velocity_u.size() + velocity_v.size()].value * dz;
          }
        }
      }
#pragma omp for schedule(static,1) nowait // pressure loop
      for (int ii = 0; ii < pressure.size(); ii++) {
        solution[velocity_u.size() + velocity_v.size() + velocity_w.size() + ii] = solution_int[nVel + ii];
      }
    }
  }
  else {
#pragma omp parallel
    {
#pragma omp for schedule(static,1) nowait
      for (int ii = 0; ii < velocity_u.size(); ii++) {
        double dx;
        // am i an interior node?
        if (interior_u[ii]) solution[ii] = solution_int[interior_u_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii].type == 1) solution[ii] = boundary[ii].value;
        // neumann bc node
        else {
          // check for a node to the right, if yes, calculate value using prescribed flux and u value to the right
          if (velocity_u[ii].neighbors[1] != -1) {
            dx = velocity_u[velocity_u[ii].neighbors[1]].coords[0] - velocity_u[ii].coords[0];
            solution[ii] = solution_int[interior_u_nums[velocity_u[ii].neighbors[1]]] - boundary[ii].value * dx;
          }
          // else left, calculate value using prescribed flux and u value to left
          else {
            dx = velocity_u[ii].coords[0] - velocity_u[velocity_u[ii].neighbors[3]].coords[0];
            solution[ii] = solution_int[interior_u_nums[velocity_u[ii].neighbors[3]]] + boundary[ii].value * dx;
          }
        }
      }
#pragma omp for schedule(static,1) nowait
      for (int ii = 0; ii < velocity_v.size(); ii++) {
        double dy;
        // am i an interior node?
        if (interior_v[ii]) solution[velocity_u.size() + ii] = solution_int[nU + interior_v_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii + velocity_u.size()].type == 1) solution[ii + velocity_u.size()] = boundary[ii + velocity_u.size()].value;
        // neumann bc node
        else {
          // check for a node above, if yes, calculate value using prescribed flux and v value above
          if (velocity_v[ii].neighbors[2] != -1) {
            dy = velocity_v[velocity_v[ii].neighbors[2]].coords[1] - velocity_v[ii].coords[1];
            solution[ii + velocity_u.size()] = solution_int[nU + interior_v_nums[velocity_v[ii].neighbors[2]]] - boundary[ii + velocity_u.size()].value * dy;
          }
          // else below, calculate value using prescribed flux and v value below
          else {
            dy = velocity_v[ii].coords[1] - velocity_v[velocity_v[ii].neighbors[0]].coords[1];
            solution[ii + velocity_u.size()] = solution_int[nU + interior_v_nums[velocity_v[ii].neighbors[0]]] + boundary[ii + velocity_u.size()].value * dy;
          }
        }
      }
#pragma omp for schedule(static,1) nowait
      for (int ii = 0; ii < pressure.size(); ii++) {
        solution[velocity_u.size() + velocity_v.size() + ii] = solution_int[nVel + ii];
      }
    }
  }
}

/** \brief hgf::models::stokes::output_vtk saves the solution to the Stokes flow to a file for VTK visualiztion.
 *
 * @param[in] par - parameters struct containing problem information, including problem directory.
 * @param[in] msh - mesh object containing a quadrilateral or hexagonal representation of geometry from problem folder addressed in parameters& par.
 * @param[in,out] file_name - string used to name the output file, which is placed in the problem directory contained in parameters& par.
 */
void
hgf::models::stokes::output_vtk(const parameters& par, const hgf::mesh& msh, std::string& file_name)
{

  if (par.dimension == 3) { // 3d output
    double uval, vval, wval;
    int pzero, uzero, vzero, wzero;
    int nNodes = (int)msh.gtlNode.size() / 8;
    int nEls = (int)msh.els.size();
    uzero = 0;
    vzero = uzero + (int)velocity_u.size();
    wzero = vzero + (int)velocity_v.size();
    pzero = wzero + (int)velocity_w.size();
    // build an exclusive nodes vector
    std::vector<double> nodes(nNodes * 3);
#pragma omp parallel for
    for (int ii = 0; ii < nNodes; ii++) {
      if (msh.gtlNode[idx2(ii, 0, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 0, 8)] - 1].vtx[5].coords[2];
      }
      else if (msh.gtlNode[idx2(ii, 1, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 1, 8)] - 1].vtx[4].coords[2];
      }
      else if (msh.gtlNode[idx2(ii, 2, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 2, 8)] - 1].vtx[7].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 3, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 3, 8)] - 1].vtx[6].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 4, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 4, 8)] - 1].vtx[1].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 5, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 5, 8)] - 1].vtx[0].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 6, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 6, 8)] - 1].vtx[3].coords[2];

      }
      else if (msh.gtlNode[idx2(ii, 7, 8)]) {
        nodes[idx2(ii, 0, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[1];
        nodes[idx2(ii, 2, 3)] = msh.els[msh.gtlNode[idx2(ii, 7, 8)] - 1].vtx[2].coords[2];

      }
    }
    // write to vtk file
    bfs::path output_path(par.problem_path / file_name.c_str());
    output_path += ".vtk";
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET UNSTRUCTURED_GRID\n";
    outstream << "POINTS " << nNodes << " double\n";
    for (int row = 0; row < nNodes; row++) {
      outstream << nodes[idx2(row, 0, 3)] << "\t";
      outstream << nodes[idx2(row, 1, 3)] << "\t";
      outstream << nodes[idx2(row, 2, 3)] << "\n";
    }
    outstream << "\n";
    outstream << "CELLS " << nEls << " " << 9 * nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 8 << "\t";
      outstream << msh.els[row].vtx[0].gnum << "\t";
      outstream << msh.els[row].vtx[1].gnum << "\t";
      outstream << msh.els[row].vtx[2].gnum << "\t";
      outstream << msh.els[row].vtx[3].gnum << "\t";
      outstream << msh.els[row].vtx[7].gnum << "\t";
      outstream << msh.els[row].vtx[6].gnum << "\t";
      outstream << msh.els[row].vtx[5].gnum << "\t";
      outstream << msh.els[row].vtx[4].gnum << "\t";
    }
    outstream << "\n";
    outstream << "CELL_TYPES " << nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 12 << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nEls << "\n";
    outstream << "SCALARS pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < nEls; row++) {
      outstream << solution[pzero + row] << "\n";
    }
    outstream << "\n";
    outstream << "VECTORS velocity double\n";
    for (int row = 0; row < nEls; row++) {
      uval = 0.5 * (solution[uzero+ptv[idx2(row, 0, 6)]] \
                    + solution[uzero + ptv[idx2(row, 1, 6)]]);
      vval = 0.5 * (solution[vzero + ptv[idx2(row, 2, 6)]] \
                    + solution[vzero + ptv[idx2(row, 3, 6)]]);
      wval = 0.5 * (solution[wzero + ptv[idx2(row, 4, 6)]] \
                    + solution[wzero + ptv[idx2(row, 5, 6)]]);
      outstream << uval << "\t" << vval << "\t" << wval << "\n";
    }
    outstream << "\n";
    outstream.close();
  }
  else { // 2d output
    double uval, vval;
    int pzero, uzero, vzero;
    int nNodes = (int)msh.gtlNode.size() / 4;
    int nEls = (int)msh.els.size();
    uzero = 0;
    vzero = uzero + (int)velocity_u.size();
    pzero = vzero + (int)velocity_v.size();
    // build an exclusive nodes vector
    std::vector<double> nodes(nNodes * 2);
#pragma omp parallel for
    for (int ii = 0; ii < nNodes; ii++) {
      if (msh.gtlNode[idx2(ii, 0, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 0, 4)] - 1].vtx[2].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 1, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 1, 4)] - 1].vtx[3].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 2, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 2, 4)] - 1].vtx[1].coords[1];
      }
      else if (msh.gtlNode[idx2(ii, 3, 4)]) {
        nodes[idx2(ii, 0, 2)] = msh.els[msh.gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[0];
        nodes[idx2(ii, 1, 2)] = msh.els[msh.gtlNode[idx2(ii, 3, 4)] - 1].vtx[0].coords[1];
      }
    }
    // write solution vtk file
    bfs::path output_path( par.problem_path / file_name.c_str() );
    output_path += ".vtk";
    std::ofstream outstream;
    outstream.open(output_path.string());
    outstream << "# vtk DataFile Version 3.0\n";
    outstream << "vtk output\n";
    outstream << "ASCII\n\n";
    outstream << "DATASET UNSTRUCTURED_GRID\n";
    outstream << "POINTS " << nNodes << " double\n";
    for (int row = 0; row < nNodes; row++) {
      outstream << nodes[idx2(row, 0, 2)] << "\t";
      outstream << nodes[idx2(row, 1, 2)] << "\t";
      outstream << 0.0 << "\n";
    }
    outstream << "\n";
    outstream << "CELLS " << nEls << " " << 5 * nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 4 << "\t";
      outstream << msh.els[row].vtx[0].gnum << "\t";
      outstream << msh.els[row].vtx[1].gnum << "\t";
      outstream << msh.els[row].vtx[2].gnum << "\t";
      outstream << msh.els[row].vtx[3].gnum << "\t";
    }
    outstream << "\n";
    outstream << "CELL_TYPES " << nEls << "\n";
    for (int row = 0; row < nEls; row++) {
      outstream << 9 << "\n";
    }
    outstream << "\n";
    outstream << "CELL_DATA " << nEls << "\n";
    outstream << "SCALARS pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < nEls; row++) {
      outstream << solution[pzero + row] << "\n";
    }
    outstream << "\n";
    outstream << "VECTORS velocity double\n";
    for (int row = 0; row < nEls; row++) {
      uval = 0.5 * (solution[uzero+ptv[idx2(row, 0, 4)]] \
                    + solution[uzero + ptv[idx2(row, 1, 4)]]);
      vval = 0.5 * (solution[vzero + ptv[idx2(row, 2, 4)]] \
                    + solution[vzero + ptv[idx2(row, 3, 4)]]);
      outstream << uval << "\t" << vval << "\t" << 0 << "\n";
    }
    outstream << "\n";
    outstream.close();
  }
}
