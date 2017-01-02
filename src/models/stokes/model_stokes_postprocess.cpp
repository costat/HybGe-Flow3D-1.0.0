/* stokes postprocess source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::solution_insert_boundaries(void)
{
  // Get sizes. In 2d, nW and velocity_w.size should both == 0
  int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int nW = std::accumulate(interior_w.begin(), interior_w.end(), 0);
  int nVel = nU + nV + nW;
  solution.resize(velocity_u.size() + velocity_v.size() + velocity_w.size() + pressure.size());
    
  if (nW) {
#pragma omp parallel
    {



    }
  }
  else {
#pragma omp parallel
    {
#pragma omp for schedule(dynamic) nowait
      for (int ii = 0; ii < velocity_u.size(); ii++) {
        double dx;
        // am i an interior node?
        if (interior_u[ii]) solution[ii] = interior[interior_u_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii].type == 1) solution[ii] = boundary[ii].value;
        // neumann bc node
        else {
          // check for a node to the right, if yes, calculate value using prescribed flux and u value to the right
          if (velocity_u[ii].neighbors[1] != -1) {
            dx = velocity_u[velocity_u[ii].neighbors[1]].coords[0] - velocity_u[ii].coords[0];
            solution[ii] = interior[interior_u_nums[velocity_u[ii].neighbors[1]]] - boundary[ii].value * dx;
          }
          // else left, calculate value using prescribed flux and u value to left
          else {
            dx = velocity_u[ii].coords[0] - velocity_u[velocity_u[ii].neighbors[3]].coords[0];
            solution[ii] = interior[interior_u_nums[velocity_u[ii].neighbors[3]]] + boundary[ii].value * dx;
          }
        }
      }
#pragma omp for schedule(dynamic) nowait
      for (int ii = 0; ii < velocity_v.size(); ii++) {
        double dy;
        // am i an interior node?
        if (interior_v[ii]) solution[velocity_u.size() + ii] = interior[nU + interior_v_nums[ii]];
        // am i a dirichlet bc node?
        else if (boundary[ii + velocity_u.size()].type == 1) solution[ii + velocity_u.size()] = boundary[ii + velocity_u.size()].value;
        // neumann bc node
        else {
          // check for a node above, if yes, calculate value using prescribed flux and v value above
          if (velocity_v[ii].neighbors[2] != -1) {
            dy = velocity_v[velocity_v[ii].neighbors[2]].coords[1] - velocity_v[ii].coords[1];
            solution[ii + velocity_u.size()] = interior[nU + interior_v_nums[velocity_v[ii].neighbors[2]]] - boundary[ii + velocity_u.size()].value * dy;
          }
          // else below, calculate value using prescribed flux and v value below
          else {
            dy = velocity_v[ii].coords[1] - velocity_v[velocity_v[ii].neighbors[0]].coords[1];
            solution[ii + velocity_u.size()] = interior[nU + interior_v_nums[velocity_v[ii].neighbors[0]]] + boundary[ii + velocity_u.size()].value * dy;
          }
        }
      }
#pragma omp for schedule(dynamic) nowait
      for (int ii = 0; ii < pressure.size(); ii++) {
        solution[velocity_u.size() + velocity_v.size() + velocity_w.size() + ii] = interior[nVel + ii];
      }
    }
  }
}

void
hgf::models::stokes::output_vtk(const parameters& par, const hgf::mesh& msh)
{

  if (par.nz) { // 3d output
    //--- not yet implemented ---//
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
    // post process the solution to get colocated values for velocities and pressure
    bfs::path output_path( par.problem_path / "Solution.vtk" );
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
#ifdef _DOF_VISUAL_DEBUG // sets solution values based on coordinates for dof sanity check
    outstream << "\n";
    outstream << "CELL_DATA " << nEls << "\n";
    outstream << "SCALARS pressure double\n";
    outstream << "LOOKUP_TABLE default\n";
    for (int row = 0; row < nEls; row++) {
      outstream << (pressure[row].coords[0] + pressure[row].coords[1]) << "\n";
    }
    outstream << "\n";
    outstream << "VECTORS velocity double\n";
    for (int row = 0; row < nEls; row++) {
      uval = 0.5 * (velocity_u[ptv[idx2(row, 0, 4)]].coords[0] \
                    + velocity_u[ptv[idx2(row, 1, 4)]].coords[0]);
      vval = 0.5 * (velocity_v[ptv[idx2(row, 2, 4)]].coords[1] \
                    + velocity_v[ptv[idx2(row, 3, 4)]].coords[1]);
      outstream << uval << "\t" << vval << "\t" << 0 << "\n";
    }
    outstream << "\n";
#else // actual output section
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
#endif
    outstream.close();
  }
}
