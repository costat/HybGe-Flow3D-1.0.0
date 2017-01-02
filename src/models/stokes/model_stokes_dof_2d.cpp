/* stokes dof 2d source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh& msh)
{
  // this functions sets degrees of freedom for the velocity components and pressure in 2d
  velocity_u.reserve((par.nx + 1)*par.ny);
  velocity_v.reserve((par.ny + 1)*par.nx);
  pressure.reserve(msh.els.size());

#pragma omp parallel
  {
#pragma omp sections
    {
#pragma omp section
      { // u section
        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          // if there's no neighbor cell to the left, then we have 2 new dofs for u
          if (msh.els[cell].edg[3].neighbor == -1) {
            //--- dof on edge 3 ---//
            // coordinates
            dof_temp.coords[0] = 0.5 * \
              (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[3].coords[0]);
            dof_temp.coords[1] = 0.5 * \
              (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[3].coords[1]);
            dof_temp.coords[2] = 0;
            // doftype
            dof_temp.doftype = 1;
            // cell_numbers
            dof_temp.cell_numbers[1] = cell;
            dof_temp.cell_numbers[0] = -1;
            // push_back
            velocity_u.push_back(dof_temp);
          }
          //--- dof on edge 1 ---//
          // coordinates
          dof_temp.coords[0] = 0.5 * \
            (msh.els[cell].vtx[1].coords[0] + msh.els[cell].vtx[2].coords[0]);
          dof_temp.coords[1] = 0.5 * \
            (msh.els[cell].vtx[1].coords[1] + msh.els[cell].vtx[2].coords[1]);
          dof_temp.coords[2] = 0;
          // doftype
          dof_temp.doftype = 1;
          // cell_numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = msh.els[cell].edg[1].neighbor;
          // push_back
          velocity_u.push_back(dof_temp);
        }
      }
#pragma omp section
      { // v section
        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          // if there's no neighbor cell below, then we have 2 new dofs for v
          if (msh.els[cell].edg[0].neighbor == -1) {
            //--- dof on edge 0 ---//
            // coordinates
            dof_temp.coords[0] = 0.5 * \
              (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[1].coords[0]);
            dof_temp.coords[1] = 0.5 * \
              (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[1].coords[1]);
            dof_temp.coords[2] = 0;
            // doftype
            dof_temp.doftype = 1;
            // cell_numbers
            dof_temp.cell_numbers[1] = cell;
            dof_temp.cell_numbers[0] = -1;
            // push_back
            velocity_v.push_back(dof_temp);
          }
          //--- dof on edge 2 ---//
          // coordinates
          dof_temp.coords[0] = 0.5 * \
            (msh.els[cell].vtx[3].coords[0] + msh.els[cell].vtx[2].coords[0]);
          dof_temp.coords[1] = 0.5 * \
            (msh.els[cell].vtx[3].coords[1] + msh.els[cell].vtx[2].coords[1]);
          dof_temp.coords[2] = 0;
          // doftype
          dof_temp.doftype = 1;
          // cell_numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = msh.els[cell].edg[2].neighbor;
          // push_back
          velocity_v.push_back(dof_temp);
        }
      }
#pragma omp section
      { // pressure section
        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          // coordinates
          dof_temp.coords[0] = 0.25 * \
            (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[1].coords[0] + \
              msh.els[cell].vtx[2].coords[0] + msh.els[cell].vtx[3].coords[0]);
          dof_temp.coords[1] = 0.25 * \
            (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[1].coords[1] + \
              msh.els[cell].vtx[2].coords[1] + msh.els[cell].vtx[3].coords[1]);
          // doftype
          dof_temp.doftype = 0;
          // cell numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = -1;
          // neighbors
          for (int nbr = 0; nbr < 4; nbr++) {
            dof_temp.neighbors[nbr] = msh.els[cell].edg[nbr].neighbor;
          }
          // push_back
          pressure.push_back(dof_temp);
        }
      }
    }
  }

  // sort the v component so that 2nd block of linear system is a standard Poisson array
  std::sort(velocity_v.begin(), velocity_v.end(), byXbyY());

#ifdef _DOF_SORT_DEBUG
  std::cout << "\nChecking velocity sort U:\n";
  for (int ii = 0; ii < velocity_u.size(); ii++) {
    std::cout << velocity_u[ii].coords[0] << "\t" << velocity_u[ii].coords[1] << "\n";
  }
  std::cout << "\nChecking velocity sort V:\n";
  for (int ii = 0; ii < velocity_v.size(); ii++) {
    std::cout << velocity_v[ii].coords[0] << "\t" << velocity_v[ii].coords[1] << "\n";
  }
#endif

  // pressure to velocity relationships
  ptv.resize(pressure.size() * 4, -1);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      // pressure cells should know who u cells are
      if (velocity_u[ii].cell_numbers[0] != -1) ptv[idx2(velocity_u[ii].cell_numbers[0], 1, 4)] = ii;
      if (velocity_u[ii].cell_numbers[1] != -1) ptv[idx2(velocity_u[ii].cell_numbers[1], 0, 4)] = ii;
    }
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      // pressure cells should know who v cells are
      if (velocity_v[ii].cell_numbers[0] != -1) ptv[idx2(velocity_v[ii].cell_numbers[0], 3, 4)] = ii;
      if (velocity_v[ii].cell_numbers[1] != -1) ptv[idx2(velocity_v[ii].cell_numbers[1], 2, 4)] = ii;
    }
  }

  // set neighbors for velocity components
  dof_neighbors_2d(par, msh);

  interior_u.resize(velocity_u.size());
  interior_v.resize(velocity_v.size());

#pragma omp parallel
  {
    // determine boundary and interior cells
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      if (velocity_u[ii].neighbors[1] != -1 && velocity_u[ii].neighbors[3] != -1) {
        interior_u[ii] = 1;
      }
    }
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      if (velocity_v[ii].neighbors[0] != -1 && velocity_v[ii].neighbors[2] != -1) {
        interior_v[ii] = 1;
      }
    }
  }

  interior_u_nums.resize(interior_u.size(), -1);
  interior_v_nums.resize(interior_v.size(), -1);
  int cell = -1;
  for (int ii = 0; ii < interior_u.size(); ii++)
    if (interior_u[ii]) {
      cell++;
      interior_u_nums[ii] = cell;
    }
  cell = -1;
  for (int ii = 0; ii < interior_v.size(); ii++)
    if (interior_v[ii]) {
      cell++;
      interior_v_nums[ii] = cell;
    }

#ifdef _DOF_NEIGHBOR_DEBUG
  std::cout << "\nVelocity numbers in each pressure cell:\n";
  for (int ii = 0; ii < (ptv.size()/4); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 4; jj++) {
      std::cout << "\t" << ptv[idx2(ii, jj, 4)];
    }
    std::cout << "\n";
  }
  std::cout << "\nU velocity neighbor lists:\n";
  for (int ii = 0; ii < velocity_u.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 4; jj++) {
      std::cout << "\t" << velocity_u[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
  std::cout << "\nV velocity neighbor lists:\n";
  for (int ii = 0; ii < velocity_v.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 4; jj++) {
      std::cout << "\t" << velocity_v[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
  std::cout << "\npressure neighbor lists:\n";
  for (int ii = 0; ii < pressure.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 4; jj++) {
      std::cout << "\t" << pressure[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
#endif
}

void
hgf::models::stokes::dof_neighbors_2d(const parameters& par, const hgf::mesh& msh)
{
 
#pragma omp parallel
  {
#pragma omp for schedule(dynamic)
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      int no_neighbor_u[4] = { 1, 1, 1, 1 };
      //-- neighbor 0 --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[0] != -1) {
          // lower neighbor exists through pressure on the left
          no_neighbor_u[0] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[0];
          velocity_u[ii].neighbors[0] = ptv[idx2(pcell, 1, 4)];
        }
      if (no_neighbor_u[0] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[0] != -1) {
          // lower neighbor exists through pressure on the right
          no_neighbor_u[0] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[0];
          velocity_u[ii].neighbors[0] = ptv[idx2(pcell, 0, 4)];
        }
      if (no_neighbor_u[0]) velocity_u[ii].neighbors[0] = -1;

      //-- neighbor 1 --//
      if (velocity_u[ii].cell_numbers[1] != -1) {
        no_neighbor_u[1] = 0;
        velocity_u[ii].neighbors[1] = ptv[idx2(velocity_u[ii].cell_numbers[1], 1, 4)];
      }
      if (no_neighbor_u[1]) velocity_u[ii].neighbors[1] = -1;

      //-- neighbor 2 --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[2] != -1) {
          // upper neighbor exists through pressure on left
          no_neighbor_u[2] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[2];
          velocity_u[ii].neighbors[2] = ptv[idx2(pcell, 1, 4)];
        }
      if (no_neighbor_u[2] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[2] != -1) {
          // upper neighbor exists through pressure on right
          no_neighbor_u[2] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[2];
          velocity_u[ii].neighbors[2] = ptv[idx2(pcell, 0, 4)];
        }
      if (no_neighbor_u[2]) velocity_u[ii].neighbors[2] = -1;

      //-- neighbor 3 --//
      if (velocity_u[ii].cell_numbers[0] != -1) {
        no_neighbor_u[3] = 0;
        velocity_u[ii].neighbors[3] = ptv[idx2(velocity_u[ii].cell_numbers[0], 0, 4)];
      }
      if (no_neighbor_u[3]) velocity_u[ii].neighbors[3] = -1;
    }
#pragma omp for schedule(dynamic) nowait 
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      int no_neighbor_v[4] = { 1, 1, 1, 1 };
      //-- neighbor 0 --//
      if (velocity_v[ii].cell_numbers[0] != -1) {
        no_neighbor_v[0] = 0;
        velocity_v[ii].neighbors[0] = ptv[idx2(velocity_v[ii].cell_numbers[0], 2, 4)];
      }
      if (no_neighbor_v[0]) velocity_v[ii].neighbors[0] = -1;

      //-- neighbor 1 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[1] != -1) {
          // right neighbor exists through pressure below
          no_neighbor_v[1] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[1];
          velocity_v[ii].neighbors[1] = ptv[idx2(pcell, 3, 4)];
        }
      if (no_neighbor_v[1] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[1] != -1) {
          // right neighbor exists through pressure above
          no_neighbor_v[1] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[1];
          velocity_v[ii].neighbors[1] = ptv[idx2(pcell, 2, 4)];
        }
      if (no_neighbor_v[1]) velocity_v[ii].neighbors[1] = -1;

      //-- neighbor 2 --//
      if (velocity_v[ii].cell_numbers[1] != -1) {
        no_neighbor_v[2] = 0;
        velocity_v[ii].neighbors[2] = ptv[idx2(velocity_v[ii].cell_numbers[1], 3, 4)];
      }
      if (no_neighbor_v[2]) velocity_v[ii].neighbors[2] = -1;

      //-- neighbor 3 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[3] != -1) {
          // left neighbor exists through pressure below
          no_neighbor_v[3] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[3];
          velocity_v[ii].neighbors[3] = ptv[idx2(pcell, 3, 4)];
        }
      if (no_neighbor_v[3] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[3] != -1) {
          // left neighbor exists through pressure above
          no_neighbor_v[3] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[3];
          velocity_v[ii].neighbors[3] = ptv[idx2(pcell, 2, 4)];
        }
      if (no_neighbor_v[3]) velocity_v[ii].neighbors[3] = -1;
    }
  }

}
