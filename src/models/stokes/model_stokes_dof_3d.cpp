/* stokes dof 3d source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh& msh)
{
  // this function sets the degrees of freedom for the velocity components and pressure in 3d
  velocity_u.reserve((par.nx + 1)*par.ny*par.nz);
  velocity_v.reserve((par.ny + 1)*par.nz*par.nx);
  velocity_w.reserve((par.nz + 1)*par.nx*par.ny);
  pressure.resize(msh.els.size());

#pragma omp parallel
  {
#pragma omp sections nowait
    {
#pragma omp section
      { // u section
        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          if (msh.els[cell].fac[3].neighbor == -1) {
            // if there's no neighbor cell backwards in x, then we have 2 new dofs for u
            dof_temp.coords[0] = 0.5 * \
              (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[3].coords[0]);
            dof_temp.coords[1] = 0.5 * \
              (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[3].coords[1]);
            dof_temp.coords[2] = 0.5 * \
              (msh.els[cell].vtx[0].coords[2] + msh.els[cell].vtx[7].coords[2]);
            // doftype
            dof_temp.doftype = 1;
            // cell numbers
            dof_temp.cell_numbers[1] = cell;
            dof_temp.cell_numbers[0] = -1;
            // push back
            velocity_u.push_back(dof_temp);
          }
          // -- dof on face 1 --//
          // coordinates
          dof_temp.coords[0] = 0.5 * \
            (msh.els[cell].vtx[1].coords[0] + msh.els[cell].vtx[2].coords[0]);
          dof_temp.coords[1] = 0.5 * \
            (msh.els[cell].vtx[1].coords[1] + msh.els[cell].vtx[2].coords[1]);
          dof_temp.coords[2] = 0.5 * \
            (msh.els[cell].vtx[1].coords[2] + msh.els[cell].vtx[6].coords[2]);
          // doftype
          dof_temp.doftype = 1;
          // cell_numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = msh.els[cell].fac[1].neighbor;
          // push_back
          velocity_u.push_back(dof_temp);
        }
      }
#pragma omp section
      { // v section

        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          if (msh.els[cell].fac[0].neighbor == -1) {
            // if there's no neighbor back in y, then we have 2 new dofs for v
            //--- dof on face 0 ---//
            // coordinates
            dof_temp.coords[0] = 0.5 * \
              (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[1].coords[0]);
            dof_temp.coords[1] = 0.5 * \
              (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[1].coords[1]);
            dof_temp.coords[2] = 0.5 * \
              (msh.els[cell].vtx[0].coords[2] + msh.els[cell].vtx[7].coords[2]);
            // doftype
            dof_temp.doftype = 1;
            // cell_numbers
            dof_temp.cell_numbers[1] = cell;
            dof_temp.cell_numbers[0] = -1;
            // push_back
            velocity_v.push_back(dof_temp);
          }
          //--- dof on face 2 ---//
          // coordinates
          dof_temp.coords[0] = 0.5 * \
            (msh.els[cell].vtx[3].coords[0] + msh.els[cell].vtx[2].coords[0]);
          dof_temp.coords[1] = 0.5 * \
            (msh.els[cell].vtx[3].coords[1] + msh.els[cell].vtx[2].coords[1]);
          dof_temp.coords[2] = 0.5 * \
            (msh.els[cell].vtx[3].coords[2] + msh.els[cell].vtx[4].coords[2]);
          // doftype
          dof_temp.doftype = 1;
          // cell_numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = msh.els[cell].fac[2].neighbor;
          // push_back
          velocity_v.push_back(dof_temp);
        }
      }
#pragma omp section
      { // w section
        degree_of_freedom dof_temp;
        for (int cell = 0; cell < msh.els.size(); cell++) {
          if (msh.els[cell].fac[4].neighbor == -1) {
            // if there's no neighbor back in z, then we have 2 new dofs for w
            //--- dof on face 4 ---//
            // coordinates
            dof_temp.coords[0] = 0.5 * \
              (msh.els[cell].vtx[0].coords[0] + msh.els[cell].vtx[1].coords[0]);
            dof_temp.coords[1] = 0.5 * \
              (msh.els[cell].vtx[0].coords[1] + msh.els[cell].vtx[3].coords[1]);
            dof_temp.coords[2] = 0.5 * \
              (msh.els[cell].vtx[0].coords[2] + msh.els[cell].vtx[3].coords[2]);
            // doftype
            dof_temp.doftype = 1;
            // cell_numbers
            dof_temp.cell_numbers[1] = cell;
            dof_temp.cell_numbers[0] = -1;
            // push back
            velocity_w.push_back(dof_temp);
          }
          //--- dof on face 5 ---//
          // coordinates
          dof_temp.coords[0] = 0.5 * \
            (msh.els[cell].vtx[7].coords[0] + msh.els[cell].vtx[6].coords[0]);
          dof_temp.coords[1] = 0.5 * \
            (msh.els[cell].vtx[7].coords[1] + msh.els[cell].vtx[4].coords[1]);
          dof_temp.coords[2] = 0.5 * \
            (msh.els[cell].vtx[7].coords[2] + msh.els[cell].vtx[4].coords[2]);
          // doftype
          dof_temp.doftype = 1;
          // cell_numbers
          dof_temp.cell_numbers[0] = cell;
          dof_temp.cell_numbers[1] = msh.els[cell].fac[5].neighbor;
          // push back
          velocity_w.push_back(dof_temp);
        }
      }
    }
#pragma omp for // pressure section can be done in parallel and doesn't need to wait for the velocities
    for (int cell = 0; cell < msh.els.size(); cell++) {
      // p section
      degree_of_freedom dof_temp;
      // coordinates
      for (int dir = 0; dir < 3; dir++) {
        dof_temp.coords[dir] = 0;
        for (int ii = 0; ii < 8; ii++) dof_temp.coords[dir] += msh.els[cell].vtx[ii].coords[dir];
        dof_temp.coords[dir] = dof_temp.coords[dir] / 8;
      }
      // doftype
      dof_temp.doftype = 0;
      // cell numbers
      dof_temp.cell_numbers[0] = cell;
      dof_temp.cell_numbers[1] = -1;
      // neighbors
      for (int nbr = 0; nbr < 6; nbr++) dof_temp.neighbors[nbr] = msh.els[cell].fac[nbr].neighbor;
      // place pressure dof
      pressure[cell] = dof_temp;
    }
  }
  // sort the v component for matrix condition #
  std::sort(velocity_v.begin(), velocity_v.end(), byZbyXbyY());

  // sort the w component for matrix condition #
  std::sort(velocity_w.begin(), velocity_w.end(), byYbyXbyZ());

#ifdef _DOF_SORT_DEBUG
  std::cout << "\nChecking velocity sort U:\n";
  for (int ii = 0; ii < velocity_u.size(); ii++) {
    std::cout << velocity_u[ii].coords[0] << "\t" << velocity_u[ii].coords[1] << "\t" << velocity_u[ii].coords[2] << "\n";
  }
  std::cout << "\nChecking velocity sort V:\n";
  for (int ii = 0; ii < velocity_v.size(); ii++) {
    std::cout << velocity_v[ii].coords[0] << "\t" << velocity_v[ii].coords[1] << "\t" << velocity_v[ii].coords[2] << "\n";
  }
  std::cout << "\nChecking velocity sort W:\n";
  for (int ii = 0; ii < velocity_w.size(); ii++) {
    std::cout << velocity_w[ii].coords[0] << "\t" << velocity_w[ii].coords[1] << "\t" << velocity_w[ii].coords[2] << "\n";
  }
#endif

#ifdef _DOF_CELL_NUMBERS_DEBUG
  std::cout << "\nChecking U Pressure Numbers:\n";
  for (int ii = 0; ii < velocity_u.size(); ii++) {
    std::cout << velocity_u[ii].cell_numbers[0] << "\t" << velocity_u[ii].cell_numbers[1] << "\n";
  }
  std::cout << "\nChecking V Pressure Numbers:\n";
  for (int ii = 0; ii < velocity_v.size(); ii++) {
    std::cout << velocity_v[ii].cell_numbers[0] << "\t" << velocity_v[ii].cell_numbers[1] << "\n";
  }
  std::cout << "\nChecking W Pressure Numbers:\n";
  for (int ii = 0; ii < velocity_w.size(); ii++) {
    std::cout << velocity_w[ii].cell_numbers[0] << "\t" << velocity_w[ii].cell_numbers[1] << "\n";
  }
#endif

  // pressure to velocity relationships
  ptv.resize(pressure.size() * 6, -1);
#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      // pressure cells should know who u cells are
      if (velocity_u[ii].cell_numbers[0] != -1) ptv[idx2(velocity_u[ii].cell_numbers[0], 1, 6)] = ii;
      if (velocity_u[ii].cell_numbers[1] != -1) ptv[idx2(velocity_u[ii].cell_numbers[1], 0, 6)] = ii;
    }
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      // pressure cells should know who v cells are
      if (velocity_v[ii].cell_numbers[0] != -1) ptv[idx2(velocity_v[ii].cell_numbers[0], 3, 6)] = ii;
      if (velocity_v[ii].cell_numbers[1] != -1) ptv[idx2(velocity_v[ii].cell_numbers[1], 2, 6)] = ii;
    }
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_w.size(); ii++) {
      // pressure cells should know who w cells are
      if (velocity_w[ii].cell_numbers[0] != -1) ptv[idx2(velocity_w[ii].cell_numbers[0], 5, 6)] = ii;
      if (velocity_w[ii].cell_numbers[1] != -1) ptv[idx2(velocity_w[ii].cell_numbers[1], 4, 6)] = ii;
    }
  }
#ifdef _DOF_PTV_DEBUG
  std::cout << "\nChecking PTV Array:\n";
  for (int ii = 0; ii < ptv.size() / 6; ii++) {
    for (int jj = 0; jj < 6; jj++) {
      std::cout << ptv[idx2(ii, jj, 6)] << "\t";
    }
    std::cout << "\n";
  }
#endif
  
  // set neighbors for velocity components
  dof_neighbors_3d(par, msh);

  // determine boundary and interior cells
  interior_u.resize(velocity_u.size());
  interior_v.resize(velocity_v.size());
  interior_w.resize(velocity_w.size());
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
#pragma omp for schedule(dynamic) nowait
    for (int ii = 0; ii < velocity_w.size(); ii++) {
      if (velocity_w[ii].neighbors[4] != -1 && velocity_w[ii].neighbors[5] != -1) {
        interior_w[ii] = 1;
      }
    }
  }

#ifdef _INTERIOR_DEBUG
  std::cout << "\nInterior U Lists:\n";
  for (int ii = 0; ii < interior_u.size(); ii++) {
    std::cout << interior_u[ii] << "\n";
  }
  std::cout << "\nInterior V Lists:\n";
  for (int ii = 0; ii < interior_v.size(); ii++) {
    std::cout << interior_v[ii] << "\n";
  }
  std::cout << "\nInterior W Lists:\n";
  for (int ii = 0; ii < interior_w.size(); ii++) {
    std::cout << interior_w[ii] << "\n";
  }
#endif

  interior_u_nums.resize(interior_u.size(), -1);
  interior_v_nums.resize(interior_v.size(), -1);
  interior_w_nums.resize(interior_w.size(), -1);
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
  cell = -1;
  for (int ii = 0; ii < interior_w.size(); ii++)
    if (interior_w[ii]) {
      cell++;
      interior_w_nums[ii] = cell;
    }

#ifdef _DOF_NEIGHBOR_DEBUG
  std::cout << "\nVelocity numbers in each pressure cell:\n";
  for (int ii = 0; ii < (ptv.size() / 6); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 6; jj++) {
      std::cout << "\t" << ptv[idx2(ii, jj, 6)];
    }
    std::cout << "\n";
  }
  std::cout << "\nU velocity neighbor lists:\n";
  for (int ii = 0; ii < velocity_u.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 6; jj++) {
      std::cout << "\t" << velocity_u[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
  std::cout << "\nV velocity neighbor lists:\n";
  for (int ii = 0; ii < velocity_v.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 6; jj++) {
      std::cout << "\t" << velocity_v[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
  std::cout << "\nW velocity neighbor lists:\n";
  for (int ii = 0; ii < velocity_w.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 6; jj++) {
      std::cout << "\t" << velocity_w[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
  std::cout << "\npressure neighbor lists:\n";
  for (int ii = 0; ii < pressure.size(); ii++) {
    std::cout << ii << "\t";
    for (int jj = 0; jj < 6; jj++) {
      std::cout << "\t" << pressure[ii].neighbors[jj];
    }
    std::cout << "\n";
  }
#endif

}

void
hgf::models::stokes::dof_neighbors_3d(const parameters& par, const hgf::mesh& msh)
{
#pragma omp parallel
  {

#pragma omp for schedule(dynamic) nowait // u loop
    for (int ii = 0; ii < velocity_u.size(); ii++) {
      int no_neighbor_u[6] = { 1, 1, 1, 1, 1, 1 };
      //-- neighbor 0 (y- direction) --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[0] != -1) {
          // y- neighbor exists through pressure in x- direction
          no_neighbor_u[0] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[0];
          velocity_u[ii].neighbors[0] = ptv[idx2(pcell, 1, 6)];
        }
      if (no_neighbor_u[0] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[0] != -1) {
          // y- neighbor exists through pressure in x+ direction
          no_neighbor_u[0] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[0];
          velocity_u[ii].neighbors[0] = ptv[idx2(pcell, 0, 6)];
        }
      if (no_neighbor_u[0]) velocity_u[ii].neighbors[0] = -1;

      //-- neighbor 1 (x+ direction) --//
      if (velocity_u[ii].cell_numbers[1] != -1) {
        // x+ neighbor exists through pressure on right
        no_neighbor_u[1] = 0;
        velocity_u[ii].neighbors[1] = ptv[idx2(velocity_u[ii].cell_numbers[1], 1, 6)];
      }
      if (no_neighbor_u[1]) velocity_u[ii].neighbors[1] = -1;

      //-- neighbor 2 (y+ direction) --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[2] != -1) {
          // y+ neighbor exists through pressure on left
          no_neighbor_u[2] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[2];
          velocity_u[ii].neighbors[2] = ptv[idx2(pcell, 1, 6)];
        }
      if (no_neighbor_u[2] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[2] != -1) {
          // y+ neighbor exists through pressure on right
          no_neighbor_u[2] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[2];
          velocity_u[ii].neighbors[2] = ptv[idx2(pcell, 0, 6)];
        }
      if (no_neighbor_u[2]) velocity_u[ii].neighbors[2] = -1;

      //-- neighbor 3 (x- direction) --//
      if (velocity_u[ii].cell_numbers[0] != -1) {
        no_neighbor_u[3] = 0;
        velocity_u[ii].neighbors[3] = ptv[idx2(velocity_u[ii].cell_numbers[0], 0, 6)];
      }
      if (no_neighbor_u[3]) velocity_u[ii].neighbors[3] = -1;

      //-- neighbor 4 (z- direction) --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[4] != -1) {
          // z- neighbor exists through pressure on left
          no_neighbor_u[4] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[4];
          velocity_u[ii].neighbors[4] = ptv[idx2(pcell, 1, 6)];
        }
      if (no_neighbor_u[4] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[4] != -1) {
          // z- neighbor exists through pressure on the right
          no_neighbor_u[4] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[4];
          velocity_u[ii].neighbors[4] = ptv[idx2(pcell, 0, 6)];
        }
      if (no_neighbor_u[4]) velocity_u[ii].neighbors[4] = -1;

      //-- neighbor 5 (z+ direction) --//
      if (velocity_u[ii].cell_numbers[0] != -1)
        if (pressure[velocity_u[ii].cell_numbers[0]].neighbors[5] != -1) {
          // z+ neighbor exists through pressure on the left
          no_neighbor_u[5] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[0]].neighbors[5];
          velocity_u[ii].neighbors[5] = ptv[idx2(pcell, 1, 6)];
        }
      if (no_neighbor_u[5] && velocity_u[ii].cell_numbers[1] != -1)
        if (pressure[velocity_u[ii].cell_numbers[1]].neighbors[5] != -1) {
          // z+ neighbor exists through pressure on the right
          no_neighbor_u[5] = 0;
          int pcell = pressure[velocity_u[ii].cell_numbers[1]].neighbors[5];
          velocity_u[ii].neighbors[5] = ptv[idx2(pcell, 0, 6)];
        }
      if (no_neighbor_u[5]) velocity_u[ii].neighbors[5] = -1;

    }

#pragma omp for schedule(dynamic) nowait // v loop
    for (int ii = 0; ii < velocity_v.size(); ii++) {
      int no_neighbor_v[6] = { 1, 1, 1, 1, 1, 1 };
      //-- neighbor 0 --//
      if (velocity_v[ii].cell_numbers[0] != -1) {
        no_neighbor_v[0] = 0;
        velocity_v[ii].neighbors[0] = ptv[idx2(velocity_v[ii].cell_numbers[0], 2, 6)];
      }
      if (no_neighbor_v[0]) velocity_v[ii].neighbors[0] = -1;

      //-- neighbor 1 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[1] != -1) {
          // x+ neighbor exists through pressure in y- direction
          no_neighbor_v[1] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[1];
          velocity_v[ii].neighbors[1] = ptv[idx2(pcell, 3, 6)];
        }
      if (no_neighbor_v[1] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[1] != -1) {
          // x+ neighbor exists through pressure in y+ direction
          no_neighbor_v[1] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[1];
          velocity_v[ii].neighbors[1] = ptv[idx2(pcell, 2, 6)];
        }
      if (no_neighbor_v[1]) velocity_v[ii].neighbors[1] = -1;

      //-- neighbor 2 --//
      if (velocity_v[ii].cell_numbers[1] != -1) {
        no_neighbor_v[2] = 0;
        velocity_v[ii].neighbors[2] = ptv[idx2(velocity_v[ii].cell_numbers[1], 3, 6)];
      }
      if (no_neighbor_v[2]) velocity_v[ii].neighbors[2] = -1;

      //-- neighbor 3 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[3] != -1) {
          // x- neighbor exists through pressure in y- direction
          no_neighbor_v[3] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[3];
          velocity_v[ii].neighbors[3] = ptv[idx2(pcell, 3, 6)];
        }
      if (no_neighbor_v[3] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[3] != -1) {
          // x- neighbor exists through pressure in y+ direction
          no_neighbor_v[3] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[3];
          velocity_v[ii].neighbors[3] = ptv[idx2(pcell, 2, 6)];
        }
      if (no_neighbor_v[3]) velocity_v[ii].neighbors[3] = -1;

      //-- neighbor 4 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[4] != -1) {
          // z- neighbor exists through pressure in y- direction
          no_neighbor_v[4] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[4];
          velocity_v[ii].neighbors[4] = ptv[idx2(pcell, 3, 6)];
        }
      if (no_neighbor_v[4] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[4] != -1) {
          // z- neighbor exists through pressure in y+ direction
          no_neighbor_v[4] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[4];
          velocity_v[ii].neighbors[4] = ptv[idx2(pcell, 2, 6)];
        }
      if (no_neighbor_v[4]) velocity_v[ii].neighbors[4] = -1;

      //-- neighbor 5 --//
      if (velocity_v[ii].cell_numbers[0] != -1)
        if (pressure[velocity_v[ii].cell_numbers[0]].neighbors[5] != -1) {
          // z+ neighbor exists through pressure in y- direction
          no_neighbor_v[4] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[0]].neighbors[5];
          velocity_v[ii].neighbors[5] = ptv[idx2(pcell, 3, 6)];
        }
      if (no_neighbor_v[5] && velocity_v[ii].cell_numbers[1] != -1)
        if (pressure[velocity_v[ii].cell_numbers[1]].neighbors[5] != -1) {
          // z+ neighbor exists through pressure in y+ direction
          no_neighbor_v[5] = 0;
          int pcell = pressure[velocity_v[ii].cell_numbers[1]].neighbors[5];
          velocity_v[ii].neighbors[5] = ptv[idx2(pcell, 2, 6)];
        }
      if (no_neighbor_v[5]) velocity_v[ii].neighbors[5] = -1;

    }

#pragma omp for schedule(dynamic) nowait // w loop
    for (int ii = 0; ii < velocity_w.size(); ii++) {
      int no_neighbor_w[6] = { 1, 1, 1, 1, 1, 1 };
      //-- neighbor 0 (y- direction) --//
      if (velocity_w[ii].cell_numbers[0] != -1)
        if (pressure[velocity_w[ii].cell_numbers[0]].neighbors[0] != -1) {
          no_neighbor_w[0] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[0]].neighbors[0];
          velocity_w[ii].neighbors[0] = ptv[idx2(pcell, 5, 6)];
        }
      if (no_neighbor_w[0] && velocity_w[ii].cell_numbers[1] != -1)
        if (pressure[velocity_w[ii].cell_numbers[1]].neighbors[0] != -1) {
          no_neighbor_w[0] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[1]].neighbors[0];
          velocity_w[ii].neighbors[0] = ptv[idx2(pcell, 4, 6)];
        }
      if (no_neighbor_w[0]) velocity_w[ii].neighbors[0] = -1;

      //-- neighbor 1 (x+ direction) --//
      if (velocity_w[ii].cell_numbers[0] != -1)
        if (pressure[velocity_w[ii].cell_numbers[0]].neighbors[1] != -1) {
          no_neighbor_w[1] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[0]].neighbors[1];
          velocity_w[ii].neighbors[1] = ptv[idx2(pcell, 5, 6)];
        }
      if (no_neighbor_w[1] && velocity_w[ii].cell_numbers[1] != -1)
        if (pressure[velocity_w[ii].cell_numbers[1]].neighbors[1] != -1) {
          no_neighbor_w[1] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[1]].neighbors[1];
          velocity_w[ii].neighbors[1] = ptv[idx2(pcell, 4, 6)];
        }
      if (no_neighbor_w[1]) velocity_w[ii].neighbors[1] = -1;

      //-- neighbor 2 (y+ direction) --//
      if (velocity_w[ii].cell_numbers[0] != -1)
        if (pressure[velocity_w[ii].cell_numbers[0]].neighbors[2] != -1) {
          no_neighbor_w[2] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[0]].neighbors[2];
          velocity_w[ii].neighbors[2] = ptv[idx2(pcell, 5, 6)];
        }
      if (no_neighbor_w[2] && velocity_w[ii].cell_numbers[1] != -1)
        if (pressure[velocity_w[ii].cell_numbers[1]].neighbors[2] != -1) {
          no_neighbor_w[2] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[1]].neighbors[2];
          velocity_w[ii].neighbors[2] = ptv[idx2(pcell, 4, 6)];
        }
      if (no_neighbor_w[2]) velocity_w[ii].neighbors[2] = -1;

      //-- neighbor 3 (x- direction) --//
      if (velocity_w[ii].cell_numbers[0] != -1)
        if (pressure[velocity_w[ii].cell_numbers[0]].neighbors[3] != -1) {
          no_neighbor_w[3] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[0]].neighbors[3];
          velocity_w[ii].neighbors[3] = ptv[idx2(pcell, 4, 6)];
        }
      if (no_neighbor_w[3] && velocity_w[ii].cell_numbers[1] != -1)
        if (pressure[velocity_w[ii].cell_numbers[1]].neighbors[3] != -1) {
          no_neighbor_w[3] = 0;
          int pcell = pressure[velocity_w[ii].cell_numbers[1]].neighbors[3];
          velocity_w[ii].neighbors[3] = ptv[idx2(pcell, 5, 6)];
        }
      if (no_neighbor_w[3]) velocity_w[ii].neighbors[3] = -1;

      //-- neighbor 4 (z- direction) --//
      if (velocity_w[ii].cell_numbers[0] != -1) {
        no_neighbor_w[4] = 0;
        velocity_w[ii].neighbors[4] = ptv[idx2(velocity_w[ii].cell_numbers[0], 4, 6)];
      }
      if (no_neighbor_w[4]) velocity_w[ii].neighbors[4] = -1;

      //-- neighbor 5 (z+ direction) --//
      if (velocity_w[ii].cell_numbers[1] != -1) {
        no_neighbor_w[5] = 0;
        velocity_w[ii].neighbors[5] = ptv[idx2(velocity_w[ii].cell_numbers[1], 5, 6)];
      }
      if (no_neighbor_w[5]) velocity_w[ii].neighbors[5] = -1;

    }
  }
}
