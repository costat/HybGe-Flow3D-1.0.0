/* stokes bc 2d source */

// hgf includes
#include "model_stokes.hpp"

// simple distance formula
#define distance(x1,y1,x2,y2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2))

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::xflow_2d(const parameters& par, const hgf::mesh& msh)
{
  
  boundary.resize(velocity_u.size() + velocity_v.size());

  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int shift_rows = shift_v + nV;

  // setup threading parameters
  int NTHREADS = omp_get_max_threads();
  int block_size_u = ((int)velocity_u.size() % NTHREADS) ? (int)((velocity_u.size() / NTHREADS) + 1) : (int)(velocity_u.size() / NTHREADS);
  int block_size_v = ((int)velocity_v.size() % NTHREADS) ? (int)((velocity_v.size() / NTHREADS) + 1) : (int)(velocity_v.size() / NTHREADS);
  int block_size_p = ((int)pressure.size() % NTHREADS) ? (int)((pressure.size() / NTHREADS) + 1) : (int)(pressure.size() / NTHREADS);

  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_u_arrays, temp_v_arrays;
  temp_u_arrays.resize(NTHREADS);
  temp_v_arrays.resize(NTHREADS);
  int maxu = block_size_u;
  int maxv = block_size_v;
  for (int ii = 0; ii < NTHREADS; ii++) { temp_u_arrays[ii].reserve(maxu); temp_v_arrays[ii].reserve(maxv); }

  // xmin, xmax
  double xmin = 0.0;
  double xmax = par.length;
  double ymin = 0.0;
  double ymax = par.width;
  double eps = 1E-12;

  double maxin = 1.0;

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {

      int nbrs[4];
      array_coo temp_coo;
      double dx, dy;

      for (int ii = kk*block_size_u; ii < std::min((kk + 1)*block_size_u, (int)interior_u_nums.size()); ii++) {
        double value = 0;
        int bc_contributor[4] = { 0, 0, 0, 0 };
        int nnbr = 0;
        if (!interior_u[ii]) goto uexit;
        for (int jj = 0; jj < 4; jj++) {
          nbrs[jj] = velocity_u[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_u[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 4) goto uexit;

        dx = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ii].coords[0] - pressure[velocity_u[ii].cell_numbers[0]].coords[0])) : \
          (2 * (pressure[velocity_u[ii].cell_numbers[1]].coords[0] - velocity_u[ii].coords[0]));
          
        dy = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[0], 3, 4)]].coords[1] - velocity_u[ii].coords[1])) : \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[1], 3, 4)]].coords[1] - velocity_u[ii].coords[1]));

        // S neighbor?
        if (bc_contributor[0]) {
          // Type Dirichlet
          value += par.viscosity * dx / (0.5*dy);
        }

        // E neighbor?
        if (bc_contributor[1]) {
          // Type Dirichlet?
          if (velocity_u[ii].coords[0] + dx <= xmax - eps) {
            value += par.viscosity * dx / dy;
            boundary[nbrs[1]].type = 1;
            boundary[nbrs[1]].value = 0.0;
          }
          // Type Neumann?
          else {// nothing to do... outflow is 0 neumann 
            boundary[nbrs[1]].type = 2;
            boundary[nbrs[1]].value = 0.0;
          }
        }

        // N neighbor?
        if (bc_contributor[2]) {
          // Type Dirichlet
          value += par.viscosity * dx / (0.5*dy);
        }

        // W neighbor?
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dx / dy;
          boundary[nbrs[3]].type = 1;
          boundary[nbrs[3]].value = 0.0;
          // is it an inflow bdr?
          if (velocity_u[ii].coords[0] - dx < xmin + eps) {
            double bvalue = (velocity_u[ii].coords[1] - ymin) * (ymax - velocity_u[ii].coords[1]);
            rhs[interior_u_nums[ii]] += bvalue * par.viscosity * dy / dx;
            boundary[nbrs[3]].value += bvalue;
          }
        }

        temp_coo.i_index = interior_u_nums[ii];
        temp_coo.j_index = interior_u_nums[ii];
        temp_coo.value = value;

        temp_u_arrays[kk].push_back(temp_coo);

      uexit:;
      }
    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {

      int nbrs[4];
      array_coo temp_coo;
      double dx, dy;

      for (int ii = kk*block_size_v; ii < std::min((kk + 1)*block_size_v, (int)interior_v_nums.size()); ii++) {
        double value = 0;
        int bc_contributor[4] = { 0, 0, 0, 0 };
        int nnbr = 0;
        if (!interior_v[ii]) goto vexit;
        for (int jj = 0; jj < 4; jj++) {
          nbrs[jj] = velocity_v[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_v[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 4) goto vexit;

        dx = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[0], 1, 4)]].coords[0] - velocity_v[ii].coords[0])) : \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[1], 1, 4)]].coords[0] - velocity_v[ii].coords[0]));
        dy = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ii].coords[1] - pressure[velocity_v[ii].cell_numbers[0]].coords[1])) : \
          (2 * (pressure[velocity_v[ii].cell_numbers[1]].coords[1] - velocity_v[ii].coords[1]));


        // S neighbor?
        if (bc_contributor[0]) {
          // Type Dirichlet
          boundary[nbrs[0] + velocity_u.size()].type = 1;
          boundary[nbrs[0] + velocity_u.size()].value = 0.0;
          value += par.viscosity * dy / dx;
        }

        // E neighbor?
        if (bc_contributor[1]) {
          // Type Dirichlet?
          if (velocity_v[ii].coords[0] + 0.5*dx <= xmax - eps) value += par.viscosity * dx / dy;
          // Type Neumann?
          else;

        }

        // N neighbor?
        if (bc_contributor[2]) {
          // Type Dirichlet
          boundary[nbrs[2] + velocity_u.size()].type = 1;
          boundary[nbrs[2] + velocity_u.size()].value = 0.0;
          value += par.viscosity * dy / dx;
        }

        // W neighbor?
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dy / (0.5*dx);
        }

        temp_coo.i_index = shift_v + interior_v_nums[ii];
        temp_coo.j_index = shift_v + interior_v_nums[ii];
        temp_coo.value = value;

        temp_v_arrays[kk].push_back(temp_coo);

      vexit:;
      }
    }
#pragma omp for schedule(dynamic) nowait // continuity equation
    for (int kk = 0; kk < NTHREADS; kk++) {
      double dxy[2], uval;
      int i_index;
      for (int ii = kk*block_size_p; ii < std::min((kk + 1)*block_size_p, (int)pressure.size()); ii++) {

        dxy[0] = distance(velocity_u[ptv[idx2(ii, 0, 4)]].coords[0], velocity_u[ptv[idx2(ii, 0, 4)]].coords[1], \
          velocity_u[ptv[idx2(ii, 1, 4)]].coords[0], velocity_u[ptv[idx2(ii, 1, 4)]].coords[1]);
        dxy[1] = distance(velocity_v[ptv[idx2(ii, 2, 4)]].coords[0], velocity_v[ptv[idx2(ii, 2, 4)]].coords[1], \
          velocity_v[ptv[idx2(ii, 3, 4)]].coords[0], velocity_v[ptv[idx2(ii, 3, 4)]].coords[1]);

        // ux
        if (interior_u_nums[ptv[idx2(ii, 0, 4)]] == -1) {
          if (pressure[ii].coords[0] - 0.5*dxy[0] < xmin + eps) {
            i_index = shift_rows + ii;
            uval = (velocity_u[ptv[idx2(ii, 0, 4)]].coords[1] - ymin) * (ymax - velocity_u[ptv[idx2(ii, 0, 4)]].coords[1]);
            rhs[i_index] -= (dxy[0] * dxy[1] / dxy[0]) * uval;
          }
        }

        if (interior_u_nums[ptv[idx2(ii, 1, 4)]] == -1) {
          if (pressure[ii].coords[0] + 0.5*dxy[0] > xmax - eps) {
            i_index = shift_rows + ii;
            uval = (velocity_u[ptv[idx2(ii, 0, 4)]].coords[1] - ymin) * (ymax - velocity_u[ptv[idx2(ii, 0, 4)]].coords[1]); // ugh convert neumann to dirichlet here
            rhs[i_index] += (dxy[0] * dxy[1] / dxy[0]) * uval;
          }
        }
      }
    }
  }
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_u_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_u_arrays[ii][jj]);
    }
  }
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_v_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_v_arrays[ii][jj]);
    }
  }
}

void
hgf::models::stokes::yflow_2d(const parameters& par, const hgf::mesh& msh)
{
  boundary.resize(velocity_u.size() + velocity_v.size());

  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int shift_rows = shift_v + nV;

  // setup threading parameters
  int NTHREADS = omp_get_max_threads();
  int block_size_u = ((int)velocity_u.size() % NTHREADS) ? (int)((velocity_u.size() / NTHREADS) + 1) : (int)(velocity_u.size() / NTHREADS);
  int block_size_v = ((int)velocity_v.size() % NTHREADS) ? (int)((velocity_v.size() / NTHREADS) + 1) : (int)(velocity_v.size() / NTHREADS);
  int block_size_p = ((int)pressure.size() % NTHREADS) ? (int)((pressure.size() / NTHREADS) + 1) : (int)(pressure.size() / NTHREADS);

  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_u_arrays, temp_v_arrays;
  temp_u_arrays.resize(NTHREADS);
  temp_v_arrays.resize(NTHREADS);
  int maxu = block_size_u;
  int maxv = block_size_v;
  for (int ii = 0; ii < NTHREADS; ii++) { temp_u_arrays[ii].reserve(maxu); temp_v_arrays[ii].reserve(maxv); }

  // xmin, xmax
  double xmin = 0.0;
  double xmax = par.length;
  double ymin = 0.0;
  double ymax = par.width;
  double eps = 1E-12;

  double maxin = 1.0;

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {

      int nbrs[4];
      array_coo temp_coo;
      double dx, dy;

      for (int ii = kk*block_size_u; ii < std::min((kk + 1)*block_size_u, (int)interior_u_nums.size()); ii++) {
        double value = 0;
        int bc_contributor[4] = { 0, 0, 0, 0 };
        int nnbr = 0;
        if (!interior_u[ii]) goto uexit;
        for (int jj = 0; jj < 4; jj++) {
          nbrs[jj] = velocity_u[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_u[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 4) goto uexit;

        dx = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ii].coords[0] - pressure[velocity_u[ii].cell_numbers[0]].coords[0])) : \
          (2 * (pressure[velocity_u[ii].cell_numbers[1]].coords[0] - velocity_u[ii].coords[0]));

        dy = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[0], 3, 4)]].coords[1] - velocity_u[ii].coords[1])) : \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[1], 3, 4)]].coords[1] - velocity_u[ii].coords[1]));

        // S neighbor?
        if (bc_contributor[0]) {
          // Type Dirichlet
          value += par.viscosity * dx / (0.5*dy);
        }

        // E neighbor?
        if (bc_contributor[1]) {
          // Type Dirichlet
          value += par.viscosity * dy / dx;
          boundary[nbrs[1]].type = 1;
          boundary[nbrs[1]].value = 0.0;
        }

        // N neighbor?
        if (bc_contributor[2]) {
          // Type Dirichlet?
          if (velocity_u[ii].coords[1] + 0.5*dy <= ymax - eps) value += par.viscosity * dx / (0.5*dy);
          // Type Neumann
          else;
        }

        // W neighbor?
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dy / dx;
          boundary[nbrs[3]].type = 1;
          boundary[nbrs[3]].value = 0.0;
        }

        temp_coo.i_index = interior_u_nums[ii];
        temp_coo.j_index = interior_u_nums[ii];
        temp_coo.value = value;

        temp_u_arrays[kk].push_back(temp_coo);

      uexit:;
      } // ii loop, u section
    } // kk loop, u section
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {

      int nbrs[4];
      array_coo temp_coo;
      double dx, dy;

      for (int ii = kk*block_size_v; ii < std::min((kk + 1)*block_size_v, (int)interior_v_nums.size()); ii++) {
        double value = 0;
        int bc_contributor[4] = { 0, 0, 0, 0 };
        int nnbr = 0;
        if (!interior_v[ii]) goto vexit;
        for (int jj = 0; jj < 4; jj++) {
          nbrs[jj] = velocity_v[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_v[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 4) goto vexit;

        dx = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[0], 1, 4)]].coords[0] - velocity_v[ii].coords[0])) : \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[1], 1, 4)]].coords[0] - velocity_v[ii].coords[0]));
        dy = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ii].coords[1] - pressure[velocity_v[ii].cell_numbers[0]].coords[1])) : \
          (2 * (pressure[velocity_v[ii].cell_numbers[1]].coords[1] - velocity_v[ii].coords[1]));


        // S neighbor?
        if (bc_contributor[0]) {
          // Type Dirichlet
          boundary[nbrs[0] + velocity_u.size()].type = 1;
          boundary[nbrs[0] + velocity_u.size()].value = 0.0;
          value += par.viscosity * dx / dy;
          // is this an inflow boundary?
          if (velocity_v[ii].coords[1] - dy < ymin + eps) {
            double bvalue = (velocity_v[ii].coords[0] - xmin) * (xmax - velocity_v[ii].coords[0]);
            rhs[interior_v_nums[ii] + shift_v] += bvalue * par.viscosity * dx / dy;
            boundary[nbrs[3] + velocity_u.size()].value += bvalue;
          }
        }

        // E neighbor?
        if (bc_contributor[1]) {
          // Type Dirichlet?
          value += par.viscosity * dy / (0.5*dx);
        }

        // N neighbor?
        if (bc_contributor[2]) {
          // Type Dirichlet?
          if (velocity_v[ii].coords[1] + dy <= ymin - eps) {
            boundary[nbrs[2] + velocity_u.size()].type = 1;
            boundary[nbrs[2] + velocity_u.size()].value = 0.0;
            value += par.viscosity * dx / dy;
          }
          else {// nothing to do... outflow is 0 neumann 
            boundary[nbrs[2] + velocity_u.size()].type = 2;
            boundary[nbrs[2] + velocity_u.size()].value = 0.0;
          }
        }

        // W neighbor?
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dy / (0.5*dx);
        }

        temp_coo.i_index = shift_v + interior_v_nums[ii];
        temp_coo.j_index = shift_v + interior_v_nums[ii];
        temp_coo.value = value;

        temp_v_arrays[kk].push_back(temp_coo);

      vexit:;
      } // ii loop, v section
    } // kk loop, v section
#pragma omp for schedule(dynamic) nowait // continuity equation
    for (int kk = 0; kk < NTHREADS; kk++) {
      double dxy[2], vval;
      int i_index;
      for (int ii = kk*block_size_p; ii < std::min((kk + 1)*block_size_p, (int)pressure.size()); ii++) {

        dxy[0] = distance(velocity_u[ptv[idx2(ii, 0, 4)]].coords[0], velocity_u[ptv[idx2(ii, 0, 4)]].coords[1], \
          velocity_u[ptv[idx2(ii, 1, 4)]].coords[0], velocity_u[ptv[idx2(ii, 1, 4)]].coords[1]);
        dxy[1] = distance(velocity_v[ptv[idx2(ii, 2, 4)]].coords[0], velocity_v[ptv[idx2(ii, 2, 4)]].coords[1], \
          velocity_v[ptv[idx2(ii, 3, 4)]].coords[0], velocity_v[ptv[idx2(ii, 3, 4)]].coords[1]);

        // vy
        if (interior_v_nums[ptv[idx2(ii, 2, 4)]] == -1) {
          if (pressure[ii].coords[1] - 0.5*dxy[1] < ymin + eps) {
            i_index = shift_rows + ii;
            vval = (velocity_v[ptv[idx2(ii, 2, 4)]].coords[0] - xmin) * (xmax - velocity_v[ptv[idx2(ii, 2, 4)]].coords[0]);
            rhs[i_index] -= (dxy[0] * dxy[1] / dxy[1]) * vval;
          }
        }

        if (interior_v_nums[ptv[idx2(ii, 3, 4)]] == -1) {
          if (pressure[ii].coords[1] + 0.5*dxy[1] > ymax - eps) {
            i_index = shift_rows + ii;
            vval = (velocity_v[ptv[idx2(ii, 2, 4)]].coords[0] - xmin) * (xmax - velocity_v[ptv[idx2(ii, 2, 4)]].coords[0]); // ugh convert neumann to dirichlet here
            rhs[i_index] += (dxy[0] * dxy[1] / dxy[1]) * vval;
          }
        }
      }
    }

  } // omp parallel
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_u_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_u_arrays[ii][jj]);
    }
  }
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_v_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_v_arrays[ii][jj]);
    }
  }
}
