/* stokes bc 2d source */

// hgf includes
#include "model_stokes.hpp"

// simple distance formula
#define distance(x1,y1,z1,x2,y2,z2) sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2))

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::xflow_3d(const parameters& par, const hgf::mesh& msh)
{
  
  boundary.resize(velocity_u.size() + velocity_v.size() + velocity_w.size());

  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int shift_w = std::accumulate(interior_v.begin(), interior_v.end(), shift_v);
  int nW = std::accumulate(interior_w.begin(), interior_w.end(), 0);
  
  int shift_rows = shift_w + nW;

  // setup threading parameters
  int NTHREADS = omp_get_max_threads();
  int block_size_u = ((int)velocity_u.size() % NTHREADS) ? (int)((velocity_u.size() / NTHREADS) + 1) : (int)(velocity_u.size() / NTHREADS);
  int block_size_v = ((int)velocity_v.size() % NTHREADS) ? (int)((velocity_v.size() / NTHREADS) + 1) : (int)(velocity_v.size() / NTHREADS);
  int block_size_w = ((int)velocity_w.size() % NTHREADS) ? (int)((velocity_w.size() / NTHREADS) + 1) : (int)(velocity_w.size() / NTHREADS);
  int block_size_p = ((int)pressure.size() % NTHREADS) ? (int)((pressure.size() / NTHREADS) + 1) : (int)(pressure.size() / NTHREADS);

  // define temp coo arays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_u_arrays, temp_v_arrays, temp_w_arrays;
  temp_u_arrays.resize(NTHREADS);
  temp_v_arrays.resize(NTHREADS);
  temp_w_arrays.resize(NTHREADS);

  int maxu = block_size_u;
  int maxv = block_size_v;
  int maxw = block_size_w;
  for (int ii = 0; ii < NTHREADS; ii++) { 
    temp_u_arrays[ii].reserve(maxu); 
    temp_v_arrays[ii].reserve(maxv); 
    temp_w_arrays[ii].reserve(maxw); 
  }

  // geometry limits
  double xmin = 0.0;
  double xmax = par.length;
  double ymin = 0.0;
  double ymax = par.width;
  double zmin = 0.0;
  double zmax = par.height;
  double eps = 1E-10;

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // u section
      
      int nbrs[6];
      array_coo temp_coo;
      double dx, dy, dz;

      for (int ii = kk*block_size_u; ii < std::min((kk + 1)*block_size_u, (int)interior_u_nums.size()); ii++) {
        double value = 0;
        if (!interior_u[ii]) goto uexit;
        int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
        int nnbr = 0;
        for (int jj = 0; jj < 6; jj++) {
          nbrs[jj] = velocity_u[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_u[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 6) goto uexit;

        dx = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ii].coords[0] - pressure[velocity_u[ii].cell_numbers[0]].coords[0])) :
          (2 * (pressure[velocity_u[ii].cell_numbers[1]].coords[0] - velocity_u[ii].coords[0]));

        dy = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[0], 3, 6)]].coords[1] - velocity_u[ii].coords[1])) : \
          (2 * (velocity_v[ptv[idx2(velocity_u[ii].cell_numbers[1], 3, 6)]].coords[1] - velocity_u[ii].coords[1]));

        dz = (velocity_u[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_w[ptv[idx2(velocity_u[ii].cell_numbers[0], 5, 6)]].coords[2] - velocity_u[ii].coords[2])) : \
          (2 * (velocity_w[ptv[idx2(velocity_u[ii].cell_numbers[1], 5, 6)]].coords[2] - velocity_u[ii].coords[2]));    

        // y- neighbor?
        if (bc_contributor[0]) {
          // Type Dirichlet
          value += par.viscosity * dx * dz / (0.5 * dy);
        }

        // x+ neighbor?
        if (bc_contributor[1]) {
          // Type Dirichlet?
          if (velocity_u[ii].coords[0] + dx <= xmax - eps) {
            value += par.viscosity * dx * dy / dx;
            boundary[nbrs[1]].type = 1;
            boundary[nbrs[1]].value = 0.0;
          }
          // Type Neumann?
          else {
            boundary[nbrs[1]].type = 2;
            boundary[nbrs[1]].value = 0.0;
          }
        }

        // y+ neighbor?
        if (bc_contributor[2]) {
          // Type Dirichlet
          value += par.viscosity * dx * dz / (0.5*dy);
        }

        // x- neighbor?
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dy * dz / dx;
          boundary[nbrs[3]].type = 1;
          boundary[nbrs[3]].value = 0.0;
          // is it an inflow boundary?
          if (velocity_u[ii].coords[0] - dx < xmin + eps) {
            rhs[interior_u_nums[ii]] += -1.0;
            boundary[nbrs[3]].value += -1.0;
          }
        }

        // z- neighbor?
        if (bc_contributor[4]) {
          // Type Dirichlet
          value += par.viscosity * dx * dy / (0.5 * dz);
        }

        // z+ neighbor?
        if (bc_contributor[5]) {
          // Type Dirichlet
          value += par.viscosity * dx * dy / (0.5 * dz);
        }

        temp_coo.i_index = interior_u_nums[ii];
        temp_coo.j_index = interior_u_nums[ii];
        temp_coo.value = value;

        temp_u_arrays[kk].push_back(temp_coo);

      uexit:;
      }

    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // v section

      int nbrs[6];
      array_coo temp_coo;
      double dx, dy, dz;

      for (int ii = kk*block_size_v; ii < std::min((kk + 1)*block_size_v, (int)interior_v_nums.size()); ii++) {
        double value = 0;
        if (!interior_v[ii]) goto vexit;
        int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
        int nnbr = 0;
        for (int jj = 0; jj < 6; jj++) {
          nbrs[jj] = velocity_v[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_v[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 6) goto vexit;

        dy = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ii].coords[1] - pressure[velocity_v[ii].cell_numbers[0]].coords[1])) :
          (2 * (pressure[velocity_v[ii].cell_numbers[1]].coords[1] - velocity_v[ii].coords[1]));

        dx = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[0], 1, 6)]].coords[1] - velocity_v[ii].coords[1])) : \
          (2 * (velocity_u[ptv[idx2(velocity_v[ii].cell_numbers[1], 1, 6)]].coords[1] - velocity_v[ii].coords[1]));

        dz = (velocity_v[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_w[ptv[idx2(velocity_v[ii].cell_numbers[0], 5, 6)]].coords[2] - velocity_v[ii].coords[2])) : \
          (2 * (velocity_w[ptv[idx2(velocity_v[ii].cell_numbers[1], 5, 6)]].coords[2] - velocity_v[ii].coords[2]));

        // y- neighbor
        if (bc_contributor[0]) {
          // Type Dirichlet
          boundary[nbrs[0] + shift_v].type = 1;
          boundary[nbrs[0] + shift_v].value = 0.0;
          value += par.viscosity * dx * dz / dy;
        }

        // x+ neighbor
        if (bc_contributor[1]) {
          // Type Dirichlet?
          if (velocity_v[ii].coords[0] + 0.5*dx <= xmax - eps) value += par.viscosity * dy * dz / (0.5 * dx);
          // Type Neumann;
          else;
        }

        // y+ neighbor
        if (bc_contributor[2]) {
          // Type Dirichlet
          boundary[nbrs[2] + shift_v].type = 1;
          boundary[nbrs[2] + shift_v].value = 0.0;
          value += par.viscosity * dx * dz / dy;
        }

        // x- neighbor
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dz * dy / (0.5 * dx);
        }

        // z- neighbor
        if (bc_contributor[4]) {
          // Type Dirichlet
          value += par.viscosity * dy * dz / (0.5 * dz);
        }

        // z+ neighbor
        if (bc_contributor[5]) {
          // Type Dirichlet
          value += par.viscosity * dy * dz / (0.5 * dz);
        }

        temp_coo.i_index = shift_v + interior_v_nums[ii];
        temp_coo.j_index = shift_v + interior_v_nums[ii];
        temp_coo.value = value;

        temp_v_arrays[kk].push_back(temp_coo);

      vexit:;
      }

    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // w section

      int nbrs[6];
      array_coo temp_coo;
      double dx, dy, dz;

      for (int ii = kk*block_size_w; ii < std::min((kk + 1)*block_size_w, (int)interior_w_nums.size()); ii++) {

        double value = 0;
        if (!interior_w[ii]) goto wexit;
        int bc_contributor[6] = { 0, 0, 0, 0, 0, 0 };
        int nnbr = 0;
        for (int jj = 0; jj < 6; jj++) {
          nbrs[jj] = velocity_w[ii].neighbors[jj];
          if (nbrs[jj] != -1 && interior_w[nbrs[jj]]) nnbr++;
          else bc_contributor[jj] = 1;
        }
        if (nnbr == 6) goto wexit;

        dz = (velocity_w[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_w[ii].coords[2] - pressure[velocity_w[ii].cell_numbers[0]].coords[2])) :
          (2 * (pressure[velocity_w[ii].cell_numbers[1]].coords[2] - velocity_w[ii].coords[2]));

        dx = (velocity_w[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_u[ptv[idx2(velocity_w[ii].cell_numbers[0], 1, 6)]].coords[0] - velocity_u[ptv[idx2(velocity_w[ii].cell_numbers[0], 0, 6)]].coords[0])) : \
          (2 * (velocity_u[ptv[idx2(velocity_w[ii].cell_numbers[1], 1, 6)]].coords[0] - velocity_u[ptv[idx2(velocity_w[ii].cell_numbers[1], 0, 6)]].coords[0]));

        dy = (velocity_w[ii].cell_numbers[0] != -1) ? \
          (2 * (velocity_v[ptv[idx2(velocity_w[ii].cell_numbers[0], 3, 6)]].coords[1] - velocity_v[ptv[idx2(velocity_w[ii].cell_numbers[0], 2, 6)]].coords[1])) : \
          (2 * (velocity_v[ptv[idx2(velocity_w[ii].cell_numbers[1], 3, 6)]].coords[1] - velocity_v[ptv[idx2(velocity_w[ii].cell_numbers[1], 2, 6)]].coords[1]));

        // y- neighbor
        if (bc_contributor[0]) {
          // Type Dirichlet
          value += par.viscosity * dx * dz / (0.5 * dy);
        }

        // x+ neighbor
        if (bc_contributor[1]) {
          // Type Dirichlet?
          if (velocity_w[ii].coords[0] + 0.5*dx <= xmax - eps) value += par.viscosity * dy * dz / (0.5 * dx);
          // Type Neumann;
        }

        // y+ neighbor
        if (bc_contributor[2]) {
          // Type Dirichlet
          value += par.viscosity * dx * dz / (0.5 * dy);
        }

        // x- neighbor
        if (bc_contributor[3]) {
          // Type Dirichlet
          value += par.viscosity * dy * dz / (0.5 * dx);
        }

        // z- neighbor
        if (bc_contributor[4]) {
          // Type Dirichlet
          boundary[nbrs[0] + shift_w].type = 1;
          boundary[nbrs[0] + shift_w].value = 0.0;
          value += par.viscosity * dx * dy / dz;
        }

        // z+ neighbor
        if (bc_contributor[5]) {
          // Type Dirichlet
          boundary[nbrs[0] + shift_w].type = 1;
          boundary[nbrs[0] + shift_w].value = 0.0;
          value += par.viscosity * dx * dy / dz;
        }

        temp_coo.i_index = shift_w + interior_w_nums[ii];
        temp_coo.j_index = shift_w + interior_w_nums[ii];
        temp_coo.value = value;

        temp_w_arrays[kk].push_back(temp_coo);

      wexit:;
      }

    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // continuity equation

      double dxyz[3], uval;
      int i_index;

      for (int ii = kk*block_size_p; ii < std::min((kk + 1)*block_size_p, (int)pressure.size()); ii++) {

        dxyz[0] = distance(velocity_u[ptv[idx2(ii, 0, 6)]].coords[0], velocity_u[ptv[idx2(ii, 0, 6)]].coords[1], velocity_u[ptv[idx2(ii, 0, 6)]].coords[2], \
                           velocity_u[ptv[idx2(ii, 1, 6)]].coords[0], velocity_u[ptv[idx2(ii, 1, 6)]].coords[1], velocity_u[ptv[idx2(ii, 1, 6)]].coords[2]);

        dxyz[1] = distance(velocity_v[ptv[idx2(ii, 2, 6)]].coords[0], velocity_v[ptv[idx2(ii, 2, 6)]].coords[1], velocity_v[ptv[idx2(ii, 2, 6)]].coords[2], \
                           velocity_v[ptv[idx2(ii, 3, 6)]].coords[0], velocity_v[ptv[idx2(ii, 3, 6)]].coords[1], velocity_v[ptv[idx2(ii, 3, 6)]].coords[2]);

        dxyz[2] = distance(velocity_w[ptv[idx2(ii, 4, 6)]].coords[0], velocity_w[ptv[idx2(ii, 4, 6)]].coords[1], velocity_w[ptv[idx2(ii, 4, 6)]].coords[2], \
                           velocity_w[ptv[idx2(ii, 5, 6)]].coords[0], velocity_w[ptv[idx2(ii, 5, 6)]].coords[1], velocity_w[ptv[idx2(ii, 5, 6)]].coords[2]);

        // ux
        if (interior_u_nums[ptv[idx2(ii, 0, 6)]] == -1) {
          if (pressure[ii].coords[0] - 0.5*dxyz[0] < xmin + eps) {
            i_index = shift_rows + ii;
            uval = 1.0;
            rhs[i_index] += (dxyz[0] * dxyz[1] * dxyz[2] / dxyz[0]) * uval;
          }
        }
        if (interior_u_nums[ptv[idx2(ii, 1, 6)]] == -1) {
          if (pressure[ii].coords[0] + 0.5*dxyz[0] > xmax - eps) {
            i_index = shift_rows + ii;
            uval = 1.0;
            rhs[i_index] -= (dxyz[0] * dxyz[1] * dxyz[2] / dxyz[0]) * uval;
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
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_w_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_w_arrays[ii][jj]);
    }
  }
}
