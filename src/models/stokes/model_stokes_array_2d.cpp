/* stokes array 2d source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// simple distance formula
#define distance(x1,y1,x2,y2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2))

void
hgf::models::stokes::build_array_2d(const parameters& par, const hgf::mesh& msh)
{

  // calls to set up 2nd order velocity terms
  momentum_2d(par.viscosity);

  // calls to set up the continuity equations
  continuity_2d();

}

void
hgf::models::stokes::momentum_2d(double visc)
{

  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int shift_p = shift_v + nV;

  // setup threading parameters
  int NTHREADS = omp_get_max_threads();
#ifdef _THREADS_DEBUG
  std::cout << "\nNTHREADS in DIFFUSION = " << NTHREADS << ".\n";
#endif
  int block_size_u = ((int)velocity_u.size() % NTHREADS) ? (int)((velocity_u.size() / NTHREADS) + 1) : (int)(velocity_u.size() / NTHREADS);
  int block_size_v = ((int)velocity_v.size() % NTHREADS) ? (int)((velocity_v.size() / NTHREADS) + 1) : (int)(velocity_v.size() / NTHREADS);

  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_u_arrays, temp_v_arrays;
  temp_u_arrays.resize(NTHREADS);
  temp_v_arrays.resize(NTHREADS);
  int maxu = block_size_u * 7;
  int maxv = block_size_v * 7;
  for (int ii = 0; ii < NTHREADS; ii++) { temp_u_arrays[ii].reserve(maxu); temp_v_arrays[ii].reserve(maxv); }

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {
      int entries = 0;
      array_coo temp_coo[7] = { 0 };
      double d_dofs[4], d_edges[4];
      int nbrs[4], pres[2];
      for (int ii = kk*block_size_u; ii < std::min((kk + 1)*block_size_u, (int)interior_u_nums.size()); ii++) {
        // grab neighbor numbers
        for (int jj = 0; jj < 4; jj++) { nbrs[jj] = velocity_u[ii].neighbors[jj]; }
        // compute cell center distances
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(velocity_u[ii].coords[0], velocity_u[ii].coords[1], \
              velocity_u[nbrs[jj]].coords[0], velocity_u[nbrs[jj]].coords[1]);
          }
        }
        // compute edge distances
        int v[4]; // v cells surrounding the u cell
        if (velocity_u[ii].cell_numbers[0] > -1) {
          v[0] = ptv[idx2(velocity_u[ii].cell_numbers[0], 2, 4)];
          v[3] = ptv[idx2(velocity_u[ii].cell_numbers[0], 3, 4)];
        }
        else goto uexit;

        if (velocity_u[ii].cell_numbers[1] > -1) {
          v[1] = ptv[idx2(velocity_u[ii].cell_numbers[1], 2, 4)];
          v[2] = ptv[idx2(velocity_u[ii].cell_numbers[1], 3, 4)];
        }
        else goto uexit;

        for (int jj = 0; jj < 4; jj++) {
          int nn = (jj < 3) ? (jj + 1) : 0;
          d_edges[jj] = distance(velocity_v[v[jj]].coords[0], velocity_v[v[jj]].coords[1], \
            velocity_v[v[nn]].coords[0], velocity_v[v[nn]].coords[1]);
        }

        // off diagonal entries
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1 && interior_u[nbrs[jj]]) {
            entries++;
            temp_coo[entries - 1].value = -visc * d_edges[jj] / d_dofs[jj];
            temp_coo[entries - 1].i_index = interior_u_nums[ii];
            temp_coo[entries - 1].j_index = interior_u_nums[nbrs[jj]];
          }
        }
        // diagonal entry
        entries++;
        temp_coo[entries - 1].value = 0;
        for (int jj = 0; jj < (entries - 1); jj++) {
          temp_coo[entries - 1].value -= temp_coo[jj].value;
        }
        temp_coo[entries - 1].i_index = interior_u_nums[ii];
        temp_coo[entries - 1].j_index = interior_u_nums[ii];

        // pressure gradient
        entries += 2;
        pres[0] = velocity_u[ii].cell_numbers[0];
        pres[1] = velocity_u[ii].cell_numbers[1];
        temp_coo[entries - 2].value = -d_edges[0] * d_edges[1] / d_edges[0];
        temp_coo[entries - 2].i_index = interior_u_nums[ii];
        temp_coo[entries - 2].j_index = pres[0] + shift_p;
        temp_coo[entries - 1].value = d_edges[0] * d_edges[1] / d_edges[0];
        temp_coo[entries - 1].i_index = interior_u_nums[ii];
        temp_coo[entries - 1].j_index = pres[1] + shift_p;

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_u_arrays[kk].push_back(temp_coo[jj]);
        }

        // node is a physical boundary node
      uexit:
        entries = 0;

      }
    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {
      int entries = 0;
      array_coo temp_coo[7] = { 0 };
      double d_dofs[4], d_edges[4];
      int nbrs[4], pres[2];
      for (int ii = kk*block_size_v; ii < std::min((kk + 1)*block_size_v, (int)interior_v_nums.size()); ii++) {
        // grab neighbor numbers
        for (int jj = 0; jj < 4; jj++) { nbrs[jj] = velocity_v[ii].neighbors[jj]; }
        // compute cell center distances
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(velocity_v[ii].coords[0], velocity_v[ii].coords[1], \
              velocity_v[nbrs[jj]].coords[0], velocity_v[nbrs[jj]].coords[1]);
          }
        }
        // compute edge distances
        int u[4];
        if (velocity_v[ii].cell_numbers[0] > -1) {
          u[0] = ptv[idx2(velocity_v[ii].cell_numbers[0], 0, 4)];
          u[1] = ptv[idx2(velocity_v[ii].cell_numbers[0], 1, 4)];
        }
        else goto vexit;
        if (velocity_v[ii].cell_numbers[1] > -1) {
          u[2] = ptv[idx2(velocity_v[ii].cell_numbers[1], 1, 4)];
          u[3] = ptv[idx2(velocity_v[ii].cell_numbers[1], 0, 4)];
        }
        else goto vexit;

        for (int jj = 0; jj < 4; jj++) {
          int nn = (jj < 3) ? (jj + 1) : 0;
          d_edges[jj] = distance(velocity_u[u[jj]].coords[0], velocity_u[u[jj]].coords[1], \
            velocity_u[u[nn]].coords[0], velocity_u[u[nn]].coords[1]);
        }

        // off diagonal entries
        for (int jj = 0; jj < 4; jj++) {
          if (nbrs[jj] > -1 && interior_v[nbrs[jj]]) {
            entries++;
            temp_coo[entries - 1].value = -visc * d_edges[jj] / d_dofs[jj];
            temp_coo[entries - 1].i_index = interior_v_nums[ii]+shift_v;
            temp_coo[entries - 1].j_index = interior_v_nums[nbrs[jj]]+shift_v;
          }
        }
        // diagonal entry
        entries++;
        temp_coo[entries - 1].value = 0;
        for (int jj = 0; jj < (entries - 1); jj++) {
          temp_coo[entries - 1].value -= temp_coo[jj].value;
        }
        temp_coo[entries - 1].i_index = interior_v_nums[ii]+shift_v;
        temp_coo[entries - 1].j_index = interior_v_nums[ii]+shift_v;

        // pressure gradient
        entries += 2;
        pres[0] = velocity_v[ii].cell_numbers[0];
        pres[1] = velocity_v[ii].cell_numbers[1];
        temp_coo[entries - 2].value = -d_edges[0] * d_edges[1] / d_edges[0];
        temp_coo[entries - 2].i_index = interior_v_nums[ii] + shift_v;
        temp_coo[entries - 2].j_index = pres[0] + shift_p;
        temp_coo[entries - 1].value = d_edges[0] * d_edges[1] / d_edges[0];
        temp_coo[entries - 1].i_index = interior_v_nums[ii] + shift_v;
        temp_coo[entries - 1].j_index = pres[1] + shift_p;

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_v_arrays[kk].push_back(temp_coo[jj]);
        }

        // node is a physical boundary node
      vexit:
        entries = 0;

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
hgf::models::stokes::continuity_2d(void)
{
  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  int shift_rows = shift_v + nV;

  // setup threading parameters
  int NTHREADS = omp_get_max_threads();
#ifdef _THREADS_DEBUG
  std::cout << "\nNTHREADS in DIFFUSION = " << NTHREADS << ".\n";
#endif
  int block_size_p = ((int)pressure.size() % NTHREADS) ? (int)((pressure.size() / NTHREADS) + 1) : (int)(pressure.size() / NTHREADS);

  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = block_size_p * 4;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {
      double dxy[2];
      array_coo temp_array[4];
      for (int ii = kk*block_size_p; ii < std::min((kk + 1)*block_size_p, (int)pressure.size()); ii++) {

        dxy[0] = distance(velocity_u[ptv[idx2(ii, 0, 4)]].coords[0], velocity_u[ptv[idx2(ii, 0, 4)]].coords[1], \
                   velocity_u[ptv[idx2(ii, 1, 4)]].coords[0], velocity_u[ptv[idx2(ii, 1, 4)]].coords[1]);
        dxy[1] = distance(velocity_v[ptv[idx2(ii, 2, 4)]].coords[0], velocity_v[ptv[idx2(ii, 2, 4)]].coords[1], \
                   velocity_v[ptv[idx2(ii, 3, 4)]].coords[0], velocity_v[ptv[idx2(ii, 3, 4)]].coords[1]);

        // ux
        temp_array[0].i_index = shift_rows + ii;
        temp_array[0].j_index = interior_u_nums[ptv[idx2(ii, 0, 4)]];
        temp_array[0].value = dxy[0] * dxy[1] / dxy[0];
        
        temp_array[1].i_index = shift_rows + ii;
        temp_array[1].j_index = interior_u_nums[ptv[idx2(ii, 1, 4)]];
        temp_array[1].value = -dxy[0] * dxy[1] / dxy[0];
        
        // vy
        temp_array[2].i_index = shift_rows + ii;
        temp_array[2].j_index = (interior_v_nums[ptv[idx2(ii, 2, 4)]] != -1) ? (shift_v + interior_v_nums[ptv[idx2(ii, 2, 4)]]) : interior_v_nums[ptv[idx2(ii, 2, 4)]];
        temp_array[2].value = dxy[0] * dxy[1] / dxy[1];

        temp_array[3].i_index = shift_rows + ii;
        temp_array[3].j_index = (interior_v_nums[ptv[idx2(ii, 3, 4)]] != -1) ? (shift_v + interior_v_nums[ptv[idx2(ii, 3, 4)]]) : interior_v_nums[ptv[idx2(ii, 3, 4)]];
        temp_array[3].value = -dxy[0] * dxy[1] / dxy[1];

        for (int jj = 0; jj < 4; jj++) {
          if (temp_array[jj].j_index != -1) temp_arrays[kk].push_back(temp_array[jj]);
        }
      }
    }
  }
  // paste
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_arrays[ii][jj]);
    }
  }
}
