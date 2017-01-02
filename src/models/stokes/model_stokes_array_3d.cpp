/* stokes array 3d source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// simple distance formula
#define distance(x1,y1,z1,x2,y2,z2) \
  sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2))

void
hgf::models::stokes::build_array_3d(const parameters& par, const hgf::mesh& msh)
{

  // momentum equation entries
  momentum_3d(par.viscosity);

  // continuity equation
  continuity_3d();

}

void
hgf::models::stokes::momentum_3d(double visc)
{

  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int shift_w = std::accumulate(interior_v.begin(), interior_v.end(), shift_v);
  int shift_p = std::accumulate(interior_w.begin(), interior_w.end(), shift_w);

  // threading parameters
  int NTHREADS = omp_get_max_threads();
  int block_size_u = ((int)velocity_u.size() % NTHREADS) ? (int)((velocity_u.size() / NTHREADS) + 1) : (int)(velocity_u.size() / NTHREADS);
  int block_size_v = ((int)velocity_v.size() % NTHREADS) ? (int)((velocity_v.size() / NTHREADS) + 1) : (int)(velocity_v.size() / NTHREADS);
  int block_size_w = ((int)velocity_w.size() % NTHREADS) ? (int)((velocity_w.size() / NTHREADS) + 1) : (int)(velocity_w.size() / NTHREADS);

  // define temp coo arrays to store results in parallel region
  std::vector< std::vector< array_coo > > temp_u_arrays, temp_v_arrays, temp_w_arrays;
  temp_u_arrays.resize(NTHREADS);
  temp_v_arrays.resize(NTHREADS);
  temp_w_arrays.resize(NTHREADS);
  int maxu = block_size_u * 9;
  int maxv = block_size_v * 9;
  int maxw = block_size_w * 9;
  for (int ii = 0; ii < NTHREADS; ii++) { temp_u_arrays[ii].reserve(maxu); temp_v_arrays[ii].reserve(maxv); temp_w_arrays[ii].reserve(maxw); }

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // u blocks
      
      int entries = 0;
      array_coo temp_coo[9] = { 0 };
      double d_dofs[6], d_faces[6];
      int nbrs[6], pres[2];
      double length, height, width, volume;

      for (int ii = kk*block_size_u; ii < std::min((kk + 1)*block_size_u, (int)interior_u_nums.size()); ii++) {

        // cell center distances
        for (int jj = 0; jj < 6; jj++) nbrs[jj] = velocity_u[ii].neighbors[jj];
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(velocity_u[ii].coords[0], velocity_u[ii].coords[1], velocity_u[ii].coords[2], \
              velocity_u[nbrs[jj]].coords[0], velocity_u[nbrs[jj]].coords[1], velocity_u[nbrs[jj]].coords[2]);
          }
        }

        // face dimensions
        // currently set to uniform griding, come back and extend later
        int v[4]; // v cells surrounding the u cell
        int w[4]; // w cells surrounding the u cell
        if (velocity_u[ii].cell_numbers[0] > -1) {
          v[0] = ptv[idx2(velocity_u[ii].cell_numbers[0], 2, 6)];
          v[3] = ptv[idx2(velocity_u[ii].cell_numbers[0], 3, 6)];
          w[0] = ptv[idx2(velocity_u[ii].cell_numbers[0], 4, 6)];
          w[3] = ptv[idx2(velocity_u[ii].cell_numbers[0], 5, 6)];
        }
        else goto uexit;
        if (velocity_u[ii].cell_numbers[1] > -1) {
          v[1] = ptv[idx2(velocity_u[ii].cell_numbers[1], 2, 6)];
          v[2] = ptv[idx2(velocity_u[ii].cell_numbers[1], 3, 6)];
          w[1] = ptv[idx2(velocity_u[ii].cell_numbers[1], 4, 6)];
          w[2] = ptv[idx2(velocity_u[ii].cell_numbers[1], 5, 6)];
        }
        else goto uexit;
        // face 0
        width = distance(velocity_v[v[0]].coords[0], velocity_v[v[0]].coords[1], velocity_v[v[0]].coords[2], \
          velocity_v[v[1]].coords[0], velocity_v[v[1]].coords[1], velocity_v[v[1]].coords[2]);
        height = distance(velocity_w[w[0]].coords[0], velocity_w[w[0]].coords[1], velocity_w[w[0]].coords[2], \
            velocity_w[w[3]].coords[0], velocity_w[w[3]].coords[1], velocity_w[w[3]].coords[2]);
        length = width;
        volume = length * width * height;
        d_faces[0] = width*height;
        for (int jj = 0; jj < 5; jj++) d_faces[jj + 1] = d_faces[0];
      
        // off diagonal entries
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1 && interior_u[nbrs[jj]]) {
            entries++;
            temp_coo[entries - 1].value = -visc * d_faces[jj] / d_dofs[jj];
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
        temp_coo[entries - 2].value = -volume / length;
        temp_coo[entries - 2].i_index = interior_u_nums[ii];
        temp_coo[entries - 2].j_index = pres[0] + shift_p;
        temp_coo[entries - 1].value = volume / length;
        temp_coo[entries - 1].i_index = interior_u_nums[ii];
        temp_coo[entries - 1].j_index = pres[1] + shift_p;

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_u_arrays[kk].push_back(temp_coo[jj]);
        }

        // exit
      uexit:
        entries = 0;
      }

    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // v blocks

      int entries = 0;
      array_coo temp_coo[9] = { 0 };
      double d_dofs[6], d_faces[6];
      int nbrs[6], pres[2];
      double length, height, width, volume;

      for (int ii = kk*block_size_v; ii < std::min((kk + 1)*block_size_v, (int)interior_v_nums.size()); ii++) {

        // cell center distances
        for (int jj = 0; jj < 6; jj++) nbrs[jj] = velocity_v[ii].neighbors[jj];
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(velocity_v[ii].coords[0], velocity_v[ii].coords[1], velocity_v[ii].coords[2], \
              velocity_v[nbrs[jj]].coords[0], velocity_v[nbrs[jj]].coords[1], velocity_v[nbrs[jj]].coords[2]);
          }
        }

        // face dimensions
        // currently set to uniform griding, come back and extend later
        int u[4]; // u cells surrounding the u cell
        int w[4]; // w cells surrounding the u cell
        if (velocity_v[ii].cell_numbers[0] > -1) {
          u[0] = ptv[idx2(velocity_v[ii].cell_numbers[0], 0, 6)];
          u[3] = ptv[idx2(velocity_v[ii].cell_numbers[0], 1, 6)];
          w[0] = ptv[idx2(velocity_v[ii].cell_numbers[0], 4, 6)];
          w[3] = ptv[idx2(velocity_v[ii].cell_numbers[0], 5, 6)];
        }
        else goto vexit;
        if (velocity_v[ii].cell_numbers[1] > -1) {
          u[1] = ptv[idx2(velocity_v[ii].cell_numbers[1], 0, 6)];
          u[2] = ptv[idx2(velocity_v[ii].cell_numbers[1], 1, 6)];
          w[1] = ptv[idx2(velocity_v[ii].cell_numbers[1], 4, 6)];
          w[2] = ptv[idx2(velocity_v[ii].cell_numbers[1], 5, 6)];
        }
        else goto vexit;
        // face 0
        length = distance(velocity_u[u[0]].coords[0], velocity_u[u[0]].coords[1], velocity_u[u[0]].coords[2], \
          velocity_u[u[3]].coords[0], velocity_u[u[3]].coords[1], velocity_u[u[3]].coords[2]);
        height = distance(velocity_w[w[0]].coords[0], velocity_w[w[0]].coords[1], velocity_w[w[0]].coords[2], \
          velocity_w[w[3]].coords[0], velocity_w[w[3]].coords[1], velocity_w[w[3]].coords[2]);
        width = length;
        volume = length * width * height;
        d_faces[0] = width*height;
        for (int jj = 0; jj < 5; jj++) d_faces[jj + 1] = d_faces[0];

        // off diagonal entries
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1 && interior_v[nbrs[jj]]) {
            entries++;
            temp_coo[entries - 1].value = -visc * d_faces[jj] / d_dofs[jj];
            temp_coo[entries - 1].i_index = interior_v_nums[ii] + shift_v;
            temp_coo[entries - 1].j_index = interior_v_nums[nbrs[jj]] + shift_v;
          }
        }

        // diagonal entry
        entries++;
        temp_coo[entries - 1].value = 0;
        for (int jj = 0; jj < (entries - 1); jj++) {
          temp_coo[entries - 1].value -= temp_coo[jj].value;
        }
        temp_coo[entries - 1].i_index = interior_v_nums[ii] + shift_v;
        temp_coo[entries - 1].j_index = interior_v_nums[ii] + shift_v;

        // pressure gradient
        entries += 2;
        pres[0] = velocity_v[ii].cell_numbers[0];
        pres[1] = velocity_v[ii].cell_numbers[1];
        temp_coo[entries - 2].value = -volume / width;
        temp_coo[entries - 2].i_index = interior_v_nums[ii] + shift_v;
        temp_coo[entries - 2].j_index = pres[0] + shift_p;
        temp_coo[entries - 1].value = volume / width;
        temp_coo[entries - 1].i_index = interior_v_nums[ii] + shift_v;
        temp_coo[entries - 1].j_index = pres[1] + shift_p;

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_v_arrays[kk].push_back(temp_coo[jj]);
        }

        // exit
      vexit:
        entries = 0;
      }

    }
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) { // w blocks
      int entries = 0;
      array_coo temp_coo[9] = { 0 };
      double d_dofs[6], d_faces[6];
      int nbrs[6], pres[2];
      double length, height, width, volume;

      for (int ii = kk*block_size_w; ii < std::min((kk + 1)*block_size_w, (int)interior_w_nums.size()); ii++) {

        // cell center distances
        for (int jj = 0; jj < 6; jj++) nbrs[jj] = velocity_w[ii].neighbors[jj];
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1) {
            d_dofs[jj] = distance(velocity_w[ii].coords[0], velocity_w[ii].coords[1], velocity_w[ii].coords[2], \
              velocity_w[nbrs[jj]].coords[0], velocity_w[nbrs[jj]].coords[1], velocity_w[nbrs[jj]].coords[2]);
          }
        }

        // face dimensions
        // currently set to uniform griding, come back and extend later
        int u[4]; // u cells surrounding the w cell
        int v[4]; // v cells surrounding the w cell
        if (velocity_w[ii].cell_numbers[0] > -1) {
          u[0] = ptv[idx2(velocity_w[ii].cell_numbers[0], 0, 6)];
          u[3] = ptv[idx2(velocity_w[ii].cell_numbers[0], 1, 6)];
          v[0] = ptv[idx2(velocity_w[ii].cell_numbers[0], 2, 6)];
          v[3] = ptv[idx2(velocity_w[ii].cell_numbers[0], 3, 6)];
        }
        else goto wexit;
        if (velocity_w[ii].cell_numbers[1] > -1) {
          u[1] = ptv[idx2(velocity_w[ii].cell_numbers[1], 0, 6)];
          u[2] = ptv[idx2(velocity_w[ii].cell_numbers[1], 1, 6)];
          v[1] = ptv[idx2(velocity_w[ii].cell_numbers[1], 2, 6)];
          v[2] = ptv[idx2(velocity_w[ii].cell_numbers[1], 3, 6)];
        }
        else goto wexit;
        // face 0
        length = distance(velocity_u[u[0]].coords[0], velocity_u[u[0]].coords[1], velocity_u[u[0]].coords[2], \
          velocity_u[u[3]].coords[0], velocity_u[u[3]].coords[1], velocity_u[u[3]].coords[2]);
        width = distance(velocity_v[v[0]].coords[0], velocity_v[v[0]].coords[1], velocity_v[v[0]].coords[2], \
          velocity_v[v[3]].coords[0], velocity_v[v[3]].coords[1], velocity_v[v[3]].coords[2]);
        height = length;
        volume = length * width * height;
        d_faces[0] = width*height;
        for (int jj = 0; jj < 5; jj++) d_faces[jj + 1] = d_faces[0];

        // off diagonal entries
        for (int jj = 0; jj < 6; jj++) {
          if (nbrs[jj] > -1 && interior_w[nbrs[jj]]) {
            entries++;
            temp_coo[entries - 1].value = -visc * d_faces[jj] / d_dofs[jj];
            temp_coo[entries - 1].i_index = interior_w_nums[ii] + shift_w;
            temp_coo[entries - 1].j_index = interior_w_nums[nbrs[jj]] + shift_w;
          }
        }

        // diagonal entry
        entries++;
        temp_coo[entries - 1].value = 0;
        for (int jj = 0; jj < (entries - 1); jj++) {
          temp_coo[entries - 1].value -= temp_coo[jj].value;
        }
        temp_coo[entries - 1].i_index = interior_w_nums[ii] + shift_w;
        temp_coo[entries - 1].j_index = interior_w_nums[ii] + shift_w;

        // pressure gradient
        entries += 2;
        pres[0] = velocity_w[ii].cell_numbers[0];
        pres[1] = velocity_w[ii].cell_numbers[1];
        temp_coo[entries - 2].value = -volume / width;
        temp_coo[entries - 2].i_index = interior_w_nums[ii] + shift_w;
        temp_coo[entries - 2].j_index = pres[0] + shift_p;
        temp_coo[entries - 1].value = volume / width;
        temp_coo[entries - 1].i_index = interior_w_nums[ii] + shift_w;
        temp_coo[entries - 1].j_index = pres[1] + shift_p;

        // place values into temporoary coo array
        for (int jj = 0; jj < entries; jj++) {
          temp_w_arrays[kk].push_back(temp_coo[jj]);
        }

        // exit
      wexit:
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
  for (int ii = 0; ii < NTHREADS; ii++) {
    for (int jj = 0; jj < temp_w_arrays[ii].size(); jj++) {
      coo_array.push_back(temp_w_arrays[ii][jj]);
    }
  }

}

void
hgf::models::stokes::continuity_3d(void)
{
  int shift_v = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int shift_w = std::accumulate(interior_v.begin(), interior_v.end(), shift_v);
  int shift_rows = std::accumulate(interior_w.begin(), interior_w.end(), shift_w);

  // threading
  int NTHREADS = omp_get_max_threads();
  int block_size_p = ((int)pressure.size() % NTHREADS) ? (int)((pressure.size() / NTHREADS) + 1) : (int)(pressure.size() / NTHREADS);

  // temp coo arrays for parallel regions
  std::vector< std::vector< array_coo > > temp_arrays;
  temp_arrays.resize(NTHREADS);
  int maxp = block_size_p * 6;
  for (int ii = 0; ii < NTHREADS; ii++) temp_arrays[ii].reserve(maxp);

#pragma omp parallel
  {
#pragma omp for schedule(dynamic) nowait
    for (int kk = 0; kk < NTHREADS; kk++) {
      double dxyz[3];
      array_coo temp_array[6];
      for (int ii = kk*block_size_p; ii < std::min((kk + 1)*block_size_p, (int)pressure.size()); ii++) {

        dxyz[0] = distance(velocity_u[ptv[idx2(ii, 0, 6)]].coords[0], velocity_u[ptv[idx2(ii, 0, 6)]].coords[1], velocity_u[ptv[idx2(ii, 0, 6)]].coords[2], \
          velocity_u[ptv[idx2(ii, 1, 6)]].coords[0], velocity_u[ptv[idx2(ii, 1, 6)]].coords[1], velocity_u[ptv[idx2(ii, 1, 6)]].coords[2]);
        dxyz[1] = distance(velocity_v[ptv[idx2(ii, 2, 6)]].coords[0], velocity_v[ptv[idx2(ii, 2, 6)]].coords[1], velocity_v[ptv[idx2(ii, 2, 6)]].coords[2], \
          velocity_v[ptv[idx2(ii, 3, 6)]].coords[0], velocity_v[ptv[idx2(ii, 3, 6)]].coords[1], velocity_v[ptv[idx2(ii, 3, 6)]].coords[2]);
        dxyz[2] = distance(velocity_w[ptv[idx2(ii, 4, 6)]].coords[0], velocity_w[ptv[idx2(ii, 4, 6)]].coords[1], velocity_w[ptv[idx2(ii, 4, 6)]].coords[2], \
          velocity_w[ptv[idx2(ii, 5, 6)]].coords[0], velocity_w[ptv[idx2(ii, 5, 6)]].coords[1], velocity_w[ptv[idx2(ii, 5, 6)]].coords[2]);

        // ux
        temp_array[0].i_index = shift_rows + ii;
        temp_array[0].j_index = interior_u_nums[ptv[idx2(ii, 0, 6)]];
        temp_array[0].value = dxyz[0] * dxyz[1] * dxyz[2] / dxyz[0];

        temp_array[1].i_index = shift_rows + ii;
        temp_array[1].j_index = interior_u_nums[ptv[idx2(ii, 1, 6)]];
        temp_array[1].value = -dxyz[0] * dxyz[1] * dxyz[2] / dxyz[0];

        // vy
        temp_array[2].i_index = shift_rows + ii;
        temp_array[2].j_index = (interior_v_nums[ptv[idx2(ii, 2, 6)]] != -1) ? (shift_v + interior_v_nums[ptv[idx2(ii, 2, 6)]]) : interior_v_nums[ptv[idx2(ii, 2, 6)]];
        temp_array[2].value = dxyz[0] * dxyz[1] * dxyz[2] / dxyz[1];

        temp_array[3].i_index = shift_rows + ii;
        temp_array[3].j_index = (interior_v_nums[ptv[idx2(ii, 3, 6)]] != -1) ? (shift_v + interior_v_nums[ptv[idx2(ii, 3, 6)]]) : interior_v_nums[ptv[idx2(ii, 3, 6)]];
        temp_array[3].value = -dxyz[0] * dxyz[1] * dxyz[2] / dxyz[1];

        // wz
        temp_array[4].i_index = shift_rows + ii;
        temp_array[4].j_index = (interior_w_nums[ptv[idx2(ii, 4, 6)]] != -1) ? (shift_w + interior_w_nums[ptv[idx2(ii, 4, 6)]]) : interior_w_nums[ptv[idx2(ii, 4, 6)]];
        temp_array[4].value = dxyz[0] * dxyz[1] * dxyz[2] / dxyz[2];

        temp_array[5].i_index = shift_rows + ii;
        temp_array[5].j_index = (interior_w_nums[ptv[idx2(ii, 5, 6)]] != -1) ? (shift_w + interior_w_nums[ptv[idx2(ii, 5, 6)]]) : interior_w_nums[ptv[idx2(ii, 5, 6)]];
        temp_array[5].value = -dxyz[0] * dxyz[1] * dxyz[2] / dxyz[2];

        for (int jj = 0; jj < 6; jj++) {
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
