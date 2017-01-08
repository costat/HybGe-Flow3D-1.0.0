#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::immersed_boundary(const parameters& par, double eta)
{

  array_coo temp_coo;
  temp_coo.value = (double)1 / eta;
  int ib;
  int dim_mult = 2 * par.dimension;
  int n_u = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int n_v = std::accumulate(interior_v.begin(), interior_v.end(), 0);

  std::vector< int > pressure_ib_list;
  int void_node = -1;
  for (int ii = 0; ii < (int)par.voxel_geometry.size(); ii++) {
    if (par.voxel_geometry[ii] != 1) void_node++;
    if (par.voxel_geometry[ii] == 2) pressure_ib_list.push_back(void_node);
  }

  for (int ii = 0; ii < (int)pressure_ib_list.size(); ii++) {
    ib = pressure_ib_list[ii];
    // u component
    if (interior_u[ptv[idx2(ib, 0, dim_mult)]]) {
      temp_coo.i_index = interior_u_nums[ptv[idx2(ib, 0, dim_mult)]];
      temp_coo.j_index = interior_u_nums[ptv[idx2(ib, 0, dim_mult)]];
      coo_array.push_back(temp_coo);
    }
    if (interior_u[ptv[idx2(ib, 1, dim_mult)]]) {
      temp_coo.i_index = interior_u_nums[ptv[idx2(ib, 1, dim_mult)]];
      temp_coo.j_index = interior_u_nums[ptv[idx2(ib, 1, dim_mult)]];
      coo_array.push_back(temp_coo);
    }

    // v component
    if (interior_v[ptv[idx2(ib, 2, dim_mult)]]) {
      temp_coo.i_index = interior_v_nums[ptv[idx2(ib, 2, dim_mult)]] + n_u;
      temp_coo.j_index = interior_v_nums[ptv[idx2(ib, 2, dim_mult)]] + n_u;
      coo_array.push_back(temp_coo);
    }
    if (interior_v[ptv[idx2(ib, 3, dim_mult)]]) {
      temp_coo.i_index = interior_v_nums[ptv[idx2(ib, 3, dim_mult)]] + n_u;
      temp_coo.j_index = interior_v_nums[ptv[idx2(ib, 3, dim_mult)]] + n_u;
      coo_array.push_back(temp_coo);
    }

    // w component
    if (par.dimension == 3) {
      if (interior_w[ptv[idx2(ii, 4, dim_mult)]]) {
        temp_coo.i_index = interior_w_nums[ptv[idx2(ii, 4, dim_mult)]] + n_u + n_v;
        temp_coo.j_index = interior_w_nums[ptv[idx2(ii, 4, dim_mult)]] + n_u + n_v;
        coo_array.push_back(temp_coo);
      }
      if (interior_w[ptv[idx2(ii, 5, dim_mult)]]) {
        temp_coo.i_index = interior_w_nums[ptv[idx2(ii, 5, dim_mult)]] + n_u + n_v;
        temp_coo.j_index = interior_w_nums[ptv[idx2(ii, 5, dim_mult)]] + n_u + n_v;
        coo_array.push_back(temp_coo);
      }
    }
  }

}
