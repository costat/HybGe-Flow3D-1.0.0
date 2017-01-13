#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

// includes IB from geo file in array
void
hgf::models::stokes::immersed_boundary(const parameters& par, double eta)
{

  array_coo temp_coo;
  temp_coo.value = (double)1 / eta;
  int ib;
  int dim_mult = 2 * par.dimension;
  int n_u = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int n_v = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  pressure_ib_list.assign(pressure.size(), 0);

  int void_node = -1;
  for (int ii = 0; ii < (int)par.voxel_geometry.size(); ii++) {
    if (par.voxel_geometry[ii] != 1) void_node++;
    if (par.voxel_geometry[ii] == 2) pressure_ib_list[void_node] = 1;
  }

  for (int ib = 0; ib < (int)pressure.size(); ib++) {
    if (pressure_ib_list[ib]) {
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
        if (interior_w[ptv[idx2(ib, 4, dim_mult)]]) {
          temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
          temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
          coo_array.push_back(temp_coo);
        }
        if (interior_w[ptv[idx2(ib, 5, dim_mult)]]) {
          temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
          temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
          coo_array.push_back(temp_coo);
        }
      }
    }
  }

}

// generates random IB cells up to vol_frac% of the flow domain and places in array
void 
hgf::models::stokes::random_immersed_boundary(const parameters& par, double eta, double vol_frac)
{
  std::vector< unsigned long > temp_list;
  int n_flow = pressure.size();

  int n_switch = (vol_frac / 100) * n_flow;
  // nothing to do if domain is too small to fill vol_frac% of flow cells with IB
  if (n_switch <= 0) return;

  int seed_cell;
  int n_switched = 0;

  // determine new ib cells
new_ib_cell:
  seed_cell = rand() % n_flow;
  if (pressure_ib_list[seed_cell]) goto new_ib_cell;
  n_switched++;
  pressure_ib_list[seed_cell] = 1; // total ibs tracker
  temp_list.push_back(seed_cell);        // new ibs tracker
  if (n_switched < n_switch) goto new_ib_cell;

  // add ibs to array
  array_coo temp_coo;
  temp_coo.value = (double)1 / eta;
  int ib;
  int dim_mult = 2 * par.dimension;
  int n_u = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int n_v = std::accumulate(interior_v.begin(), interior_v.end(), 0);
  for (int ii = 0; ii < (int)temp_list.size(); ii++) {
    ib = temp_list[ii];
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
      if (interior_w[ptv[idx2(ib, 4, dim_mult)]]) {
        temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
        temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
        coo_array.push_back(temp_coo);
      }
      if (interior_w[ptv[idx2(ib, 5, dim_mult)]]) {
        temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
        temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
        coo_array.push_back(temp_coo);
      }
    }
  }

}

// imports ib to array from an input vector of 1s and 0s. size of input must match size of pressure array
void
hgf::models::stokes::import_immersed_boundary(parameters& par, std::vector< int >& input_ib)
{
  pressure_ib_list = input_ib;
  array_coo temp_coo;
  int dim_mult = par.dimension * 2;
  for (int ib = 0; ib < (int)pressure.size(); ib++) {
    if (pressure_ib_list[ib]) {
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
        if (interior_w[ptv[idx2(ib, 4, dim_mult)]]) {
          temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
          temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 4, dim_mult)]] + n_u + n_v;
          coo_array.push_back(temp_coo);
        }
        if (interior_w[ptv[idx2(ib, 5, dim_mult)]]) {
          temp_coo.i_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
          temp_coo.j_index = interior_w_nums[ptv[idx2(ib, 5, dim_mult)]] + n_u + n_v;
          coo_array.push_back(temp_coo);
        }
      }
    }
  }
}

