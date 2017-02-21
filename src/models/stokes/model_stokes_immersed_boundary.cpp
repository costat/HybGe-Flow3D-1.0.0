#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::stokes::immersed_boundary applies the immersed boundary given in the input geometry file to the Stokes linear system.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] eta - penalization parameter for the immersed boundary. 
 */
void
hgf::models::stokes::immersed_boundary(const parameters& par, double eta)
{

  array_coo temp_coo;
  temp_coo.value = (double)1 / eta;
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

/** \brief hgf::models::stokes::random_immersed_boundary generates random immersed boundary cells until vol_frac% of the domain is covered, and applies the associated immersed boundary to the Stokes linear system.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] eta - penalization parameter for the immersed boundary. 
 * @param[in] vol_frac - between 0 and 1, gives the percentage of the void space to be converted to immersed boundary cells.
 */
void
hgf::models::stokes::random_immersed_boundary(const parameters& par, double eta, double vol_frac)
{
  std::vector< unsigned long > temp_list;
  int n_flow = (int)pressure.size();

  int n_switch = (int)((vol_frac / 100) * n_flow);
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

/** \brief hgf::models::stokes::random_immersed_boundary_clump generates random immersed boundary cells with a given affinity to clump together, until vol_frac% of the domain is covered, and applies the associated immersed boundary to the Stokes linear system.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] eta - penalization parameter for the immersed boundary. 
 * @param[in] vol_frac - between 0 and 1, gives the percentage of the void space to be converted to immersed boundary cells.
 * @param[in] likelihood - gives the likelihood that a newly generated immersed boundary cell must be connected to an existing clump of immersed boundary cells.
 */
int
hgf::models::stokes::random_immersed_boundary_clump(const parameters& par, double eta, double vol_frac, double likelihood)
{

  std::vector< unsigned long > temp_list;
  int n_flow = (int)pressure.size();
  int FAIL_MAX = n_flow;
  int nfails = 0; // track failed attempts to plant ib cells. exits with failure after FAIL_MAX faklures

  for (int i = 0; i < n_flow; i++) if (pressure_ib_list[i]) temp_list.push_back(i);

  int n_switch = (int)((vol_frac / 100) * n_flow);
  // nothing to do if domain is too small to fill vol_frac% of flow cells with IB
  if (n_switch <= 0) return 1;

  int seed_cell, seed_nbr, nnbr;
  int n_switched = 0;

  // determine new ib cells, no clump requirement
new_ib_cell:
  seed_cell = rand() % n_flow;
  if (pressure_ib_list[seed_cell]) {
    nfails++;
    goto toss_coin; // already an IB cell
  }
  n_switched++;
  pressure_ib_list[seed_cell] = 1; // total ibs tracker
  temp_list.push_back(seed_cell);  // new ibs tracker
  if (n_switched < n_switch) goto toss_coin;
  else goto set_ib;

  // find a new ib cell, require clump
ib_clump_cell:
  seed_cell = rand() % temp_list.size(); // pick a current ib cell
  seed_nbr = rand() % par.dimension*2;
  if (pressure[seed_cell].neighbors[seed_nbr] == -1) {
    goto toss_coin; // the seed cell doesn't have the appropriate neighbor
  }
  if (pressure_ib_list[pressure[seed_cell].neighbors[seed_nbr]]) {
    nfails++;
    goto toss_coin; // already an IB cell
  }
  n_switched++;
  pressure_ib_list[pressure[seed_cell].neighbors[seed_nbr]] = 1;
  if (n_switched < n_switch) goto toss_coin;
  else goto set_ib;

  // decide if new ib cell will require clumping
toss_coin:
  if (nfails > FAIL_MAX) return 1;
  if (((double) rand() / (RAND_MAX)) > likelihood) goto new_ib_cell;
  else goto ib_clump_cell;

set_ib:
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

  return 0;
}

/** \brief hgf::models::stokes::import_immersed_boundary applies an immersed boundary to the Stokes linear system according to an input vector identifying immersed boundary cells.
 *
 * @param[in] par - parameters struct containing problem information.
 * @param[in] input_ib - if input_ib[i] == 1, then pressure[i] and the associated staggered velocity cells are immersed boundary cells.
 * @param[in] eta - penalization parameter for the immersed boundary. 
 */
void
hgf::models::stokes::import_immersed_boundary(parameters& par, std::vector< int >& input_ib, double eta)
{
  pressure_ib_list = input_ib;
  array_coo temp_coo;
  temp_coo.value = (double)1 / eta;
  int dim_mult = par.dimension * 2;
  int n_u = std::accumulate(interior_u.begin(), interior_u.end(), 0);
  int n_v = std::accumulate(interior_v.begin(), interior_v.end(), 0);

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

