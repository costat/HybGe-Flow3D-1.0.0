/* porenetwork boundary source */

// hgf includes
#include "model_porenetwork.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::porenetwork::setup_xflow_bc(const parameters& par)
{
  array_coo temp_coo[7];
  double pn_epsilon = 1E-12;
  int idx_swap[6] = { 2, 3, 0, 1, 5, 4 };
  int entries;
  for (int ii = 0; ii < (int)pressure.size(); ii++) {
    entries = 0;
    if (pressure[ii].doftype) {
      // inflow
      if (pressure[ii].coords[0] <= pn_epsilon) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
        rhs[ii] += 1;
      }
      // outflow
      else if (pressure[ii].coords[0] >= (par.length - pn_epsilon)) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
      }
      // other
      else {
        for (int dir = 0; dir < (par.dimension * 2); dir++) {
          if (pressure[ii].neighbors[dir] != -1) {
            entries++;
            temp_coo[entries - 1].i_index = ii;
            temp_coo[entries - 1].j_index = pressure[ii].neighbors[dir];
            temp_coo[entries - 1].value = -0.5 * (permeability[idx2(temp_coo[entries - 1].j_index, idx_swap[dir], (par.dimension * 2))] + permeability[idx2(ii, dir, (par.dimension * 2))]);
          }
        }
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 0;
        for (int dir = 0; dir < (entries-1); dir++) {
          temp_coo[entries - 1].value -= temp_coo[dir].value;
        }
      }
      // place in array
      for (int jj = 0; jj < entries; jj++) {
        coo_array.push_back(temp_coo[jj]);
      }

    }
  }
}

void
hgf::models::porenetwork::setup_yflow_bc(const parameters& par)
{
  array_coo temp_coo[7];
  double pn_epsilon = 1E-12;
  int idx_swap[6] = { 2, 3, 0, 1, 5, 4 };
  int entries;
  for (int ii = 0; ii < (int)pressure.size(); ii++) {
    entries = 0;
    if (pressure[ii].doftype) {
      // inflow
      if (pressure[ii].coords[1] <= pn_epsilon) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
        rhs[ii] += 1;
      }
      // outflow
      else if (pressure[ii].coords[1] >= (par.width - pn_epsilon)) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
      }
      // other
      else {
        for (int dir = 0; dir < (par.dimension * 2); dir++) {
          if (pressure[ii].neighbors[dir] != -1) {
            entries++;
            temp_coo[entries - 1].i_index = ii;
            temp_coo[entries - 1].j_index = pressure[ii].neighbors[dir];
            temp_coo[entries - 1].value = -0.5 * (permeability[idx2(temp_coo[entries - 1].j_index, idx_swap[dir], (par.dimension * 2))] + permeability[idx2(ii, dir, (par.dimension * 2))]);
          }
        }
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 0;
        for (int dir = 0; dir < (entries - 1); dir++) {
          temp_coo[entries - 1].value -= temp_coo[dir].value;
        }
      }
      // place in array
      for (int jj = 0; jj < entries; jj++) {
        coo_array.push_back(temp_coo[jj]);
      }

    }
  }
}

void
hgf::models::porenetwork::setup_zflow_bc(const parameters& par)
{
  // quick exit
  if (par.dimension == 2) {
    std::cout << "\nError: zflow boundary conditions are not compatible with 2-dimensional problem. Exiting\n";
    exit(1);
  }
  array_coo temp_coo[7];
  double pn_epsilon = 1E-12;
  int idx_swap[6] = { 2, 3, 0, 1, 5, 4 };
  int entries;
  for (int ii = 0; ii < (int)pressure.size(); ii++) {
    entries = 0;
    if (pressure[ii].doftype) {
      // inflow
      if (pressure[ii].coords[2] <= pn_epsilon) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
        rhs[ii] += 1;
      }
      // outflow
      else if (pressure[ii].coords[2] >= (par.height - pn_epsilon)) {
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 1;
      }
      // other
      else {
        for (int dir = 0; dir < (par.dimension * 2); dir++) {
          if (pressure[ii].neighbors[dir] != -1) {
            entries++;
            temp_coo[entries - 1].i_index = ii;
            temp_coo[entries - 1].j_index = pressure[ii].neighbors[dir];
            temp_coo[entries - 1].value = -0.5 * (permeability[idx2(temp_coo[entries - 1].j_index, idx_swap[dir], (par.dimension * 2))] + permeability[idx2(ii, dir, (par.dimension * 2))]);
          }
        }
        entries++;
        temp_coo[entries - 1].i_index = ii;
        temp_coo[entries - 1].j_index = ii;
        temp_coo[entries - 1].value = 0;
        for (int dir = 0; dir < (entries - 1); dir++) {
          temp_coo[entries - 1].value -= temp_coo[dir].value;
        }
      }
      // place in array
      for (int jj = 0; jj < entries; jj++) {
        coo_array.push_back(temp_coo[jj]);
      }

    }
  }
}
