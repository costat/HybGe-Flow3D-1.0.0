/* porenetwork array source */

// hgf includes
#include "model_porenetwork.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

/** \brief hgf::models::uniform_porenetwork::build builds the linear system for a uniform porenetwork.
 *
 * It is assumed that the P-N DOFs and permeabilities have already been set.
 *
 * @param[in] par - parameters struct containing problem information.
 */
void
hgf::models::uniform_porenetwork::build(const parameters& par)
{
  rhs.resize(pressure.size());
  solution.resize(pressure.size());
  array_coo temp_coo[7];
  if (par.dimension == 2) {
    for (int ii = 0; ii < (int)pressure.size(); ii++) {
      if (!pressure[ii].doftype) {
        // indexes
        for (int dir = 0; dir < 4; dir++) {
          temp_coo[dir].i_index = ii;
          temp_coo[dir].j_index = pressure[ii].neighbors[dir];
        }
        temp_coo[4].i_index = ii;
        temp_coo[4].j_index = ii;
        // values
        temp_coo[0].value = -0.5 * (permeability[idx2(temp_coo[0].j_index, 2, 4)] + permeability[idx2(ii, 0, 4)]);
        temp_coo[1].value = -0.5 * (permeability[idx2(temp_coo[1].j_index, 3, 4)] + permeability[idx2(ii, 1, 4)]);
        temp_coo[2].value = -0.5 * (permeability[idx2(temp_coo[2].j_index, 0, 4)] + permeability[idx2(ii, 2, 4)]);
        temp_coo[3].value = -0.5 * (permeability[idx2(temp_coo[3].j_index, 1, 4)] + permeability[idx2(ii, 3, 4)]);
        temp_coo[4].value = -(temp_coo[0].value + temp_coo[1].value + temp_coo[2].value + temp_coo[3].value);
        for (int dir = 0; dir < 5; dir++) {
          coo_array.push_back(temp_coo[dir]);
        }
      }
    }
  }
  else {
    for (int ii = 0; ii < (int)pressure.size(); ii++) {
      if (!pressure[ii].doftype) {
        // indexes
        for (int dir = 0; dir < 6; dir++) {
          temp_coo[dir].i_index = ii;
          temp_coo[dir].j_index = pressure[ii].neighbors[dir];
        }
        temp_coo[6].i_index = ii;
        temp_coo[6].j_index = ii;
        // values
        temp_coo[0].value = -0.5 * (permeability[idx2(temp_coo[0].j_index, 2, 6)] + permeability[idx2(ii, 0, 6)]);
        temp_coo[1].value = -0.5 * (permeability[idx2(temp_coo[1].j_index, 3, 6)] + permeability[idx2(ii, 1, 6)]);
        temp_coo[2].value = -0.5 * (permeability[idx2(temp_coo[2].j_index, 0, 6)] + permeability[idx2(ii, 2, 6)]);
        temp_coo[3].value = -0.5 * (permeability[idx2(temp_coo[3].j_index, 1, 6)] + permeability[idx2(ii, 3, 6)]);
        temp_coo[4].value = -0.5 * (permeability[idx2(temp_coo[4].j_index, 5, 6)] + permeability[idx2(ii, 4, 6)]);
        temp_coo[5].value = -0.5 * (permeability[idx2(temp_coo[5].j_index, 4, 6)] + permeability[idx2(ii, 5, 6)]);
        temp_coo[6].value = -(temp_coo[0].value + temp_coo[1].value + temp_coo[2].value + \
                              temp_coo[3].value + temp_coo[4].value + temp_coo[5].value);
        for (int dir = 0; dir < 7; dir++) {
          coo_array.push_back(temp_coo[dir]);
        }
      }
    }
  }
}
