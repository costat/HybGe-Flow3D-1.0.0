/* porenetwork permeability functions source */

// hgf includes
#include "model_porenetwork.hpp"

#include <random>

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::porenetwork::init_permeability_one(const parameters& par)
{
  permeability.assign(pressure.size() * par.dimension * 2, 0.0);
  for (int ii = 0; ii < pressure.size(); ii++) {
    for (int jj = 0; jj < (par.dimension * 2); jj++) {
      if (pressure[ii].neighbors[jj] != -1) {
        permeability[idx2(ii, jj, (par.dimension * 2))] = 1.0;
      }
    }
  }
}

void 
hgf::models::porenetwork::init_permeability_random(const parameters& par, double min, double max)
{
  permeability.assign(pressure.size() * par.dimension * 2, 0.0);
  std::uniform_real_distribution<double> unif(min, max);
  std::default_random_engine re;

  for (int ii = 0; ii < pressure.size(); ii++) {
    for (int jj = 0; jj < (par.dimension * 2); jj++) {
      if (pressure[ii].neighbors[jj] != -1) {
        permeability[idx2(ii, jj, (par.dimension * 2))] = unif(re);
      }
    }
  }
}
