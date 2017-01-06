#ifndef _MULTISCALE_FLOW_H
#define _MULTISCALE_FLOW_H

#include "../../../include/hgflow.hpp"

namespace hgf
{
  namespace multiscale
  {
    namespace flow
    {
      double compute_permeability_x(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
      double compute_permeability_y(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
      double compute_permeability_z(const parameters& par, const std::vector< degree_of_freedom >& velocity_u, \
                                                           const std::vector< degree_of_freedom >& velocity_v, \
                                                           const std::vector< degree_of_freedom >& velocity_w, \
                                                           const std::vector< double > solution);
    }
  }
}

#endif
