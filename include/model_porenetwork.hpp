#ifndef _MODELS_PORENETWORK_H
#define _MODELS_PORENETWORK_H

#include "hgflow.hpp"

// system includes
#include <vector>

namespace hgf
{
  /** \brief Contains functions and classes that setup and post-process flow models.
   *
   */
  namespace models
  {
    /** \brief Contains functionality for setup and post-processessing the solution of a simple uniform porenetwork model.
     * 
     */
    class uniform_porenetwork
    {

      public:

        std::vector< degree_of_freedom > pressure;             /**< Degrees of freedom associated with the fluid pressure. */
        std::vector< double > permeability;                    /**< Vector of permeabilities in P-N model throats */
        std::vector< array_coo > coo_array;                    /**< Linear system associated to the P-N problem stored in COO (coordinate) sparse format */
        std::vector< double > rhs;                             /**< Right-hand side vector (force). */
        std::vector< double > solution;                        /**< Vector for P-N solution */
        void build_uniform_network(const parameters& par, int n_pores_x, int n_pores_y, int n_pores_z);
        void init_permeability_one(const parameters& par);
        void init_permeability_random(const parameters& par, double min, double max);
        void build(const parameters& par);
        void setup_xflow_bc(const parameters& par);
        void setup_yflow_bc(const parameters& par);
        void setup_zflow_bc(const parameters& par);
        void output_vtk(const parameters& par);

    };
  }
}

#endif
