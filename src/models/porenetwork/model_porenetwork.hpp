#ifndef _MODELS_PORENETWORK_H
#define _MODELS_PORENETWORK_H

#include "../../../include/hgflow.hpp"

// system includes
#include <vector>

namespace hgf
{
  namespace models
  {
    class porenetwork
    {

      public:

        std::vector< degree_of_freedom > pressure;
        std::vector< double > permeability;
        std::vector< array_coo > coo_array;
        std::vector< double > rhs, solution;
        void build_uniform_network(const parameters& par, int n_pores_x, int n_pores_y, int n_pores_z);
        void import_network(const parameters& par);
        void init_permeability_one(const parameters& par);
        void init_permeability_random(const parameters& par, double min, double max);
        void build(const parameters& par);
        void setup_xflow_bc(const parameters& par);
        void setup_yflow_bc(const parameters& par);
        void setup_zflow_bc(const parameters& par);
        void output_vtk(const parameters& par);

      private:

        //void build_degrees_of_freedom_2d(const parameters& par);
        //void build_array_2d(const parameters& par);
        //void xflow_2d(const parameters& par);

        //void build_degrees_of_freedom_3d(const parameters& par);
        //void build_array_3d(const parameters& par);
        //void xflow_3d(const parameters& par);

    };
  }
}

#endif
