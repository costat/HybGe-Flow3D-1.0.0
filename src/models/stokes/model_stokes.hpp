#ifndef _MODELS_STOKES_H
#define _MODELS_STOKES_H

#include "../../../include/hgflow.hpp"

// system includes
#include <vector>
#include <array>
#include <iostream>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <boost/filesystem.hpp>

namespace hgf
{
  namespace models
  {
    // Stokes
    class stokes
    {

      public:
      
        std::vector< degree_of_freedom > velocity_u, velocity_v, velocity_w, pressure;
        std::vector< int > ptv;       /* pressure.size() x (DIM x 2) size array containing
                                         IDs of velocity cells related to this pressure cell */
        std::vector< int > interior_u, interior_v, interior_w; 
        std::vector< int > interior_u_nums, interior_v_nums, interior_w_nums;
        std::vector< boundary_nodes > boundary;
        std::vector< array_coo > coo_array;
        std::vector< double > rhs, solution, interior;
        void build(const parameters& par, const hgf::mesh& msh);
        void solution_insert_boundaries(void);
        void output_vtk(const parameters& par, const hgf::mesh& msh);
        void setup_xflow_bc(const parameters& par, const hgf::mesh& msh);
      
      private:
        
        void build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh& msh);
        void dof_neighbors_2d(const parameters& par, const hgf::mesh& msh);
        void build_array_2d(const parameters& par, const hgf::mesh& msh);
        void momentum_2d(double visc);
        void continuity_2d(void);
        void xflow_2d(const parameters& par, const hgf::mesh& msh);

        void build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh& msh);
        void dof_neighbors_3d(const parameters& par, const hgf::mesh& msh);
        void build_array_3d(const parameters& par, const hgf::mesh& msh);
        void momentum_3d(double visc);
        void continuity_3d(void);
        void xflow_3d(const parameters& par, const hgf::mesh& msh);

    };
  }
}

#endif
