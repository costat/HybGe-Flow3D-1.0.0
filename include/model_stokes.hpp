#ifndef _MODELS_STOKES_H
#define _MODELS_STOKES_H

#include "hgflow.hpp"

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
  /**
   * \brief Contains functions and classes that setup and post-process flow models.
   *
   */
  namespace models
  {
    /** \brief Contains functionality for setup and post-processessing the solution of Stokes fluid flow models in 2d or 3d.
     * 
     */
    class stokes
    {

      public:
      
        std::vector< degree_of_freedom > velocity_u;                  /**< Degrees of freedom associated with the x-component of the fluid velocity. */
        std::vector< degree_of_freedom > velocity_v;                  /**< Degrees of freedom associated with the y-component of the fluid velocity. */
        std::vector< degree_of_freedom > velocity_w;                  /**< Degrees of freedom associated with the z-component of the fluid velocity. */
        std::vector< degree_of_freedom > pressure;                    /**< Degrees of freedom associated with the fluid pressure. */
        std::vector< int > interior_u;                                /**< If interior_u[i] == 1, then velocity_u[i] is an internal degree of freedom. */
        std::vector< int > interior_v;                                /**< If interior_v[i] == 1, then velocity_v[i] is an internal degree of freedom. */
        std::vector< int > interior_w;                                /**< If interior_w[i] == 1, then velocity_w[i] is an internal degree of freedom. */
        std::vector< int > pressure_ib_list;                          /**< If pressure_ib_list[i] == 1, then pressure[i] and it's associated staggered velocity components are immersed boundary cells */
        std::vector< array_coo > coo_array;                           /**< Linear system associated to the Stokes' problem stored in COO (coordinate) sparse format */
        std::vector< double > rhs;                                    /**< Right-hand side vector (force). */
        std::vector< double > solution;                               /**< Vector for storing full solution, including interior and boundary DOFs. */
        std::vector< double > solution_int;                           /**< Vector for storing solution for interior DOFs. Corresponds to produced coo_array and RHS, which are built with boundary DOFs eliminated. */
        void build(const parameters& par, const hgf::mesh& msh);
        void solution_build(void);
        void output_vtk(const parameters& par, const hgf::mesh& msh, std::string& file_name);
        void setup_xflow_bc(const parameters& par, const hgf::mesh& msh);
        void setup_yflow_bc(const parameters& par, const hgf::mesh& msh);
        void setup_zflow_bc(const parameters& par, const hgf::mesh& msh);
        void random_immersed_boundary(const parameters& par, double eta, double vol_frac);
        int random_immersed_boundary_clump(const parameters& par, double eta, double vol_frac, double likelihood);
        void immersed_boundary(const parameters& par, double eta);
        void import_immersed_boundary(parameters& par, std::vector< int >& input_ib, double eta);
    
      private:
        
        std::vector< boundary_nodes > boundary;                    
        std::vector< int > interior_u_nums, interior_v_nums, interior_w_nums;
        std::vector< int > ptv;
        void build_degrees_of_freedom_2d(const parameters& par, const hgf::mesh& msh);
        void dof_neighbors_2d(const parameters& par, const hgf::mesh& msh);
        void build_array_2d(const parameters& par, const hgf::mesh& msh);
        void momentum_2d(double visc);
        void continuity_2d(void);
        void xflow_2d(const parameters& par, const hgf::mesh& msh);
        void yflow_2d(const parameters& par, const hgf::mesh& msh);

        void build_degrees_of_freedom_3d(const parameters& par, const hgf::mesh& msh);
        void dof_neighbors_3d(const parameters& par, const hgf::mesh& msh);
        void build_array_3d(const parameters& par, const hgf::mesh& msh);
        void momentum_3d(double visc);
        void continuity_3d(void);
        void xflow_3d(const parameters& par, const hgf::mesh& msh);
        void yflow_3d(const parameters& par, const hgf::mesh& msh);
        void zflow_3d(const parameters& par, const hgf::mesh& msh);

    };
  }
}

#endif
