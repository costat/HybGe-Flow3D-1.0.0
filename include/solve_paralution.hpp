#ifndef _SOLVE_PARALUTION_H
#define _SOLVE_PARALUTION_H

#include "hgflow.hpp"

// system includes
#include <vector>
#include <paralution.hpp>
#include <omp.h>

namespace hgf
{
  /** \brief Contains functions that solve hgf models.
   *
   */
  namespace solve
  {
    /** \brief Contains functions that use paralution preconditions and solvers to solve hgf models.
     * 
     */
    namespace paralution
    {
      void init_solver(void);
      void finalize_solver(void);
      void solve(const parameters& par, \
        const std::vector< array_coo >& array, \
        const std::vector< double >& rhs, \
        std::vector<double>& solution);
      void solve_ps_flow(const parameters& par, \
        const std::vector< array_coo >& array, \
        const std::vector< double >& rhs, \
        std::vector<double>& solution, \
        int n_u, int n_v, int n_w, int n_p);
    }
  }
}

#endif
