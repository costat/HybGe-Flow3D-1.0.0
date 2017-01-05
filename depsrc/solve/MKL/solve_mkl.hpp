#ifndef _SOLVE_MKL_H
#define _SOLVE_MKL_H

#include "../../../include/hgflow.hpp"

// system includes
#include <vector>
#include <omp.h>
#include <mkl.h>

namespace hgf
{
  namespace solve
  {
    namespace mkl
    {
      void init_solver(void);
      void finalize_solver(void);
      void solve(const parameters& par, \
        const std::vector< array_coo >& array, \
        const std::vector< double >& rhs, \
        std::vector<double>& solution);
    }
  }
}

#endif
