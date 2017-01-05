/* solvers header */
#ifndef _SOLVE_H
#define _SOLVE_H

/* To add a solver include the entry point header here. */

// hooks into MKL
#ifdef _MKL
  #include "../depsrc/solve/MKL/solve_mkl.hpp"
#endif

// hooks into Paralution
#ifdef _PARALUTION
  #include "../depsrc/solve/paralution/solve_paralution.hpp"
#endif

#endif
