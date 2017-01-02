/* solvers header */
#ifndef _SOLVE_H
#define _SOLVE_H

/* To add a solver include the entry point header here. 
   A solver should take a COO array and a pointer to an empty
   solution vector as inputs and output the solution to that pointer. 
   If the solver src is placed within HGFROOT/src/models/ cmake will include all src 
   and headers.*/

// HGF ships publically with hooks into Paralution
#ifdef _PARALUTION
  #include "../depsrc/solve/paralution/solve_paralution.hpp"
#endif

#endif
