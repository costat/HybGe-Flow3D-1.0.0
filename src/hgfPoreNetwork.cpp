#include <vector>
#include <omp.h>
#include <math.h>
#include <paralution.hpp>

// hgf includes
#ifndef CUDA_BUILD
# define CUDA_BUILD 0
#endif

#if CUDA_BUILD
#include "hgfMeshCu.cuh"
#else
#include "hgfMesh.hpp"
#endif

#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfPoreNetwork.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

void
PoreNetworkSolveDirect( const PoreNetwork& pn, const std::vector<double>& Ks, \
                        std::vector<double>& Solution, int direction )
{
  // delcarations
  std::vector<int> matIs, matJs;
  std::vector<double> matVals, force;
  matIs.reserve( (pn.nPores * pn.DIM * 2 + 1) );
  matJs.reserve( (pn.nPores * pn.DIM * 2 + 1) );
  matVals.reserve( (pn.nPores * pn.DIM * 2 + 1) );
  force.resize( pn.nPores );

  // interior pores
  PoreNetworkArray( pn, matIs, matJs, matVals, Ks );

  // boundary pores
  PoreNetworkBoundary( pn, matIs, matJs, matVals, force, Ks, direction );

  // paralution arrays
  LocalVector<double> sol;
  LocalVector<double> forceP;
  LocalMatrix<double> mat;

  // initialize force and solution vectors
  forceP.Allocate("force vector", pn.nPores);
  forceP.Zeros();
  sol.Allocate("solution", pn.nPores);
  sol.Zeros();

    // assemble paralution arrays from COO data
  mat.Assemble( &matIs[0], &matJs[0], &matVals[0], matIs.size(), \
                "operator", pn.nPores, pn.nPores );

  for (int cl = 0; cl < pn.nPores; cl++) {
    forceP[cl] = force[cl];
  }
  std::cout << "\nCheck\n";

  // GMRES object
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
  ls.SetOperator(mat);

  // build
  ls.Build();

  // solve the system
  ls.Solve(forceP, &sol);

  // pass solution to input vector
  for (int cl = 0; cl < pn.nPores; cl++) {
    Solution[cl] = sol[cl];
  }

  // clear paralution objects
  ls.Clear();
  mat.Clear();
  forceP.Clear();
  sol.Clear();

}
