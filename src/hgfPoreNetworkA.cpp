#include <vector>
#include <omp.h>
#include <math.h>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfPoreNetworkA.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

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

  // magma solve

}
