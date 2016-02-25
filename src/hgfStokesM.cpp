#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <magma.h>
#include <magmasparse.h>
#include <magma_lapack.h>
#include <stdio.h>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

// solves the stokes system directly as a single linear system
void
StokesSolveDirect( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par )
{

  // delcarations
  std::vector<int> matIs, matJs, rowPTR;;
  std::vector<double> matVals, force;
  matIs.reserve(Mesh.maxNNZ);
  matJs.reserve(Mesh.maxNNZ);
  matVals.reserve(Mesh.maxNNZ);
  force.resize(Mesh.dofTotal);

  // interior cells
  StokesArray( Mesh, Par.visc, matIs, matJs, matVals );

  // boundary conditions
  AxisFlowDrive( Mesh, matIs, matJs, matVals, force, Par.visc, Par.direction );

  // immersed boundary
  immersedBoundary( Mesh, matIs, matJs, matVals );

  // build the rowPTR vector for CSR rep of the array. first we sort the COO vecs
  sortCOO( matIs, matJs, matVals );
  buildCSR( matIs, matJs, matVals, rowPTR );

  // magma declarations
  magma_int_t info = 0;
  magma_init();
  magma_dopts opts;
  magma_queue_t queue = NULL;
  magma_queue_create( &queue );

  magma_d_matrix A={Magma_CSR}, d_A={Magma_CSR};
  magma_d_vector b, d_b, d_x;

  // Set A and b from CSR vectors
  magma_dcsrset( Mesh.dofTotal, Mesh.dofTotal, &rowPTR[0], &matJs[0], &matVals[0], &A, queue );
  magma_dvset( Mesh.dofTotal, 1, &force[0], &b, queue );

  printf( "\n%% matrix info: %d-by-%d with %d nonzeros\n\n",
            (int) A.num_rows,(int) A.num_cols,(int) A.nnz );

  // magma setup
  magma_dsolverinfo_init( &opts.solver_par, &opts.precond_par, queue );
  opts.solver_par.solver = Magma_GMRES;
  opts.solver_par.maxiter = Par.maxIt;
  opts.solver_par.rtol = Par.tolRel;
  opts.solver_par.atol = Par.tolAbs;
  opts.solver_par.restart = 30;
  opts.precond_par.solver = Magma_NONE;
  opts.blocksize = 32;
  opts.alignment = 1;
  A.blocksize = opts.blocksize;
  A.alignment = opts.alignment;
  magma_d_precondsetup( A, b, &opts.solver_par, &opts.precond_par, queue );
  magma_dmtransfer( A, &d_A, Magma_CPU, Magma_DEV, queue );
  magma_dmtransfer( b, &d_b, Magma_CPU, Magma_DEV, queue );
  magma_dvinit( &d_x, Magma_DEV, Solution.size(), 1, Solution[0], queue );

  // solve the system
  magma_d_solver( d_A, d_b, &d_x, &opts, queue );
  magma_dsolverinfo( &opts.solver_par, &opts.precond_par, queue );

  magma_int_t mout,nout;
  double * valout;
  // bring the solution back to host, fill std::vector solution
  magma_dvget( d_x, &mout, &nout, &valout, queue );
  for (int ii = 0; ii < Solution.size(); ii++) {
    Solution[ii] = valout[ii];
  }

  // magma frees
  magma_dsolverinfo_free( &opts.solver_par, &opts.precond_par, queue );
  magma_d_mfree( &d_A, queue );
  magma_d_vfree( &d_b, queue );
  magma_d_vfree( &d_x, queue );
  magma_queue_destroy( queue );
  magma_finalize();
}
