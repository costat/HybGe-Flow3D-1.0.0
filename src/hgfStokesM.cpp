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
StokesSolveDirect( const FluidMesh& Mesh, double visc, int direction, \
                   std::vector<double>& Solution, double tolAbs, double tolRel, \
                   int maxIt, int nThreads, int prec )
{

  // delcarations
  std::vector<int> matIs, matJs, rowPTR;;
  std::vector<double> matVals, force;
  matIs.reserve(Mesh.maxNNZ);
  matJs.reserve(Mesh.maxNNZ);
  matVals.reserve(Mesh.maxNNZ);
  force.resize(Mesh.dofTotal);

  // interior cells
  StokesArray( Mesh, visc, matIs, matJs, matVals );

  // boundary conditions
  AxisFlowDrive( Mesh, matIs, matJs, matVals, force, visc, direction );

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
  opts.solver_par.maxiter = 10000;
  opts.solver_par.rtol = 1e-8;
  opts.solver_par.atol = 1e-8;
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
// initializes the pressure solution to a linear. for use in decoupled iterative schemes, e.g. the urzawa schemes
void
initPressure( const FluidMesh& Mesh, std::vector<double>& Solution, int direction )
{
  int shift, cl2;
  if (Mesh.DIM == 2) shift = Mesh.DOF[1] + Mesh.DOF[2];
  else shift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
  switch (direction)
  {
    case 0 :
    {
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        cl2 = cl + shift;
        Solution[ cl2 ] = 1 - (Mesh.xLim[0] - Mesh.PCellCenters[ idx2( cl, 0, Mesh.DIM ) ]) / (Mesh.xLim[1] - Mesh.xLim[0]);
      }
      break;
    }
    case 1 :
    {
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        cl2 = cl + shift;
        Solution[ cl2 ] = 1 - (Mesh.yLim[0] - Mesh.PCellCenters[ idx2( cl, 1, Mesh.DIM ) ]) / (Mesh.yLim[1] - Mesh.yLim[0]);
      }
      break;
    }
    case 2 :
    {
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        cl2 = cl + shift;
        Solution[ cl2 ] = 1 - (Mesh.zLim[0] - Mesh.PCellCenters[ idx2( cl, 2, Mesh.DIM ) ]) / (Mesh.zLim[1] - Mesh.zLim[0]);
      }
      break;
    }
  }
}

void
setForceRich( const FluidMesh& Mesh, const std::vector<double>& Solution, std::vector<double>& b, \
                                                   const std::vector<double>& force, int component )
{
  int velShift;
  if ( Mesh.DIM == 3 ) velShift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
  else velShift = Mesh.DOF[1] + Mesh.DOF[2];
  switch (component)
  {
    case 0 :
    {
      for (int cl = 0; cl < Mesh.DOF[1]; cl++) {
        if ( Mesh.UCellPressureNeighbor[ idx2( cl, 0, 2 ) ] && Mesh.UCellPressureNeighbor[ idx2( cl, 1, 2 ) ] )
        { // has two pressure neighbors, treated as interior cell with regards to pressure in force
          b[cl] = force[cl] - ( Solution[ Mesh.UCellPressureNeighbor[ idx2( cl, 1, 2 ) ] + velShift - 1 ] \
                               -Solution[ Mesh.UCellPressureNeighbor[ idx2( cl, 0, 2 ) ] + velShift - 1 ] ) \
                              / Mesh.UCellWidths[ idx2( cl, 0, Mesh.DIM ) ];
        }
        else // this cell center is on the true geo boundary - no adjacent pressure term in 1 direction, no pressure term in force
        {
          b[cl] = force[cl];
        }
      }
      break; // case 0 break
    }
    case 1 :
    {
      for (int cl = 0; cl < Mesh.DOF[2]; cl++) {
        if ( Mesh.VCellPressureNeighbor[ idx2( cl, 0, 2 ) ] && Mesh.VCellPressureNeighbor[ idx2( cl, 1, 2 ) ] )
        { // has two pressure neighbors, treated as interior cell with regards to pressure in force
          b[cl] = force[cl] - ( Solution[ Mesh.VCellPressureNeighbor[ idx2( cl, 1, 2 ) ] + velShift - 1 ] \
                               -Solution[ Mesh.VCellPressureNeighbor[ idx2( cl, 0, 2 ) ] + velShift - 1 ] ) \
                              / Mesh.VCellWidths[ idx2( cl, 1, Mesh.DIM ) ];
        }
        else // this cell center is on the true geo boundary - no adjacent pressure term in 1 direction, no pressure term in force
        {
          b[cl] = force[cl];
        }
      }
      break; // case 1 break
    }
    case 2 :
    {
      for (int cl = 0; cl < Mesh.DOF[3]; cl++) {
        if ( Mesh.WCellPressureNeighbor[ idx2( cl, 0, 2 ) ] && Mesh.WCellPressureNeighbor[ idx2( cl, 1, 2 ) ] )
        { // has two pressure neighbors, treated as interior cell with regards to pressure in force
          b[cl] = force[cl] - ( Solution[ Mesh.WCellPressureNeighbor[ idx2( cl, 1, 2 ) ] + velShift - 1 ] \
                               -Solution[ Mesh.WCellPressureNeighbor[ idx2( cl, 0, 2 ) ] + velShift - 1 ] ) \
                              / Mesh.WCellWidths[ idx2( cl, 1, Mesh.DIM ) ];
        }
        else // this cell center is on the true geo boundary - no adjacent pressure term in 1 direction, no pressure term in force
        {
          b[cl] = force[cl];
        }
      }
      break; // case 2 break
    }
  }
}
void
updatePressureRich( const FluidMesh& Mesh, std::vector<double>& Solution, \
                    const std::vector<double>& solU, const std::vector<double>& solV, \
                    const std::vector<double>& solW, double& residual, double relax )
{
  int velShift;
  std::vector<double> pressureHold;
  pressureHold.resize( Mesh.DOF[0] );

  switch ( Mesh.DIM )
  {
    case 2 :
    {
      velShift = Mesh.DOF[1] + Mesh.DOF[2];
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        pressureHold[cl] = Solution[ (velShift + cl) ] - relax * ( \
                           ( solU[ Mesh.PressureCellUNeighbor[ idx2( cl, 1, 2 ) ]-1 ] - solU[ Mesh.PressureCellUNeighbor[ idx2( cl, 0, 2 ) ]-1 ] ) \
                           / Mesh.PCellWidths[ idx2( cl, 0, Mesh.DIM ) ] + \
                           ( solV[ Mesh.PressureCellVNeighbor[ idx2( cl, 1, 2 ) ]-1 ] - solV[ Mesh.PressureCellVNeighbor[ idx2( cl, 0, 2 ) ]-1 ] ) \
                           / Mesh.PCellWidths[ idx2( cl, 1, Mesh.DIM ) ] );
      }
      break;
    }
    default :
    {
      velShift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
      for (int cl = 0; cl < Mesh.DOF[0]; cl++) {
        pressureHold[cl] = Solution[ (velShift + cl) ] - relax * ( \
                           ( solU[ Mesh.PressureCellUNeighbor[ idx2( cl, 1, 2 ) ]-1 ] - solU[ Mesh.PressureCellUNeighbor[ idx2( cl, 0, 2 ) ]-1 ] ) \
                           / Mesh.PCellWidths[ idx2( cl, 0, Mesh.DIM ) ] + \
                           ( solV[ Mesh.PressureCellVNeighbor[ idx2( cl, 1, 2 ) ]-1 ] - solV[ Mesh.PressureCellVNeighbor[ idx2( cl, 0, 2 ) ]-1 ] ) \
                           / Mesh.PCellWidths[ idx2( cl, 1, Mesh.DIM ) ] + \
                           ( solW[ Mesh.PressureCellWNeighbor[ idx2( cl, 1, 2 ) ]-1 ] - solW[ Mesh.PressureCellWNeighbor[ idx2( cl, 0, 2 ) ]-1 ] ) \
                           / Mesh.PCellWidths[ idx2( cl, 2, Mesh.DIM ) ] );
      }
      break;
    }
  }
  // loop over vector components, calc residual contribution then set new pressure
  residual = 0;
  for (int cl = 0; cl < Mesh.DOF[0]; cl++)
  {
    residual = residual + pow( (pressureHold[ cl ] - Solution[ (velShift + cl) ]), 2.0 );
    Solution[ (velShift + cl) ] = pressureHold[ cl ];
  }
  residual = sqrt( residual );
}

/* solves the stokes system with a first order richardson scheme for the schur complement.
   velocity solve step broken into components */
void
StokesSolveRich( const FluidMesh& Mesh, double visc, int direction, \
                 std::vector<double>& Solution, double tolAbs, double tolRel, \
                 int maxIt, int nThreads, int prec, double relax )
{

  int richIt = 0;
  int richItMax = 1000;
  double residual = 1.0;
  double tolerance = 1e-10;

  // delcarations
  std::vector<int> matIsU, matJsU, matIsV, matJsV, matIsW, matJsW;
  std::vector<double> matValsU, matValsV, matValsW, forceU, forceV, forceW, bU, bV, bW;

  matIsU.reserve( Mesh.maxNNZ );
  matJsU.reserve( Mesh.maxNNZ );
  matValsU.reserve( Mesh.maxNNZ );
  forceU.resize( Mesh.DOF[1] );
  bU.resize( Mesh.DOF[1] );

  matIsV.reserve( Mesh.maxNNZ );
  matJsV.reserve( Mesh.maxNNZ );
  matValsV.reserve( Mesh.maxNNZ );
  forceV.resize( Mesh.DOF[2] );
  bV.resize( Mesh.DOF[2] );

  if (Mesh.DIM == 3)
  {
    matIsW.reserve( Mesh.maxNNZ );
    matJsW.reserve( Mesh.maxNNZ );
    matValsW.reserve( Mesh.maxNNZ );
    forceW.resize( Mesh.DOF[3] );
    bW.resize( Mesh.DOF[3] );
  }

  // interior cells
  VelocityArray( Mesh, visc, matIsU, matJsU, matValsU, 0 );
  VelocityArray( Mesh, visc, matIsV, matJsV, matValsV, 1 );
  if (Mesh.DIM == 3) VelocityArray( Mesh, visc, matIsW, matJsW, matValsW, 2 );

  // boundary conditions
  AxisFlowSingleComponent( Mesh, matIsU, matJsU, matValsU, forceU, visc, direction, 0 );
  AxisFlowSingleComponent( Mesh, matIsV, matJsV, matValsV, forceV, visc, direction, 1 );
  if (Mesh.DIM == 3) AxisFlowSingleComponent( Mesh, matIsW, matJsW, matValsW, forceW, visc, direction, 2 );

  // immersed boundary
  immersedBoundarySingleComponent( Mesh, matIsU, matJsU, matValsU, 0 );
  immersedBoundarySingleComponent( Mesh, matIsV, matJsV, matValsV, 1 );
  if (Mesh.DIM == 3) immersedBoundarySingleComponent( Mesh, matIsW, matJsW, matValsW, 2 );

  // initialize pressure solution
  initPressure( Mesh, Solution, direction );

  // magma setup

  // richardson loop
  while (1)
  {
    richIt++;
    // compute new forces
    setForceRich( Mesh, Solution, bU, forceU, 0 );
    setForceRich( Mesh, Solution, bV, forceV, 1 );
    if (Mesh.DIM == 3) setForceRich( Mesh, Solution, bW, forceW, 2 );

    // magma solves

    // update the pressure vector and compute residual

    std::cout << "\nRichardson Iteration " << richIt << " \t Residual: " << residual << "\n";

    if ( residual < tolerance )
    {
      std::cout << "\n Residual tolerance criteria satisfied, exiting Richardson loop.";
      break;
    }
    else if ( richIt >= richItMax )
    {
      std::cout << "\n Maximum iterations reached, exiting Richardson loop.";
      break;
    }

  }

  // post process??

}
