#include <vector>
#include <omp.h>

#include <paralution.hpp>

#include "hgfMesh.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;
// solves the stokes system directly as a single linear system
void
StokesSolveDirect( const FluidMesh& Mesh, double visc, int direction, \
                   std::vector<double>& Solution, double tolAbs, double tolRel, \
                   int maxIt, int nThreads, int prec )
{

  // delcarations
  std::vector<int> matIs, matJs;
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

  set_omp_threads_paralution(nThreads);

  // paralution arrays
  LocalVector<double> sol;
  LocalVector<double> forceP;
  LocalMatrix<double> mat;

  // initialize force and solution vectors
  forceP.Allocate("force vector", Mesh.dofTotal);
  forceP.Zeros();
  sol.Allocate("solution", Mesh.dofTotal);
  sol.Zeros();

  // assemble paralution arrays from COO data
  mat.Assemble( &matIs[0], &matJs[0], &matVals[0], matIs.size(), \
                "operator", Mesh.dofTotal, Mesh.dofTotal);

  for (int cl = 0; cl < Mesh.dofTotal; cl++) {
    forceP[cl] = force[cl];
  }

  // GMRES object
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
  ls.Init(tolAbs, tolRel, 1e8, maxIt);
  ls.SetOperator(mat);
  ls.Verbose(2);
  ls.SetBasisSize(100);

  // preconditioning
  ILU<LocalMatrix<double>, LocalVector<double>, double> p;
  p.Set(prec);
  ls.SetPreconditioner(p);

  // build
  ls.Build();

  // solve the system
  ls.Solve(forceP, &sol);

  // pass solution from paralution object to std vector input
  for (int cl = 0; cl < Mesh.dofTotal; cl++)
  {
    Solution[cl] = sol[cl];
  }

  // clear paralution objects
  ls.Clear();
  mat.Clear();
  p.Clear();
  forceP.Clear();
  sol.Clear();

}
/* solves the stokes system with a first order richardson scheme for the schur complement.
   velocity solve step broken into components */
void
StokesSolveRich( const FluidMesh& Mesh, double visc, int direction, \
                   std::vector<double>& Solution, double tolAbs, double tolRel, \
                   int maxIt, int nThreads, int prec )
{

  // delcarations
  std::vector<int> matIsU, matJsU, matIsV, matJsV, matIsW, matJsW;
  std::vector<double> matValsU, matValsV, matValsW, forceU, forceV, forceW, bV, bW, bU;
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
  if ( Mesh.DIM == 3 )
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
  if ( Mesh.DIM == 3 ) VelocityArray( Mesh, visc, matIsW, matJsW, matValsW, 2 );

  // boundary conditions
  AxisFlowSingleComponent( Mesh, matIsU, matJsU, matValsU, forceU, visc, direction, 0 );
  AxisFlowSingleComponent( Mesh, matIsV, matJsV, matValsV, forceV, visc, direction, 1 );
  if ( Mesh.DIM == 3 ) AxisFlowSingleComponent( Mesh, matIsW, matJsW, matValsW, forceW, visc, direction, 2 );

  // immersed boundary
  immersedBoundarySingleComponent( Mesh, matIsU, matJsU, matValsU, 0 );
  immersedBoundarySingleComponent( Mesh, matIsV, matJsV, matValsV, 1 );
  if ( Mesh.DIM == 3 ) immersedBounarySingleComponent( Mesh, matIsW, matJsW, matValsW, 2 );

  set_omp_threads_paralution(nThreads);

  // paralution arrays
  LocalVector<double> solU, solV, solW;
  LocalVector<double> forceUP, forceVP, forceWP;
  LocalMatrix<double> matU, matV, matW;

  // initialize solution and force, define initial guess for pressure
  forceUP.Allocate("force vector U", Mesh.DOF[1]);
  forceUP.Zeros();
  solU.Allocate("solution U", Mesh.DOF[1]);
  solU.Zeros();
  forceVP.Allocate("force vector V", Mesh.DOF[2]);
  forceVP.Zeros();
  solV.Allocate("solution V", Mesh.DOF[2]);
  solV.Zeros();
  if ( Mesh.DIM == 3 )
  {
    forceWP.Allocate("force vector W", Mesh.DOF[2]);
    forceWP.Zeros();
    solW.Allocate("solution W", Mesh.DOF[2]);
    solW.Zeros();
  }
  initPressure( );

  // assemble the matrix from COO data

  // GMRES object

  // preconditioning

  // build

  // richardson loop

  // pass solution from paralution object to std vector input

  // clear paralution objects

}

