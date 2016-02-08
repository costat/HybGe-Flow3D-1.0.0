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
  initPressure( Mesh, Solution );

  // assemble the matrices from COO data
  matU.Assemble( &matIsU[0], &matJsU[0], &matValsU[0], matIsU.size(), \
                 "U operator", Mesh.DOF[1], Mesh.DOF[1]);
  matV.Assemble( &matIsV[0], &matJsV[0], &matValsV[0], matIsV.size(), \
                 "V operator", Mesh.DOF[2], Mesh.DOF[2]);
  if (Mesh.DIM == 3) matW.Assemble( &matIsW[0], &matJsW[0], &matValsW[0], matIsW.size(), \
                     "W operator", Mesh.DOF[3], Mesh.DOF[3]);

  // GMRES object
  GMRES<LocalMatrix<double>, LocalVector<double>, double> lsU;
  lsU.Init(tolAbs, tolRel, 1e8, maxIt);
  lsU.SetOperator(matU);
  lsU.Verbose(1);
  lsU.SetBasisSize(100);
  GMRES<LocalMatrix<double>, LocalVector<double>, double> lsV;
  lsV.Init(tolAbs, tolRel, 1e8, maxIt);
  lsV.SetOperator(matV);
  lsV.Verbose(1);
  lsV.SetBasisSize(100);
  GMRES<LocalMatrix<double>, LocalVector<double>, double> lsW;
  if ( Mesh.DIM == 3) {
    lsW.Init(tolAbs, tolRel, 1e8, maxIt);
    lsW.SetOperator(matW);
    lsW.Verbose(1);
    lsW.SetBasisSize(100);
  }

  // preconditioning
  ILU<LocalMatrix<double>, LocalVector<double>, double> pU;
  pU.Set(prec);
  lsU.SetPreconditioner(pU);
  ILU<LocalMatrix<double>, LocalVector<double>, double> pV;
  pV.Set(prec);
  lsV.SetPreconditioner(pV);
  ILU<LocalMatrix<double>, LocalVector<double>, double> pW;
  if ( Mesh.DIM == 3)
  {
    pW.Set(prec);
    lsW.SetPreconditioner(pW);
  }

  // build
  lsU.Build();
  lsV.Build();
  if ( Mesh.DIM == 3 ) lsW.Build();

  // richardson loop
  while (1)
  {

    // compute new forces
    setForceRich( Mesh, Solution, bU, forceU, 0 );
    setForceRich( Mesh, Solution, bV, forceV, 1 );
    if (Mesh.DIM == 3) setForceRich( Mesh, Solution, bW, forceW, 2 );
    // set paralution force object
    for (int cl = 0; cl < Mesh.DOF[1]; cl++) {
      forceUP[cl] = bU[cl];
    }
    for (int cl = 0; cl < Mesh.DOF[2]; cl++) {
      forceVP[cl] = bV[cl];
    }
    if ( Mesh.DIM == 3 ) {
      for (int cl = 0; cl < Mesh.DOF[3]; cl++) {
        forceWP[cl] = bW[cl];
      }
    }

    // solve systems
    lsU.Solve(forceUP, &solU);
    lsV.Solve(forceVP, &solV);
    if (Mesh.DIM == 3) lsW.Solve(forceWP, &solW);

    // update the pressure vector and compute residual
    updatePressureRich( Mesh, Solution, solU, solV, solW, residual );

    if ( (residual < tolerance) || (richIt > richItMax) ) break;

  }
  // pass solution from paralution object to std vector input
  for (int cl = 0; cl < Mesh.DOF[1]; cl++) {
    Solution[cl] = solU[cl];
  }
  for (int cl = Mesh.DOF[1]; cl < (Mesh.DOF[1]+Mesh.DOF[2]); cl++)
  {
    Solution[cl] = solV[cl];
  }
  if (Mesh.DIM == 3) {
    for (int cl = (Mesh.DOF[1]+Mesh.DOF[2]); cl < (Mesh.DOF[1]+Mesh.DOF[2]+Mesh.DOF[3]); cl++) {
      Solution[cl] = solW[cl];
    }
  }

  // clear paralution objects
  lsU.Clear();
  matU.Clear();
  pU.Clear();
  forcePU.Clear();
  solU.Clear();

  lsV.Clear();
  matV.Clear();
  pV.Clear();
  forcePV.Clear();
  solV.Clear();

  lsW.Clear();
  matW.Clear();
  pW.Clear();
  forcePW.Clear();
  solW.Clear();

}

void
initPressure( const FluidMesh& Mesh, std::vector<double>& Solution )
{
  int shift, cl2;
  if (Mesh.DIM == 2) shift = Mesh.DOF[1] + Mesh.DOF[2];
  else shift = Mesh.DOF[1] + Mesh.DOF[2] + Mehs.DOF[3];
  for (int cl = 0; cl < Mesh.DOF[0]; cl++)
  {
    cl2 = cl + shift;
    Solution[ cl2 ] = 1 - (Mesh.xLim[0] - Mesh.PCellCenters[ idx( cl, 0, Mesh.DIM ) ]) / (Mesh.xLim[1] - Mesh.xLim[0]);
  }
}

void
setForceRich( const FluidMesh& Mesh, const std::vector<double>& Solution, std::vector<double>& b, \
                                                         std::vector<double>& force, int component )
{

}

void
updatePressureRich( const FluidMesh& Mesh, const std::vector<double>& Solution, \
                    paralution::LocalVector<double>& solU, paralution::LocalVector<double>& solV, \
                    paralution::LocalVector<double>& solW )
{

}

