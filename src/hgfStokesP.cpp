#include <vector>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <paralution.hpp>
#include <stdio.h>

// hgf includes
#include "hgfMeshCu.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

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

  set_omp_threads_paralution( Par.nThreads );

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
  ls.Init(Par.tolAbs, Par.tolRel, 1e8, Par.maxIt);
  ls.SetOperator(mat);
  ls.Verbose(2);
  ls.SetBasisSize(100);

  // preconditioning
/*  BlockPreconditioner<LocalMatrix<double>, LocalVector<double>, double> p;
  Solver<LocalMatrix<double>, LocalVector<double>, double> **p2;
  int n = Mesh.DIM;
  int *size;
  size = new int[n];
  size[0] = Mesh.DOF[1];
  size[1] = Mesh.DOF[2];
  if (Mesh.DIM == 3) {
    size[2] = Mesh.DOF[3];
  }

  p2 = new Solver<LocalMatrix<double>, LocalVector<double>, double> *[n];

  for (int i = 0; i < n; ++i) {
    ILU<LocalMatrix<double>, LocalVector<double>, double> *mc;
    mc = new ILU<LocalMatrix<double>, LocalVector<double>, double>;
    mc->Set(Par.prec);
    p2[i] = mc;
  }

  p.Set(n, size, p2);
*/
  ILU<LocalMatrix<double>, LocalVector<double>, double> p;
  p.Set(Par.prec);
  ls.SetPreconditioner(p);

  // build
  ls.Build();

  // solve the system
  ls.Solve(forceP, &sol);

  // pass solution from paralution object to std vector input
  for (int cl = 0; cl < Mesh.dofTotal; cl++) {
    Solution[cl] = sol[cl];
  }

  // clear paralution objects
  ls.Clear();
  mat.Clear();
  //p.Clear();
  forceP.Clear();
  sol.Clear();

}

void InitPressure( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par )
{
  int shift, cl2;
  if (Mesh.DIM == 2) shift = Mesh.DOF[1] + Mesh.DOF[2];
  else shift = Mesh.DOF[1] + Mesh.DOF[2] + Mesh.DOF[3];
  switch (Par.direction)
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
SetForceUZCG( const FluidMesh& Mesh, const std::vector<double>& Solution, LocalVector<double>& b, \
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

/* This function solves the stokes ib system with a krylov accelerated uzawa iteration scheme */
void StokesSolveUZCG( const FluidMesh& Mesh, std::vector<double>& Solution, const ProbParam& Par )
{

  //========= problems setup ============//

  // declarations
  std::vector<int> i1, i2, i3, j1, j2, j3;
  std::vector<double> val1, val2, val3, f1, f2, f3;
  int continueKrylov = 1;
  double res = Par.tolAbs + 1;

  i1.reserve( Mesh.DOF[1] * 7 );
  j1.reserve( Mesh.DOF[1] * 7 );
  val1.reserve( Mesh.DOF[1] * 7 );
  f1.resize( Mesh.DOF[1] );

  i2.reserve( Mesh.DOF[2] * 7 );
  j2.reserve( Mesh.DOF[2] * 7 );
  val2.reserve( Mesh.DOF[2] * 7 );
  f2.resize( Mesh.DOF[2] );

  if (Mesh.DIM == 3) {
    i3.reserve( Mesh.DOF[3] * 7 );
    j3.reserve( Mesh.DOF[3] * 7 );
    val3.reserve( Mesh.DOF[3] * 7 );
    f3.resize( Mesh.DOF[3] );
  }

  // build arrays
  // interior cells
  VelocityArray( Mesh, Par.visc, i1, j1, val1, 0 );
  VelocityArray( Mesh, Par.visc, i2, j2, val2, 1 );
  if (Mesh.DIM == 3) VelocityArray( Mesh, Par.visc, i3, j3, val3, 2 );

  // boundary conditions
  AxisFlowSingleComponent( Mesh, i1, j1, val1, f1, Par.visc, Par.direction, 0 );
  AxisFlowSingleComponent( Mesh, i2, j2, val2, f2, Par.visc, Par.direction, 1 );
  if (Mesh.DIM == 3) AxisFlowSingleComponent( Mesh, i3, j3, val3, f3, Par.visc, Par.direction, 2 );

  // immersed boundary
  immersedBoundarySingleComponent( Mesh, i1, j1, val1, 0 );
  immersedBoundarySingleComponent( Mesh, i2, j2, val2, 1 );
  if (Mesh.DIM == 3) immersedBoundarySingleComponent( Mesh, i3, j3, val3, 2 );

  //============ Paralution setup =============//

  set_omp_threads_paralution( Par.nThreads );

  // delcare paralution objects
  LocalVector<double> x1, x2, x3, y;
  LocalVector<double> b1, b2, b3, by;
  LocalMatrix<double> mat1, mat2, mat3;

  // initialize force and solution vectors
  b1.Allocate("force vector u", Mesh.DOF[1]);
  b1.Zeros();
  b2.Allocate("force vector v", Mesh.DOF[2]);
  b2.Zeros();
  if (Mesh.DIM == 3)
  x1.Allocate("solution u", Mesh.DOF[1]);
  x1.Zeros();
  x2.Allocate("solution v", Mesh.DOF[2]);
  x3.Zeros();
  if (Mesh.DIM == 3) {
    b3.Allocate("force vector w", Mesh.DOF[3]);
    b3.Zeros();
    x3.Allocate("solution w", Mesh.DOF[3]);
    x3.Zeros();
  }

  // assemble the paralution matrices from coo data
  mat1.Assemble( i1.data(), j1.data(), val1.data(), i1.size(), "u operator", Mesh.DOF[1], Mesh.DOF[1] );
  mat2.Assemble( i2.data(), j2.data(), val2.data(), i2.size(), "v operator", Mesh.DOF[2], Mesh.DOF[2] );
  if (Mesh.DIM == 3) mat3.Assemble( i3.data(), j3.data(), val3.data(), i3.size(), "w operator", Mesh.DOF[3], Mesh.DOF[3] );

  // solver objects
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls1;
  ls1.Init( Par.tolAbs, Par.tolRel, 1e8, Par.maxIt );
  ls1.SetOperator( mat1 );
  ls1.Verbose( 1 );
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls2;
  ls2.Init( Par.tolAbs, Par.tolRel, 1e8, Par.maxIt );
  ls2.SetOperator( mat2 );
  ls2.Verbose( 1 );
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls3;
  if (Mesh.DIM == 3) {
    ls3.Init( Par.tolAbs, Par.tolRel, 1e8, Par.maxIt );
    ls3.SetOperator( mat3 );
    ls3.Verbose( 1 );
  }

  // preconditioning
  ILU<LocalMatrix<double>, LocalVector<double>, double> p1;
  p1.Set(Par.prec);
  ls1.SetPreconditioner( p1 );

  ILU<LocalMatrix<double>, LocalVector<double>, double> p2;
  p2.Set(Par.prec);
  ls2.SetPreconditioner( p2 );

  ILU<LocalMatrix<double>, LocalVector<double>, double> p3;
  if (Mesh.DIM == 3) {
    p3.Set(Par.prec);
    ls3.SetPreconditioner( p3 );
  }

  // build
  ls1.Build();
  ls2.Build();
  if (Mesh.DIM == 3) ls3.Build();

  // set the initial pressure
  InitPressure( Mesh, Solution, Par );

  //========= krylov solve section =========//

  SolveMomentum :
    // set the force in the momentum equations
    SetForceUZCG( Mesh, Solution, b1, f1, 0 );
    SetForceUZCG( Mesh, Solution, b2, f2, 1 );
    if (Mesh.DIM == 3) SetForceUZCG( Mesh, Solution, b3, f3, 2 );

    // solve momentum equations



//    if (continueKrylov) goto KrylovDirection;
//    else goto cleanup;

//  KrylovDirection :

    // compute residual r

    // set p2

    //

//  NewPressure :


//  cleanup :

}

void SolverInit( void )
{
  init_paralution();
}

void SolverFinalize( void )
{
  stop_paralution();
}
