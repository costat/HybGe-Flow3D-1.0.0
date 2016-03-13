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
  FGMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
  ls.Init(Par.tolAbs, Par.tolRel, 1e8, Par.maxIt);
  ls.SetOperator(mat);
  ls.Verbose(1);
  ls.SetBasisSize(100);

  // preconditioning
  DiagJacobiSaddlePointPrecond<LocalMatrix<double>, LocalVector<double>, double> p;
  // Upper preconditioner is a block preconditioner broken up for velocity components
  BlockPreconditioner<LocalMatrix<double>, LocalVector<double>, double> p_k;
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
    MultiColoredILU<LocalMatrix<double>, LocalVector<double>, double> *mc;
    mc = new MultiColoredILU<LocalMatrix<double>, LocalVector<double>, double>;
    mc->Set(Par.prec);
    p2[i] = mc;
  }
  p_k.Set(n, size, p2);

  // lower preconditioner
  Jacobi<LocalMatrix<double>, LocalVector<double>, double> p_s;
  p.Set(p_k, p_s);

//  ILU<LocalMatrix<double>, LocalVector<double>, double> p;
//  p.Set(Par.prec);
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
  p.Clear();
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

  // declarations for momentum solves
  std::vector<int> i1, i2, i3, j1, j2, j3;
  std::vector<double> val1, val2, val3, f1, f2, f3;

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

  // build momentum arrays
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

  // declarations for SpMVs; \grad^T matrix and grad_x/y/z matrices
  std::vector<int> bti, btj, bi1, bj1, bi2, bj2, bi3, bj3;
  std::vector<double> btval, bval1, bval2, bval3;

  bti.reserve( 6 * Mesh.DOF[0] );
  btj.reserve( 6 * Mesh.DOF[0] );
  btval.reserve( 6 * Mesh.DOF[0] );

  bi1.reserve( 2 * Mesh.DOF[1] );
  bj1.reserve( 2 * Mesh.DOF[1] );
  bval1.reserve( 2 * Mesh.DOF[1] );

  bi2.reserve( 2 * Mesh.DOF[2] );
  bj2.reserve( 2 * Mesh.DOF[2] );
  bval2.reserve( 2 * Mesh.DOF[2] );

  if (Mesh.DIM == 3) {
    bi3.reserve( 2 * Mesh.DOF[3] );
    bj3.reserve( 2 * Mesh.DOF[3] );
    bval3.reserve( 2 * Mesh.DOF[3] );
  }

  GradientTranspose( Mesh, bti, btj, btval );
  Gradient( Mesh, bi1, bj1, bval1, 0 );
  Gradient( Mesh, bi2, bj2, bval2, 1 );
  if (Mesh.DIM == 3) Gradient( Mesh, bi3, bj3, bval3, 2 );

  //============ Paralution setup =============//

  set_omp_threads_paralution( Par.nThreads );

  // paralution objects for momentum equations
  LocalVector<double> x1, x2, x3, y;
  LocalVector<double> b1, b2, b3;
  LocalMatrix<double> mat1, mat2, mat3;

  // initialize force and solution vectors
  b1.Allocate("force vector u", Mesh.DOF[1]);
  b1.Zeros();
  b2.Allocate("force vector v", Mesh.DOF[2]);
  b2.Zeros();
  x1.Allocate("solution u", Mesh.DOF[1]);
  x1.Zeros();
  x2.Allocate("solution v", Mesh.DOF[2]);
  x2.Zeros();
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
  ls1.Verbose( 2 );
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls2;
  ls2.Init( Par.tolAbs, Par.tolRel, 1e8, Par.maxIt );
  ls2.SetOperator( mat2 );
  ls2.Verbose( 2 );
  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls3;
  if (Mesh.DIM == 3) {
    ls3.Init( Par.tolAbs, Par.tolRel, 1e8, Par.maxIt );
    ls3.SetOperator( mat3 );
    ls3.Verbose( 2 );
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

  // paralution objects for Krylov steps
  LocalVector<double> r, pone1, pone2, pone3, poneall, ptwo, vall, ptp1, ptp2, ptp3, atwo, resvec;
  LocalMatrix<double> gmat1, gmat2, gmat3, gtmat;

  int vdof = Mesh.DOF[1] + Mesh.DOF[2];
  if (Mesh.DIM == 3) vdof += Mesh.DOF[3];

  r.Allocate("residual vector", Mesh.DOF[0]);
  r.Zeros();
  pone1.Allocate("p upper vector 1 in cg", Mesh.DOF[1]);
  pone1.Zeros();
  pone2.Allocate("p upper vector 2 in cg", Mesh.DOF[2]);
  pone2.Zeros();
  if (Mesh.DIM == 3) {
    pone3.Allocate("p upper vector 3 in cg", Mesh.DOF[3]);
    pone3.Zeros();
  }
  poneall.Allocate("p upper full", vdof);
  ptwo.Allocate("p vector in cg", Mesh.DOF[0]);
  ptwo.Zeros();
  vall.Allocate("full velocity", vdof);
  vall.Zeros();
  ptp1.Allocate("from p2 to p1", Mesh.DOF[1]);
  ptp1.Zeros();
  ptp2.Allocate("from p2 to p1", Mesh.DOF[2]);
  ptp2.Zeros();
  if (Mesh.DIM == 3) {
    ptp3.Allocate("from p2 to p1", Mesh.DOF[3]);
    ptp3.Zeros();
  }
  atwo.Allocate("atwo", Mesh.DOF[0]);
  atwo.Zeros();

  gmat1.Assemble( bi1.data(), bj1.data(), bval1.data(), bi1.size(), "grad x operator", Mesh.DOF[1], Mesh.DOF[0] );
  gmat2.Assemble( bi2.data(), bj2.data(), bval2.data(), bi2.size(), "grad y operator", Mesh.DOF[2], Mesh.DOF[0] );
  if (Mesh.DIM == 3) gmat3.Assemble( bi3.data(), bj3.data(), bval3.data(), bi3.size(), "grad z operator", Mesh.DOF[3], Mesh.DOF[0] );
  gtmat.Assemble( bti.data(), btj.data(), btval.data(), bti.size(), "grad transpose operator", Mesh.DOF[0], vdof );

  // set the initial pressure
  InitPressure( Mesh, Solution, Par );
  y.Allocate("pressure para", Mesh.DOF[0]);
  y.CopyFromData( &Solution[vdof] );
  resvec.Allocate("residual vector", Mesh.DOF[0]);
  resvec.Zeros();

  std::cout << "\nSolving model with CG-Uzawa Scheme...\n";

  //=========================================//
  //========= krylov uzawa section =========//
  //=========================================//
  int continueKrylov = 1;
  int kit = 1;
  double res = Par.tolAbs + 1;
  double alpha, beta;

  SolveMomentum :
  {
    // set the force in the momentum equations
    SetForceUZCG( Mesh, Solution, b1, f1, 0 );
    SetForceUZCG( Mesh, Solution, b2, f2, 1 );
    if (Mesh.DIM == 3) SetForceUZCG( Mesh, Solution, b3, f3, 2 );

    // solve momentum equations
    ls1.Solve( b1, &x1 );
    ls2.Solve( b2, &x2 );
    if (Mesh.DIM == 3) ls3.Solve( b3, &x3 );
    goto KryOne;
  }

  KryOne :
  {
    vall.CopyFrom( x1, 0, 0, Mesh.DOF[1] );
    vall.CopyFrom( x2, 0, Mesh.DOF[1], Mesh.DOF[2] );
    if (Mesh.DIM == 3) vall.CopyFrom( x3, 0, (Mesh.DOF[1] + Mesh.DOF[2]), Mesh.DOF[3] );
    // compute residual r = \grad^T * [u,v,w]^T
    gtmat.Apply( vall, &r );
    ptwo.CopyFrom( r );
    goto KryIt;
  }

  KryNotOne :
  {
    kit++;
    beta = r.Dot( atwo )/ptwo.Dot( atwo );
    ptwo.ScaleAddScale( (-beta), r, 1 );
    goto KryIt;
  }

  KryIt :
  {
    gmat1.Apply( ptwo, &ptp1 );
    gmat2.Apply( ptwo, &ptp2 );
    if (Mesh.DIM == 3 ) gmat3.Apply( ptwo, &ptp3 );
    ls1.Solve( ptp1, &pone1 );
    ls2.Solve( ptp2, &pone2 );
    if (Mesh.DIM == 3) ls3.Solve( ptp3, &pone3 );
    poneall.CopyFrom( pone1, 0, 0, Mesh.DOF[1] );
    poneall.CopyFrom( pone2, 0, Mesh.DOF[1], Mesh.DOF[2] );
    if (Mesh.DIM == 3) poneall.CopyFrom( pone3, 0, (Mesh.DOF[1]+Mesh.DOF[2]), Mesh.DOF[3] );
    gtmat.Apply( poneall, &atwo );
    alpha = ptwo.Dot( r ) / ptwo.Dot( atwo );

    // updates
    y.AddScale( ptwo, alpha );
    r.AddScale( atwo, (-alpha) );
    vall.AddScale( poneall, (-alpha) );

    // check continuity eq residual
    gtmat.Apply( vall, &resvec );
    res = resvec.Norm();
    std::cout << "\nKrylov-Uzawa Residual " << res << " at iteration " << kit << "\n";
    if (res < Par.tolAbs) continueKrylov = 0;

    // end or advance an iteration
    if (continueKrylov) goto KryNotOne;
    else goto cleanup;
  }

  cleanup :
    // get solution data to Solution vector
    for (int ii = 0; ii < vdof; ii++) {
      Solution[ii] = vall[ii];
    }
    for (int ii = 0; ii < Mesh.DOF[0]; ii++) {
      Solution[vdof+ii] = y[ii];
    }
    // paralution clears
    // momentum solve objects
    ls1.Clear();
    mat1.Clear();
    p1.Clear();
    b1.Clear();
    x1.Clear();
    ls2.Clear();
    mat2.Clear();
    p2.Clear();
    b2.Clear();
    x2.Clear();
    ls3.Clear();
    mat3.Clear();
    p3.Clear();
    b3.Clear();
    x3.Clear();
    // krylov it objects
    gtmat.Clear();
    r.Clear();
    pone1.Clear();
    pone2.Clear();
    pone3.Clear();
    poneall.Clear();
    ptwo.Clear();
    vall.Clear();
    atwo.Clear();
    resvec.Clear();
    gmat1.Clear();
    gmat2.Clear();
    gmat3.Clear();
}

void SolverInit( void )
{
  init_paralution();
}

void SolverFinalize( void )
{
  stop_paralution();
}
