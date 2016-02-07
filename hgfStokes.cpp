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

void
StokesSolve( const FluidMesh& Mesh, double visc, int direction, \
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

