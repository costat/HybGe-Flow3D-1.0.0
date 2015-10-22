#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include <paralution.hpp>

#include "hgf.hpp"
#include "hgfMesh.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index,
   uses row major indexing */
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

/* hgfStokesDrive solves the Stokes equations on a cartesian grid
   with an immersed boundary. The model is solved in 2d or 3d, either
   once for upscaling a constant conductivity, or 2/3 times for
   upscaling a conductivity tensor. */
void
hgfStokesDrive( unsigned long *gridin, int size1, int ldi1, int ldi2, \
                int nx, int ny, int nz, \
                double length, double width, double height, int direction, \
                double visc, int nThreads, int prec, int numSims, int simNum )
{

  if (simNum == 1) init_paralution();

  switch ( direction )
  {
    case 3 : // Solve all 3 flow directions for upscaled tensor
    {
      std::string outNameX = "flowrunX.dat";
      std::string outNameY = "flowrunY.dat";
      std::string outNameZ = "flowrunZ.dat";
      int dofTotal, maxNZ;
      double mesh_duration, array_duration, solve_duration, total_duration;
      double start, array_start, solve_start;

      start = omp_get_wtime();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );

      // Mesh Timer
      mesh_duration = ( omp_get_wtime() - start );
      array_start = omp_get_wtime();

      // Compute total degrees of freedom and maximum # of nonzero values in matrix.
      dofTotal = Mesh.TotalDOF();
      maxNZ = Mesh.MaxNonZero();

      /* The set matIs, matJs, and matVals are
         used as a COO sparse representation of the linear system. */
      std::vector<int> matIsX, matJsX;
      std::vector<double> matValsX, forceX;
      matIsX.reserve(maxNZ);
      matJsX.reserve(maxNZ);
      matValsX.reserve(maxNZ);
      forceX.resize(dofTotal);

      /* Compute vectors matVals, matIs, matJs
         describing values and i/j locations for matrix assembly.
         This call does not include boundary conditions or
         immersed boundary, thus, this step can be performed once
         for any number of realizations or boundary setups.  */
      StokesArray( Mesh, visc, matIsX, matJsX, matValsX );

      // Immersed Boundary
      immersedBoundary ( Mesh, matIsX, matJsX, matValsX );

      switch ( Mesh.DIM )
      {
        case 3 :
        {
          // Copy pre-boundary condition vectors to Y flow and Z flow versions
          std::vector<int> matIsY = matIsX;
          std::vector<int> matIsZ = matIsX;
          std::vector<int> matJsY = matJsX;
          std::vector<int> matJsZ = matJsX;
          std::vector<double> matValsY = matValsX;
          std::vector<double> matValsZ = matValsX;
          std::vector<double> forceY, forceZ;
          forceY.resize( dofTotal );
          forceZ.resize( dofTotal );

          // Boundary Conditions & Force
          BoundaryConditions( \
            Mesh, visc, matIsX, matJsX, matValsX, forceX, 0 );
          BoundaryConditions( \
            Mesh, visc, matIsY, matJsY, matValsY, forceY, 1 );
          BoundaryConditions( \
            Mesh, visc, matIsZ, matJsZ, matValsZ, forceZ, 2 );

          // Array timer
          array_duration = ( omp_get_wtime() - array_start );
          solve_start = omp_get_wtime();

          set_omp_threads_paralution(nThreads);

          // Declare paralution vector and matrix objects
          LocalVector<double> solX;
          LocalVector<double> forcePX;
          LocalMatrix<double> matX;
          LocalVector<double> solY;
          LocalVector<double> forcePY;
          LocalMatrix<double> matY;
          LocalVector<double> solZ;
          LocalVector<double> forcePZ;
          LocalMatrix<double> matZ;

          // Initialize force and solution vectors
          forcePX.Allocate("force vector", dofTotal);
          forcePX.Zeros();
          solX.Allocate("solution", dofTotal);
          solX.Zeros();
          forcePY.Allocate("force vector", dofTotal);
          forcePY.Zeros();
          solY.Allocate("solution", dofTotal);
          solY.Zeros();
          forcePZ.Allocate("force vector", dofTotal);
          forcePZ.Zeros();
          solZ.Allocate("solution", dofTotal);
          solZ.Zeros();

          // Assemble paralution arrays from previously built COO arrays.
          matX.Assemble(&matIsX[0], &matJsX[0], &matValsX[0], matIsX.size(), \
                        "operator", dofTotal, dofTotal);
          for (int cl = 0; cl < dofTotal; cl++) {
            forcePX[cl] = forceX[cl];
          }
          matY.Assemble(&matIsY[0], &matJsY[0], &matValsY[0], matIsY.size(), \
                        "operator", dofTotal, dofTotal);
          for (int cl = 0; cl < dofTotal; cl++) {
            forcePY[cl] = forceY[cl];
          }
          matZ.Assemble(&matIsZ[0], &matJsZ[0], &matValsZ[0], matIsZ.size(), \
                        "operator", dofTotal, dofTotal);
          for (int cl = 0; cl < dofTotal; cl++) {
            forcePZ[cl] = forceZ[cl];
          }

          // Setup a GMRES solver object and an ILU preconditioner.
          GMRES<LocalMatrix<double>, LocalVector<double>, double > lsX;
          ILU<LocalMatrix<double>, LocalVector<double>, double> pX;
          GMRES<LocalMatrix<double>, LocalVector<double>, double > lsY;
          ILU<LocalMatrix<double>, LocalVector<double>, double> pY;
          GMRES<LocalMatrix<double>, LocalVector<double>, double > lsZ;
          ILU<LocalMatrix<double>, LocalVector<double>, double> pZ;

          // ILU Level
          pX.Set(prec);
          pY.Set(prec);
          pZ.Set(prec);

          // Pass the matrix and preconditioner to the solver.
          lsX.SetOperator(matX);
          lsX.SetPreconditioner(pX);
          lsY.SetOperator(matY);
          lsY.SetPreconditioner(pY);
          lsZ.SetOperator(matZ);
          lsZ.SetPreconditioner(pZ);
          lsX.Verbose(2);
          lsY.Verbose(2);
          lsZ.Verbose(2);

          // Build the solver.
          lsX.Build();
          lsY.Build();
          lsZ.Build();

          // Solve the problem
          lsX.Solve(forcePX, &solX);
          lsY.Solve(forcePY, &solY);
          lsZ.Solve(forcePZ, &solZ);

          // Linear solve timer
          solve_duration = ( omp_get_wtime() - solve_start );

          // Write solution to file
          writeSolutionL ( Mesh, solX, outNameX );
          writeSolutionL ( Mesh, solY, outNameY );
          writeSolutionL ( Mesh, solZ, outNameZ );
          computeKTensorL ( Mesh, solX, solY, solZ );

          // Clear arrays no longer in use.
          matX.Clear();
          forcePX.Clear();
          solX.Clear();
          matY.Clear();
          forcePY.Clear();
          solY.Clear();
          matZ.Clear();
          forcePZ.Clear();
          solZ.Clear();

          std::cout << "Arrays constructed in " << array_duration << "seconds\n";
          std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
          std::cout << "Linear system solved in ";
          std::cout  << solve_duration << "seconds\n";
          // Total timers
          total_duration = ( omp_get_wtime() - start );
          std::cout << "Total time: " << total_duration << "seconds\n";
          break;
        }
        case 2 :
        {
          // Copy pre-boundary condition vectors to Y flow and Z flow versions
          std::vector<int> matIsY = matIsX;
          std::vector<int> matJsY = matJsX;
          std::vector<double> matValsY = matValsX;
          std::vector<double> forceY;
          forceY.resize( dofTotal );
          // Boundary Conditions & Force
          BoundaryConditions( \
            Mesh, visc, matIsX, matJsX, matValsX, forceX, 0 );
          BoundaryConditions( \
            Mesh, visc, matIsY, matJsY, matValsY, forceY, 1 );

          // Array timer
          array_duration = ( omp_get_wtime() - array_start );
          solve_start = omp_get_wtime();

          set_omp_threads_paralution(nThreads);

          // Declare paralution vector and matrix objects
          LocalVector<double> solX;
          LocalVector<double> forcePX;
          LocalMatrix<double> matX;
          LocalVector<double> solY;
          LocalVector<double> forcePY;
          LocalMatrix<double> matY;
          LocalVector<double> solZ;

          // Initialize force and solution vectors
          forcePX.Allocate("force vector", dofTotal);
          forcePX.Zeros();
          solX.Allocate("solution", dofTotal);
          solX.Zeros();
          forcePY.Allocate("force vector", dofTotal);
          forcePY.Zeros();
          solY.Allocate("solution", dofTotal);
          solY.Zeros();
          solZ.Allocate("solutin", dofTotal);
          solZ.Zeros();

          // Assemble paralution arrays from previously built COO arrays.
          matX.Assemble(&matIsX[0], &matJsX[0], &matValsX[0], matIsX.size(), \
                        "operator", dofTotal, dofTotal);
          for (int cl = 0; cl < dofTotal; cl++) {
            forcePX[cl] = forceX[cl];
          }
          matY.Assemble(&matIsY[0], &matJsY[0], &matValsY[0], matIsY.size(), \
                        "operator", dofTotal, dofTotal);
          for (int cl = 0; cl < dofTotal; cl++) {
            forcePY[cl] = forceY[cl];
          }

          // Setup a GMRES solver object and an ILU preconditioner.
          GMRES<LocalMatrix<double>, LocalVector<double>, double > lsX;
          ILU<LocalMatrix<double>, LocalVector<double>, double> pX;
          GMRES<LocalMatrix<double>, LocalVector<double>, double > lsY;
          ILU<LocalMatrix<double>, LocalVector<double>, double> pY;

          // ILU Level
          pX.Set(prec);
          pY.Set(prec);

          // Pass the matrix and preconditioner to the solver.
          lsX.SetOperator(matX);
          lsX.SetPreconditioner(pX);
          lsY.SetOperator(matY);
          lsY.SetPreconditioner(pY);
          lsX.Verbose(2);
          lsY.Verbose(2);

          // Build the solver.
          lsX.Build();
          lsY.Build();

          // Solve the problem
          lsX.Solve(forcePX, &solX);
          lsY.Solve(forcePY, &solY);

          // Linear solve timer
          solve_duration = ( omp_get_wtime() - solve_start );

          // Write solution to file
          writeSolutionL ( Mesh, solX, outNameX );
          writeSolutionL ( Mesh, solY, outNameY );
          computeKTensorL ( Mesh, solX, solY, solZ );

          // Clear arrays no longer in use.
          matX.Clear();
          forcePX.Clear();
          solX.Clear();
          matY.Clear();
          forcePY.Clear();
          solY.Clear();
          solZ.Clear();

          std::cout << "Arrays constructed in " << array_duration << "seconds\n";
          std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
          std::cout << "Linear system solved in " << solve_duration;
          std::cout  << "seconds\n";
          // Total timers
          total_duration = ( omp_get_wtime() - start );
          std::cout << "Total time: " << total_duration << "seconds\n";
          break;
        }
      }
      break;
    }
    default : // Solve a single flow direction, upscale constant conductivity
    {
      std::string outName = "flowrun.dat";
      int dofTotal, maxNZ;
      double mesh_duration, array_duration, solve_duration, total_duration;
      double start, array_start, solve_start;

      start = omp_get_wtime();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );
      // Mesh Timer
      mesh_duration = ( omp_get_wtime() - start );
      array_start = omp_get_wtime();

      // Compute total degrees of freedom and maximum # of nonzero values in matrix.
      dofTotal = Mesh.TotalDOF();
      maxNZ = Mesh.MaxNonZero();

      /* The set matIs, matJs, and matVals are
         used as a COO sparse representation of the linear system. */
      std::vector<int> matIs, matJs;
      std::vector<double> matVals, force;
      matIs.reserve(maxNZ);
      matJs.reserve(maxNZ);
      matVals.reserve(maxNZ);
      force.resize(dofTotal);

      /* Compute vectors matVals, matIs, matJs
         describing values and i/j locations for matrix assembly.
         This call does not include boundary conditions or
         immersed boundary, thus, this step can be performed once
         for any number of realizations or boundary setups.  */
      StokesArray( Mesh, visc, matIs, matJs, matVals );

      // Boundary Conditions & Force
      BoundaryConditions( Mesh, visc, matIs, matJs, matVals, \
                                  force, direction );

      // Immersed Boundary
      immersedBoundary ( Mesh, matIs, matJs, matVals );

      // Array timer
      array_duration = ( omp_get_wtime() - array_start );
      solve_start = omp_get_wtime();

      set_omp_threads_paralution(nThreads);

      // Declare paralution vector and matrix objects
      LocalVector<double> sol;
      LocalVector<double> forceP;
      LocalMatrix<double> mat;

      // Initialize force and solution vectors
      forceP.Allocate("force vector", dofTotal);
      forceP.Zeros();
      sol.Allocate("solution", dofTotal);
      sol.Zeros();

      // Assemble paralution arrays from previously built COO arrays.
      mat.Assemble(&matIs[0], &matJs[0], &matVals[0], matIs.size(), \
                   "operator", dofTotal, dofTotal);

      for (int cl = 0; cl < dofTotal; cl++) {
        forceP[cl] = force[cl];
      }

      // Setup a GMRES solver object.
      GMRES<LocalMatrix<double>, LocalVector<double>, double > ls;
      ls.Init(1e-6, 1e-6, 1e8, 1500);
      ls.SetOperator(mat);
      ls.Verbose(2);

      ILU<LocalMatrix<double>, LocalVector<double>, double> p;
      p.Set(prec);
      ls.SetPreconditioner(p);

      // Build the solver.
      ls.Build();

      // Solve the problem
      ls.Solve(forceP, &sol);

      // Linear solve timer
      solve_duration = ( omp_get_wtime() - solve_start );

      // Write solution to file
      writeSolutionL ( Mesh, sol, outName );
      computeKConstantDrive ( Mesh, sol, direction );

      // Clear arrays no longer in use.
      ls.Clear();
      mat.Clear();
      p.Clear();
      forceP.Clear();
      sol.Clear();

      std::cout << "Mesh constructed in " << mesh_duration << "seconds\n";
      std::cout << "Arrays constructed in " << array_duration << "seconds\n";
      std::cout << "Linear system solved in " << solve_duration << "seconds\n";
      // Total timers
      total_duration = ( omp_get_wtime() - start );
      std::cout << "Total time: " << total_duration << "seconds\n";
      break;
    }
  }

  if (simNum == numSims) stop_paralution();
}
