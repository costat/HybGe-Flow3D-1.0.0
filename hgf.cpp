#include <vector>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>

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
   with an immersed boundary. */
void
hgfStokesDrive( unsigned long *gridin, int size1, int ldi1, int ldi2, \
                int nx, int ny, int nz, \
                double length, double width, double height, int direction, \
                double visc, int nThreads )
{
  switch ( direction )
  {
    case 3 : // Solve all 3 flow directions for upscaled tensor
    {
      std::string outNameX = "flowrunX.dat";
      std::string outNameY = "flowrunY.dat";
      std::string outNameZ = "flowrunZ.dat";
      int dofTotal, maxNZ;
      double mesh_duration, array_duration, solve_duration, total_duration;
      std::clock_t start, array_start, solve_start;

      start = std::clock();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );

      // Mesh Timer
      mesh_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout << "\n Mesh constructed in " << mesh_duration << "seconds\n\n";
      array_start = std::clock();

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
          array_duration = ( std::clock() - array_start ) / (double) CLOCKS_PER_SEC;
          std::cout << " Arrays constructed in " << array_duration << "seconds\n\n";
          solve_start = std::clock();

          init_paralution();

          set_omp_threads_paralution(nThreads);

          // Prints threads and accelerator info to console
          info_paralution();

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

          // Pass the matrix and preconditioner to the solver.
          lsX.SetOperator(matX);
          lsX.SetPreconditioner(pX);
          lsY.SetOperator(matY);
          lsY.SetPreconditioner(pY);
          lsZ.SetOperator(matZ);
          lsZ.SetPreconditioner(pZ);

          // Build the solver.
          lsX.Build();
          lsY.Build();
          lsZ.Build();

          // Solve the problem
          lsX.Solve(forcePX, &solX);
          lsY.Solve(forcePY, &solY);
          lsZ.Solve(forcePZ, &solZ);

          // Linear solve timer
          solve_duration = ( std::clock() - solve_start ) / (double) CLOCKS_PER_SEC;
          std::cout << "\n Linear system solved in ";
          std::cout  << solve_duration << "seconds\n\n";

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

          stop_paralution();

          // Total timers
          total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
          std::cout << " Total time: " << total_duration << "seconds\n\n";
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
          array_duration = ( std::clock() - array_start ) / (double) CLOCKS_PER_SEC;
          std::cout << " Arrays constructed in " << array_duration << "seconds\n\n";
          solve_start = std::clock();

          init_paralution();

          set_omp_threads_paralution(nThreads);

          // Prints threads and accelerator info to console
          info_paralution();

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

          // Pass the matrix and preconditioner to the solver.
          lsX.SetOperator(matX);
          lsX.SetPreconditioner(pX);
          lsY.SetOperator(matY);
          lsY.SetPreconditioner(pY);

          // Build the solver.
          lsX.Build();
          lsY.Build();

          // Solve the problem
          lsX.Solve(forcePX, &solX);
          lsY.Solve(forcePY, &solY);

          // Linear solve timer
          solve_duration = ( std::clock() - solve_start ) / (double) CLOCKS_PER_SEC;
          std::cout << "\n Linear system solved in " << solve_duration;
          std::cout  << "seconds\n\n";

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

          stop_paralution();

          // Total timers
          total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
          std::cout << " Total time: " << total_duration << "seconds\n\n";
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
      std::clock_t start, array_start, solve_start;

      start = std::clock();

      // Build mesh object
      FluidMesh Mesh;
      Mesh.BuildUniformMesh( gridin, ldi1, ldi2, nx, ny, nz, length, width, height );

      // Mesh Timer
      mesh_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout << "\n Mesh constructed in " << mesh_duration << "seconds\n\n";
      array_start = std::clock();

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
      array_duration = ( std::clock() - array_start ) / (double) CLOCKS_PER_SEC;
      std::cout << " Arrays constructed in " << array_duration << "seconds\n\n";
      solve_start = std::clock();

      init_paralution();

      set_omp_threads_paralution(nThreads);

      // Prints threads and accelerator info to console
      info_paralution();

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

      // Setup a GMRES solver object and an ILU preconditioner.
      GMRES<LocalMatrix<double>, LocalVector<double>, double > ls;
      ILU<LocalMatrix<double>, LocalVector<double>, double> p;

      // Pass the matrix and preconditioner to the solver.
      ls.SetOperator(mat);
      ls.SetPreconditioner(p);

      // Build the solver.
      ls.Build();

      // Solve the problem
      ls.Solve(forceP, &sol);

      // Linear solve timer
      solve_duration = ( std::clock() - solve_start ) / (double) CLOCKS_PER_SEC;
      std::cout << "\n Linear system solved in " << solve_duration << "seconds\n\n";

      // Write solution to file
      writeSolutionL ( Mesh, sol, outName );
      computeKConstantDrive ( Mesh, sol, direction );

      // Clear arrays no longer in use.
      mat.Clear();
      forceP.Clear();
      sol.Clear();

      stop_paralution();

      // Total timers
      total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      std::cout << " Total time: " << total_duration << "seconds\n\n";
      break;
    }
  }
}
