#include "solve_paralution.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

using namespace paralution;

void
hgf::solve::paralution::init_solver(void)
{
  init_paralution();
}

void
hgf::solve::paralution::solve(const parameters& par, \
  const std::vector< array_coo >& array, \
  const std::vector< double >& rhs, std::vector< double >& solution)
{
  set_omp_threads_paralution(omp_get_max_threads());

  int *i_index, *j_index;
  double *value;

  i_index = (int *)malloc(array.size() * sizeof(int));
  j_index = (int *)malloc(array.size() * sizeof(int));
  value = (double *)malloc(array.size() * sizeof(double));

  for (int ii = 0; ii < array.size(); ii++) {
    i_index[ii] = array[ii].i_index;
    j_index[ii] = array[ii].j_index;
    value[ii] = array[ii].value;
  }

  LocalVector<double> sol;
  LocalVector<double> force;
  LocalMatrix<double> mat;

  force.Allocate("force vector", (int)rhs.size());
  for (int ii = 0; ii < rhs.size(); ii++) {
    force[ii] = rhs[ii];
  }

  sol.Allocate("solution", (int)rhs.size());
  sol.Zeros();

  mat.Assemble(i_index, j_index, value, (int)array.size(), \
    "operator", (int)rhs.size(), (int)rhs.size());

#ifdef _PARALUTION_MATRIX_DEBUG
  mat.WriteFileMTX("MatrixCheck.dat");
  force.WriteFileASCII("RHS.dat");
#endif

  GMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
  ls.Init(1e-8, 1e-6, 1e8, 1000);
  ls.SetOperator(mat);
  ls.Verbose(2);

  ILU<LocalMatrix<double>, LocalVector<double>, double> p;
  p.Set(2);
  ls.SetPreconditioner(p);

  ls.Build();
  ls.Solve(force, &sol);

#ifdef _PARALUTION_SOLUTION_DEBUG
  sol.WriteFileASCII("SOL.dat");
#endif

  for (int ii = 0; ii < solution.size(); ii++) {
    solution[ii] = sol[ii];
  }

  ls.Clear();
  mat.Clear();
  force.Clear();
  sol.Clear();

  free(i_index);
  free(j_index);
  free(value);

}

void
hgf::solve::paralution::finalize_solver(void)
{
  stop_paralution();
}
