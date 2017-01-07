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
hgf::solve::paralution::solve_ps_flow(const parameters& par, \
  const std::vector< array_coo >& array, \
  const std::vector< double >& rhs, \
  std::vector<double>& solution, \
  int n_u, int n_v, int n_w, int n_p)
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

  // GMRES object
  FGMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
  ls.Init(1e-8, 1e-6, 1e8, 500);
  ls.SetOperator(mat);
  ls.Verbose(2);
  ls.SetBasisSize(100);

  // preconditioning
  DiagJacobiSaddlePointPrecond<LocalMatrix<double>, LocalVector<double>, double> p;
  // Upper preconditioner is a block preconditioner broken up for velocity components
  BlockPreconditioner<LocalMatrix<double>, LocalVector<double>, double> p_k;
  Solver<LocalMatrix<double>, LocalVector<double>, double> **p2;
  int n = par.dimension;
  int *size;
  size = new int[n];
  size[0] = n_u;
  size[1] = n_v;
  if (par.dimension == 3) {
    size[2] = n_w;
  }
  p2 = new Solver<LocalMatrix<double>, LocalVector<double>, double> *[n];
  for (int i = 0; i < n; ++i) {
    MultiColoredILU<LocalMatrix<double>, LocalVector<double>, double> *mc;
    mc = new MultiColoredILU<LocalMatrix<double>, LocalVector<double>, double>;
    mc->Set(3);
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
  ls.Solve(force, &sol);

  int dof_total = (par.dimension == 3) ? n_u + n_v + n_w + n_p : n_u + n_v + n_p;
  // pass solution from paralution object to std vector input
  for (int cl = 0; cl < dof_total; cl++) {
    solution[cl] = sol[cl];
  }

  // clear paralution objects
  ls.Clear();
  mat.Clear();
  p.Clear();
  force.Clear();
  sol.Clear();

}

void
hgf::solve::paralution::finalize_solver(void)
{
  stop_paralution();
}
