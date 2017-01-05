#include "solve_mkl.hpp"

// 1d ->2d index
#define idx2(i, j, li) ((i * ldi) + j)

void
hgf::solve::mkl::init_solver(void)
{

}

void
hgf::solve::mkl::finalize_solver(void)
{

}

void
hgf::solve::mkl::solve(const parameters& par, \
  const std::vector< array_coo >& array, \
  const std::vector< double >& rhs, \
  std::vector<double>& solution)
{

  // Matrix data
  int n = (int)rhs.size();
  int *i_csr, *i_coo, *j_arr;
  double *a_arr;
  // values
  a_arr = (double *)malloc(array.size() * sizeof(double));
  i_coo = (int *)malloc(array.size() * sizeof(int));
  j_arr = (int *)malloc(array.size() * sizeof(int));
  i_csr = (int *)malloc((n + 1) * sizeof(int));

  for (int ii = 0; ii < array.size(); ii++) {
    i_coo[ii] = array[ii].i_index;
    j_arr[ii] = array[ii].j_index;
    a_arr[ii] = array[ii].value;
  }
  // create csr array
  i_csr[0] = 0;
  int rcount = 0;
  for (int ii = 1; ii < array.size(); ii++) {
    if (j_arr[ii] > rcount) {
      rcount++;
      i_csr[rcount] = ii;
    }
  }
  i_csr[n] = (int)array.size();

  // pardiso arguments

}

