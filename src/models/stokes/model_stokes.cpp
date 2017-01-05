/* stokes main source */

// hgf includes
#include "model_stokes.hpp"

// 1d->2d index
#define idx2(i, j, ldi) ((i * ldi) + j)

void
hgf::models::stokes::build(const parameters& par, const hgf::mesh& msh)
{
  if (par.dimension == 2) {

    // setup the degrees of freedom
    build_degrees_of_freedom_2d(par, msh);

    // initialize solution and rhs
    int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
    int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
    int nP = (int)pressure.size();
    interior.resize(nU + nV + nP);
    rhs.resize(nU + nV + nP);

    // setup the linear system
    build_array_2d(par, msh);

  }
  else {

    // setup the degrees of freedom
    build_degrees_of_freedom_3d(par, msh);

    // initialize solution and rhs
    int nU = std::accumulate(interior_u.begin(), interior_u.end(), 0);
    int nV = std::accumulate(interior_v.begin(), interior_v.end(), 0);
    int nW = std::accumulate(interior_w.begin(), interior_w.end(), 0);
    int nP = (int)pressure.size();
    interior.resize(nU + nV + nW + nP);
    rhs.resize(nU + nV + nW + nP);

    // setup the linear system
    build_array_3d(par, msh);

  }
#ifdef _ARRAY_DEBUG
  std::cout << "\nArray size = " << coo_array.size() << "\n";
  std::cout << "\nnU = " << std::accumulate(interior_u.begin(), interior_u.end(), 0) << ",\tnV = " << std::accumulate(interior_v.begin(), interior_v.end(), 0);
  if (par.dimension == 3) {
    std::cout << ",\tNW = " << std::accumulate(interior_w.begin(), interior_w.end(), 0);
  }
  std::cout << ",\tnP = " << (int)pressure.size() << "\n";
  for (int ii = 0; ii < coo_array.size(); ii++) {
    std::cout << coo_array[ii].i_index << "\t" << coo_array[ii].j_index << "\t" << coo_array[ii].value << "\n";
  }
  std::cout << "\n";
#endif
}

void
hgf::models::stokes::setup_xflow_bc(const parameters& par, const hgf::mesh& msh)
{

  if (par.dimension == 2) xflow_2d(par, msh);
  else xflow_3d(par, msh);

}


