#ifndef _TYPES_H
#define _TYPES_H

#include <boost/filesystem.hpp>

// problem parameters
struct parameters
{
  double viscosity;
  int dimension;                                 // dimension
  int nx;                                        // x direction voxels in voxel geometry
  int ny;                                        // y direction voxels in voxel geometry
  int nz;                                        // z direction voxels in voxel geometry
  double length;
  double width;
  double height;
  std::vector< unsigned long > voxel_geometry;   // voxel geometry
  boost::filesystem::path problem_path;
  int verbose;                                   // defaults to 0
};

// vertexs
struct vertex
{
  double coords[3];                       // coordinates of the vertex
  int gnum;                               // global node number
};

// edges
struct edge
{
  int vns[2];                             // local vectex numbers
  int gnum;								                // global edge number
  int neighbor;                           // for 2d, who's across the face? -1 for a boundary face
  int bctype;                             // for 2d, boundary condition type; 0 - interior, 1 - dirichlet, 2 - neumann
};

struct face
{
  int vns[4];                             // local vertex numbers
  int gnum;                               // global face number
  int neighbor;                           // who's across the face? -1 for a boundary
  int bctype;                             // boundary condition type; 0 - interior, 1 - dirichlet, 2 - neumann
};

// hex (or quad) cells
struct qcell
{
  struct edge edg[12];                         // cell edge information
  struct face fac[6];                          // cell face information
  struct vertex vtx[8];                        // cell vertices information
  double dx;                                   // length in x direction
  double dy;                                   // length in y direction
  double dz;
};

struct degree_of_freedom
{
  int doftype;                   // 0 for interior, 1 for edge, 2 for face
  double coords[3];              // coordinates of the degree of freedom
  int cell_numbers[2];           // list of mesh cell numbers containing dof
  int neighbors[6];              // list of neighbors for the dof
};

struct dof_unstructured
{
  int doftype;
  double coords[3];
  std::vector<double> normals;
  std::vector<int> cell_numbers;
  std::vector<int> neighbors;
};

struct array_coo
{
  int i_index;
  int j_index;
  double value;
};

// sort struct for v velocity component in stokes / navier-stokes
struct byXbyY
{
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[0] < two.coords[0] || \
      (one.coords[0] == two.coords[0] && one.coords[1] < two.coords[1]));
  }
};

// sort struct for v velocity component in 3d stokes / navier-stokes
struct byZbyXbyY
{
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[2] < two.coords[2] || (one.coords[2] == two.coords[2] && one.coords[0] < two.coords[0]) \
      || (one.coords[2] == two.coords[2] && one.coords[0] == two.coords[0] && one.coords[1] < two.coords[1]));
  }
};

// sort struct for w velocity component in 3d stokes / navier-stokes
struct byYbyXbyZ
{
  bool operator()(degree_of_freedom const &one, degree_of_freedom const &two)
  {
    return (one.coords[1] < two.coords[1] || (one.coords[1] == two.coords[1] && one.coords[0] < two.coords[0]) \
      || (one.coords[1] == two.coords[1] && one.coords[0] == two.coords[0] && one.coords[2] < two.coords[2]));
  }
};

// sort struct for coo array
struct byIbyJ
{
  bool operator()(array_coo const &one, array_coo const &two)
  {
    return (one.i_index < two.i_index || \
      (one.i_index == two.i_index && one.j_index < two.j_index));
  }
};

// sort struct for coo array
struct byJbyI
{
  bool operator()(array_coo const &one, array_coo const &two)
  {
    return (one.j_index < two.j_index || \
      (one.j_index == two.j_index && one.i_index < two.i_index));
  }
};

// bc types and values
struct boundary_nodes
{
  int type;
  double value;
};

#endif
