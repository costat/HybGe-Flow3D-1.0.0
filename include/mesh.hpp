#ifndef _MESH_H
#define _MESH_H

// system includes
#include <vector>

namespace hgf
{
  /** \brief Contains functions and classes related to meshing multiscale flow problems.
   *
   */
  namespace mesh 
  {
    /** \brief Class creates quadrilateral and hexagonal meshes from voxel input files.
     *
     */
    class voxel
    {
      public:
        std::vector< qcell > els;                           /**< Vector of quadrilateral or hexagonal cells in the mesh. */
        std::vector< int > gtlNode;                         /**< Global to local node map. */
        void build( parameters& par);
        void printVTK(const parameters& par);
      private:
        void build_from_voxel_quad( parameters& par);
        void build_from_voxel_hex(parameters& par);
    };
  }
}

#endif
