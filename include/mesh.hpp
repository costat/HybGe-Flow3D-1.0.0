/* mesh header */
#ifndef _MESH_H
#define _MESH_H

// system includes
#include <vector>

namespace hgf
{
  class mesh
  {
    public:
      std::vector< qcell > els;
      std::vector< int > gtlNode;
      int order;
      void build( parameters& par);
    private:
      void build_from_voxel_quad( parameters& par);
      void build_from_voxel_hex(parameters& par);
      void printVTK(const parameters& par);
  };
}

#endif
