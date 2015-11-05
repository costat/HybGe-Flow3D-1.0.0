// hgfOPCA.hpp
#ifndef OptPCA_H
#define OptPCA_H

#include<vector>

class OptPCA
{
  public:
    // Public data
    std::vector<unsigned long> newSnaps;
    int snapShotLength;
    // Public functions
    void genSnap( unsigned long* oldSnaps );
  private:
    int surfaceAreaFunc( std::vector<unsigned long> candidateSnap, const FluidMesh& sampleMesh );
};

#endif
