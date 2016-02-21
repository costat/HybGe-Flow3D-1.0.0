// hgfPoreNetwork.hpp
#ifndef HGFPORENETWORK_H
#define HGFPORENETWORK_H

#include <vector>
#include <paralution.hpp>

// hgf includes
#ifndef CUDA_BUILD
# define CUDA_BUILD 0
#endif

#if CUDA_BUILD
#include "hgfMeshCu.cuh"
#else
#include "hgfMesh.hpp"
#endif

void
PoreNetworkSolveDirect( const PoreNetwork& pn, const std::vector<double>& Ks, \
                        std::vector<double>& Solution, int direction );

#endif
