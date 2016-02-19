// hgfPoreNetwork.hpp
#ifndef HGFPORENETWORK_H
#define HGFPORENETWORK_H

#include <vector>
#include <paralution.hpp>

#include "hgfMesh.hpp"

void
PoreNetworkSolveDirect( const PoreNetwork& pn, const std::vector<double>& Ks, \
                        std::vector<double>& Solution, int direction );

#endif