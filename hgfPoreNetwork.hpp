#include <vector>
#include <paralution.hpp>

#include "hgfMesh.hpp"

void
PoreNetworkSolveDirect( const PoreNetwork& pn, const std::vector<double>& Ks, \
                        double& KPN, int direction );

