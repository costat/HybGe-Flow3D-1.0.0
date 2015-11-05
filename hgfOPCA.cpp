#include <iostream>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <cmath>
#include <omp.h>
#include "hgfMesh.hpp"
#include "hgf.hpp"
#include "hgfArrays.hpp"
#include "hgfBC.hpp"
#include "hgfIB.hpp"
#include "hgfPP.hpp"

/* Define a 2d -> 1d array index, uses row major ordering */
#define idx2(i, j, ldi) ((i * ldi) + j)
/* Define a 3d -> 1d array index, uses row major ordering */
#define idx3(i, j, k, ldi1, ldi2) (k + (ldi2 * (j + ldi1 * i)))
