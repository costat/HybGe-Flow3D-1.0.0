/* models header */
#ifndef _MODELS_H
#define _MODELS_H

/* To add a model include the entry point header here.
   A model should take an hgf parameter file and hgf msh class as inputs
   and handle the construction of a COO array describing the discrete system. 
   The model should also contain a pointer to a solution output
   to be handed to a solver class and be able to appropriately postprocess the solution. 
   If the model src is placed within HGFROOT/src/models/ cmake will include all src 
   and headers.*/

// HGF ships with stationary Stokes and a simple PoreNetwork model
#include "../src/models/stokes/model_stokes.hpp"
#include "../src/models/porenetwork/model_porenetwork.hpp"

#endif
