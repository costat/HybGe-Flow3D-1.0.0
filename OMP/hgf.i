%module hgf

%{
  #define SWIG_FILE_WITH_INIT
  #include "../src/hgf.hpp"
  #include "../src/hgfMesh.hpp"
  #include "../src/hgfStokes.hpp"
  #include "../src/hgfPoreNetwork.hpp"
  #include "../src/hgfArrays.hpp"
  #include "../src/hgfBC.hpp"
  #include "../src/hgfIB.hpp"
  #include "../src/hgfPP.hpp"
%}

%include "../src/numpy.i"

%init %{
  import_array();
%}

%apply (unsigned long *IN_ARRAY1, int DIM1) {(unsigned long *gridin, int size1)}

%include "../src/hgf.hpp"
%include "../src/hgfMesh.hpp"
%include "../src/hgfStokes.hpp"
%include "../src/hgfPoreNetwork.hpp"
%include "../src/hgfArrays.hpp"
%include "../src/hgfBC.hpp"
%include "../src/hgfIB.hpp"
%include "../src/hgfPP.hpp"
