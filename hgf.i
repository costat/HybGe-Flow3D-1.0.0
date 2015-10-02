%module hgf

%{
  #define SWIG_FILE_WITH_INIT
  #include "hgf.hpp"
  #include "hgfMesh.hpp"
  #include "hgfArrays.hpp"
  #include "hgfBC.hpp"
  #include "hgfIB.hpp"
  #include "hgfPP.hpp"
%}

%include "numpy.i"

%init %{
  import_array();
%}

%apply (unsigned long *IN_ARRAY1, int DIM1) {(unsigned long *gridin, int size1)}

%include "hgf.hpp"
%include "hgfMesh.hpp"
%include "hgfArrays.hpp"
%include "hgfBC.hpp"
%include "hgfIB.hpp"
%include "hgfPP.hpp"

