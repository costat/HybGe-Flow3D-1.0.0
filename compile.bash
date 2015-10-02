# build swig template and compile C++ modules
swig -c++ -python hgf.i
python hgfSetup.py build_ext --inplace
