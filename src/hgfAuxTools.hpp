// hgfAuxTools.hpp
#ifndef HGFAUXTOOLS_H
#define HGFAUXTOOLS_H

#include <vector>

struct arrayCOO
{
  int I;
  int J;
  double Val;
};

struct byIbyJ
{
  bool operator()(arrayCOO const &one, arrayCOO const &two)
  {
    return ( one.I < two.I || (one.I == two.I && one.J < two.J ) );
  }
};

struct sortStruc2
{
  double xx;
  double yy;
  unsigned long ind;
};

struct sortStruc3
{
  double xx;
  double yy;
  double zz;
  unsigned long ind;
};

struct byXbyY
{
  bool operator()(sortStruc2 const &one, sortStruc2 const &two)
  {
    return ( one.xx < two.xx || (one.xx == two.xx && one.yy < two.yy) );
  }
};

struct byZbyXbyY
{
  bool operator()(sortStruc3 const &one, sortStruc3 const &two)
  {
    return ( one.zz < two.zz || (one.zz == two.zz && one.xx < two.xx) \
           || (one.zz == two.zz && one.xx == two.xx && one.yy < two.yy));
  }
};

struct byYbyXbyZ
{
  bool operator()(sortStruc3 const &one, sortStruc3 const &two)
  {
    return ( one.yy < two.yy || (one.yy == two.yy && one.xx < two.xx) \
           || (one.yy == two.yy && one.xx == two.xx && one.zz < two.zz) );
  }
};

struct ProbParam
{
  int nx, ny, nz, nThreads, prec, direction, solver, output;
  double length, width, height, visc, relax;
  bool isMesh;
  std::vector<unsigned long> gridin;
};

#endif
