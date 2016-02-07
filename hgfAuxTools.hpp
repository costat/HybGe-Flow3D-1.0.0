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

