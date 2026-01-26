 // ********************



// geometry.hh



// ********************



#ifndef __GEOMETRY_HH__

#define __GEOMETRY_HH__


// ********************
extern bool open_boundary_conditions;


inline int get_index(const int t, const int x, const int y, const int z, const int T, const int L) 
{
  int tt;
  if (!open_boundary_conditions) tt= (t+T)%T;  
  if (open_boundary_conditions) {
      if (t>=0 && t<T) tt= t; 
      else return -1;
  }
  int xx = (x+L)%L;
  int yy = (y+L)%L;
  int zz = (z+L)%L;

  return ((tt*L+xx)*L+yy)*L+zz;
}



inline int get_index_timeslice(const int x, const int y, const int z, const int T, const int L)
{
  int xx = (x+L)%L;
  int yy = (y+L)%L;
  int zz = (z+L)%L;

  return (xx*L+yy)*L+zz;
}

inline int get_index_timeslice_t(const int t, const int x, const int y, const int z, const int T, const int L)
{
  int xx = (x+L)%L;
  int yy = (y+L)%L;
  int zz = (z+L)%L;

  return (xx*L+yy)*L+zz;
}


inline int get_index_II(const int t, const int x, const int y, const int z, const int T, const int L, const int mu) 
{
  int tt;
  if (!open_boundary_conditions) tt= (t+T)%T;  
  if (open_boundary_conditions) {
      if (t>=0 && t<T) tt= t; 
      else return -1;
      if (t==T && mu==0) return -1; 
  }
  int xx = (x+L)%L;
  int yy = (y+L)%L;
  int zz = (z+L)%L;

  

  int ix = ((tt*L+xx)*L+yy)*L+zz;
  return (4*ix+mu) * 8;
}



inline int ggi(const int ix, const int mu)
{
  if (ix==-1) return -1;
  //if ()
  else return (4*ix+mu) * 8;
}



// ********************



#endif



// ********************
