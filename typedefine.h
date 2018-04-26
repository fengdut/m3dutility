/* Type definitions */
typedef struct H_datatype {
  int    ncells, npoints, planes[2], vpp;
  int    *neighbor_list, *conn_list;
  float  *x, *y, *z, *Bx, *By, *Bz, *BR, *Bphi, *psi;
  double *R, *ctab, *stab;
  double dphi;
  float tim;
  int nframes;

} hdatatype;


typedef struct D_edge{
  int  el0, v;
  char side;
} d_edge;

typedef struct Edge{
  d_edge o[6];
  char   n;
} edge;

typedef enum Output_type
{
	ASCII, HDF5
} 
output_type;

typedef struct Settype {
  char outname[128];  /* Output filename, if specified. */
  double convtol, aguess, bguess, theta_sym, wingspan;
  int  plane0;  /* Index of plane for Poincare plotting */
  int  frame;
  int  poindir; /* Direction for Poincare plotting (default=forward) */
  int  numsurf;  /* Number of surfaces to trace out */
  int  outerdots; /* Number of dots on outermost surface */
  int  mindots; /* Number of dots on innermost surface */
  int  qsurf, qtrans, maxtrans;
  int  oflag;       /* Output filename specified? */
  int  filtermode; /* Toroidal mode number for filtering */
  int  maxpass;
  int planes_per_cpu;
  output_type filetype;
} settype;
