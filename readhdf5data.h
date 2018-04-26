
#define H5_USE_16_API
#include <hdf5.h>
#include"typedefine.h"
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

/*============================================================*/
void readHDF5scalar(int idnum, int nvars,char *name, float *buf, hid_t groupid);
/* Read time stamp */
float readHDF5time(hid_t h5fid, int frame, hdatatype *hptr);

void read_plane_frame(hid_t h5fid,int frame,hdatatype *hdata);

void read_connection(hid_t h5fid,hdatatype *hdata);

void read_xyz(settype  settings,hid_t h5fid,int frame,hdatatype *hdata,double *parity ,hsize_t  *npoints);

void read_B(settype  settings,hid_t h5fid,int frame,hsize_t npoints,double parity,int eqflag,hdatatype *hdata);


