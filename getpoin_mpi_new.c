/* To compile on Sunfire: */
/* module load lff95; module load hdf5 */
/* gcc -O3 -DH5_USE_16_API -DQPROFILE ng12counter.c -L${HDF5_HOME}/lib -lhdf5 -lm -lz -lpthread */
/*  or  */
/* gcc -O3 -DH5_USE_16_API -DPOINCARE ng12counter.c -L${HDF5_HOME}/lib -lhdf5 -lm -lz -lpthread */

#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define H5_USE_16_API
#include <hdf5.h>

#include<time.h>
#include"typedefine.h"


#define SAFETY 0.9
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define PGROW -0.2
//#define USE_MPI

#ifdef USE_MPI
#include"mpi.h"
#endif

/* Prototypes */
void   poincare(hdatatype *, settype *), qprofile(hdatatype *, settype *);
int    nrrkqs(double *, double *, double, double, double *, double *, double *,hdatatype *);
int    nrrkck(double *, double *, double, double *, double *, hdatatype *);
int    getBfield(double *, double *, hdatatype *);
int    get_tri_weights(double *, int, int *, double *, hdatatype *);
void   findElementNeighbors(hdatatype *);
void   add_edge(edge *, int *, char, int, int *);
double find_axis(hdatatype *, settype *);
double get_excursion(double, hdatatype *, settype *, double *, double *);
void   readHDF5scalar(int,int, char *, float *, hid_t);
void   filter(float *, int, hdatatype *);
float  readHDF5time(hid_t, int, hdatatype *);
settype parameter_input();
void printhelp(char *argv[],int frame,settype settings);
void getoption(int argc, char *argv[],char *fname,int *eqflag,int * frame,settype * settings);
void read_plane_frame(hid_t h5fid,int frame,hdatatype *hdata);
void read_connection(hid_t h5fid,hdatatype *hdata);
void read_xyz(settype  settings,hid_t h5fid,int frame,hdatatype *hdata,double *parity ,hsize_t  *npoints);
void read_B(settype  settings,hid_t h5fid,int frame,hsize_t npoints,double parity,int eqflag,hdatatype *hdata);
/*============================================================*/
int main(int argc, char *argv[])
{
	hid_t h5fid;
	hsize_t       npoints;

	double        parity = 1.0;
	int            frame=0,  eqflag=0;
	char          fname[64];
	hdatatype     hdata;
	settype       settings;
	clock_t start, end;
	double cpu_time_used;
	int myid,nump;
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nump);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
	myid=0;
	nump=1;
#endif

	settings=parameter_input();

	if (argc < 2) 
	{ 
		if(myid==0)
			printhelp(argv,frame,settings);
		exit(0);
	}

	getoption(argc,argv,fname,&eqflag,&frame,&settings);

	/* Open the HDF5 file, read-only */
	if ((h5fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) 
	{
		fprintf(stderr, "my id %d, Could not access HDF5 file %s.\n",myid, fname);
		return 1;
	}
	/* Read number of planes, planes/processor */
	/* Read number of time slices in file */
	read_plane_frame(h5fid,frame,&hdata);
	/* Read connectivity data */
	read_connection(h5fid,&hdata);
	read_xyz(settings,h5fid,frame,&hdata,&parity ,&npoints);

	read_B(settings,h5fid,frame,npoints,parity,eqflag,&hdata);

	start	= clock();
	poincare(&hdata, &settings);
	end		= clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

	if(myid==0)
		printf("poincare time:%f \n",cpu_time_used);


	H5Fclose(h5fid);

	/* Deallocate the dynamic arrays */
	free(hdata.Bz);
	free(hdata.BR);
	free(hdata.Bphi);

	free(hdata.stab-1); 
	free(hdata.ctab-1);
	free(hdata.R); 
	free(hdata.z);
	free(hdata.conn_list);
#ifdef USE_MPI
	MPI_Finalize();
#endif
  return 0;
}


/*============================================================*/
void poincare(hdatatype *hptr, settype *setptr)
{
	double stepfrac = 0.001, eps;
	FILE   *fp;
	int    fillfactor = setptr->outerdots/setptr->numsurf;
	int    dot, rad, oldplane, dotspersurf, newplane, vv[3];
	float  Rmaj, Zmid, Redg, Rinn, minibuf[2];
	float	 *buffer,*buffersum;
	int	 *bufindex;
	double rmstep, rmstep2, h,th, hdid,thdid, hnext,thnext, hmax, ofac;
	double x0[3], dx[3], txscal[3],xscal[3], xold[3];
	int i=0;

	hid_t    h5file, h5fspace, h5dataset;
	hsize_t  cursize[2];
	hssize_t start[2];
	int myid,nump;
#ifdef USE_MPI
	MPI_Comm_size(MPI_COMM_WORLD,&nump);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
	myid=0;
	nump=1;
#endif
	if(myid==0)
		fprintf(stderr,"%d surfaces; %d dots on outermost surface; %d on innermost.\n", setptr->numsurf, setptr->outerdots, setptr->mindots);

	/* Find R at edge on midplane in target poloidal plane */
	eps = setptr->convtol;
	Rmaj = hptr->R[setptr->plane0*hptr->vpp];
	Zmid = hptr->z[setptr->plane0*hptr->vpp];
	for (dot=1, Redg=0.0, Rinn=2.0*Rmaj; dot<hptr->vpp; dot++) 
	{
		if(hptr->R[dot + setptr->plane0*hptr->vpp] > Redg)
			Redg = hptr->R[dot + setptr->plane0*hptr->vpp];
		if(hptr->R[dot + setptr->plane0*hptr->vpp] < Rinn)
			Rinn = hptr->R[dot + setptr->plane0*hptr->vpp];
	}
	fprintf(stderr, "R range is %f to %f.\n", Rinn, Redg);
	Rmaj = find_axis(hptr, setptr);
	setptr->theta_sym=0;
	Zmid += (Rmaj - hptr->R[setptr->plane0*hptr->vpp])*tan(setptr->theta_sym);
	fprintf(stderr, "Major radius = %f\n", Rmaj);
	fprintf(stderr, "Minor radius = %f\n", Redg - Rmaj);
	fprintf(stderr, "z_ax = %f\n", Zmid);
	h = stepfrac * (Redg - Rmaj);   /* Step distance for line trace */
	fprintf(stderr, "h_search = %le\n", h);

	/* Find Redg in theta direction */
	for(x0[0]=Rmaj;; x0[0]+=h) 
	{
		x0[2] = Zmid + (x0[0] - Rmaj)*tan(setptr->theta_sym);
		if ((dot = get_tri_weights(x0, 0, vv, dx, hptr)) < 0.0) 
			break;
		oldplane = dot;
	}
	fprintf(stderr, "Exit at element %d\n", oldplane);
	fprintf(stderr, "v = %d, %d, %d\n", hptr->conn_list[3*oldplane],hptr->conn_list[3*oldplane+1], hptr->conn_list[3*oldplane+2]);
	fprintf(stderr, "wts = %le, %le, %le\n", dx[0], dx[1], dx[2]);
	hmax = dx[0];
	vv[0]=0;
	if(dx[1] < hmax) 
	{
		hmax=dx[1]; 
		vv[0]=1;
	}
	if(dx[2] < hmax)
		vv[0]=2;
	fprintf(stderr, "min weight at %d\n", vv[0]);
	vv[1] = (vv[0] + 1) % 3;  vv[2] = (vv[0] + 2) % 3;
	hdid  = (hptr->z[hptr->conn_list[3*oldplane + vv[2]]]
			- hptr->z[hptr->conn_list[3*oldplane + vv[1]]])/
			(hptr->R[hptr->conn_list[3*oldplane + vv[2]]]
			- hptr->R[hptr->conn_list[3*oldplane + vv[1]]]);
	fprintf(stderr, "hdid = %le\n", hdid);
	if(hdid == tan(setptr->theta_sym)) 
	{
		if(dx[vv[1]] > dx[vv[2]])
			Redg = hptr->R[hptr->conn_list[3*oldplane + vv[1]]];
		else
			Redg = hptr->R[hptr->conn_list[3*oldplane + vv[2]]];
	} 
	else 
	{
		Redg = (hptr->z[hptr->conn_list[3*oldplane + vv[1]]] - Zmid + Rmaj*tan(setptr->theta_sym)
			 - hptr->R[hptr->conn_list[3*oldplane + vv[1]]]*hdid)/(tan(setptr->theta_sym) - hdid);
	}
	fprintf(stderr, "Redg ?= %lf\n", Redg);

	/* Find Rinn in theta direction */
	for(x0[0]=Rmaj;; x0[0]-=h) 
	{
		x0[2] = Zmid + (x0[0] - Rmaj)*tan(setptr->theta_sym);
		if((dot = get_tri_weights(x0, 0, vv, dx, hptr)) < 0.0)
			break;
		oldplane = dot;
	}
	fprintf(stderr, "Exit at element %d\n", oldplane);
	fprintf(stderr, "v = %d, %d, %d\n", hptr->conn_list[3*oldplane],hptr->conn_list[3*oldplane+1], hptr->conn_list[3*oldplane+2]);
	fprintf(stderr, "wts = %le, %le, %le\n", dx[0], dx[1], dx[2]);
	hmax = dx[0];
	vv[0]=0;
	if(dx[1] < hmax)
	{
		hmax=dx[1];
		vv[0]=1;
	}
	if(dx[2] < hmax)
		vv[0]=2;
	fprintf(stderr, "min weight at %d\n", vv[0]);
	vv[1]	= (vv[0] + 1) % 3;  vv[2] = (vv[0] + 2) % 3;
	hdid	= (hptr->z[hptr->conn_list[3*oldplane + vv[2]]]
				- hptr->z[hptr->conn_list[3*oldplane + vv[1]]])/
				(hptr->R[hptr->conn_list[3*oldplane + vv[2]]]
				- hptr->R[hptr->conn_list[3*oldplane + vv[1]]]);
	Rinn	= (hptr->z[hptr->conn_list[3*oldplane + vv[1]]] - Zmid + Rmaj*tan(setptr->theta_sym)
				- hptr->R[hptr->conn_list[3*oldplane + vv[1]]]*hdid)/(tan(setptr->theta_sym) - hdid);
	fprintf(stderr, "Rinn ?= %lf\n", Rinn);

	rmstep	= (Redg - Rmaj)/(float)setptr->numsurf;  /* Dist. btwn surfs */
	rmstep2 = (Rmaj - Rinn)/(float)setptr->numsurf;
	hmax = 0.25*Rinn*hptr->dphi;
	fprintf(stderr, "Increments: fw=%le, rv=%le\n", rmstep, rmstep2);
	fprintf(stderr, "Initial step size = %le; max = %le\n", h, hmax);

	/* Set scales for error tolerance checking */
	txscal[1]=xscal[1] = 6.283185307;
	txscal[2]=xscal[2] = fabs((double)(Redg - Rinn));

	bufindex=malloc(sizeof(int)*(setptr->numsurf+1));
	bufindex[0]=0;
	for(rad=1;rad<setptr->numsurf+1;rad++)
	{
		bufindex[rad]=bufindex[rad-1]+	(int)(setptr->mindots +((rad-1)/(setptr->numsurf-1.0))*(setptr->outerdots - setptr->mindots));
	}
	
	buffer=malloc(sizeof(float)*bufindex[setptr->numsurf]*3);
	buffersum=malloc(sizeof(float)*bufindex[setptr->numsurf]*3);
	memset(buffer,0,sizeof(float)*bufindex[setptr->numsurf]*3);
	memset(buffersum,0,sizeof(float)*bufindex[setptr->numsurf]*3);

	start[0]=0;
	h =th= stepfrac * (Redg - Rmaj); 
	hdid=thdid=h;
	hnext=thnext=h;


	/* Outer loop (over flux surfaces) */
	//for (rad=0, *start=0; rad<setptr->numsurf; rad++) 

	for(rad=myid; rad<setptr->numsurf; rad=rad+nump) 
	{
		xscal[1] = txscal[1];
		xscal[2] = txscal[1];
		h = th; 
		hdid=thdid;
		hnext=thnext;

		/* Initialize 1st pt. on this flux surface */
		switch (setptr->poindir) 
		{
			case 1: /* Backward */
				x0[0] = Rmaj - rad*rmstep2;
				break;
			case 2: /* Both */
				x0[0] = (rad % 2) ? Rmaj + rad*rmstep : Rmaj - rad*rmstep2;
				break;
			default: /* Forward */
			x0[0] = Rmaj + rad*rmstep;
		}
		xold[0] = x0[0];
		xold[1] = x0[1] = hptr->dphi*(double)setptr->plane0;
		xold[2] = x0[2] = Zmid + (x0[0] - Rmaj)*tan(setptr->theta_sym);
		oldplane = setptr->plane0;

		buffer[bufindex[rad]*3]=(float)rad;
		buffer[bufindex[rad]*3+1]=(float)x0[0];
		buffer[bufindex[rad]*3+2]=(float)x0[2];

		dotspersurf = (int)(setptr->mindots +(rad/(setptr->numsurf-1.0))*(setptr->outerdots - setptr->mindots));

		/* Inner loop (trace the current flux surface) */
		for (dot=1,i=1; dot<dotspersurf;) 
		{
			i++;
		
			if(getBfield(x0, dx, hptr))
				break;      /* Evaluate derivatives */
			xscal[0] = fabs(x0[0]) + 1.0e-30;             /* Set tolerances */
			if(nrrkqs(x0, dx, h, eps, xscal, &hdid, &hnext, hptr)) 
				break; /* stepper */
			if(x0[1] > 2.0*M_PI)
				x0[1] -= 2.0*M_PI;
			else if (x0[1] < 0.0)
				x0[1] += 2.0*M_PI;

			newplane = (int)(x0[1]/hptr->dphi);

			/* Intersection test */
			if((newplane == setptr->plane0) && (newplane != oldplane)) 
			{
				if(setptr->plane0 == 0) 
					xold[1] -= 2.0*M_PI;
				ofac = (x0[1] - newplane*hptr->dphi)/(x0[1] - xold[1]);
				minibuf[0] = (float)(ofac*xold[0] + (1.0-ofac)*x0[0]);
				minibuf[1] = (float)(ofac*xold[2] + (1.0-ofac)*x0[2]);

				buffer[bufindex[rad]*3+(dot)*3]=(float)rad;
				buffer[bufindex[rad]*3+(dot)*3+1]=minibuf[0];
				buffer[bufindex[rad]*3+(dot)*3+2]=minibuf[1];
				dot++;
			}

			h = (hnext > hmax) ? hmax : hnext;
			oldplane = newplane;
			xold[0] = x0[0]; 
			xold[1] = x0[1]; 
			xold[2] = x0[2];
		}
		txscal[1] = xscal[1];
		txscal[2] = xscal[1];
		th = th; 
		thdid=thdid;
		thnext=thnext;
		printf("Surface %d:\t %d intersections.\n", rad,dot);
	} /* end outer loop (rad) */

#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(buffer, buffersum, bufindex[setptr->numsurf]*3, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	memcpy(buffersum,buffer,sizeof(float)*bufindex[setptr->numsurf]*3);
#endif
	if(myid==0)
	{
		if (setptr->filetype == ASCII) 
		{
			fputs("Creating ascii output file.\n", stderr);
			if(!setptr->oflag) 
			{
				sprintf(setptr->outname,"poin%04.0f_%d.dat",hptr->tim,setptr->plane0);
			}
			if((fp = fopen(setptr->outname, "w")) == NULL) 
			{
				fprintf(stderr, "Unable to open file \"%s\" for writing.\n", setptr->outname);
				return;
			}
			for(i=0;i<bufindex[setptr->numsurf];i++)
			{
				fprintf(fp, "%f\t%e\t%e\n", buffersum[i*3],buffer[i*3+1],buffer[i*3+2] );
			}
			fclose(fp);
		}
		else 
		{
			fputs("Creating HDF5 output file.\n", stderr);
			if(!setptr->oflag) 
				sprintf(setptr->outname,"poin%08.2f_%d.h5",hptr->tim,setptr->plane0);
				h5file = H5Fcreate(setptr->outname, H5F_ACC_TRUNC,H5P_DEFAULT, H5P_DEFAULT);
			if(h5file < 0) 
				fputs("Error creating HDF5 file.\n", stderr);

			cursize[0] = bufindex[setptr->numsurf]; 
			cursize[1] = 3;
			h5fspace = H5Screate_simple(2, cursize, NULL);
			if (h5fspace < 0)
				fputs("Error creating file dataspace.\n", stderr);
			h5dataset = H5Dcreate(h5file, "/dset",H5T_NATIVE_FLOAT,h5fspace,H5P_DEFAULT);
			H5Dwrite(h5dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buffersum);
			H5Dclose(h5dataset);
			H5Sclose(h5fspace);
			H5Fclose(h5file);
		}
	}
	free(bufindex); 
	free(buffer);
	free(buffersum);
}


/*============================================================*/
/* stepper */
/* Adapted from routine rkqs(), Numerical Recipes in C, 2nd ed., p. 719 */
int nrrkqs(double *x, double *dxdt, double htry, double eps, double *xscal,
	    double *hdid, double *hnext, hdatatype *hptr)
{
  double errmax, etmp, h, htmp;
  double xtemp[3], xerr[3];

  for (h=htry;;) {
    /* Take a step */
    if (nrrkck(x, dxdt, h, xtemp, xerr, hptr)) return 1;

    /* Evaluate accuracy */
    errmax = fabs(*xerr/(*xscal));
    if ((etmp = fabs(xerr[1]/xscal[1])) > errmax) errmax = etmp;
    if ((etmp = fabs(xerr[2]/xscal[2])) > errmax) errmax = etmp;
    errmax /= eps;
    if (errmax <= 1.0) break;  /* step succeeded; keep the results */

    /* Error too large. Reduce step size and try again. */
    htmp = SAFETY*h*pow(errmax, PSHRNK);
    if ((etmp = 0.1*h) > htmp) h=etmp; else h=htmp;

    /* Check for underflow condition */
    if ((1.0 + h) == 1.0) {
      fputs("Underflow in RK stepper!\n", stderr);
      return 2;
    }
  }

  /* See if the step size can safely be increased */
  if (errmax > ERRCON) *hnext = SAFETY*h*pow(errmax,PGROW);
  else *hnext = 5.0*h;

  /* Advise driver of step size just used */
  *hdid = h;

  /* Set new values */
  x[0] = xtemp[0];  x[1] = xtemp[1];  x[2] = xtemp[2];

  return 0;
}

/*============================================================*/
/* algorithm */
/* Adapted from routine rkck(), Numerical Recipes in C, 2nd ed., pp. 719-720 */
/* Advance vector x a single 5th-order Runge-Kutta step of size h,       */
/* using a different combination of the same intermediates to provide an */
/* estimate of the highest-order truncation error.                       */
int nrrkck(double *x, double *dxdt, double h, double *xout, double *xerr,
	   hdatatype *hptr)
{
  int ierr;
  double xtemp[3], ak2[3], ak3[3], ak4[3], ak5[3], ak6[3];
  static double b21=0.2;
  static double b31=3.0/40.0, b32=9.0/40.0;
  static double b41=0.3, b42=-0.9, b43=1.2;
  static double b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0;
  static double b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0,
    b64=44275.0/110592.0, b65=253.0/4096.0;
  static double c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0;
  static double dc1=37.0/378.0-2825.0/27648.0, dc3=250.0/621.0-18575.0/48384.0,
    dc4=125.0/594.0-13525.0/55296.0, dc5=-277.0/14336.0, dc6=512.0/1771.0-0.25;

  /* First step */
  xtemp[0] = x[0] + b21*h*dxdt[0];
  xtemp[1] = x[1] + b21*h*dxdt[1];
  xtemp[2] = x[2] + b21*h*dxdt[2];

  /* Second step */
  if (ierr = getBfield(xtemp, ak2, hptr)) return ierr;;
  xtemp[0] = x[0] + h*(b31*dxdt[0] + b32*ak2[0]);
  xtemp[1] = x[1] + h*(b31*dxdt[1] + b32*ak2[1]);
  xtemp[2] = x[2] + h*(b31*dxdt[2] + b32*ak2[2]);

  /* Third step */
  if (ierr = getBfield(xtemp, ak3, hptr)) return ierr;
  xtemp[0] = x[0] + h*(b41*dxdt[0] + b42*ak2[0] + b43*ak3[0]);
  xtemp[1] = x[1] + h*(b41*dxdt[1] + b42*ak2[1] + b43*ak3[1]);
  xtemp[2] = x[2] + h*(b41*dxdt[2] + b42*ak2[2] + b43*ak3[2]);

  /* Fourth step */
  if (ierr = getBfield(xtemp, ak4, hptr)) return ierr;
  xtemp[0] = x[0] + h*(b51*dxdt[0] + b52*ak2[0] + b53*ak3[0] + b54*ak4[0]);
  xtemp[1] = x[1] + h*(b51*dxdt[1] + b52*ak2[1] + b53*ak3[1] + b54*ak4[1]);
  xtemp[2] = x[2] + h*(b51*dxdt[2] + b52*ak2[2] + b53*ak3[2] + b54*ak4[2]);

  /* Fifth step */
  if (ierr = getBfield(xtemp, ak5, hptr)) return ierr;
  xtemp[0] = x[0] + h*(b61*dxdt[0] + b62*ak2[0] + b63*ak3[0] + b64*ak4[0] +
		       b65*ak5[0]);
  xtemp[1] = x[1] + h*(b61*dxdt[1] + b62*ak2[1] + b63*ak3[1] + b64*ak4[1] +
		       b65*ak5[1]);
  xtemp[2] = x[2] + h*(b61*dxdt[2] + b62*ak2[2] + b63*ak3[2] + b64*ak4[2] +
		       b65*ak5[2]);

  /* Sixth step */
  if (ierr = getBfield(xtemp, ak6, hptr)) return ierr;

  /* Best guess of new value */
  xout[0] = x[0] + h*(c1*dxdt[0] + c3*ak3[0] + c4*ak4[0] + c6*ak6[0]);
  xout[1] = x[1] + h*(c1*dxdt[1] + c3*ak3[1] + c4*ak4[1] + c6*ak6[1]);
  xout[2] = x[2] + h*(c1*dxdt[2] + c3*ak3[2] + c4*ak4[2] + c6*ak6[2]);

  /* Best estimate of new highest-order truncation error */
  xerr[0] = h*(dc1*dxdt[0] + dc3*ak3[0] + dc4*ak4[0] + dc5*ak5[0] +
	       dc6*ak6[0]);
  xerr[1] = h*(dc1*dxdt[1] + dc3*ak3[1] + dc4*ak4[1] + dc5*ak5[1] +
	       dc6*ak6[1]);
  xerr[2] = h*(dc1*dxdt[2] + dc3*ak3[2] + dc4*ak4[2] + dc5*ak5[2] +
	       dc6*ak6[2]);

  return 0;
}

/*====================================================================*/
int interp2d(float *var, int plane, double R, double z, hdatatype *hptr,
	     float *val)
{
  double x[3], w[3];
  int    v[3];

  x[0] = R;  x[2] = z;
  if (get_tri_weights(x, plane*hptr->vpp, v, w, hptr) < 0) {
    *val = -1.0;
    return 1;
  }

  *val = w[0]*var[v[0]] + w[1]*var[v[1]] + w[2]*var[v[2]];

  return 0;
}

/*============================================================*/
int getBfield(double *x1, double *dx, hdatatype *hptr)
{
  double Rinv = 1.0/(*x1), relphi = x1[1]/hptr->dphi;
  int    nplanes = hptr->planes[0];
  int    baseplane = (int)(relphi + 0.5);
  int    lastplane = baseplane-1, nextplane = baseplane+1;
  int    v[3];
  double tw1, tw2, tw3, w[3], factor;
  double Bzlast, Bzthis, Bznext;
  double BRlast, BRthis, BRnext, Bphilast, Bphithis, Bphinext;

  /* Compute weights for cross-plane quadratic interpolation */
  tw2 = (relphi - lastplane)*(nextplane - relphi);
  tw3 = 0.5*((relphi - lastplane) - tw2);
  tw1 = 1.0 - (tw2 + tw3);

  /* Compute offsets of other planes' vertices */
  baseplane = baseplane % nplanes;
  nextplane = hptr->vpp*((nextplane % nplanes) - baseplane);
  lastplane = hptr->vpp*(((lastplane + nplanes) % nplanes) - baseplane);

  /* Find the triangle in the base plane containing the given point */
  if (get_tri_weights(x1, baseplane*hptr->vpp, v, w, hptr) < 0) 
	  return 1; 

  BRlast	= w[0]*hptr->BR[v[0]  + lastplane] + w[1]*hptr->BR[v[1] + lastplane] +
			w[2]*hptr->BR[v[2]	  + lastplane];
  Bphilast	= (w[0]*hptr->R[v[0]  + lastplane]*hptr->Bphi[v[0] + lastplane] +
				w[1]*hptr->R[v[1] + lastplane]*hptr->Bphi[v[1] + lastplane] +
				w[2]*hptr->R[v[2] + lastplane]*hptr->Bphi[v[2] + lastplane])*Rinv;
  Bzlast	= w[0]*hptr->Bz[v[0]  + lastplane] + w[1]*hptr->Bz[v[1] + lastplane] +
			w[2]*hptr->Bz[v[2]    + lastplane];

  BRthis	= w[0]*hptr->BR[v[0]] + w[1]*hptr->BR[v[1]] + w[2]*hptr->BR[v[2]];
  Bphithis	= (w[0]*hptr->R[v[0]]*hptr->Bphi[v[0]] +
				w[1]*hptr->R[v[1]]*hptr->Bphi[v[1]] +
				w[2]*hptr->R[v[2]]*hptr->Bphi[v[2]])*Rinv;
  Bzthis	= w[0]*hptr->Bz[v[0]] + w[1]*hptr->Bz[v[1]] + w[2]*hptr->Bz[v[2]];

  BRnext	= w[0]*hptr->BR[v[0] + nextplane] + w[1]*hptr->BR[v[1] + nextplane] +
			w[2]*hptr->BR[v[2] + nextplane];
  Bphinext	= (w[0]*hptr->R[v[0] + nextplane]*hptr->Bphi[v[0] + nextplane] +
			w[1]*hptr->R[v[1] + nextplane]*hptr->Bphi[v[1] + nextplane] +
			w[2]*hptr->R[v[2] + nextplane]*hptr->Bphi[v[2] + nextplane])*Rinv;
  Bznext	= w[0]*hptr->Bz[v[0] + nextplane] + w[1]*hptr->Bz[v[1] + nextplane] +
			w[2]*hptr->Bz[v[2] + nextplane];

  /* Interpolate across planes */
  dx[0] = tw1*BRlast   + tw2*BRthis   + tw3*BRnext;
  dx[1] = tw1*Bphilast + tw2*Bphithis + tw3*Bphinext;
  dx[2] = tw1*Bzlast   + tw2*Bzthis   + tw3*Bznext;

  /* Normalize */
  factor = 1.0/sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
  dx[0]*=factor;  dx[1]*=factor*Rinv;  dx[2]*=factor;

  return 0;
}

/*============================================================*/
int get_tri_weights(double *x1, int node_offset, int *v, double *w,
		    hdatatype *hptr)
{
  int        *tri;
  double     f, dx, dy, dR1, dR2, dz1, dz2;
  float      xA, xB, xC, yA, yB, yC;
  static int el=0;
  int        next, inflag;
  int i =0;
  while (1) 
  {
	i++;
	if(i>10000)
	{
		fprintf(stderr,"getweight err i=:%d\n",i);
		return -1;
	}
    inflag = 1;
    tri = hptr->conn_list + 3*el;
    v[0] = tri[0] + node_offset;
    v[1] = tri[1] + node_offset;
    v[2] = tri[2] + node_offset;

    xA = hptr->R[v[0]]; 
	yA = hptr->z[v[0]];
    xB = hptr->R[v[1]]; 
	yB = hptr->z[v[1]];
    if ((xB - xA)*(x1[2] - yA) < (yB - yA)*(*x1 - xA)) 
		goto lout0;

  lin1:
    xC = hptr->R[v[2]];  yC = hptr->z[v[2]];
    dR1 = xC - xB;  dz1 = yC - yB;
    if (dR1*(x1[2] - yB) < dz1*(*x1 - xB)) 
		goto lout1;

  lin2:
    dx = *x1 - xC;  dy = x1[2] - yC;
    dR2 = xA - xC;  dz2 = yA - yC;
    if (dR2*dy < dz2*dx) 
		goto lout2;

    /* side == 3 */
    if (inflag) 
	{
      /* Compute weights */
      f = 1.0/(dR1*dz2 - dR2*dz1);
      w[0] = f*(dy*dR1 - dx*dz1);
      w[1] = f*(dy*dR2 - dx*dz2);
      w[2] = 1.0 - (w[0] + w[1]);

      return el;
    } 
	return -1;

    /* side < 3 */
  lout0:
    if ((next = hptr->neighbor_list[3*el]) >= 0) 
	{
      el = next;        /* Try the nearest neighbor */
      continue;
    } inflag=0; 
	goto lin1;

  lout1:
    if ((next = hptr->neighbor_list[3*el + 1]) >= 0) 
	{
      el = next;        /* Try the nearest neighbor */
      continue;
    } inflag=0; 
	goto lin2;

  lout2:
    if ((next = hptr->neighbor_list[3*el + 2]) >= 0) 
	{
      el = next;        /* Try the nearest neighbor */
      continue;
    } 
	return -1;
  } /* end while */
}

/*============================================================*/
void findElementNeighbors(hdatatype *hptr)
{
  int  v, el, nel=hptr->ncells, vpp=hptr->vpp, *nlist;
  int  *tri, *conn_list=hptr->conn_list;
  edge *edge_list;

  /* Allocate space for neighbor list, initialize */
  if ((nlist = (int *) malloc(3*nel*sizeof(int))) == NULL) {
    fputs("Not enough memory in findElementNeighbors.\n", stderr);
    return;
  }
  for (el=0; el<3*nel; el++) nlist[el] = -1;

  /* Allocate and initialize edge list hash table */
  edge_list = (edge *) malloc((vpp-1) * sizeof(edge));
  for (v=0; v<vpp-1; v++) edge_list[v].n = 0;

  /* Build tables */
  for (el=0, tri=conn_list; el<nel; el++, tri+=3) {
    add_edge(edge_list, tri, 0, el, nlist);
    add_edge(edge_list, tri, 1, el, nlist);
    add_edge(edge_list, tri, 2, el, nlist);
  } /* end loop el */
  free(edge_list);

  hptr->neighbor_list = nlist;
}

/*============================================================*/
void add_edge(edge *list, int *tri, char side, int el,
	      int *nlist)
{
  int  i, v1, v2, vo;
  edge *ed;

  /* Sort the vertices */
  v1 = tri[side];  v2 = tri[(side+1)%3];
  if (v1 < v2) { ed = list+v1;  vo = v2; }
  else         { ed = list+v2;  vo = v1; }

  /* See if this edge is already present */
  for (i=0; i<ed->n; i++)
    if (ed->o[i].v == vo) {           /* It is! Update the neighbor table. */
      nlist[3*el + side] = ed->o[i].el0;
      nlist[3*ed->o[i].el0 + ed->o[i].side] = el;
      return;
    }

  /* The edge was not present; add it. */
  ed->o[ed->n].v = vo;
  ed->o[ed->n].el0 = el;
  ed->o[ed->n].side = side;
  ed->n++;
}

/*============================================================*/
/* Find the real magnetic axis by minimizing excursions       */
double find_axis(hdatatype *hptr, settype *setptr)
{
  int    vert, pass;
  double Rmaj, Redg, rmin, p0, p00, pm, pp, pn, f0, fn, newtheta, newR;
  double globmin=1.0e+7, Rglob, thetaglob;

  /* Find R at edge on midplane in target poloidal plane */
  Rmaj = hptr->R[setptr->plane0*hptr->vpp];
  for (vert=1, Redg=0.0; vert<hptr->vpp; vert++)
    if (hptr->R[vert + setptr->plane0*hptr->vpp] > Redg)
      Redg = hptr->R[vert + setptr->plane0*hptr->vpp];
  rmin = Redg - Rmaj;

  if (setptr->bguess != -1.0) 
  {
    setptr->theta_sym = atan2(setptr->bguess - hptr->z[setptr->plane0 * hptr->vpp],
			      setptr->aguess - hptr->R[setptr->plane0 * hptr->vpp]);
    fprintf(stderr, "Using theta0 = %lf\n", setptr->theta_sym);
  }
  thetaglob = setptr->theta_sym;

  if (setptr->aguess < 0.0)
    Rglob = p0 = Rmaj;
  else
    Rglob = p0 = setptr->aguess;

  for (pass=0; pass<setptr->maxpass; pass++) 
  { 
	/* Find the right theta_sym */
    fprintf(stderr, "Pass %d:\n", pass+1);
    p00 = p0;
    if (setptr->aguess < 0.0) 
	{
      pm = p0 - 0.55*rmin;
      pp = p0 + 0.55*rmin;
    } else 
	{
      pm = p0 - setptr->wingspan*rmin;
      pp = p0 + setptr->wingspan*rmin;
    }

    f0 = get_excursion(p0, hptr, setptr, &newtheta, &newR);
    printf("Initial bracket: %lf, %lf, %lf (exc0=%le)\n", pm,p0,pp,f0);
    while (pp - pm > 1.0e-6*rmin) 
	{
      if ((pp - p0) > (p0 - pm))
	  {
		pn = p0 + 0.381966*(pp - p0);
		fn = get_excursion(pn, hptr, setptr, &newtheta, &newR);
      } 
	  else 
	  {
		pn = p0; fn = f0;
		p0 -= 0.381966*(p0 - pm);
		f0 = get_excursion(p0, hptr, setptr, &newtheta, &newR);
      }

      if (fn < f0)
	  {
		pm = p0;
		p0 = pn;
		f0 = fn;
      } 
	  else
	  {
		pp = pn;
      }
      //fprintf(stderr, "%lf, %lf, %lf (exc0=%le)\n", pm,p0,pp,f0);
    }
    printf("Best guess: %lf vs %lf (exc=%le)\n", p0, p00, f0);
		
    if (f0 < globmin)
	{
      globmin = f0;
      Rglob = p0;
      thetaglob = setptr->theta_sym;
    }
	
    //printf("New theta = %lf\n", newtheta);
    //printf("New R = %lf\n", newR);
    setptr->theta_sym = newtheta;
    p0 = newR;
	if(globmin<1e-4)
		break;
  }

  printf("Global best: R=%lf, theta=%lf: exc=%le\n", Rglob, thetaglob, globmin);
  setptr->theta_sym = thetaglob;
  return Rglob;
}

/*============================================================*/
double get_excursion(double R0, hdatatype *hptr, settype *setptr, double *newangle,
		     double *newR)
{
  long int  N=1L;
  double    h, hdid, hnext, hmax;
  double    xsum, xssum, ysum, yssum;
  double    x0[3], dx[3], xscal[3], phi0;
  int       transit, last=0;

  x0[0] = xsum = R0;
  xssum = R0*R0;
  x0[1] = phi0 = setptr->plane0*hptr->dphi;
  x0[2] = ysum = hptr->z[setptr->plane0*hptr->vpp] +
    (R0 - hptr->R[setptr->plane0*hptr->vpp])*tan(setptr->theta_sym);
  yssum = ysum*ysum;

  xscal[0] = xscal[2] = R0 + 1.0e-30;
  xscal[1] = 6.283;
  h = 0.001*xscal[0];
  hmax = 0.01*x0[0]*hptr->dphi;

  do {
    if (getBfield(x0, dx, hptr)) break;           /* Evaluate derivatives */
    if (nrrkqs(x0, dx, h, 1.0e-9, xscal, &hdid, &hnext, hptr)) break; /* stepper */
    h = (hnext > hmax) ? hmax : hnext;

    transit = (int)((x0[1] - phi0)/(2.0*M_PI));
    if (transit != last) {
      xsum += x0[0];
      xssum += x0[0]*x0[0];
      ysum += x0[2];
      yssum += x0[2]*x0[2];
      N++;
    }
    last = transit;

  } while (N < setptr->maxtrans);

  if (N < setptr->maxtrans) return 1.0e+6;  /* R0 is out of range. */

  xsum /= N;  xssum /= N;
  ysum /= N;  yssum /= N;
  *newangle = atan2(ysum - hptr->z[0], xsum - hptr->R[0]);
  *newR = xsum;
  /* *newangle = setptr->theta_sym; */
  return sqrt(xssum - xsum*xsum + yssum - ysum*ysum);
}



/*============================================================*/
void filter(float *buf, int mode, hdatatype *hptr)
{
  int    plane, i;
  float  *ptr;
  double factor, factor0, fcmphi, fsmphi;
  double *ccoef, *scoef;

  ccoef = (double *)malloc(hptr->vpp*sizeof(double));
  scoef = (double *)malloc(hptr->vpp*sizeof(double));

  /* Error check */
  if (mode > hptr->planes[0]/2)
    fputs("Warning: mode requested too high; aliasing may occur.\n", stderr);

  /* Initialize factor */
  factor0 = 1.0/hptr->planes[0];
  factor = (mode) ? 2.0*factor0 : 0.0;

  /* Zeroth plane */
  for (i=0; i<hptr->vpp; i++) {
    ccoef[i] = buf[i] * factor;
    scoef[i] = 0.0;
  }

  /* Loop over other planes to complete the DFT */
  for (plane=1, ptr=buf+hptr->vpp; plane<hptr->planes[0];
       plane++, ptr+=hptr->vpp) {

    /* Update the coefficients */
    fcmphi = factor*hptr->ctab[(mode*plane) % hptr->planes[0]];
    fsmphi = factor*hptr->stab[(mode*plane) % hptr->planes[0]];
    for (i=0; i<hptr->vpp; i++) {
      ccoef[i] += ptr[i] * fcmphi;
      scoef[i] += ptr[i] * fsmphi;
      buf[i] += ptr[i];
    } /* end loop i */
  } /* end loop plane */

  /* Loop over planes again to reconstruct requested mode */
  for (plane=1, ptr=buf+hptr->vpp; plane<hptr->planes[0];
       plane++, ptr+=hptr->vpp)
    for (i=0; i<hptr->vpp; i++)
      ptr[i] = factor0*buf[i] +
	ccoef[i]*hptr->ctab[(mode*plane) % hptr->planes[0]] +
	scoef[i]*hptr->stab[(mode*plane) % hptr->planes[0]];

  for (i=0; i<hptr->vpp; i++)
    buf[i] = factor0*buf[i] + ccoef[i];

  free(scoef); free(ccoef);
}

