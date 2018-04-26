/* module load hdf5 */
/* cc  -DH5_USE_16_API -DUSE_MPI getfw_f.c -o getfw ${HDF5_INCLUDE_OPTS} ${HDF5_POST_LINK_OPTS} */

#define H5_USE_16_API
#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <hdf5.h>
#include "readhdf5data_new.h"

#ifdef USE_MPI
#include<mpi.h>
#endif

#define SAFETY 0.9
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define PGROW -0.2

/* Type definitions */

typedef struct {
	int  el0, v;
	char side;
} d_edge;

typedef struct{
	d_edge o[6];
	char   n;
} edge;

typedef struct {
	int  plane0;  /* Index of plane for Poincare plotting */
	int  poindir; /* Direction for Poincare plotting (default=forward) */
	double convtol, aguess, qseek, epsilon;
	int  qsurf, qtrans, frame;
	
	int  oflag;       /* Output filename specified? */
	char outname[128];  /* Output filename, if specified. */
	int  filtermode; /* Toroidal mode number for filtering */
	int  jump; /* First surface to trace */
	double bguess, theta_sym, wingspan;
	int maxtrans;
} settype;

typedef struct
{
	float planeno;
	float time;
	float islandwidth;
	float Rw;
}islanddata;

/* Prototypes */
int    nrrkqs(double *, double *, double, double, double *, double *, double *,
			  hdatatype *);
int    nrrkck(double *, double *, double, double *, double *, hdatatype *);
int    getBfield(double *, double *, hdatatype *);
int    get_tri_weights(double *, int, int *, double *, hdatatype *);
void   findElementNeighbors(hdatatype *);
void   add_edge(edge *, int *, char, int, int *);
double find_axis(hdatatype *, settype *);
double get_excursion(double, hdatatype *, settype *);

void   filter(float *, int, hdatatype *);
int    myround(double);
double getislandwidth(hdatatype *hptr, settype *setptr, double* Rw);

double find_O_point(hdatatype *hptr, settype *setptr);
void setdefault(settype * settings);
void printhelp(settype  settings);
void getoption(int argc,char *argv[],settype * settings,int *eqflag);
void alloc_pointmem(hid_t h5fid,hdatatype* hdata);
void free_pointmem(hdatatype* hdata);

/*============================================================*/
int main(int argc, char *argv[])
{
	hid_t h5fid, groupid,   dataid;
            hsize_t       npoints;
	double        parity = 1.0;
	int            frame=0, eqflag=0, first=1;
	char          fname[64],fout[100];
	FILE		* fhand;
	hdatatype     hdata;
	settype       settings;
	
	/*total island width need to get*/
	int index;
	int totalpt,nopt,remainpt,idpt;
	int frame0,frame1;
	int plane,plane0,plane1;
	int myid,numprocs;
	islanddata* islandwidthtotal;
	islanddata* islandwidth;
	islanddata* maxislandwidth;
	int * noptarray, * disp;
	double Rw;
	
#ifdef USE_MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
	numprocs=1;
	myid=0;
#endif
	setdefault(&settings);
	
	if (argc < 2)
	{ /* Print usage message if no arguments are given */
		if(0==myid)
			printhelp(settings);
		exit(0);
	}
	/* Parse command-line options */
	strcpy(fname, argv[1]);
	getoption(argc,argv,&settings,&eqflag);
	
	/* Open the HDF5 file, read-only */
	if ((h5fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0)
	{
		fprintf(stderr, "myid =%d,Could not access HDF5 file %s.\n",myid, fname);
		return 1;
	}
	
	/* Read number of time slices in file */
	groupid = H5Gopen(h5fid, "/");
	dataid      = H5Aopen_name(groupid, "nsteps");
	H5Aread(dataid, H5T_NATIVE_INT_g, &hdata.nframes);
	H5Aclose(dataid);  H5Gclose(groupid);
	if(0==myid)
		fprintf(stderr, "%d time slices in file.\n", hdata.nframes);
	
	if (settings.frame >= hdata.nframes)
	{
		fputs("Frame requested is out of bounds.\n", stderr);
		H5Fclose(h5fid);
		return 2;
	}
        read_plane_frame(h5fid,0,&hdata);
	read_connection(h5fid,&hdata);
	
	totalpt     =hdata.nframes*hdata.planes[0];
	nopt  =(int)floor((float)totalpt/(float)numprocs);
	remainpt=totalpt-nopt*numprocs;
	
	islandwidthtotal  =     malloc(totalpt*sizeof(islanddata));
	maxislandwidth	=	malloc(hdata.nframes*sizeof(islanddata));
	
	noptarray		=	malloc(sizeof(int)*numprocs);
	disp			=	malloc(sizeof(int)*numprocs);
	
	for(index=0;index<numprocs;index++)
	{
		noptarray[index]=nopt;
	}
	for(index=0;index<remainpt;index++)
	{
		noptarray[index]++;
	}
	for(index=0;index<numprocs;index++)
	{
		noptarray[index]=noptarray[index]*sizeof(islanddata)/sizeof(float);
	}
	disp[0]=0;
	for (index=1; index<numprocs; index++) 
	{
		disp[index]=disp[index-1]+noptarray[index-1];
	}
	
	if(myid<remainpt)
	{
		nopt++;     
		frame0      =nopt*myid/hdata.planes[0];
		frame1      =nopt*(myid+1)/hdata.planes[0];
		plane0	    =(nopt*(myid))%hdata.planes[0];
		plane1      =(nopt*(myid+1))%hdata.planes[0];
	}
	else
	{
		frame0      =((nopt+1)*remainpt+(myid-remainpt)*nopt)/hdata.planes[0];
		frame1      =((nopt+1)*remainpt+(myid+1-remainpt)*nopt)/hdata.planes[0];
		plane0      =((nopt+1)*remainpt+(myid-remainpt)*nopt)%hdata.planes[0];
		plane1      =((nopt+1)*remainpt+(myid+1-remainpt)*nopt)%hdata.planes[0];
	}
	islandwidth=malloc(nopt*sizeof(islanddata));
	
	
	
	fprintf(stderr,"myid=%d,frame0=%d, plane0=%d, frame1=%d,plane1=%d \n",myid,frame0,plane0,frame1,plane1);
	if(frame0>(hdata.nframes-1)||frame1>(hdata.nframes))
	{
		fprintf(stderr,"parallel alloc error\n");
		exit(3);
	} 
	alloc_pointmem(h5fid,&hdata);

        read_xyz(h5fid,&hdata,&parity ,&npoints);
	for(idpt=0;idpt<nopt;idpt++)
	{
		
		frame=frame0+(idpt+plane0)/hdata.planes[0];
		plane=(idpt+plane0)%hdata.planes[0];

		settings.plane0=plane;
		settings.frame =frame;
		islandwidth[idpt].planeno=(float)plane;
		printf("begin myid=%d,frame=%d,plane=%d\n",myid,frame,plane);
		hdata.tim=islandwidth[idpt].time=readHDF5time(h5fid,frame, &hdata);
	        read_B(h5fid,frame,npoints,parity,eqflag,&hdata);	
		islandwidth[idpt].islandwidth=(float)getislandwidth(&hdata, &settings,&Rw);
		islandwidth[idpt].Rw=Rw;	
		printf("finish myid=%d,frame=%d,plane=%d\n",myid,frame,plane);
	}
	
	//printf("myid=%d finish \n",myid);
#ifdef USE_MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(islandwidth, noptarray[myid],MPI_FLOAT,islandwidthtotal, noptarray, disp, MPI_FLOAT,0, MPI_COMM_WORLD);
#else
//	memcpy(islandwidth,islandwidthtotal,sizeof(islanddata)*nopt);
        for(idpt=0;idpt<nopt;idpt++)
	{
		islandwidthtotal[idpt].time=islandwidth[idpt].time;
		islandwidthtotal[idpt].islandwidth=islandwidth[idpt].islandwidth;
	}
#endif

	
	if(0==myid)
	{
		sprintf(fout,"%s.island",fname);
		fhand=fopen(fout,"w");
		for (index=0;index<totalpt; index++)
		{
			fprintf(fhand, "%f\t%f\t%f %f \n",islandwidthtotal[index].time,islandwidthtotal[index].planeno,islandwidthtotal[index].islandwidth,islandwidthtotal[index].Rw);
		}
		fclose(fhand);				
	}
	H5Fclose(h5fid);
	
	/* Deallocate the dynamic arrays */
	free(islandwidthtotal);
	free(islandwidth);
	free(maxislandwidth);
	free(noptarray);
	free(disp);
												 
	free_pointmem(&hdata);
	free(hdata.conn_list);
	
	if(0==myid)
	{

		printf("finish work \n");
	}
#ifdef USE_MPI
	MPI_Finalize();
#endif
	return 0;
}

void printhelp(settype settings)
{
	fprintf(stderr, "Usage: <HDF5 filename> [<options>]\n");
	fputs("\nOptions:\n", stderr);
	fputs(" -tol[erance] <tolerance>  -  convergence tol. for integration (default 1e-9).\n", stderr);
	fputs(" -eq  -  assume file is equilibrium (mag. axis at vert. 0).\n",stderr);
	fputs(" -g[uess] <radius>  -  initial guess for axis position.\n", stderr);
	fputs(" -m[ode] <n> - restrict B pert to toroidal mode # n.\n", stderr);
	fprintf(stderr," -j[ump] <index of innermost surface> (default %d).\n",settings.jump);
	fputs(" -ba[ck]  - starting pts of surfaces go inboard from axis.\n",stderr);
	fprintf(stderr," -qs[urfaces] <# of pts in q profile> (default %d).\n",settings.qsurf);
	fprintf(stderr," -qt[ransits] <# of toroidal transits for q> (default %d).\n",settings.qtrans);
	fprintf(stderr," -qk <target q val> (default %le).\n", settings.qseek);
	fprintf(stderr," -e[psilon] <q target tolerance> (default %le).\n", settings.epsilon);
	fputs(" -o[utput] <output filename> (default \"q.dat\").\n",stderr);
}
void setdefault(settype * settings)
{
	/* Default settings */
	settings->plane0 = 0;         /* Index of plane for Poincare plotting */
	settings->poindir=0;                /* Direction for Poincare plotting (default=forward) */
	settings->convtol=1.0e-9;
	settings->aguess=-1.0;
	settings->qsurf=1128;
	settings->qtrans=1200;
	settings->oflag=0;                  /* Output filename specified? */
	settings->filtermode=-1;            /* Toroidal mode number for filtering */
	settings->jump		= 640;
	settings->qseek		= 2.0;
	settings->frame		= -1;
	settings->epsilon	= 2.0e-3;
	settings->aguess	=-1.0;
    settings->bguess	=-1.0;
    settings->theta_sym = 0.0;
	settings->wingspan=0.05;
	settings->maxtrans = 10;                        /* Dots to determine axis position */
}

void getoption(int argc,char *argv[],settype * settings, int *eqflag)
{
	int arg;
	for (arg=2; arg<argc; arg++)
	{
		if (!strncmp(argv[arg], "-tol", 4))
		{
			settings->convtol = atof(argv[++arg]);
			fprintf(stderr, "Resetting tolerance to %le.\n", settings->convtol);
			continue;
		}
		if (!strncmp(argv[arg], "-eq", 3))
		{
			*eqflag = 1;
			continue;
		}
		if (!strncmp(argv[arg], "-o", 2))
		{
			settings->oflag = 1;
			strcpy(settings->outname, argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-g", 2))
		{
			settings->aguess = atof(argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-m", 2))
		{
			settings->filtermode = atoi(argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-ba", 3))
		{
			settings->poindir = 1;
			fputs("Stepping inward.\n", stderr);
			continue;
		}
		if (!strncmp(argv[arg], "-qs", 3))
		{
			settings->qsurf = atoi(argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-qt", 3))
		{
			settings->qtrans = atoi(argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-qk", 3))
		{
			settings->qseek = atof(argv[++arg]);
			continue;
		}
		if (!strncmp(argv[arg], "-e", 2))
		{
			settings->epsilon = atof(argv[++arg]);
			continue;
		}
		fprintf(stderr, "Unrecognized option %s.\n", argv[arg]);
	} /* end loop arg */
}
void alloc_pointmem(hid_t h5fid,hdatatype* hdata)
{
	hid_t  groupid,  dataid, dataspace;
	hssize_t      npoints;
	char          node_name[64];
	
	/* Determine group name for coords. in this time slice, open the group */
	sprintf(node_name, "/time_coordinates[%d]/coordinates", 0);
	if ((groupid = H5Gopen(h5fid, node_name)) < 0)
		exit(2);
	/* Open the dataset for reading */
	dataid            = H5Dopen(groupid, "values");
	/* Determine the size of the data */
	dataspace   = H5Dget_space(dataid);
	npoints     = H5Sget_simple_extent_npoints(dataspace);
	H5Gclose(groupid);
	
	
	/* Allocate memory for data */
	npoints /= 3;
	hdata->npoints=(int)npoints;
	
	
	hdata->x = (float *)malloc((int)npoints * sizeof(float));
	hdata->y = (float *)malloc((int)npoints * sizeof(float));
	hdata->z = (float *)malloc((int)npoints * sizeof(float));
	hdata->R = (double*)malloc((int)npoints * sizeof(double));
	
	
	hdata->psi = (float *)malloc((int)npoints * sizeof(float));
	
	
	hdata->Bx = (float *)malloc((int)npoints * sizeof(float));
	hdata->By = (float *)malloc((int)npoints * sizeof(float));
	hdata->Bz = (float *)malloc((int)npoints * sizeof(float));
	
	
	hdata->ctab = (double *)malloc((hdata->planes[0]+2) * sizeof(double))+1;
	hdata->stab = (double *)malloc((hdata->planes[0]+2) * sizeof(double))+1;
	
	
	hdata->BR   = (float *)malloc((int)npoints * sizeof(float));
	hdata->Bphi = (float *)malloc((int)npoints * sizeof(float));
}
void free_pointmem(hdatatype* hdata)
{
	free(hdata->x);
	free(hdata->y);   
	free(hdata->z);
	free(hdata->R);
	free(hdata->psi);
	
	free(hdata->Bz);
	free(hdata->BR);
	free(hdata->Bphi);
	
	free(hdata->stab-1);
	free(hdata->ctab-1);
}

double getq(int plane0,hdatatype *hptr,double *zax,double *xax,int qtrans,double h,double hmax,double *xscal,double *x0, double *dx)
{
	int		ttransit, ptransit,dum1;
	int		iplane;
	double	dphi=0;
	double	dtheta=0;
	double	oldplane=0, newplane,theta;
	const int naxplane = 360;
	double  phi_old = x0[1];
	double	thold = 0.0;
	double	eps=1.0e-9;
	double	hdid,hnext;
	int		index=0;
	/*Inner loop*/
	for (ttransit=ptransit=0;(ttransit < qtrans) && (ptransit < qtrans);)
	{
		if (getBfield(x0, dx, hptr))
		{
			fprintf(stderr,"error: Get B field out of region\n" );
			break;            
		}
		/* Evaluate derivatives */
		xscal[0] = fabs(x0[0]) + 1.0e-30;                                               /* Set tolerances */
		if (nrrkqs(x0, dx, h, eps, xscal, &hdid, &hnext, hptr))
		{
			fprintf(stderr,"error rk fail \n" );
			break;																			/* stepper */
		}
		dphi += x0[1] - phi_old;
		if (x0[1] < 0.0)
			x0[1] += xscal[1];
		else if (x0[1] > xscal[1])
		{
			x0[1] -= xscal[1];
		}
		phi_old = x0[1];
		
		iplane	= (int)(x0[1]/(xscal[1]/naxplane));
		dum1	= (int)(x0[1]/(xscal[1]/naxplane) - iplane);
		theta	= atan2(x0[2]-((1.0-dum1)*zax[iplane] +dum1*zax[(iplane+1)%naxplane]), x0[0]-((1.0-dum1)*xax[iplane] +dum1*xax[(iplane+1)%naxplane]));
		
		dtheta += theta - thold;
		if ((theta - thold) > 1.0)
		{
			dtheta -= 2.0*M_PI;
			ptransit++;
		}
		else if ((theta - thold) < -1.0)
		{
			dtheta += 2.0*M_PI;
			ptransit++;
		}
		thold = theta;
		newplane = (int)(x0[1]/hptr->dphi);
		if ((newplane == plane0) && (newplane != oldplane)) /* Intersection */
			ttransit++;
		
		h = (hnext > hmax) ? hmax : hnext;
		oldplane = newplane;
		
		index++;
		if(index>1085800)
		{
			fprintf(stderr,"error max index up flow: %d dphi=%f,dtheta=%f,ttransit=%d,ptransit=%d\n",index,dphi,dtheta,ttransit,ptransit);
			return fabs(dphi/dtheta);
		}
		
	} /* end inner loop (ttransit) */
	
	//fprintf(stderr,"index: %d dphi=%f,dtheta=%f \n",index,dphi,dtheta);
	return fabs(dphi/dtheta);
}



/*============================================================*/
double getislandwidth(hdatatype *hptr, settype *setptr, double *Rw)
{
	double stepfrac=0.001, eps;
	double rmin, rmstep, h, hdid, hnext, hmax, dum1;
	double epsilon = setptr->epsilon;
	double  phi_old, psi, R_old, z_old;
	double shere, slast, lts, psi0, qold=setptr->qseek,alpha, beta;
        double Ra,Rb;
	
	double x0[3], dx[3], xscal[3], *xax, *zax;
	float  Rmaj, Zmid, Redg, Ztop, Rinn;
	int    ttransit,  rad,  newplane, tempo[3];
	int    started=0;
	const int naxplane = 360;
	
	/* data for get the jump point*/
	int max_step;
	double q0, q1,qmid;
	double R0,R1;
	double RL,RR,Rmid;
	int myid;
   	*Rw=0;

#ifdef USE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#else
	myid=0;
#endif
	
	psi0 = hptr->psi[setptr->plane0 * hptr->vpp];
	
	xax = (double *)malloc(naxplane * sizeof(double));
	zax = (double *)malloc(naxplane * sizeof(double));
	
	/* Find R at edge on midplane in target poloidal plane */
	eps = setptr->convtol;
	Zmid = hptr->z[0];
	Rinn = (float)hptr->R[0];
	for (ttransit=1, Redg=Ztop=0.0; ttransit<hptr->vpp; ttransit++)
	{
		if (hptr->R[ttransit] > Redg)  Redg = (float)hptr->R[ttransit];
		if (hptr->R[ttransit] < Rinn)  Rinn = (float)hptr->R[ttransit];
		if (hptr->z[ttransit] > Ztop)  Ztop = (float)hptr->z[ttransit];
	}
	
	/* Find magnetic axis by minimizing excursions */
	Rmaj = (float)find_axis(hptr, setptr);
	
	/* Find psi_0 */
	x0[0] = Rmaj; 
	x0[1] = setptr->plane0*hptr->dphi; 
	x0[2] = Zmid;
	get_tri_weights(x0, 0, tempo, dx, hptr);
	psi = dx[0]*hptr->psi[tempo[0]] + dx[1]*hptr->psi[tempo[1]] +dx[2]*hptr->psi[tempo[2]];

	
	
	/* kappa = (Ztop - Zmid)/(Redg - Rmaj);*/  /* Elongation */
	Rmaj += (float)hptr->R[0];
	if (setptr->poindir)
		rmstep = (Rinn - Rmaj)/(float)(setptr->qsurf+1);  /* Dist btwn surfs */
	else
		rmstep = (Redg - Rmaj)/(float)(setptr->qsurf+1);  /* Dist btwn surfs */
	
	h = stepfrac * (Redg - Rmaj);   /* Step distance for line trace */
	hmax = 0.05*Rmaj*hptr->dphi;
	//fprintf(stderr, "step size = %le; max = %le\n", h, hmax);
	
	/* Set scales for error tolerance checking */
	xscal[1] = hptr->planes[0] * hptr->dphi;
	xscal[2] = fabs((double)(Ztop - Zmid));
	
	//if(myid==10)
	//	printf("begin to find 3d axis position\n");
	
	/* Find 3D axis position */
	xax[0] = x0[0] = Rmaj - hptr->R[0];  x0[1] = 0;  zax[0] = x0[2] = Zmid;
	for (newplane=1; newplane<naxplane; newplane++)
	{
		while (x0[1] < newplane*xscal[1]/naxplane)
		{
			getBfield(x0, dx, hptr);
			xscal[0]    = fabs(x0[0]) + 1.0e-30;
			R_old       = x0[0];
			phi_old     = x0[1];
			z_old       = x0[2];
			nrrkqs(x0, dx, h, eps, xscal, &hdid, &hnext, hptr);
			h           = (hnext > hmax) ? hmax : hnext;
		}
		dum1 = (newplane*xscal[1]/naxplane - phi_old)/(x0[1] - phi_old);
		xax[newplane] = dum1*x0[0] + (1.0 - dum1)*R_old;
		zax[newplane] = dum1*x0[2] + (1.0 - dum1)*z_old;
	}
	
	
	/* Get the jump point*/
	q0=q1=0;
	R0=Redg;
	R1=Rmaj;
	RL=0;
	RR=R0-R1;
	max_step=0;
	h = stepfrac*(Redg - Rmaj);

	rmin=RL+0.2;
	x0[0] = R1 + rmin- hptr->R[0];
	x0[1] = setptr->plane0*hptr->dphi;
	x0[2] = Zmid;
	
	q0=getq(setptr->plane0,hptr,zax,xax,setptr->qtrans,h,hmax,xscal,x0,dx);

	rmin=RR-0.01;
	x0[0] = R1 + rmin- hptr->R[0];
	x0[1] = setptr->plane0*hptr->dphi;
	x0[2] = Zmid;

	q1=getq(setptr->plane0,hptr,zax,xax,setptr->qtrans,h,hmax,xscal,x0,dx);

	do
	{
		Rmid	=	(RL+RR)*0.5;
		x0[0] = R1 + Rmid- hptr->R[0];
		x0[1] = setptr->plane0*hptr->dphi;
		x0[2] = Zmid;
		qmid=getq(setptr->plane0,hptr,zax,xax,setptr->qtrans,h,hmax,xscal,x0,dx);
		
		if(qmid>=setptr->qseek-0.05)
		{
			RR=Rmid;
			q1=qmid;
		}
		else
		{
			RL=Rmid;
			q0=qmid;
		}
		max_step++;
		if(max_step>30)
		{
			fprintf(stderr,"error: max_step =%d, find jump point err \n",max_step);
			fprintf(stderr,"frame %d, plane %d \n",setptr->frame,setptr->plane0);
			fprintf(stderr,"q0=%f, q1=%f \n",q0,q1);
			break;
		}
		if(q0>=setptr->qseek-0.05)
		{
			fprintf(stderr,"error: max_step =%d, find jump point err \n",max_step);
			fprintf(stderr,"frame %d, plane %d \n",setptr->frame,setptr->plane0);
			fprintf(stderr,"q0=%f, q1=%f \n",q0,q1);
			break;
		}
		
	}while(setptr->qseek-q0>0.10);


	setptr->jump=(int)((RL)/rmstep);
		//printf("finish jump frame=%d plane=%d, jump to %d at R=%f q=%f \n",setptr->frame,setptr->plane0,setptr->jump,RL+R1,q0);

	/* Outer loop (over flux surfaces) */
	for (rad=setptr->jump, rmin=rmstep*(setptr->jump);rad<setptr->qsurf; rad++, rmin+=rmstep)
	{
		/* Initialize 1st pt. on this flux surface */
		x0[0] = Rmaj + rmin - hptr->R[0];
		x0[1] = setptr->plane0*hptr->dphi;
		x0[2] = Zmid;
		
		/* Find psi here */
		get_tri_weights(x0, 0, tempo, dx, hptr);
		psi = dx[0]*hptr->psi[tempo[0]] + dx[1]*hptr->psi[tempo[1]] +dx[2]*hptr->psi[tempo[2]];
		shere = sqrt(1.0 - psi/psi0);
		Rb =x0[0];
		q1=getq(setptr->plane0,hptr,zax,xax,setptr->qtrans,h,hmax,xscal,x0,dx);

		alpha = setptr->qseek - qold;
		qold = q1;
		beta = qold - setptr->qseek;
		if (alpha*beta > 0.0)
		{
			fprintf(stderr, "q=%.5le at R = %lf +/- %.1le\n", q1,Rmaj + rmin - rmstep + alpha*rmstep/(alpha + beta),0.1*fabs(rmstep));
		}
		//fprintf(stderr, "%d: R=%lf, s=%lf, q=%le\n", rad, Rmaj+rmin, shere,qold);
		if (started)
		{
			if (fabs(qold- setptr->qseek)> epsilon)
			{
				fprintf(stderr, "last s = %le, rad=%d\n", shere, rad);
				printf("%d\t%d\t%.6f\t%.3le\n", setptr->frame,setptr->plane0,hptr->tim, shere - lts);
		                *Rw=Rb-Ra;		
                                printf(" Rw= %f\n ",*Rw);
				return shere-lts;
			}
		} 
		else if (fabs(qold - setptr->qseek) <= epsilon)
		{
			started = 1;
			fprintf(stderr, "1st s = %le, rad=%d\n", shere, rad);
			lts = shere;
                        Ra=Rb;
			if (rad == setptr->jump)
				fprintf(stderr,"error: *** WARNING: jump TOO LARGE,p=%d,time=%f",setptr->plane0,hptr->tim);
		} 
		else if (qold- setptr->qseek > epsilon*50)
		{
			printf("%.6f\t%.3le\n", hptr->tim, 0.00000);
			return 0;
		}
		slast = shere;
	} /* end outer loop (rad) */

	free(xax);
	free(zax);
	printf("%.6f\t%.3le\n", hptr->tim, 0.00000);
	return 0;

}


/*============================================================*/
/* stepper */
/* Adapted from routine rkqs(), Numerical Recipes in C, 2nd ed., p. 719 */
int nrrkqs(double *x, double *dxdt, double htry, double eps, double *xscal,
		   double *hdid, double *hnext, hdatatype *hptr)
{
	double errmax, etmp, h, htmp;
	double xtemp[3], xerr[3];
	int i=0;
	int err=0;
	for (h=htry;;) 
	{		
		/* Take a step */
		if (nrrkck(x, dxdt, h, xtemp, xerr, hptr)) 
			return 1;

		/* Evaluate accuracy */
		errmax = fabs(*xerr/(*xscal));
		if ((etmp = fabs(xerr[1]/xscal[1])) > errmax) 
			errmax = etmp;
		if ((etmp = fabs(xerr[2]/xscal[2])) > errmax) 
			errmax = etmp;
		errmax /= eps;
		if (errmax <= 1.0) 
			break;  /* step succeeded; keep the results */
		
		/* Error too large. Reduce step size and try again. */
		htmp = SAFETY*h*pow(errmax, PSHRNK);
		if ((etmp = 0.1*h) > htmp) 
			h=etmp; 
		else
			h=htmp;
		
		/* Check for underflow condition */
		if ((1.0 + h) == 1.0) 
		{
			fputs("Underflow in RK stepper!\n", stderr);
			return 2;
		}
		//if (h<1e-8) 
		//{
		//	fprintf(stderr,"Underflow in RK stepper! i=%d,h=%e\n",i,h);
			//return 2;
		//	err=1;
		//}
		//if(err==1)
		//{
		//	fprintf(stderr,"err next i=%d,h=%e\n",i,h);

		//}
		i++;
		//if(i>14)
		{
			//fprintf(stderr,"Underflow in RK stepper! i=%d,h=%e\n",i,h);
			//return 2;
		}
	}
	
	/* See if the step size can safely be increased */
	if (errmax > ERRCON) 
		*hnext = SAFETY*h*pow(errmax,PGROW);
	else 
		*hnext = 5.0*h;
	
	/* Advise driver of step size just used */
	*hdid = h;
	
	/* Set new values */
	x[0] = xtemp[0]; 
	x[1] = xtemp[1];
	x[2] = xtemp[2];
	
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
	if (ierr = getBfield(xtemp, ak2, hptr)) 
		return ierr;
	xtemp[0] = x[0] + h*(b31*dxdt[0] + b32*ak2[0]);
	xtemp[1] = x[1] + h*(b31*dxdt[1] + b32*ak2[1]);
	xtemp[2] = x[2] + h*(b31*dxdt[2] + b32*ak2[2]);
	
	/* Third step */
	if (ierr = getBfield(xtemp, ak3, hptr)) 
		return ierr;
	xtemp[0] = x[0] + h*(b41*dxdt[0] + b42*ak2[0] + b43*ak3[0]);
	xtemp[1] = x[1] + h*(b41*dxdt[1] + b42*ak2[1] + b43*ak3[1]);
	xtemp[2] = x[2] + h*(b41*dxdt[2] + b42*ak2[2] + b43*ak3[2]);
	
	/* Fourth step */
	if (ierr = getBfield(xtemp, ak4, hptr)) 
		return ierr;
	xtemp[0] = x[0] + h*(b51*dxdt[0] + b52*ak2[0] + b53*ak3[0] + b54*ak4[0]);
	xtemp[1] = x[1] + h*(b51*dxdt[1] + b52*ak2[1] + b53*ak3[1] + b54*ak4[1]);
	xtemp[2] = x[2] + h*(b51*dxdt[2] + b52*ak2[2] + b53*ak3[2] + b54*ak4[2]);
	
	/* Fifth step */
	if (ierr = getBfield(xtemp, ak5, hptr))
		return ierr;
	xtemp[0] = x[0] + h*(b61*dxdt[0] + b62*ak2[0] + b63*ak3[0] + b64*ak4[0] +
						 b65*ak5[0]);
	xtemp[1] = x[1] + h*(b61*dxdt[1] + b62*ak2[1] + b63*ak3[1] + b64*ak4[1] +
						 b65*ak5[1]);
	xtemp[2] = x[2] + h*(b61*dxdt[2] + b62*ak2[2] + b63*ak3[2] + b64*ak4[2] +
						 b65*ak5[2]);
	
	/* Sixth step */
	if (ierr = getBfield(xtemp, ak6, hptr)) 
		return ierr;
	
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
	
	*val = (float)(w[0]*var[v[0]] + w[1]*var[v[1]] + w[2]*var[v[2]]);
	
	return 0;
}


/*============================================================*/
int getBfield(double *x1, double *dx, hdatatype *hptr)
{
	double Rinv, relphi = x1[1]/hptr->dphi;
	int    nplanes = hptr->planes[0];
	int    baseplane = (int)(relphi + 0.5);
	int    lastplane = baseplane-1, nextplane = baseplane+1;
	int    v[3];
	double tw1, tw2, tw3, w[3], factor;
	double Bzlast, Bzthis, Bznext;
	double BRlast, BRthis, BRnext, Bphilast, Bphithis, Bphinext;
	
	Rinv = 1.0/(*x1 + hptr->R[0]);
	
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
	
	BRlast = w[0]*hptr->BR[v[0] + lastplane] + w[1]*hptr->BR[v[1] + lastplane] +
    w[2]*hptr->BR[v[2] + lastplane];
	Bphilast = (w[0]*hptr->R[v[0] + lastplane]*hptr->Bphi[v[0] + lastplane] +
				w[1]*hptr->R[v[1] + lastplane]*hptr->Bphi[v[1] + lastplane] +
				w[2]*hptr->R[v[2] + lastplane]*hptr->Bphi[v[2] + lastplane])*Rinv;
	Bzlast = w[0]*hptr->Bz[v[0] + lastplane] + w[1]*hptr->Bz[v[1] + lastplane] +
    w[2]*hptr->Bz[v[2] + lastplane];
	
	BRthis = w[0]*hptr->BR[v[0]] + w[1]*hptr->BR[v[1]] + w[2]*hptr->BR[v[2]];
	Bphithis = (w[0]*hptr->R[v[0]]*hptr->Bphi[v[0]] +
				w[1]*hptr->R[v[1]]*hptr->Bphi[v[1]] +
				w[2]*hptr->R[v[2]]*hptr->Bphi[v[2]])*Rinv;
	Bzthis = w[0]*hptr->Bz[v[0]] + w[1]*hptr->Bz[v[1]] + w[2]*hptr->Bz[v[2]];
	
	BRnext = w[0]*hptr->BR[v[0] + nextplane] + w[1]*hptr->BR[v[1] + nextplane] +
    w[2]*hptr->BR[v[2] + nextplane];

	Bphinext = (w[0]*hptr->R[v[0] + nextplane]*hptr->Bphi[v[0] + nextplane] +
				w[1]*hptr->R[v[1] + nextplane]*hptr->Bphi[v[1] + nextplane] +
				w[2]*hptr->R[v[2] + nextplane]*hptr->Bphi[v[2] + nextplane])*Rinv;
	Bznext = w[0]*hptr->Bz[v[0] + nextplane] + w[1]*hptr->Bz[v[1] + nextplane] +
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
	double     f, dx, dy, dR1, dR2, dz1, dz2, R;
	float      xA, xB, xC, yA, yB, yC;
	static int el=0;
	int        next, inflag;
	
	R = *x1 + hptr->R[0];
	int i=0;
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
		
		xA =	(float)hptr->R[v[0]];
		yA =	hptr->z[v[0]];
		xB =	(float)hptr->R[v[1]];
		yB =	hptr->z[v[1]];
		
		if ((xB - xA)*(x1[2] - yA) < (yB - yA)*(R - xA)) 
			goto lout0;
		
	lin1:
		xC = (float)hptr->R[v[2]];
		yC = hptr->z[v[2]];
		dR1 = xC - xB;  dz1 = yC - yB;
		if (dR1*(x1[2] - yB) < dz1*(R - xB)) 
			goto lout1;
		
	lin2:
		dx = R - xC; 
		dy = x1[2] - yC;
		dR2 = xA - xC;
		dz2 = yA - yC;
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
		} 
		inflag=0; 
		goto lin1;
		
	lout1:
		if ((next = hptr->neighbor_list[3*el + 1]) >= 0) 
		{
			el = next;        /* Try the nearest neighbor */
			continue;
		} 
		inflag=0; 
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
	edge *edge_list;
	int  *tri, *nlist, *conn_list=hptr->conn_list;
	int  v, el, nel=hptr->ncells, vpp=hptr->vpp;
	
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
	int    vert;
	double Rmaj, Redg, rmin, p0, p00, pm, pp, pn, f0, fn;
	
	/* Find R at edge on midplane in target poloidal plane */
	Rmaj = hptr->R[setptr->plane0*hptr->vpp];
	
	for (vert=1, Redg=0.0; vert<hptr->vpp; vert++)
	{
		if (hptr->R[vert + setptr->plane0*hptr->vpp] > Redg)
			Redg = hptr->R[vert + setptr->plane0*hptr->vpp];
	}
	rmin = Redg - Rmaj;
	if (setptr->aguess < 0.0)
	{
		p00 = p0 = Rmaj - hptr->R[0];
		pm = p0 - 0.55*rmin;
		pp = p0 + 0.55*rmin;
	}
	else
	{
		p00 = p0 = setptr->aguess - hptr->R[0];
		pm = p0 - 0.05*rmin;
		pp = p0 + 0.05*rmin;
	}
	
	f0 = get_excursion(p0, hptr, setptr);

	/* for (vert=0; vert<41; vert++) {
     pn = pm + (vert/40.0)*(pp - pm);
     fn = get_excursion(pn, hptr, setptr);
     printf("f(%lf) = %le\n", pn, fn);
     }*/
	while (pp - pm > 1.0e-6*rmin)
	{
		/* printf("%lf, %lf, %lf (exc0=%le)\n", pm,p0,pp,f0); */
		if ((pp - p0) > (p0 - pm))
		{
			pn = p0 + 0.381966*(pp - p0);
			fn = get_excursion(pn, hptr, setptr);
		}
		else
		{
			pn = p0; fn = f0;
			p0 -= 0.381966*(p0 - pm);
			f0 = get_excursion(p0, hptr, setptr);
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
	}
	//fprintf(stderr, "Best guess: %lf vs %lf (exc=%le)\n", p0, p00, f0);
	
	return p0;
}


/*============================================================*/
double get_excursion(double R0, hdatatype *hptr, settype *setptr)
{
	long int  N=1L;
	double    h, hdid, hnext, hmax;
	double    xsum, xssum, ysum, yssum;
	double    x0[3], dx[3], xscal[3], phi0, xrel;
	int       transit, last=0;
	const int maxtrans = 10;
	
	/* fprintf(stderr, "Getting excursion for R=%lf...\n", R0);*/
	x0[0] = R0;
	xsum = xssum = 0.0;
	x0[1] = phi0 = setptr->plane0*hptr->dphi;
	x0[2] = ysum = yssum = 0.0;
	
	xscal[0] = xscal[2] = 1.0;
	xscal[1] = hptr->dphi * hptr->planes[0];
	h = 0.001*xscal[0];
	hmax = 0.01;
	
	do {
		if (getBfield(x0, dx, hptr)) break;           /* Evaluate derivatives */
		if (nrrkqs(x0, dx, h, 1.0e-9, xscal, &hdid, &hnext, hptr)) break; /* stepper */
		h = (hnext > hmax) ? hmax : hnext;
		
		transit = (int)((x0[1] - phi0)/xscal[1]);
		if (transit != last) {
			xrel = x0[0] - R0;
			xsum += xrel;
			xssum += xrel*xrel;
			ysum += x0[2];
			yssum += x0[2]*x0[2];
			N++;
		}
		last = transit;
		
	} while (N < maxtrans);
	
	if (N < maxtrans) return 1.0e+6;  /* R0 is out of range. */
	
	xsum /= N;  xssum /= N;
	ysum /= N;  yssum /= N;
	/* fprintf(stderr, "x_av = %lf; variance = %le\n", xsum + R0, xssum-xsum*xsum);
     fprintf(stderr, "y_av = %lf; variance = %le\n", ysum, yssum-ysum*ysum); */
	return sqrt(fabs(xssum - xsum*xsum + yssum - ysum*ysum));
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
			ptr[i] = (float)(factor0*buf[i] +
							 ccoef[i]*hptr->ctab[(mode*plane) % hptr->planes[0]] +
							 scoef[i]*hptr->stab[(mode*plane) % hptr->planes[0]]);
	
	for (i=0; i<hptr->vpp; i++)
		buf[i] = (float)(factor0*buf[i] + ccoef[i]);
	
	free(scoef); free(ccoef);
}
/*============================================================*/
int myround(double x)
{
	int y = (int)x;
	
	if (x - y < 0.5) 
		return y;
	else 
 		return y+1;
}

