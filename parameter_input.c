#include"typedefine.h"
#include"stdio.h"
settype parameter_input()
{

	settype       settings;
	  /* Default settings */
	settings.plane0 = 0;		/* Index of plane for Poincare plotting */
	settings.poindir=0;		/* Direction for Poincare plotting (default=forward) */
	settings.numsurf=65;		/* Number of surfaces to trace out */
	settings.outerdots=9000;	/* Number of dots on outermost surface */
	settings.mindots=30;		/* Number of dots on innermost surface */
	settings.convtol=1.0e-9;
	settings.aguess=-1.0;
	settings.bguess=-1.0;
	settings.theta_sym = 0.0;
	settings.qsurf=128; 
	settings.qtrans=1024;
	settings.maxtrans = 20; /* Dots to determine axis position */
	settings.oflag=0;       /* Output filename specified? */
	settings.filtermode=-1; /* Toroidal mode number for filtering */
	settings.filetype=HDF5;
	settings.maxpass=8;
	settings.wingspan=0.1;
	settings.frame=0;
	return settings;
}

void printhelp(char *argv[],int frame,settype settings)
{
	
	/* Print usage message if no arguments are given */
    fprintf(stderr, "Usage: %s <HDF5 filename> [<options>]\n", *argv);
    fputs("\nOptions:\n", stderr);
    fprintf(stderr," -f[rame] <frame #>  - select timeslice from file (default %d).\n",frame);
    fputs(" -p[lane] <plane index>  - select plane for plot or q profile.\n",stderr);
    fputs(" -tol[erance] <tolerance>  -  convergence tol. for integration (default 1e-9).\n", stderr);
    fputs(" -eq  -  assume file is equilibrium (mag. axis at vert. 0).\n",stderr);
    fputs(" -g[uess] <radius>  -  initial guess for axis position.\n", stderr);
    fputs(" -m[ode] <n> - restrict B pert to toroidal mode # n.\n", stderr);
    fputs(" -th[eta] <angle> - angle of symmetry axis to midplane.\n", stderr);
    fprintf(stderr, " -a[xis] <n> - number of points for axis search (default %d).\n", settings.maxtrans);
    fprintf(stderr, " -n <n> - maximum axis search passes (default %d).\n", settings.maxpass);
    fprintf(stderr, " -w <width> - half-width of search radius for axis (default %le).\n",settings.wingspan);
    fprintf(stderr,"-np <number of planes per cpu>");

    /*fputs(" -H[DF5] - specifies HDF5 output instead of ascii.\n", stderr);*/
    fprintf(stderr," -s[urfaces] <# of Poincare surfaces to plot> (default %d).\n",settings.numsurf);
    fprintf(stderr," -d[ots] <# of dots on outermost surface> (default %d).\n",settings.outerdots);
    fprintf(stderr," -i[nner] <# of dots on innermost surface> (default %d).\n",settings.mindots);
    fputs(" -ba[ck]  - starting pts of surfaces go inboard from axis.\n",stderr);
    fputs(" -bo[th]  - starting pts of surfaces go in & out from axis.\n",stderr);
    fputs(" -o[utput] <output filename> (default \"poincare.dat\").\n",stderr);

}

void getoption(int argc, char *argv[],char *fname,int *eqflag,int * frame,settype * settings)
{
	int arg;
/* Parse command-line options */
  strcpy(fname, argv[1]);
  for (arg=2; arg<argc; arg++) 
  {
    if (!strncmp(argv[arg], "-f", 2)) {
      (*frame) = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-p", 2)) {
      (*settings).plane0 = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-tol", 4)) {
      (*settings).convtol = atof(argv[++arg]);
      fprintf(stderr, "Resetting tolerance to %lf.\n", (*settings).convtol);
      continue;
    }
    if (!strncmp(argv[arg], "-eq", 3)) {
      (*eqflag) = 1;
      continue;
    }
    if (!strncmp(argv[arg], "-o", 2)) {
      (*settings).oflag = 1;
      strcpy((*settings).outname, argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-g", 2)) {
      (*settings).aguess = atof(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-m", 2)) {
      (*settings).filtermode = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-th", 3)) {
      (*settings).theta_sym = atof(argv[++arg]);
      fprintf(stderr, "Resetting theta_sym to %lf.\n", (*settings).theta_sym);
      continue;
    }
    if (!strncmp(argv[arg], "-z", 2)) {
      (*settings).bguess = atof(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-a", 2)) {
      (*settings).maxtrans = atoi(argv[++arg]);
      fprintf(stderr, "Resetting axis search points to %d.\n", (*settings).maxtrans);
      continue;
    }
    if (!strncmp(argv[arg], "-n", 2)) {
      (*settings).maxpass = atoi(argv[++arg]);
      continue;
    }
    if(!strncmp(argv[arg], "-np", 4)) {
      (*settings).numsurf = atoi(argv[++arg]);
      continue;
    }

    if (!strncmp(argv[arg], "-w", 2)) {
      (*settings).wingspan = atof(argv[++arg]);
      continue;
    }

    /* if (!strncmp(argv[arg], "-H", 2) || !strncmp(argv[arg], "-h", 2)) {
       settings.filetype = HDF5;
       continue;
       }
    */
    if (!strncmp(argv[arg], "-s", 2)) {
      (*settings).numsurf = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-d", 2)) {
      (*settings).outerdots = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-i", 2)) {
      (*settings).mindots = atoi(argv[++arg]);
      continue;
    }
    if (!strncmp(argv[arg], "-ba", 3)) {
      (*settings).poindir = 1;
      fputs("Stepping inward.\n", stderr);
      continue;
    }
    if (!strncmp(argv[arg], "-bo", 3)) {
      (*settings).poindir = 2;
      fputs("Stepping inward and outward.\n", stderr);
      continue;
    }
    fprintf(stderr, "Unrecognized option %s.\n", argv[arg]);
  } /* end loop arg */
}
