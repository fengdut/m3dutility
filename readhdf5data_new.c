
#define H5_USE_16_API
#include <hdf5.h>
#include"typedefine.h"
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

/*============================================================*/
void readHDF5scalar(int idnum, int nvars,char *name, float *buf, hid_t groupid)
{
  char     node_name[64], vname[64], flag;
  int      var;
  hid_t    atid, varid, datatype, dataid, dataspace, mem_space;
  hsize_t  npoints, count[2];
  hssize_t offset[2];

  /* Find number of variables */

  if (idnum >= nvars) {  /* Error check */
    fputs("Error: Variable index out of range in readHDF5scalar.\n", stderr);
    return;
  }

  if (idnum < 0) {  /* Load by name; search all vars to get index */
    flag = 0;
    for (var=0; var<nvars; var++) {
      /* Find the name of this variable number */
      sprintf(node_name, "node_data[%d]", var);
      varid = H5Gopen(groupid, node_name);
      atid = H5Aopen_name(varid, "labels");
      datatype = H5Aget_type(atid);
      H5Aread(atid, datatype, vname);
      H5Tclose(datatype); H5Aclose(atid);

      /* Compare it to the target */
      if (!strncmp(vname, name, strlen(name))) {  /* Match found */
	flag = 1; break;
      } /* end if !strn... */
      H5Gclose(varid);
    } /* end loop var */
    if (!flag) {
      fprintf(stderr, "Error: Variable \"%s\" not found.\n", name);
      return;
    }
    fprintf(stderr, "\nReading variable \"%s\"...\n", name);
  } else {          /* Load by number */
    sprintf(node_name, "node_data[%d]", idnum);
    varid = H5Gopen(groupid, node_name);
    fprintf(stderr, "\nReading variable %d...\n", idnum);
  } /* endif idnum... */

  /* Open the dataset for reading */
  dataid = H5Dopen(varid, "values");

  /* Determine the size of the data */
  dataspace = H5Dget_space(dataid);
  npoints = H5Sget_simple_extent_npoints(dataspace);
  fprintf(stderr, "%ld points in data set.\n", (long)npoints);

  /* Set up memory dataspace for reading */
  offset[0] = offset[1] = 0;
  count[0] = npoints;   count[1] = 1;
  mem_space = H5Screate_simple(2, count, count);
  H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset, NULL, count, NULL);

  /* Set up file dataspace for reading */
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

  /* Read data */
  H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,
	  (void *)buf);

  fputs("Scalar data read.\n", stderr);
  fprintf(stderr, "Value 0 = %e\n", *buf);

  /* Close the memory and file dataspaces, dataset and variable group */
  H5Sclose(mem_space); H5Sclose(dataspace);
  H5Dclose(dataid);    H5Gclose(varid);
}
/* Read time stamp */
float readHDF5time(hid_t h5fid, int frame, hdatatype *hptr)
{
  hid_t groupid, atid;
  float *times, t;
  char node_name[50];
   sprintf(node_name, "/time_node_data[%d]", frame);
  groupid = H5Gopen(h5fid, node_name);  /* Open the root group, where "time" lives */
  if ((atid = H5Aopen_name(groupid, "time")) == -1) { /* Attribute not found */
    H5Gclose(groupid);
    return -1.0;
  }
  if (frame < hptr->nframes) {
    times = (float *)malloc(hptr->nframes * sizeof(float)); /* Allocate */
    H5Aread(atid, H5T_NATIVE_FLOAT, (void *)times); /* Read time info */
    t = times[0];
    printf("time= %f\n",times[0]);
    free(times);                                  /* Deallocate */
  } else t = -1.0;
  H5Aclose(atid); H5Gclose(groupid);           /* Clean up */
  return t;
}

void read_plane_frame(hid_t h5fid,int frame,hdatatype *hdata)
{
	hid_t  groupid,  dataid, dataspace;

	groupid = H5Gopen(h5fid, "/");
	dataid  = H5Aopen_name(groupid, "nsteps");
	H5Aread(dataid, H5T_NATIVE_INT_g, &((*hdata).nframes));
	H5Aclose(dataid);

	(*hdata).tim = readHDF5time(h5fid,frame,hdata);
 groupid = H5Gopen(h5fid, "/planes");
        dataid = H5Dopen(groupid, "values");
        dataspace = H5Dget_space(dataid);
        H5Sselect_all(dataspace);
        H5Dread(dataid, H5T_NATIVE_INT_g, dataspace, dataspace, H5P_DEFAULT,(*hdata).planes);
        H5Sclose(dataspace);
        H5Dclose(dataid);
        H5Gclose(groupid);

}

void read_connection(hid_t h5fid,hdatatype *hdata)
{
	hid_t  groupid,  atid,  dataid, dataspace, mem_space;
	hsize_t       count[2], stride[2];
	hssize_t      offset[2];
	groupid = H5Gopen(h5fid, "/cell_set[0]");
	atid = H5Aopen_name(groupid, "ncells");
	H5Aread(atid, H5T_NATIVE_INT_g, (void *)&(hdata->ncells));
	H5Aclose(atid);
	hdata->ncells /= hdata->planes[0];
	fprintf(stderr, "%d triangles/plane.\n", hdata->ncells);
	hdata->conn_list = (int *) malloc(3 * hdata->ncells * sizeof(int));
	dataid = H5Dopen(groupid, "node_connect_list");
	dataspace = H5Dget_space(dataid);
	offset[0] = 0;
	offset[1] = 0;
	stride[0] = hdata->planes[1];
	stride[1] = 1;
	count[0]  = hdata->ncells;
	count[1] = 3;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, stride, count, NULL);
	mem_space = H5Screate_simple(2, count, count);
	H5Sselect_all(mem_space);
	H5Dread(dataid, H5T_NATIVE_INT_g, mem_space, dataspace, H5P_DEFAULT,hdata->conn_list);
	fputs("Connectivity data read.\n", stderr);
	H5Sclose(mem_space);  
	H5Sclose(dataspace);
	H5Dclose(dataid); 
	H5Gclose(groupid);
}

void read_xyz(hid_t h5fid,hdatatype *hdata,double *parity ,hsize_t  *npoints)
{
	hid_t  groupid, dataid, dataspace, mem_space;
	char   node_name[64];
	hsize_t  count[2];
	hssize_t      offset[2];
	int i;
	/* Determine group name for coords. in this time slice, open the group */
	sprintf(node_name, "/time_coordinates[%d]/coordinates",0);
	if ((groupid = H5Gopen(h5fid, node_name)) < 0)
	{
	  exit(2);
	}

	/* Open the dataset for reading */
	dataid = H5Dopen(groupid, "values");

	/* Determine the size of the data */
	dataspace = H5Dget_space(dataid);
	(*npoints) = H5Sget_simple_extent_npoints(dataspace);

	/* Allocate memory for data */
	(*npoints) /= 3;
/*	hdata->x = (float *)malloc((*npoints) * sizeof(float));
	hdata->y = (float *)malloc((*npoints) * sizeof(float));
	hdata->z = (float *)malloc((*npoints) * sizeof(float));
	hdata->R = (double *)malloc((*npoints) * sizeof(double));
*/
	/* Set up memory dataspace for reading */
	offset[0] = 0;
	offset[1] = 0;
	count[0] = (*npoints);
	count[1] = 1;
	mem_space = H5Screate_simple(2, count, count);
	H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset, NULL, count, NULL);

	/* Set up file dataspace for reading */
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/* Read data */
	H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->x);
	offset[1] = 1;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->y);
	offset[1] = 2;
	H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
	H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->z);

	fprintf(stderr, "Coordinate data read for %d planes, %d vertices/plane.\n", hdata->planes[0], (int)(*npoints)/hdata->planes[0]);
	//fprintf(stderr, "x[0] = (%e, %e, %e)\n", hdata->x[0], hdata->y[0], hdata->z[0]);

	/* Test parity */
	if (hdata->y[(*npoints)/hdata->planes[0]] < 0.0) 
	{
	fputs("Switching to M3D coordinates (R,z,phi)\n", stderr);
	(*parity) = -1.0;
	}

	for (i=0; i<(*npoints); i++)
	hdata->R[i] = sqrt(hdata->x[i]*hdata->x[i] + hdata->y[i]*hdata->y[i]);
//	free(hdata->x);
//	free(hdata->y);

	/* Close the memory and file dataspaces, the dataset, and the group */
	H5Sclose(mem_space); 
	H5Sclose(dataspace);
	H5Dclose(dataid); 
	H5Gclose(groupid);

}

void read_B(hid_t h5fid,int frame,hsize_t npoints,double parity,int eqflag,hdatatype *hdata)
{
	char   node_name[64],fname[64];
	hid_t  groupid, dataid, varid,atid,dataspace, mem_space,datatype;
	hsize_t  count[2];
	hssize_t      offset[2];
	int i,nvars,var,veclen;
	/* Determine the group name for data in this time slice, open the group */
	groupid = H5Gopen(h5fid, "/");
	atid = H5Aopen_name(groupid, "nnode_data");
	H5Aread(atid, H5T_NATIVE_INT_g, (void *)&nvars);
	H5Aclose(atid);
        H5Gclose(groupid);

	sprintf(node_name, "/time_node_data[%d]", frame);
	groupid = H5Gopen(h5fid, node_name);

	//printf("Size of psi = %d\n", (int)npoints);
	hdata->psi = (float *)malloc((int)npoints * sizeof(float));
	readHDF5scalar(-1, nvars,"psi", hdata->psi, groupid);

	/* Find number of variables */
	fprintf(stderr, "%d variables in this time slice.\n", nvars);

	/* Search for vectors... */
	for (var=0; var<nvars; var++) 
	{
		/* Determine the group name for this variable, open the group */
		sprintf(node_name, "node_data[%d]", var);
		varid = H5Gopen(groupid, node_name);

		/* Check vector length */
		atid = H5Aopen_name(varid, "veclen");
		H5Aread(atid, H5T_NATIVE_INT_g, (void *)&veclen);
		H5Aclose(atid);
		if (veclen == 3) 
		{
			//fprintf(stderr, "Subgroup %s: veclen = %d.\n", node_name, veclen);

			/* Check name; make sure it's the B-field... */
			atid = H5Aopen_name(varid, "labels");
			datatype = H5Aget_type(atid);
			H5Aread(atid, datatype, fname);
			H5Tclose(datatype);
			H5Aclose(atid);
			//fprintf(stderr, "Data label is \"%s\".\n", fname);
			if (*fname == 'B') 
			{
				/* Open the dataset for reading */
				dataid = H5Dopen(varid, "values");

				/* Determine the size of the data */
				dataspace = H5Dget_space(dataid);
				hdata->npoints = (int)(npoints = H5Sget_simple_extent_npoints(dataspace) / 3);

				/* Allocate memory for data */
				hdata->Bx = (float *)malloc(npoints * sizeof(float));
				hdata->By = (float *)malloc(npoints * sizeof(float));
				hdata->Bz = (float *)malloc(npoints * sizeof(float));
				//fprintf(stderr, "%d floats allocated for each component.\n",(int)npoints);

				/* Set up memory dataspace for reading */
				offset[0] = 0; 
				offset[1] = 0;
				count[0] = npoints; 
				count[1] = 1;
				mem_space = H5Screate_simple(2, count, count);
				H5Sselect_hyperslab(mem_space, H5S_SELECT_SET, offset, NULL,count, NULL);

				/* Set up file dataspace for reading */
				H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,count, NULL);

				/* Read data */
				H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->Bx);
				offset[1] = 1;
				H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,count, NULL);
				H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->By);
				offset[1] = 2;
				H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,count, NULL);
				H5Dread(dataid, H5T_NATIVE_FLOAT, mem_space, dataspace, H5P_DEFAULT,hdata->Bz);

				fputs("Field data read.\n", stderr);
				//fprintf(stderr, "B[0] = %e x  +  %e y  +  %e z\n",hdata->Bx[0], hdata->By[0], hdata->Bz[0]);

				/* Close the memory and file dataspaces, the dataset, and the group */
				H5Sclose(mem_space);
				H5Sclose(dataspace);
				H5Dclose(dataid);
				H5Gclose(varid);

				/* Calculate data-based quantities */
				hdata->vpp = hdata->npoints / hdata->planes[0];  /* vertices per plane */
				hdata->dphi = 2.0*M_PI/(double)hdata->planes[0];  /* dphi btwn planes */

				/* Build trig tables */
				hdata->ctab = (double *)malloc((hdata->planes[0]+2) * sizeof(double))+1;
				hdata->stab = (double *)malloc((hdata->planes[0]+2) * sizeof(double))+1;
				for (i=-1; i <= hdata->planes[0]; i++) 
				{
					hdata->ctab[i] = cos(i*hdata->dphi);
					hdata->stab[i] = sin(i*hdata->dphi);
				}

				/* Convert field values to cylindrical coordinates */
				hdata->BR   = (float *)malloc(npoints * sizeof(float));
				hdata->Bphi = (float *)malloc(npoints * sizeof(float));
				for (i=0; i<hdata->planes[0]; i++) 
				{
					for (var=0; var<hdata->vpp; var++) 
					{
						hdata->BR[i*hdata->vpp + var]	=	hdata->ctab[i]*hdata->Bx[i*hdata->vpp + var] +
															parity*hdata->stab[i]*hdata->By[i*hdata->vpp + var];
						hdata->Bphi[i*hdata->vpp + var] =	parity*hdata->ctab[i]*hdata->By[i*hdata->vpp + var] -
															hdata->stab[i]*hdata->Bx[i*hdata->vpp + var];
					}
				}
				free(hdata->Bx);
				free(hdata->By);


				/* Zero Bp on axis if equilibrium */
			if (eqflag) 
			{
				fputs("Setting Bp to zero at vertex zero.\n", stderr);
				for (i=0; i<hdata->planes[0]; i++) 
				{
					hdata->BR[i*hdata->vpp] = 0.0;
					hdata->Bz[i*hdata->vpp] = 0.0;
				}
			}
		}

			/* Find element neighbors */
		findElementNeighbors(hdata);
		
		break;  /* ignore the rest of the variables in this timeslice */
	}
    }

    /* Close the variable group */

	 H5Gclose(groupid);
  } /* end loop var */


