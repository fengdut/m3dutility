/*feng wang for m3dout wxt??? particle information to distribution function */

#include<stdio.h>
#include<iostream>
#include<fstream>
#include"AllocArray.h"
#include"netcdfcpp.h"
static const int NC_ERR =3;

int main(int argc, char** argv)
{
	using namespace std;
	cout<<"distribution function"<<endl;
		
	float ***distf;
	float ***distf0;
	double pphi0=0,pphi1=0,pangle0=0,pangle1=0,e0=0,e1=0;
	double dpphi=0,dpangle=0,de=0;
	
	int n_pphi=30,n_pangle=30,n_e=30;
	
	char ifilename[100]="fort.45";
	char ofilename[100]="fort.45_30_30_30.nc";
	int  iofile=0;
	//sprintf(ofilename,"dist3D_%d_%d_%d.nc",n_pphi,n_pangle,n_e);
	cout<<"Usage: <particle data: fort.45/wxt??? [<options>]>"<<endl;
	cout<<"Options:"<<endl;
	cout<<"-i[nputfile] particle data filename (default fort.45)"<<endl;
	cout<<"-o[utputfile] 3D distribution (default: "<<ofilename<<")"<<endl;
	cout<<"-p[npphi] pphi grid No. (default 30)"<<endl;
	cout<<"-a[npangle] pitch angle grid No. (default 30)"<<endl;
	cout<<"-e[nenergy] energy angle grid No. (default 30)"<<endl;

	for(int arg=1;arg<argc;arg++)
	{
		if(!strncmp(argv[arg],"-i",2))
		{
			strcpy(ifilename,argv[++arg]);
			continue;
		}
		if(!strncmp(argv[arg],"-o",2))
		{
			iofile=1;
			strcpy(ofilename,argv[++arg]);
			continue;
		}		
		if(!strncmp(argv[arg],"-p",2))
		{
			n_pphi=atoi(argv[++arg]);
			continue;
		}		
		if(!strncmp(argv[arg],"-a",2))
		{
			n_pangle=atoi(argv[++arg]);
			continue;
		}
		if(!strncmp(argv[arg],"-e",2))
		{
			n_e=atoi(argv[++arg]);
			continue;
		}
		
	}
	if(!iofile)
		sprintf(ofilename,"%s_%d_%d_%d.nc",ifilename,n_pphi,n_pangle,n_e);

	int n_particle=0;
	double **particle_info;
	int** particle_id;
	char buf[200];
	ifstream pfile(ifilename);
 	if(pfile.is_open())
	{	
		while(!pfile.eof())
		{
			pfile.getline(buf,200);
			n_particle++;
	//		cout<<"n_particle"<<n_particle<<buf<<endl;
		}
		n_particle--;
	cout<<"Particle No: "<<n_particle<<endl;
	}
	else
	{
		cout<<"Error opening file: "<<ifilename<<endl;
		return 1;
	
	}
	pfile.clear();
	pfile.seekg(0,ios::beg);
		
	Alloc2D(particle_id,n_particle,2);
	Alloc2D(particle_info,n_particle,6);
        
	for(int i=0;i<n_particle;i++) 
	{
       		/*if(pfile.eof())
		{
			cout<<"end of the file"<<endl;
			return 2;
		}*/
		pfile.getline(buf,200);
		sscanf(buf,"%d %*s %lf %lf %lf %lf %lf %lf",&particle_id[i][0]
							,&particle_info[i][0],&particle_info[i][1],&particle_info[i][2],
							&particle_info[i][3],&particle_info[i][4],&particle_info[i][5]);	
		particle_info[i][1]=particle_info[i][1]/particle_info[i][2];
		if(0==i)
		{
			pphi0=pphi1=particle_info[i][0];
			pangle0=pangle1=particle_info[i][1];
			e0=e1=particle_info[i][2];
		}
		pphi0=pphi0<particle_info[i][0] ? pphi0:particle_info[i][0];
		pphi1=pphi1>particle_info[i][0] ? pphi1:particle_info[i][0];
		pangle0  =pangle0<particle_info[i][1] ? pangle0:particle_info[i][1];
		pangle1  =pangle1>particle_info[i][1] ? pangle1:particle_info[i][1];
		e0   =e0<particle_info[i][2] ? e0:particle_info[i][2];
		e1   =e1>particle_info[i][2] ? e1:particle_info[i][2];
	
	}
	double t=0.05*(pphi1-pphi0);
	pphi0	-=	t;
	pphi1	+=	t;
//	pphi1=0.0;
//	pphi0=-1.1;
	
	t=0.05*(pangle1-pangle0);
	pangle0-=t;
	pangle1+=t;
//	pangle1=1.5;

	t=0.05*(e1-e0)+0.1;
	
	e0-=t;
	e1+=t;
//	e1=20;

	int	i=n_particle-1;
	e0=0;
	pangle0=0;
	
	cout<<"finish read particle information:\t"<<ifilename<<endl;
	cout<<"Max pphi: "<<pphi1<<"\t Min pphi:\t"<<pphi0<<endl;
	cout<<"Max pangle: "<<pangle1<<"\t Min pangle:\t"<<pangle0<<endl;
	cout<<"Max Energy: "<<e1<<"\t Min Energy:\t"<<e0<<endl;
	
	assert(n_pphi>0&&n_pangle>0&&n_e>0);	
	Alloc3D(distf,n_pphi+1,n_pangle+1,n_e+1);
	Alloc3D(distf0,n_pphi+1,n_pangle+1,n_e+1);
	
	dpphi=(pphi1-pphi0)/(double)(n_pphi);
	float *pphi_coor=new float[n_pphi+1];
	for(i=0;i<=n_pphi;i++)
	{
		pphi_coor[i]=pphi0+dpphi*i;
	}
	dpangle=(pangle1-pangle0)/(double)(n_pangle);	
	float *pangle_coor=new float[n_pangle+1];
	for(i=0;i<=n_pangle;i++)
	{
		pangle_coor[i]=pangle0+dpangle*i;
	}
	de	=(e1-e0)/(double)(n_e);
	float *e_coor=new float[n_e+1];
	for(i=0;i<=n_e;i++)
	{
		e_coor[i]=e0+de*i;
	}
	
	int nx,ny,nz;
	double x0,x1,y0,y1,z0,z1;
	double *tp;
	double tv=dpphi*dpangle*de;
	double tw=0;
	for(i=0;i<n_particle;i++)
	{	
		
		tp=particle_info[i];
		nx=(int)((tp[0]-pphi0)/dpphi);
		x0=tp[0]-((double)nx*dpphi+pphi0);
		x1=dpphi-x0;
	
		ny=(int)((tp[1]-pangle0)/dpangle);
		y0=tp[1]-((double)ny*dpangle+pangle0);
		y1=dpangle-y0;

		nz=(double)((tp[2]-e0)/de);
		z0=tp[2]-((double)nz*de+e0);
		z1=de-z0;
		if(nx>=n_pphi||ny>=n_pangle||nz>=n_e)
			cout<<nx<<ny<<nz<<endl;
		assert(nx<n_pphi&&ny<n_pangle&&nz<n_e);
		assert(nx>=0&&ny>=0&&nz>=0);
		distf[nx][ny][nz]	+=(x1*y1*z1/tv/tv)*tp[4];
		distf[nx][ny][nz+1]	+=(x1*y1*z0/tv/tv)*tp[4];
		distf[nx][ny+1][nz]	+=(x1*y0*z1/tv/tv)*tp[4];
		distf[nx][ny+1][nz+1]	+=(x1*y0*z0/tv/tv)*tp[4];
		distf[nx+1][ny][nz]	+=(x0*y1*z1/tv/tv)*tp[4];
		distf[nx+1][ny][nz+1]	+=(x0*y1*z0/tv/tv)*tp[4];
		distf[nx+1][ny+1][nz]	+=(x0*y0*z1/tv/tv)*tp[4];
		distf[nx+1][ny+1][nz+1]	+=(x0*y0*z0/tv/tv)*tp[4];

		distf0[nx][ny][nz]	+=(x1*y1*z1/tv/tv);
                distf0[nx][ny][nz+1]	+=(x1*y1*z0/tv/tv);
                distf0[nx][ny+1][nz]	+=(x1*y0*z1/tv/tv);
                distf0[nx][ny+1][nz+1]	+=(x1*y0*z0/tv/tv);
                distf0[nx+1][ny][nz]	+=(x0*y1*z1/tv/tv);
                distf0[nx+1][ny][nz+1]	+=(x0*y1*z0/tv/tv);
                distf0[nx+1][ny+1][nz]	+=(x0*y0*z1/tv/tv);
		distf0[nx+1][ny+1][nz+1]+=(x0*y0*z0/tv/tv);

		tw	+= tp[4];
	}


	for(i=0;i<(n_pphi+1)*(n_pangle)*(n_e);i++)
	{
		if(distf0[0][0][0]>0)
		distf[0][0][i]=distf[0][0][i]/distf0[0][0][i];
	}

	cout<<"total weight is:"<<tw<<endl;
	
	NcError err(NcError::verbose_nonfatal);
	
	NcFile distfile(ofilename,NcFile::Replace);
	if(!distfile.is_valid())
	{
		cout<<"creat netcdf file error"<<endl;
		return 2;
	}
		
	NcDim *pphiDim,*pangleDim,*eDim;
	if(!(pphiDim=distfile.add_dim("pphi",n_pphi+1)))
		return NC_ERR;
	if(!(pangleDim=distfile.add_dim("pangle",n_pangle+1)))
		return NC_ERR;
	if(!(eDim=distfile.add_dim("energy",n_e+1)))
		return NC_ERR;

	NcVar *pphi_c,*pangle_c,*e_c;
	if(!(pphi_c=distfile.add_var("pphi",ncFloat,pphiDim)))
		return NC_ERR;
	if(!(pangle_c=distfile.add_var("pangle",ncFloat,pangleDim)))
		return NC_ERR;
	if(!(e_c=distfile.add_var("e",ncFloat,eDim)))	
		return NC_ERR;	
	
	NcVar *distVar;
	if(!(distVar=distfile.add_var("f0",ncFloat,pphiDim,pangleDim,eDim)))
		return NC_ERR;
	NcVar *distVar0;
	
	if(!(distVar0=distfile.add_var("f00",ncFloat,pphiDim,pangleDim,eDim)))
		return NC_ERR;


	if(!pphi_c->put(pphi_coor,n_pphi+1))
		return NC_ERR;
	if(!pangle_c->put(pangle_coor,n_pangle+1))
		return NC_ERR;
	if(!e_c->put(e_coor,n_e+1))
		return NC_ERR;
		
	if(!distVar->put(distf[0][0],n_pphi+1,n_pangle+1,n_e+1))
		return NC_ERR;

	if(!distVar0->put(distf0[0][0],n_pphi+1,n_pangle+1,n_e+1))
		return NC_ERR;


	cout<<"write date to:\t"<<ofilename<<endl;
	cout<<"finish"<<endl;		
	
	delete[] e_coor;
	delete[] pphi_coor;
	delete[] pangle_coor;
	Free2D(particle_id);
	Free2D(particle_info);
	Free3D(distf);
	Free3D(distf0);
	return 0;
}

