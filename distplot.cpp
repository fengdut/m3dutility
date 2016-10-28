/*feng wang for m3dout wxt??? particle information to distribution function */

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include"AllocArray.h"
#include"netcdf.h"
#include <iomanip>
static const int NC_ERR =3;

#pragma pack(4)
struct particledata
{
	int begin;
	int ID;
	double psibar2,mu,epatcal,wp4,pswp,wp4x;
	int end;
};

void printhelp()
{
using namespace std;
        cout<<"example:distplot -i dist.dat -0i dist.0.dat  [<Options>]>"<<endl;
        cout<<"Options:"<<endl;
        cout<<"-i[nputfile] particle data filename (default dist.dat)"<<endl;
        cout<<"-0i[nputfile] init particle data filename (default dist.0.dat)"<<endl;
        cout<<"-o[utputfile] 3D distribution (default: \"dist.82_80_81.nc\")"<<endl;
        cout<<"-p[npphi] pphi grid No. (default 30)"<<endl;
        cout<<"-0p [minpphi] (default automatic)"<<endl;
        cout<<"-1p [maxpphi] (default automatic)"<<endl;
        cout<<"-a[npangle] pitch angle grid No. (default 30)"<<endl;
        cout<<"-0a [minmu] (default automatic)"<<endl;
        cout<<"-1a [maxmu] (default automatic)"<<endl;
        cout<<"-e[nenergy] energy angle grid No. (default 30)"<<endl;
        cout<<"-0e [minenergy] (default automatic)"<<endl;
        cout<<"-1e [maxenergy] (default automatic)"<<endl;
}


int main(int argc, char** argv)
{
	using namespace std;
	cout<<"distribution function"<<endl;
		
	float ***distfpass;
	float ***distf0pass;
	float ***distftrap;
	float ***distf0trap;
	double pphi0=0,pphi1=0,pangle0=0,pangle1=0,e0=0,e1=0;
	int ip0=0,ip1=0,ia0=0,ia1=0,ie0=0,ie1=0;
	double dpphi=0,dpangle=0,de=0;
	
        double maxwp4,minwp4;
         
	int n_pphi=82,n_pangle=80,n_e=81;
	
	char ifilename[100]="dist.dat";
        char ifilename0[100]="dist.0.dat";
	char ofilename[100]="dist.82_80_81.nc";
	int  iofile=0;

        printhelp();
	for(int arg=1;arg<argc;arg++)
	{
		if(!strncmp(argv[arg],"-i",2))
		{
			strcpy(ifilename,argv[++arg]);
			continue;
		}
		if(!strncmp(argv[arg],"-0i",3))
		{
			strcpy(ifilename0,argv[++arg]);
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
		if(!strncmp(argv[arg],"-0p",3))
                {
                        pphi0=atof(argv[++arg]);
			ip0=1;
                        continue;
                }
                if(!strncmp(argv[arg],"-1p",3))
                {
                        pphi1=atof(argv[++arg]);
			ip1=1;
                        continue;
                }
                if(!strncmp(argv[arg],"-0e",3))
                {
                        e0=atof(argv[++arg]);
			ie0=1;
                        continue;
                }
                if(!strncmp(argv[arg],"-1e",3))
                {
                        e1=atof(argv[++arg]);
			ie1=1;
                        continue;
                }
                if(!strncmp(argv[arg],"-0a",3))
                {
                        pangle0=atof(argv[++arg]);
			ia0=1;
                        continue;
                }
	        if(!strncmp(argv[arg],"-1a",3))
                {
                        pangle1=atof(argv[++arg]);
			ia1=1;
                        continue;
                }

	}

	if(!iofile)
		sprintf(ofilename,"%s_%d_%d_%d.nc",ifilename,n_pphi,n_pangle,n_e);
        cout<<"dist:"<<ifilename<<endl;
        cout<<"dist0:"<<ifilename0<<endl;


	int n_particle=0;

	int particleinfosize=sizeof(particledata);
	particledata *p_info;
	cout<<"size of data:"<<particleinfosize<<endl;
	int filebegin=0;
	int fileend=0;
	ifstream pfile(ifilename,ios::binary);
 	if(pfile.is_open())
	{	
		filebegin=pfile.tellg();
		pfile.seekg(0,ios::end);
		fileend=pfile.tellg();
		n_particle=(fileend-filebegin)/particleinfosize;

	cout<<"Particle No: "<<n_particle<<endl;
	}
	else
	{
		cout<<"Error opening file: "<<ifilename<<endl;
		return 1;
	
	}


	if(n_particle<1)
	{
		cout<<"particle == 0"<<endl;	
		return 1;
	}
	p_info=new particledata[n_particle];
	pfile.clear();
	pfile.seekg(0,ios::beg);
        pfile.read((char*)p_info,fileend-filebegin);


         maxwp4=minwp4=p_info[0].wp4;
 	int maxid=0;
	int minid=abs(p_info[0].ID);
	double maxdwp4=0;
	double tpphi0=0,tpphi1=0,tpangle0=0,tpangle1=0,te0=0,te1=0;
        for(int i=0;i<n_particle;i++)
        {
                p_info[i].epatcal =p_info[i].epatcal*p_info[i].epatcal*0.5;
		 if(0==i)
                {
                        tpphi0=tpphi1=p_info[i].psibar2;
                        tpangle0=tpangle1=p_info[i].mu;
                        te0=te1=p_info[i].epatcal;
                }

                tpphi0=tpphi0<p_info[i].psibar2 ? tpphi0:p_info[i].psibar2;
                tpphi1=tpphi1>p_info[i].psibar2 ? tpphi1:p_info[i].psibar2;
                tpangle0  =tpangle0<p_info[i].mu ? tpangle0:p_info[i].mu;
                tpangle1  =tpangle1>p_info[i].mu ? tpangle1:p_info[i].mu;
                te0   =te0<p_info[i].epatcal ? te0:p_info[i].epatcal;
                te1   =te1>p_info[i].epatcal ? te1:p_info[i].epatcal;

                maxwp4 = maxwp4 >p_info[i].wp4 ? maxwp4:p_info[i].wp4;
                minwp4 = minwp4 <p_info[i].wp4 ? minwp4:p_info[i].wp4;
		maxid= maxid > abs(p_info[i].ID) ? maxid:abs(p_info[i].ID);

        }
        cout<<"maxwp4 "<<maxwp4<<endl;
	cout<<"minwp4 "<<minwp4<<endl;
	cout<<"max ID "<<maxid<<endl;
	

	double 	t=0.005*(tpphi1-tpphi0);
	tpphi0	-=	t;
	tpphi1	+=	t;
	
	t=0.005*(tpangle1-tpangle0);
	tpangle0-=t;
	if(tpangle0<0)
		tpangle0=0;
	tpangle1+=t;

	t=0.00001*(te1-te0)+0.001;
	
	te0-=t;
	if(te0<0)
		te0=0;
	te1+=t;

	int	i=n_particle-1;
	if(!ip0)
		pphi0=tpphi0;
	if(!ip1)
		pphi1=tpphi1;
	if(!ia0)
		pangle0=tpangle0;
	if(!ia1)
		pangle1=tpangle1;
	if(!ie0)
		e0=te0;
	if(!ie1)
		e1=te1;	


	cout<<"finish read particle information:\t"<<ifilename<<endl;
	cout<<"real space boundary:"<<endl;
	cout<<"Max pphi: "<<tpphi1<<"\t Min pphi:\t"<<tpphi0<<endl;
	cout<<"Max pangle: "<<tpangle1<<"\t Min pangle:\t"<<tpangle0<<endl;
	cout<<"Max Energy: "<<te1<<"\t Min Energy:\t"<<te0<<endl;
	
	cout<<"space boundary for distribution function"<<endl;
	cout<<"Max pphi: "<<pphi1<<"\t Min pphi:\t"<<pphi0<<endl;
        cout<<"Max pangle: "<<pangle1<<"\t Min pangle:\t"<<pangle0<<endl;
        cout<<"Max Energy: "<<e1<<"\t Min Energy:\t"<<e0<<endl;
	
	assert(n_pphi>=0&&n_pangle>=0&&n_e>=0);	
	Alloc3D(distfpass,n_pphi+1,n_pangle+1,n_e+1);
	Alloc3D(distf0pass,n_pphi+1,n_pangle+1,n_e+1);
	Alloc3D(distftrap,n_pphi+1,n_pangle+1,n_e+1);
	Alloc3D(distf0trap,n_pphi+1,n_pangle+1,n_e+1);

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
	particledata  tp;
//	double tv=dpphi*dpangle*de;
	double tv=1.0;

	double tw=0;
	ofstream foutp("resonce_p.dat");
	for(i=0;i<n_particle;i++)
	{	
		
		tp=p_info[i];
		nx=(int)((tp.psibar2-pphi0)/dpphi);
		x0=tp.psibar2-((double)nx*dpphi+pphi0);
		x1=dpphi-x0;
		ny=(int)((tp.mu-pangle0)/dpangle);
		y0=tp.mu-((double)ny*dpangle+pangle0);
		y1=dpangle-y0;
		nz=(double)((tp.epatcal-e0)/de);
		z0=tp.epatcal-((double)nz*de+e0);
		z1=de-z0;
		int ox,oy,oz;
		ox=oy=oz=1;
		if(n_pphi==0)
		{
			ox=0;
			x0=1;
			x1=1;
		}
		if(nx>=n_pphi)
			continue;
		if(n_e==0)
		{
			oz=0;
			z0=1;
			z1=1;
		}
		if(ny>=n_pangle)
			continue;
		if(n_pangle==0)
		{
			oy=0;
			y0=1;
			y1=1;
		}
		if(nz>=n_e)
			continue;
		if(nx>(n_pphi)||ny>n_pangle||nz>n_e)
			continue;
		if(nx<0||ny<0||nz<0)
			continue;
			
//		if(p_info[i].ID>0)
		{
		distfpass[nx][ny][nz]	+=(x1*y1*z1/tv/tv)*tp.wp4;
		distfpass[nx][ny][nz+oz]	+=(x1*y1*z0/tv/tv)*tp.wp4;
		distfpass[nx][ny+oy][nz]	+=(x1*y0*z1/tv/tv)*tp.wp4;
		distfpass[nx][ny+oy][nz+oz]	+=(x1*y0*z0/tv/tv)*tp.wp4;
		distfpass[nx+ox][ny][nz]	+=(x0*y1*z1/tv/tv)*tp.wp4;
		distfpass[nx+ox][ny][nz+oz]	+=(x0*y1*z0/tv/tv)*tp.wp4;
		distfpass[nx+ox][ny+oy][nz]	+=(x0*y0*z1/tv/tv)*tp.wp4;
		distfpass[nx+ox][ny+oy][nz+oz]	+=(x0*y0*z0/tv/tv)*tp.wp4;

		distf0pass[nx][ny][nz]       +=(x1*y1*z1/tv/tv);
                distf0pass[nx][ny][nz+oz]    +=(x1*y1*z0/tv/tv);
                distf0pass[nx][ny+oy][nz]    +=(x1*y0*z1/tv/tv);
                distf0pass[nx][ny+oy][nz+oz] +=(x1*y0*z0/tv/tv);
                distf0pass[nx+ox][ny][nz]    +=(x0*y1*z1/tv/tv);
                distf0pass[nx+ox][ny][nz+oz] +=(x0*y1*z0/tv/tv);
                distf0pass[nx+ox][ny+oy][nz] +=(x0*y0*z1/tv/tv);
                distf0pass[nx+ox][ny+oy][nz+oz]      +=(x0*y0*z0/tv/tv);

		}
/*		else 
		{
		distftrap[nx][ny][nz]       +=(x1*y1*z1/tv/tv)*tp.wp4;
                distftrap[nx][ny][nz+oz]    +=(x1*y1*z0/tv/tv)*tp.wp4;
                distftrap[nx][ny+oy][nz]    +=(x1*y0*z1/tv/tv)*tp.wp4;
                distftrap[nx][ny+oy][nz+oz] +=(x1*y0*z0/tv/tv)*tp.wp4;
                distftrap[nx+ox][ny][nz]    +=(x0*y1*z1/tv/tv)*tp.wp4;
                distftrap[nx+ox][ny][nz+oz] +=(x0*y1*z0/tv/tv)*tp.wp4;
                distftrap[nx+ox][ny+oy][nz] +=(x0*y0*z1/tv/tv)*tp.wp4;
                distftrap[nx+ox][ny+oy][nz+oz]      +=(x0*y0*z0/tv/tv)*tp.wp4;

                distf0trap[nx][ny][nz]       +=(x1*y1*z1/tv/tv);
                distf0trap[nx][ny][nz+oz]    +=(x1*y1*z0/tv/tv);
                distf0trap[nx][ny+oy][nz]    +=(x1*y0*z1/tv/tv);
                distf0trap[nx][ny+oy][nz+oz] +=(x1*y0*z0/tv/tv);
                distf0trap[nx+ox][ny][nz]    +=(x0*y1*z1/tv/tv);
                distf0trap[nx+ox][ny][nz+oz] +=(x0*y1*z0/tv/tv);
                distf0trap[nx+ox][ny+oy][nz] +=(x0*y0*z1/tv/tv);
                distf0trap[nx+ox][ny+oy][nz+oz]      +=(x0*y0*z0/tv/tv);
		}
*/
		
	}


	for(i=0;i<(n_pphi+1)*(n_pangle+1)*(n_e+1);i++)
	{
		if(distf0pass[0][0][i]>0.0)
			distfpass[0][0][i]=distfpass[0][0][i]/distf0pass[0][0][i];
		if(distf0trap[0][0][i]>0.0)
			distftrap[0][0][i]=distftrap[0][0][i]/distf0trap[0][0][i];
	}

	
	int ncid;
	int err = nc_create(ofilename,NC_CLOBBER,&ncid);	
		

        int pphi_dimid;
	nc_def_dim(ncid,"npphi",n_pphi+1,&pphi_dimid);
       
        int pangle_dimid;
        nc_def_dim(ncid,"npangle",n_pangle+1,&pangle_dimid);

        int e_dimid;
        nc_def_dim(ncid,"nenergy",n_e+1,&e_dimid);
          
	int dmins[3];
	dmins[0]=pphi_dimid;
	dmins[1]=pangle_dimid;
	dmins[2]=e_dimid;
	int pphi_id,e_id,pangle_id;

	int f0_id,f00_id;
	nc_def_var(ncid,"fpass",NC_FLOAT,3,dmins,&f0_id);
	nc_def_var(ncid,"ftrap",NC_FLOAT,3,dmins,&f00_id);
err=    nc_def_var(ncid,"pangle",NC_FLOAT,1,  &pangle_dimid,&pangle_id);
        nc_def_var(ncid,"pphi",     NC_FLOAT,1,&pphi_dimid,&pphi_id);
        nc_def_var(ncid,"e",     NC_FLOAT,1,&e_dimid,&e_id);

	
	nc_enddef(ncid);
	

	nc_put_var_float(ncid,pangle_id,pangle_coor);
err=	nc_put_var_float(ncid,pphi_id,  pphi_coor);

        printf("pphi[0] %f \n",pphi_coor[0]);
	nc_put_var_float(ncid,e_id,     e_coor);

	nc_put_var_float(ncid,f0_id,distfpass[0][0]);
	nc_put_var_float(ncid,f00_id,distftrap[0][0]);
       
        nc_close(ncid);
 
	cout<<"write date to:\t"<<ofilename<<endl;
	cout<<"finish"<<endl;		
	
	delete[] e_coor;
	delete[] pphi_coor;
	delete[] pangle_coor;

	delete[] p_info;
	Free3D(distfpass);
	Free3D(distf0pass);
	Free3D(distftrap);
	Free3D(distf0trap);
	return 0;
}


