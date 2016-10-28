#! /usr/bin/env python
#fengwang: python script for get m3d config CPU number
# 2014/09/12
from sys import argv

def get_config_cpu_no(filename):

	configdata=open(filename);
	
	A=configdata.readline();
	B=configdata.readline();
	C=configdata.readline();
	D=configdata.readline();
	E=configdata.readline();
	F=configdata.readline();
	G=configdata.readline();
	
	cpu_toroidal=int(B.split()[0]);
	cpu_r=int(D.split()[0]);
	cpu_theta=int(F.split()[0]);

	if cpu_r==1:
		total_cpu=cpu_toroidal*cpu_r*cpu_theta;
	else:
		total_cpu=cpu_toroidal*(cpu_r-1)*cpu_theta+cpu_toroidal;

	configdata.close;

	return total_cpu;




def get_pbs_cpu_no(filename):
	
	pbsfile=open(filename);
	for each_line in pbsfile:
		if each_line.find('mppwidth')!=-1:
			(Astr,Bstr)=each_line.split('=');
			pbs_cpu=int(Bstr);

		if each_line.find('aprun')!=-1:
			aprun_cpu_no=int(each_line.split()[2]);
			if(pbs_cpu != aprun_cpu_no):
				print('aprun cpu number is:', aprun_cpu_no);
				print('pbs cpu number is:',pbs_cpu);
				print('pbs cpu number does not equal to aprun cpu number!!!!');
				return False;
	return aprun_cpu_no;



import f90nml
import os

def write_wxy(wxy,filename):
	e = os.path.isfile(filename)
	if e:
        	print(filename, "exist")
        	os.remove(filename)
        	print(filename,"deleted")
	f90nml.write(wxy,filename);
	print("write a new wxy:",filename)



