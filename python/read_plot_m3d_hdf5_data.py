#!/usr/bin/env python
#FengWang: read m3d output HDF5 files 
import h5py
import numpy
import matplotlib.tri as triangle
import matplotlib.figure

def print_hdf5_base_info(filename):
	f = h5py.File(filename,'r')
	print('reading M3D output file: \t%s'%filename)
	origin_path=f.attrs['working directory']
	print('path: \t%s'%origin_path)

	plane_No=f['planes/values'][0]		
	print('planes:\t%s'%plane_No)
	nsteps=f.attrs['nsteps']	
	timeframe=f.attrs['time']
	print('steps: \t%d'%nsteps)
	str_stamp='time stamp: \t'
	for tstamp in timeframe:
		str_stamp+=format('%f, \t'%tstamp)
	        	
	print(str_stamp)
	nnodes=f.attrs['nnodes']
	print('nodes: \t%d'%nnodes)	
	cell_set=f['cell_set[0]/node_connect_list']		
	print('cells: \t%d'%len(cell_set))	

	ndatas=f.attrs['nnode_data']
	print('datas: \t%d'%ndatas)
	datagroup=f['time_node_data[0]'].items()
	print('data list: ')
	for data_id in range(0,len(datagroup)):
		datastr=format('time_node_data[0]/node_data[%d]'%data_id)
		labels=f[datastr].attrs['labels']
		print('ID \t%d, \tlabels: \t%s '%(data_id,labels))
	f.close()	

def readhdf5_3d_times(filename):
	f = h5py.File(filename,'r')

	plane_No=f['planes/values'][0]		
	
	nsteps=f.attrs['nsteps']	
	timeframes=f.attrs['time']
	nnodes=f.attrs['nnodes']
	ndatas=f.attrs['nnode_data']
	XYZ=f['time_coordinates[0]/coordinates'].values()[0]	
	X=XYZ[:,0]
	Y=XYZ[:,1]
	Z=XYZ[:,2]
	R=(X**2+Y**2)**0.5
	cell_set=f['cell_set[0]/node_connect_list']		
	trianlist=cell_set[:,0:3]
	datagroup=f['time_node_data[0]'].items()
	data_array= []	
	str_array= []
	for time_id in range(0,nsteps):
		for data_id in range(0,len(datagroup)):
			datastr=format('time_node_data[%d]/node_data[%d]'%(time_id,data_id))
			labels=f[datastr].attrs['labels']
			str_array.append(labels)
	
			datastr=format('time_node_data[%d]/node_data[%d]/values'%(time_id,data_id))
			datas=f[datastr]
			datas=datas[:,0]	
			data_array.append(datas)	
	return (plane_No,nsteps,timeframes,nnodes,ndatas,R,Z,trianlist,str_array,data_array)
	
def readhdf5_3d_time_data(filename,time_id,data_id):
	(plane_No,nsteps,timeframes,nnodes,ndatas,R,Z,trianlist,str_array,data_array)=readhdf5_3d_times(filename)
	timeframe=timeframes[time_id]	
	i_begin=time_id*ndatas
	i_end=i_begin+ndatas
	str_array=str_array[i_begin:i_end]
	data_array=data_array[i_begin:i_end]
	data=data_array[data_id]
	data_name=str_array[data_id]
	return (plane_No,timeframe,nnodes,ndatas,R,Z,trianlist,data_name,data)
		
def readhdf5_2d_time_data(filename,time_id,data_id,plane_id):
	(plane_No,timeframe,nnode,ndatas,R,Z,trianlist,data_name,data_3D)=readhdf5_3d_time_data(filename,time_id,data_id)	
	nnode_1_plane=nnode/plane_No
	i_begin=plane_id*nnode_1_plane
	i_end=i_begin+nnode_1_plane
	R=R[i_begin:i_end]
	Z=Z[i_begin:i_end]
	data_2D=data_3D[i_begin:i_end]
	trianlist=trianlist[0:len(trianlist)/plane_No]
	return (timeframe,nnode_1_plane,R,Z,trianlist,data_name,data_2D)
	
def readhdf5_f_time_data(filename,time_id,data_id,plane_id):
	(plane_No,timeframe,nnode,ndatas,R,Z,trianlist,data_name,data_3D)=readhdf5_3d_time_data(filename,time_id,data_id)	
	nnode_1_plane=nnode/plane_No
	i_begin=plane_id*nnode_1_plane
	i_end=i_begin+nnode_1_plane
	R=R[i_begin:i_end]
	Z=Z[i_begin:i_end]
	data_2D=data_3D.reshape((plane_No,nnode_1_plane))
	
	fdata_2D=numpy.fft.fft(numpy.transpose(data_2D))	
	data_2D=fdata_2D.real[:,1]
	trianlist=trianlist[0:len(trianlist)/plane_No]
	return (timeframe,nnode_1_plane,R,Z,trianlist,data_name,data_2D)
	

def plot_m3d_hdf5_2D(filename,time_id,data_id,plane_id):
	import matplotlib.pyplot as plt
	import matplotlib.tri as triangle
	from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
		grid, savefig, show
	(timeframe,nnode_1_plane,R,Z,trianlist,data_name,data_2D)=readhdf5_2d_time_data(filename,time_id,data_id,plane_id)
	#plt.gca().set_aspect('equal')
	tri=triangle.Triangulation(R,Z)
	fig = figure(figsize=(6,8))
		
	plt.tripcolor(tri,data_2D)
	plt.xlabel('R')
	plt.ylabel('Z');
	title_str=format('%s, t=%f'%(data_name,timeframe))
	plt.title(title_str)
	plt.colorbar()
	plt.show()
		
def plot_f_m3d_hdf5_2D(filename,time_id,data_id,plane_id):
	import matplotlib.pyplot as plt
	import matplotlib.tri as triangle
	from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
		grid, savefig, show
	(timeframe,nnode_1_plane,R,Z,trianlist,data_name,data_2D)=readhdf5_f_time_data(filename,time_id,data_id,plane_id)
	#plt.gca().set_aspect('equal')
	tri=triangle.Triangulation(R,Z)
	fig = figure(figsize=(6,8))
	plt.tripcolor(tri,data_2D)
	plt.xlabel('R')
	plt.ylabel('Z');
	title_str=format('%s, t=%f'%(data_name,timeframe))
	plt.title(title_str)
	plt.colorbar()
	plt.show()

		


