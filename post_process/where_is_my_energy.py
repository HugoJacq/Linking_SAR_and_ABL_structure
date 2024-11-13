import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

"""
In the context of gridnesting, there might be a possibility to reduce the distance needed for the turbulence to develop at the son domain scale. 

This discussion has been triggered by R.Manceau during the LEFE 'Atelier sur la représentation des fines échelles océaniques dans les Simulations Numériques' 9-11 october 2024 at IFREMER.
https://atelier-fines-ech.sciencesconf.org/

The added peturbation term in the momentum equation is described in details here : https://inria.hal.science/hal-04287062/document#&sid=oa&domaine=insu&doi=10.1063/5.0174381

Questions that remain:
- How this can be applied to stratified rotating turbulence ?
- R. Manceau told me that adding a perturbation term in the energy conservation equation was not needed in his tests.
- How to implement this in the grid nesing context in MNH ?

The present script will:
- look at how energy is distributed between resolved/subgrid scale
- look at kinetic energy (both mean flow and turbulent) vs intern/potential energy repartion

Notes: test are done with homogeneous SST and a rather small domain. Next step would be to do the same analysis on a bigger domain.
"""

basepath = '/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/'

path_dad = basepath+'test_HR_nesting_NoSSTgradient/06_mesonh_2models/FICHIERS_OUT/TEST2.1.001.OUT.007.nc'
path_son = [basepath+'test_HR_nesting_NoSSTgradient/06_mesonh_2models/FICHIERS_OUT/TEST2.2.001.OUT.005.nc',
			basepath+'test_HR_nesting_NoSSTgradient/06_mesonh_2models/FICHIERS_OUT/TEST2.2.001.OUT.006.nc',
			basepath+'test_HR_nesting_NoSSTgradient/06_mesonh_2models/FICHIERS_OUT/TEST2.2.001.OUT.007.nc']

ds = {'dad':xr.open_dataset(path_dad),
		'son':xr.open_mfdataset(path_son)}
nhalo = 1
dimMean = {'dad':['time','ni','nj'],
			'son':['time','nj']}
def MEAN(array,dim):
	return array.mean(dim)
	
print('* Computing turbulent kinetic energy ...')	
X,X_u,Y,Y_v,Z,Z_w = {},{},{},{},{},{}
U,V,W = {},{},{}
ET,Er,TKE = {},{},{}
meanU,meanV,meanW = {},{},{}
meanET, meanEr, meanTKE = {},{},{}
u_f,v_f,w_f = {},{},{}
for name in ds.keys():
	X[name],X_u[name] = ds[name].ni[nhalo:-nhalo],    ds[name].ni_u[nhalo:-nhalo]
	Y[name],Y_v[name] = ds[name].nj[nhalo:-nhalo],    ds[name].nj_v[nhalo:-nhalo]
	Z[name],Z_w[name] = ds[name].level[nhalo:-nhalo], ds[name].level_w[nhalo:-nhalo]
	
	U[name] = ds[name].UT.interp({'ni_u':ds[name].ni})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
	U[name] = U[name].rename(new_name_or_name_dict={'nj_u':'nj'})
	V[name] = ds[name].VT.interp({'nj_v':ds[name].nj})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
	V[name] = V[name].rename(new_name_or_name_dict={'ni_v':'ni'})
	W[name] = ds[name].WT.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
	TKE[name] = ds[name].TKET[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]

	meanU[name] = MEAN(U[name],dimMean[name])
	meanV[name] = MEAN(V[name],dimMean[name])
	meanW[name] = MEAN(W[name],dimMean[name])
	meanTKE[name] = MEAN(TKE[name],dimMean[name])

	u_f[name] = U[name] - meanU[name]
	v_f[name] = V[name] - meanV[name]
	w_f[name] = W[name] - meanW[name]
	Er[name] = 0.5*(u_f[name]**2+v_f[name]**2+w_f[name]**2)
	ET[name] = Er[name] + TKE[name]

	meanEr[name] = MEAN(Er[name],dimMean[name])
	meanET[name] = MEAN(ET[name],dimMean[name])
    

if True:
	print('* Plotting profiles of turbulent kinetic energy')
	liste_kx = [0,20,40,60,80,-1]
	figsize = (20,5)
	fig, ax = plt.subplots(1,len(liste_kx)+1,figsize = figsize,constrained_layout=True,dpi=100)
	ax[0].plot(meanET['dad'],Z['dad'],c='k',label='total')
	ax[0].plot(meanEr['dad'],Z['dad'],c='r',label='res')
	ax[0].plot(meanTKE['dad'],Z['dad'],c='b',label='sgs')
	ax[0].set_xlabel('m2/s2')
	ax[0].set_ylabel('z m')
	ax[0].legend()
	ax[0].set_title('DAD ref')

	for i,kx in enumerate(liste_kx):
		ax[1+i].plot(meanET['son'][:,kx],Z['son'],c='k',label='total')
		ax[1+i].plot(meanEr['son'][:,kx],Z['son'],c='r',label='res')
		ax[1+i].plot(meanTKE['son'][:,kx],Z['son'],c='b',label='sgs')
		ax[1+i].set_xlabel('X = '+str(X['son'][kx].values)+' m')
		ax[1+i].tick_params(axis='both',labelleft=False)
		
	for axe in ax:
		axe.set_ylim([0,1000])	
		axe.set_xlim([0,0.5])
		axe.grid()
	ax[-1].set_title('SON ref')
	fig.suptitle('Turbulent kinetic energy (m2/s2)')
	fig.savefig('PNGs_verifTurb/Z_evolution_ETvsERvsTKE.png')

if True:
	print('* Plotting X evolution of turbulent kinetic energy')
	liste_z = [10,100,400,600]
	figsize = (8,8)
	fig, ax = plt.subplots(len(liste_z),1,figsize = figsize,constrained_layout=True,dpi=100)

	for i,kz in enumerate(liste_z[::-1]):
		# dad
		ax[i].scatter(X['son'][0]/1000 - 0.2,meanET['dad'].sel(level=kz,method='nearest'),marker='x',c='k')
		ax[i].scatter(X['son'][0]/1000 - 0.2,meanEr['dad'].sel(level=kz,method='nearest'),marker='x',c='r')
		ax[i].scatter(X['son'][0]/1000 - 0.2,meanTKE['dad'].sel(level=kz,method='nearest'),marker='x',c='b')
		# son
		ax[i].plot(X['son']/1000,meanET['son'].sel(level=kz,method='nearest'),c='k',label='total')
		ax[i].plot(X['son']/1000,meanEr['son'].sel(level=kz,method='nearest'),c='r',label='res')
		ax[i].plot(X['son']/1000,meanTKE['son'].sel(level=kz,method='nearest'),c='b',label='sgs')
		ax[i].set_title('z='+str(kz)+'m',loc='right')
		ax[i].grid()
		ax[i].set_yscale('log')
		#ax[i].set_ylim([0.01,0.6])
		
		if i<len(liste_z)-1:
			ax[i].tick_params(axis='both',labelbottom=False)
		if i==len(liste_z)-1:
			ax[i].set_xlabel('X (km)')
	fig.suptitle('Turbulent kinetic energy (m2/s2)')
	fig.savefig('PNGs_verifTurb/X_evolution_ETvsERvsTKE.png')
plt.show()
    
    
    
    
    
    
    
    



