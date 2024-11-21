def return_txt_05(UNIFORM_SST,IDEALIZED,FUNC,DIR,T0,deltaT,atpos,Lx,
                  XSST_CRIT,dico_origin,path_sst,
                  IXOR,IYOR,NI_s,NJ_s,
                  nhalo,dic_nodes,dic_cpu,dic_timelimit):
    center2_replacesst = """
path = '../04_prep_real_son/'
name_in = 'INIT_SON'
name_out = 'INIT_SON_SST'
"""
    part3_replacesstbis = """
else:
    CHOICE = "SON" """

    part1_replacesst = """import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from function_replaceSST import *
import os

## INPUT ===================================
# MODEL ----------------------------------
nhalo="""+str(nhalo)+"""
# TYPE OF SST ----------------------------
IDEALIZED = """+str(IDEALIZED)+""" # 1 front, NS or EW, if False: Agulhas front
FUNC = '"""+FUNC+"""' 	# if IDEALIZED: tanh or linear
DIR = '"""+DIR+"""' # direction of change of SST
T0="""+str(T0)+""" 			# K, if IDEALIZED
deltaT="""+str(deltaT)+"""		# K, if IDEALIZED
atpos="""+str(atpos)+""" 		# m, if IDEALIZED: position of the front
Lx="""+str(Lx)+""" 		# m, if IDEALIZED: width of the front
# Name of files --------------------------"""

    center1_replacesst = """
path = '../00_prep_ideal_dad/'
name_in = 'INIT_DAD'
name_out = 'INIT_DAD_SST'
"""
    if UNIFORM_SST:
        string_except = """raise Exception("You have chosen a uniform SST no you don't need to run this step! Aborting ...")"""
    else:
        string_except = ''

    part2_replacesst = """SAVING = True

#===========================================
print('Currently replacing SST for '+path+name_in)

ds = xr.open_dataset(path+name_in+'.nc') # here xarray throw warning on duplicate 'ni' dimension but its ok
X = ds.ni.values
Y = ds.nj.values
res = X[1] - X[0]
print('resolution is (m):')
print(res)

if IDEALIZED:
    # First let us have a look at the SST prescribed
    print('ni is :')
    print(X)
    print('nj is :')
    print(Y)
    if DIR=='X':
        dim=X
    else:
        dim=Y 
    SST1D = One_Step_SST(dim,atpos,Lx,T0,deltaT,FUNC)

    fig, ax = plt.subplots(1,1,figsize = (10,5),constrained_layout=True,dpi=100) 
    ax.plot(dim/1000,SST1D-273,color='k')

    ax.set_xlabel(DIR+' (km)')
    ax.set_ylabel('SST (°C)')
    plt.show()

    ## Then if its ok, we can build the SST field 
    SST = buildSST(DIR,SST1D,len(X),len(Y))
    print('SST shape is :')
    print(SST.shape)
    fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=100) 
    s = ax.pcolormesh(X/1000,Y/1000,SST-273,cmap='rainbow')
    plt.colorbar(s,ax=ax,location='bottom')
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_title('SST (°C)')
    ax.set_aspect('equal')
    plt.show()

    ## Replacing the SST from the initial file 
    print('SST old shape :',ds.SST.shape)
    print('SST new shape :',SST.shape,' (should be = to old shape)')
    SSTold = ds.SST[:,:]
    attrs = ds.SST.attrs 
    ds['SST'] = (['nj', 'ni'], SST)
    ds.SST.attrs = attrs
    fig, ax = plt.subplots(2,1,figsize = (7,7),constrained_layout=True,dpi=100) 
    vmax=np.amax(SST)
    vmin=np.amin(SST)
    ax[0].pcolormesh(X/1000,Y/1000,SSTold,cmap='rainbow',vmin=vmin,vmax=vmax)
    s = ax[1].pcolormesh(X/1000,Y/1000,SST,cmap='rainbow',vmin=vmin,vmax=vmax)
    ax[1].set_xlabel('X (km)')
    ax[1].set_ylabel('Y (km)')
    ax[0].set_title(r'$T_0$='+str(T0)+r'K $\Delta T$='+str(deltaT)+r' atpos '+DIR+'='+str(atpos/1000)+'km Lx='+str(Lx/1000)+'km')
    plt.colorbar(s,ax=ax,location='bottom')
    fig.savefig(path+'check_SST_'+FUNC+'.png')
    plt.show()
    if SAVING:
        print('saving ...')
        ds.to_netcdf(path+name_out+'.nc')
        os.system("cp "+path+name_in+".des "+path+name_out+".des")
        print('done !')
    """

    part4_replacesst = """
    if CHOICE=='SON':
        os.system("echo 'No need to replace SST for Son domain, it is done by MNH, aborting ...'")
        raise Exception('No need to replace SST for Son domain, it is done by MNH')
    
    XSST_CRIT = """+str(XSST_CRIT)+""" # K
    
    dico_origin = """+str(dico_origin)+"""
    ds_sst = xr.open_dataset('"""+path_sst+"""').sel(lat=slice(-40,-30),lon=slice(20,30)).isel(time=0)
    
    ORIGIN = dico_origin[CHOICE]   # lower left corner coordinates of LES (°E,°N)
    figsize1 = (10,8)
    # resolution is 0.02 degree.
    DegLat = 111.2  #  =1° lat in km
    DegLon = 91.6   #  =1° lon in km at ~ -35°N
    ds_sst['lonKM'] = (ds_sst['lon']-ORIGIN[0])*DegLon
    ds_sst['latKM'] = (ds_sst['lat']-ORIGIN[1])*DegLat
    ds_sst = ds_sst.set_coords(['lonKM','latKM'])
    ds_sst['sst_on_km_grid'] = ds_sst['analysed_sst'].swap_dims({'lon':'lonKM','lat':'latKM'})
    
    ds_sst['sst_LES_grid'] = ds_sst['sst_on_km_grid'].interp(lonKM=ds.ni/1000,latKM=ds.nj/1000)
    Lx = X[-1]/1000 # km
    Ly = Y[-1]/1000 # km
    dx_d = X[1] - X[0]
    eps_x = """+str(IXOR)+"""
    eps_y = """+str(IYOR)+"""
    NI_s = """+str(NI_s)+"""
    NJ_s = """+str(NJ_s)+"""
    
    Lx_s = NI_s*dx_d
    Ly_s = NJ_s*dx_d
    
    print('lower left corner of simulation is (°E,°N) :',ORIGIN)
    print('width is (km,degLon) :',Lx,Lx/DegLon)
    print('height is (km,degLat) :',Ly,Ly/DegLat)

    vmin,vmax = 295,298

    fig, ax = plt.subplots(1,1,figsize = figsize1,constrained_layout=True,dpi=100) 
    s = ax.pcolormesh(ds_sst.lon,ds_sst.lat,ds_sst['analysed_sst'],cmap='jet',vmin=vmin,vmax=vmax)
    ax.scatter(ORIGIN[0],ORIGIN[1],marker='x',color='r')
    ax.scatter(dico_origin['SON'][0],dico_origin['SON'][1],marker='x',color='k')
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50)
    ax.set_aspect(DegLat/DegLon)
    ax.add_patch(mpl.patches.Rectangle(ORIGIN, Lx/DegLon, Ly/DegLat,
                edgecolor='red',
                fill=False,
                lw=2))
    ax.add_patch(mpl.patches.Rectangle((ORIGIN[0]+eps_x*dx_d/1000/DegLon,ORIGIN[1]+eps_y*dx_d/1000/DegLat), Lx_s/1000/DegLon, Ly_s/1000/DegLat,
                edgecolor='k',
                fill=False,
                lw=2))
    ax.set_title('SST (K) on lat/lon grid, with origin at '+str(ORIGIN))
    ax.set_xlabel('Lon (°E)')
    ax.set_ylabel('Lat (°N)')
    ax.set_xlim([23,26])
    ax.set_ylim([-37,-35])
    fig.savefig('SST_on_LAT_LON_grid_wide_view.png')
    
    fig, ax = plt.subplots(1,1,figsize = figsize1,constrained_layout=True,dpi=100) 
    s = ax.pcolormesh(ds_sst.lonKM,ds_sst.latKM,ds_sst['sst_on_km_grid'],cmap='jet',vmin=vmin,vmax=vmax)
    ax.scatter(0,0,marker='x',color='r')
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50)
    ax.set_aspect(1)
    ax.add_patch(mpl.patches.Rectangle((0, 0), Lx, Ly,
                edgecolor='r',
                fill=False,
                lw=2))
    ax.add_patch(mpl.patches.Rectangle((eps_x*dx_d/1000, eps_y*dx_d/1000), Lx_s/1000, Ly_s/1000,
                edgecolor='k',
                fill=False,
                lw=2))
    ax.set_title('SST (K) on km grid, with origin at '+str(ORIGIN))
    ax.set_xlabel('km away from '+str(ORIGIN[0])+'°E')
    ax.set_ylabel('km away from '+str(ORIGIN[1])+'°N')
    ax.set_xlim([(23-ORIGIN[0])*DegLon,(26-ORIGIN[0])*DegLon])
    ax.set_ylim([(-37-ORIGIN[1])*DegLat,(-35-ORIGIN[1])*DegLat])
    fig.savefig('SST_on_km_grid_wide_view.png')
    
    fig, ax = plt.subplots(1,1,figsize = figsize1,constrained_layout=True,dpi=100) 
    s = ax.pcolormesh(ds_sst.lonKM,ds_sst.latKM,ds_sst['sst_on_km_grid'],cmap='jet',vmin=vmin,vmax=vmax)
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50)
    ax.contour(ds_sst.lonKM,ds_sst.latKM,ds_sst['sst_on_km_grid'],levels=[XSST_CRIT],colors='k')
    ax.set_aspect(1)
    ax.add_patch(mpl.patches.Rectangle((eps_x*dx_d/1000, eps_y*dx_d/1000), Lx_s/1000, Ly_s/1000,
                edgecolor='k',
                fill=False,
                lw=2))
    ax.set_xlim([ds.ni[0]/1000,ds.ni[-1]/1000,])
    ax.set_ylim([ds.nj[0]/1000,ds.nj[-1]/1000,])
    ax.set_title('SST (K) on km grid, with origin at '+str(ORIGIN))
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.savefig('SST_on_LES_nointerp_grid.png')
    
    fig, ax = plt.subplots(1,1,figsize = figsize1,constrained_layout=True,dpi=100) 
    s = ax.pcolormesh(ds.ni/1000,ds.nj/1000,ds_sst['sst_LES_grid'],cmap='jet',vmin=vmin,vmax=vmax)
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50)
    ax.contour(ds.ni/1000,ds.nj/1000,ds_sst['sst_LES_grid'],levels=[XSST_CRIT],colors='k')
    ax.set_aspect(1)
    ax.add_patch(mpl.patches.Rectangle((eps_x*dx_d/1000, eps_y*dx_d/1000), Lx_s/1000, Ly_s/1000,
                edgecolor='k',
                fill=False,
                lw=2))
    ax.set_title('SST (K) on LES grid, with origin at '+str(ORIGIN))
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.savefig('SST_on_LES_interp_grid.png')
    
    if SAVING:
        print('saving ...')
        ds['SST'].data = ds_sst['sst_LES_grid'].data
        ds.to_netcdf(path+name_out+'.nc')
        os.system("cp "+path+name_in+".des "+path+name_out+".des")
        print('done !')
plt.show()
    """

    TXT_01_run = """#!/bin/bash
    #SBATCH -J MNH01
    #SBATCH -N """+str(dic_nodes['01'])+"""           # nodes number (=NBP)   
    #SBATCH -n """+str(dic_cpu['01'])+"""       # CPUs number (on all nodes) (=NBP*TPN) 
    #SBATCH -o MNH01.eo%j   #
    #SBATCH -e MNH01.eo%j   #
    #SBATCH -t """+ dic_timelimit['01']+"""   # time limit
    # Echo des commandes
    set -x
    python replaceSST.py | tee 'replaceSST.log'
    """

    TXT_05_replaceSST = part1_replacesst + string_except + center2_replacesst + part2_replacesst + part3_replacesstbis + part4_replacesst

    TXT_05_run = TXT_01_run

    return TXT_05_replaceSST,TXT_05_run