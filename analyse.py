import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import dask
import pathlib
from dask.distributed import Client,LocalCluster
import time
# customized modules
from mod_build_mean import build_mean_file,mass_interpolator,Reynolds_avg
from mod_build_CS import build_CS
from module_cst import *
from module_tools import *
from mod_turb_stats import *
from mod_first_look import *
from mod_spatial_stats import  *
from mod_CS_analysis import *
"""
This script is analysing mesoNH simulations to compare them to a SAR imagette

PreProcess
- define boxes
- build mean file (prognostic var and fluxes)
- build conditional sampling file

PostProcess 
- Compute L Obukhov for both LES and SAR
- Compute geometric statistics for both LES and SAR

run with : python analayse.py

"""
start = time.time()

# INPUT ZONE ============================================================================================
#
CASE = '1' # 1 is ERA5-like conditions
SAR_SIZE = 'minimal' # amount of vignette. 'all' or 'small' or 'small_only' or 'minimal'
# First look -----------------------------------------------------------
PLOT_10m_WIND = False # 10m wind and SAR roughness for first 3 boxes
WHERE_ARE_THE_BOXES = True # both SAR and LES
# Geometrical analysis -------------------------------------------------
PLOT_2D_COVARIANCE = False 
PLOT_2D_SKEWNESS = False
S2_ANALYSIS = False
# Turbulence convergence -----------------------------------------------
VERIFY_TURB_CONVERGENCE = False # plot spectrum at inflow
B_KPSD = True                   # plot k*PSD(k) ?
altZ_convergence = 200          # m 
liste_X_turb_convergence = [0,2,4,6,8,10] # km, distance from East border of son domain
A = 0.1           # y=log(x)**(-coeff)+log(A) for Kolmogorov Law in inertial subrange
Kcoeff = -5/3
# Coeherent structure analysis -----------------------------------------
PLOT_MEAN_PROGVAR = False
PLOT_MEAN_FLX = False
PLOT_TopView_with_CS = False
VAR_TOPVIEW = 'W'

# folder organisation --------------------------------------------------
workdir_path = '/home/jacqhugo/WORKDIR/MNH570/'
path_here = '/home/jacqhugo/scripts/simu_nest_aggv3/'
#path_SAR = path_here+'SAR_from_OVL/SAR_roughness_20151210t170827-20151210t170921.nc'
path_SAR = path_here+'SAR_from_IFREMER/S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_9AE1.nc'
path_SST_ODYSEA = path_here+'SST_from_Ifremer/20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc'
path_data_turb = path_here + 'DATA_TURB/'
path_save_turb_convergence = 'PNGs_turb_convergence/'
path_save_First_look = 'PNGs_First_look/'
path_save_geometric = 'PNGs_geometric/'
path_save_CS = 'PNGs_CS_analysis/'
# creating folders
os.system('mkdir -p '+path_data_turb)
os.system('mkdir -p '+path_save_turb_convergence)
os.system('mkdir -p '+path_save_First_look)
os.system('mkdir -p '+path_save_geometric)
os.system('mkdir -p '+path_save_CS)

sim_dir = {'0':{'backup_name':'', # test folder
                'out_name':[workdir_path+'simu_nest_Agulhas_0/06_mesonh_2models/FICHIERS_OUT/SON_600s_OUT_scale_offset_7.nc',
                            workdir_path+'simu_nest_Agulhas_0/06_mesonh_2models/FICHIERS_OUT/SON_600s_OUT_scale_offset_7bis.nc']},
           '1':{'backup_name':'', # ERA5 like 
                'out_name':[workdir_path+'simu_nest_Agulhas_1/06_mesonh_2models/FICHIERS_OUT/AGULA.2.001.OUT.006_comp.nc',
                            workdir_path+'simu_nest_Agulhas_1/06_mesonh_2models/FICHIERS_OUT/AGULA.2.001.OUT.007_comp.nc']},}
# the compressed files are produced with the command:
# ncks -O -7 -L 1 --ppc default=4#PABST=7#RVT=5#THT=6#PABST=7 in.nc out.nc

dsB_path = sim_dir[CASE]['backup_name']
dsO_path = sim_dir[CASE]['out_name']
#
# Global constant -------------------------------------------------------

# SST related
crit_value = 295.65 # K, = XSST_CRIT in NAM_CONDSAMP

# Boxe definition
#   SAR : According to prepare_obs.py
#   LES : not colocated because SST analysis is not colocated with SAR, 
#          so to be set to approximated position
#           Lx and Ly are the same for all boxes
d_boxes = {'SAR':{ 
                'Lx':20, # km, across track, sample
                'Ly':20, # km, along track, line
                'boxes':{
                    '1':{'O':(24.75,-35.5)}, # 'O' in °lon,°lat. [LEFT,BOTT] corner
					'2':{'O':(24.5,-35.8)}, # 
					'3':{'O':(24.75,-36)}, # 
                        }
                    },
            'LES':{ # all in meters
                'Lx':15000,
                'Ly':15000,
                'boxes':{
                    '1':{'O':(75000,85000)}, # top, highly convective
                    '2':{'O':(80000,52000)}, # middle, near neutral
                    '3':{'O':(50000,30000)}  # bott, averaged convective
                         }
                    }
            }
        
#
# Dask related ----------------------------------------------------------
DASHBOARD = False
ChunksZ = 17
NIchunks = 129
NJchunks = 201
chunksOUT = {'time':-1,
            'level':ChunksZ,
            'level_w':ChunksZ,
            'nj':NJchunks,
            'nj_u':NJchunks,
            'nj_v':NJchunks,
            'ni':NIchunks,
            'ni_u':NIchunks,
            'ni_v':NIchunks}	
# print( give_NETCDF4_chunking(temp_path) ) # = [1, 17, 201, 129]
#   if chunksOUT is different from what is saved in the netcdf4, 
#   a warning about performance is raised.
# END INPUT ZONE ========================================================================================


if __name__ == "__main__":  # This avoids infinite subprocess creation
    
    # DASK related ------------------------------------------------
   
    global client
    client = None

    if DASHBOARD:
        # sometimes dask cluster can cause problems "memoryview is too large"
        # (writing a big netcdf file for eg, hbudget_file)
        cluster = LocalCluster(threads_per_worker=32,n_workers=8)
        
        client = Client(cluster)
        print("Dashboard at :",client.dashboard_link)
    
    

    # BUILDING FILES ----------------------------------------------
    # opening OUTPUT files
    # dsB = 
    dsO = xr.open_mfdataset(dsO_path,chunks=chunksOUT) 
    dsO_i = mass_interpolator(dsO)
    dsO_i = dsO_i.isel(level=slice(nhalo,-nhalo),
                           ni=slice(nhalo,-nhalo),
                           nj=slice(nhalo,-nhalo))
    
    # ** MEAN FILE 
    #   -> mean profiles of prognostic variables
    #   -> mean profiles of resolved/subgrid fluxes
    path_origin = dsO_path
    path_out = path_data_turb+'mean_'+CASE+'.nc' 
    path_out_flx = path_data_turb+'mean_flux_'+CASE+'.nc' 
    build_mean_file(path_origin,path_out,path_out_flx,nhalo,dsO_i,d_boxes['LES'])
    dsmean = xr.open_dataset(path_out)
    dsflx = xr.open_dataset(path_out_flx)
    # ** Conditional sampling files
    #   -> 3D structures identified
    #   -> mean profiles decomposed.
    #   -> 3D volumes with structures identified
    path_out = path_data_turb+'Cond_sampling_'+CASE+'.nc' 
    gamma = 0.005
    mCS = 1.0
    build_CS(path_origin,path_out,nhalo,dsO_i,dsmean,d_boxes['LES'],gamma=gamma,mCS=mCS)
    dsCS = xr.open_dataset(path_out)


    # Opening SAR file
    dsSAR = xr.open_dataset(path_SAR).sel(pol='VV')
    # Opening SST file
    dsSST = xr.open_dataset(path_SST_ODYSEA)
    X,Y,Z = dsO_i.ni.values,dsO.nj.values,dsO.level.values

    # automatic boxe locator
    d_boxes['SAR'] = SAR_boxes_generator(dsSAR,d_boxes['SAR'],SAR_SIZE)

    # Building structure functions
    N = 2
    save_S_n_SAR(dsSAR,2,d_boxes['SAR'],path_data_turb)
    dsS2_SAR = xr.open_dataset(path_data_turb+'S_2_sig0.nc')
    save_S_n_LES(dsCS,'M',10,2,path_data_turb)
    dsS2_LES_M10 = xr.open_dataset(path_data_turb+'S_2_M10.nc')
    # N = 3
    # save_S_n_SAR(dsSAR,3,d_boxes['SAR'],path_data_turb)
    # dsS2_SAR = xr.open_dataset(path_data_turb+'S_3_sig0.nc')
    # save_S_n_LES(dsCS,'M',10,3,path_data_turb)
    # dsS2_LES_M10 = xr.open_dataset(path_data_turb+'S_3_M10.nc')

    ###----------###
    # POST-PROCESS #
    ###----------###

    # TURB_STATS
    if VERIFY_TURB_CONVERGENCE:
        print('* Verification of turbulence at inflow')
        plot_energy_spectrum(dsO,altZ_convergence,liste_X_turb_convergence,B_KPSD,A,Kcoeff,path_save_turb_convergence)
        
    # FIRST LOOK
    if PLOT_10m_WIND:
        print('* 10m wind')
        Wind_10m_in_boxes_vs_SAR(dsCS,dsSAR,indt=-1,d_boxes=d_boxes,path_save=path_save_First_look)
    if WHERE_ARE_THE_BOXES:
        print('* Plotting location of boxes on LES and SAR')
        Where_boxes(dsO_i,dsSAR,dsSST,indt=-1,d_boxes=d_boxes,SAR_SIZE=SAR_SIZE,path_save=path_save_First_look)

    # GEOMETRICS
    if PLOT_2D_COVARIANCE:
        print('* Plotting 2D covariance')
        N = 2
        print(dsS2_LES_M10)
        print(dsS2_SAR)
        Plot_S_n(dsS2_SAR,'sig0',N,path_save_geometric)
        Plot_S_n(dsS2_LES_M10,'M10',N,path_save_geometric)
    if PLOT_2D_SKEWNESS:
        print('* Plotting 2D skewness')
        N = 3
        print('TO BE DONE')

    if S2_ANALYSIS:
        print('* tests')
        # compute_integral_scale_at_tht is ok      
        # test S2_analysis ?
        path_out = 'DATA_TURB/'
        S2_analysis('SAR',dsCS,d_boxes,path_out)

    # comparer échelles caractéristiques SAR vs LES.
    # > SAR
    #       utiliser code PE Brilouet pour avoir une estimation de l'ellipse comme dans son papier.
    #       utiliser fonction de structure en polaire pour récupérer un R.
    #         
    # > LES
    #       récupérer distribution de la taille des updrafts en fonction de l'altitude( 2 méthodes, ellipse ou rayon eq.)
    # ajuster colorbar des fonctions de structure
    # paralléliser le calcul de fonction de structure

    # Régler problème de résolution du SAR
    #   check conversion deg vs km -> ok 
    #   récupérer dataset depuis copernicus puis ouverture avec le package S1 xarray -> non trop long
    #   récupérer données Alexis Mouche ou PE Brilouet ? -> j'ai, à voir les diff avec les données OVL

    # Régler problème de profiles de flux qui ne sont pas bien calculés.

    # altZ = 100
    # is_up1 = xr.where( dsCS.isel(nboxe=0).global_mask==1,1,0 )
    # is_ss1 = xr.where( dsCS.isel(nboxe=0).global_mask==2,1,0 )
    # is_up2 = xr.where( dsCS.isel(nboxe=0).global_mask==3,1,0 )
    # is_ss2 = xr.where( dsCS.isel(nboxe=0).global_mask==4,1,0 )
    # is_down = xr.where( dsCS.isel(nboxe=0).global_mask==5,1,0 )

    # B = is_up2.sel(level=altZ,method='nearest').isel(time=-1)

    # print(B.shape)
    # print(B.sum().values,
    #       400*400,
    #       B.sum().values  / (400*400),
    #     B.mean().values)

    # Lcolors = {'up1':'r',
    #         'ss1':'purple',
    #         'up2':'orange',
    #         'ss2':'pink',
    #         'down':'green'}
    
    # A = dsCS.isel(nboxe=0,time=-1).sel(level=altZ,method='nearest').global_mask

    # fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=200)
    # s = ax.pcolormesh(dsCS.X.isel(nboxe=0)/1000,dsCS.Y.isel(nboxe=0)/1000,
    #                A,cmap='jet')
    # plt.colorbar(s,ax=ax,label='objects')
    # ax.set_xlabel('X')
    # ax.set_ylabel('Y')
    # ax.set_title('objects at Z = '+str(altZ)+'m')
    # ax.set_aspect(1)

    # COHERENT STRUCTURE ANALYSIS
    if PLOT_MEAN_PROGVAR:
        print('* Plotting profiles of prognostic variables, with object decomposition')
        #Plot_mean_progvar_allboxes(dsCS,dsmean,path_save_CS,Reynolds_avg)
        
    if PLOT_MEAN_FLX:
        print('* Plotting profiles of fluxes, with object decomposition')
        #Plot_mean_flux_allboxes(dsCS,dsmean,dsflx,path_save_CS,Reynolds_avg)
    
    if PLOT_TopView_with_CS:
        print('* Plotting top view of '+VAR_TOPVIEW+' with coherent structures')
        Plot_top_view_var(dsCS,atZ=250,path_save=path_save_CS) # this is not finished

    # - Compute L Obukhov for both LES and SAR
    # LES : use L_Obukhov(tht0,u_star,surf_wtht_flx) from module_tool
    # SAR : from O'Driscoll paper (10.1029/2023GL104228)

    # - Compute geometric statistics for both LES and SAR
    # L and K statistics

    end = time.time()
    print('Total runtime of analyse.py :'+sec2hms(end - start))
    plt.show()
    # avoid memory leaks
    dsO.close()
    #dsB.close()