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
# =======================================================================================================
# INPUT ZONE                                                                                            |
# =======================================================================================================

# Selector of simulation (as in the namlist injector setup_big.py)
#   - 0 is less convective, weaker inversion than ERA5
#   - 1 is ERA5-like conditions
CASE = '1'         
 
# Amount of vignette. 'all' or 'small' or 'small_only' or 'minimal'
#   'all' is looking for a maximum of vignette on the full swath.
#   'small' is looking for a maximum of vignette but on the equivalent LES domain (~50 boxes) 
#       AND it also takes into account what asked the user (in 'd_boxes') 
#   'small_only' is same as above, without the user boxes
#   'minimal' is only the boxes asked by the user in 'd_boxes'
SAR_SIZE = 'small' 

# First look -----------------------------------------------------------
PLOT_10m_WIND       = False # 10m wind and SAR roughness for first 3 boxes
WHERE_ARE_THE_BOXES = False # both SAR and LES
# Geometrical analysis -------------------------------------------------
PLOT_2D_COVARIANCE  = True 
PLOT_2D_SKEWNESS    = False
S2_ANALYSIS         = False # this tries to fit ellipse on the 2nd order structure function (Brilouet et al. 2024)
# PLOT_ELLIPSE ?
# Turbulence convergence -----------------------------------------------
VERIFY_TURB_CONVERGENCE = False # plot spectrum at inflow
B_KPSD                  = True  # plot k*PSD(k) ?
altZ_convergence        = 200   # in meters 
liste_X_turb_convergence = [0,2,4,6,8,10] # km, distance from East border of son domain
A                       = 0.1   # y=log(x)**(-coeff)+log(A), constant for Kolmogorov Law in inertial subrange
Kcoeff                  = -5/3  # inertial subrange
# Coeherent structure analysis -----------------------------------------
PLOT_MEAN_PROGVAR       = False # mean profile in each structures
PLOT_MEAN_FLX           = False # mean flux profile with contribution from each structures
PLOT_TopView_with_CS    = False # top view with coherent srtuctures
CS_LENGTH_WITH_LABEL    = False
CS_R_EQUIVALENT         = False  # Computes equivalent radius for each structures.

# folder organisation --------------------------------------------------
workdir_path = '/home/jacqhugo/WORKDIR/MNH570/'
path_here = '/home/jacqhugo/scripts/Linking_SAR_and_ABL_structure/post_process/'
#path_SAR = path_here+'SAR_from_OVL/SAR_roughness_20151210t170827-20151210t170921.nc'
path_SAR = path_here+'SAR_from_IFREMER/S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_9AE1.nc'
path_SST_ODYSEA = path_here + 'SST_from_Ifremer/20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc'
path_data_turb = path_here + 'DATA_TURB/CASE_'+CASE+'/'
path_save_turb_convergence = 'CASE_'+CASE+'/PNGs_turb_convergence/'
path_save_First_look = 'CASE_'+CASE+'/PNGs_First_look/'
path_save_geometric = 'CASE_'+CASE+'/PNGs_geometric/'
path_save_CS = 'CASE_'+CASE+'/PNGs_CS_analysis/'


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
chunksNOHALO_interp = {'time':-1,
                    'level':16,
                    'nj':200,
                    'ni':130}
# print( give_NETCDF4_chunking(temp_path) ) # = [1, 17, 201, 129]
#   if chunksOUT is different from what is saved in the netcdf4, 
#   a warning about performance is raised.
# END INPUT ZONE ========================================================================================


if __name__ == "__main__":  # This avoids infinite subprocess creation
    
    # creating folders
    os.system('mkdir -p '+path_data_turb)
    os.system('mkdir -p '+path_save_turb_convergence)
    os.system('mkdir -p '+path_save_First_look)
    os.system('mkdir -p '+path_save_geometric)
    os.system('mkdir -p '+path_save_CS)

    # DASK related ------------------------------------------------
   
    global client
    client = None

    if DASHBOARD:
        # sometimes dask cluster can cause problems "memoryview is too large"
        # (writing a big netcdf file for eg, hbudget_file)
        cluster = LocalCluster(threads_per_worker=32,n_workers=8)
        
        client = Client(cluster)
        print("Dashboard at :",client.dashboard_link)
    
    # =======================================================================================================
    # BUILDING FILES                                                                                        |
    # =======================================================================================================

    # opening OUTPUT files
    # dsB = 
    dsO = xr.open_mfdataset(dsO_path,chunks=chunksOUT) 
    dsO_i = mass_interpolator(dsO,chunksNOHALO_interp)
    #dsO_i = dsO_i.isel(level=slice(nhalo,-nhalo),ni=slice(nhalo,-nhalo),nj=slice(nhalo,-nhalo))
    
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
    build_CS(path_origin,path_out,dsO_i,dsmean,d_boxes['LES'],gamma=gamma,mCS=mCS)
    dsCS = xr.open_dataset(path_out)


    # Opening SAR file
    dsSAR = xr.open_dataset(path_SAR).sel(pol='VV')
    # Opening SST file
    dsSST = xr.open_dataset(path_SST_ODYSEA)
    X,Y,Z = dsO_i.ni.values,dsO.nj.values,dsO.level.values

    # automatic boxe locator
    d_boxes['SAR'] = SAR_boxes_generator(dsSAR,d_boxes['SAR'],SAR_SIZE)

    # Building structure functions
    #   this process is long. Much longer if you have many boxes and big domains
    N = 2
    save_S_n_SAR(dsSAR,2,d_boxes['SAR'],client,path_data_turb)
    dsS2_SAR = xr.open_dataset(path_data_turb+'S_2_sig0.nc')
    save_S_n_LES(dsCS,'M',10,2,client,path_data_turb)
    dsS2_LES_M10 = xr.open_dataset(path_data_turb+'S_2_M10.nc')
    # N = 3
    # save_S_n_SAR(dsSAR,3,d_boxes['SAR'],path_data_turb)
    # dsS2_SAR = xr.open_dataset(path_data_turb+'S_3_sig0.nc')
    # save_S_n_LES(dsCS,'M',10,3,path_data_turb)
    # dsS2_LES_M10 = xr.open_dataset(path_data_turb+'S_3_M10.nc')

    # =======================================================================================================
    # POST-PROCESS                                                                                          |
    # =======================================================================================================

    # TURB_STATS
    if VERIFY_TURB_CONVERGENCE:
        print('\n* Verification of turbulence at inflow')
        plot_energy_spectrum(dsO,altZ_convergence,liste_X_turb_convergence,B_KPSD,A,Kcoeff,path_save_turb_convergence)
        
    # FIRST LOOK --------------------
    if PLOT_10m_WIND:
        print('\n* 10m wind')
        Wind_10m_in_boxes_vs_SAR(dsCS,dsSAR,indt=-1,d_boxes=d_boxes,path_save=path_save_First_look)
    if WHERE_ARE_THE_BOXES:
        print('\n* Plotting location of boxes on LES and SAR')
        Where_boxes(dsO_i,dsSAR,dsSST,indt=-1,d_boxes=d_boxes,SAR_SIZE=SAR_SIZE,path_save=path_save_First_look)

    # GEOMETRICS --------------------
    if PLOT_2D_COVARIANCE:
        print('\n* Plotting 2D covariance')
        N = 2
        Plot_S_n(dsS2_SAR,'sig0',N,path_txt=path_data_turb,path_save=path_save_geometric)
        Plot_S_n(dsS2_LES_M10,'M10',N,path_txt=path_data_turb,path_save=path_save_geometric)
    if PLOT_2D_SKEWNESS:
        print('\n* Plotting 2D skewness')
        N = 3
        print('TO BE DONE')

    if S2_ANALYSIS:
        print('\n* tests')
        # I need to find a test case for S2_analysis ?
        path_out = 'DATA_TURB/'
        S2_analysis('SAR','sigma0_detrend',dsS2_SAR,d_boxes,path_out)
        S2_analysis('LES','M10',dsS2_LES_M10,d_boxes,path_out)
    # To be done:
    #   - plot the fitted ellipse on S2 for LES and SAR
   
    # COHERENT STRUCTURE ANALYSIS ---
    if PLOT_MEAN_PROGVAR:
        print('\n* Plotting profiles of prognostic variables, with object decomposition')
        Plot_mean_progvar_allboxes(dsCS,dsmean,path_save_CS,Reynolds_avg)
        
    if PLOT_MEAN_FLX:
        print('\n* Plotting profiles of fluxes, with object decomposition')
        Plot_mean_flux_allboxes(dsCS,dsmean,dsflx,path_save_CS,Reynolds_avg)
    
    if PLOT_TopView_with_CS:
        VAR_TOPVIEW = 'W'               # background of the top view
        print('\n* Plotting top view of '+VAR_TOPVIEW+' with coherent structures')
        Plot_top_view_var(dsCS,atZ=250,path_save=path_save_CS) # this is not finished

    if CS_LENGTH_WITH_LABEL:
        print('\n* Find characteristic length scale of coherent structures')

    if CS_R_EQUIVALENT:
        print('\n* Computes an equivalent radius for each coherent structure, for each boxes of the LES.')
        # 1 profile or r_eq for each boxe and for each structure.
        R_equivalent_for_all_objects(dsCS,dsmean,LES_res,path_save_CS)
        


    end = time.time()
    print('Total runtime of analyse.py: '+sec2hms(end - start))
    plt.show()
    
    dsO.close() # avoid memory leaks
    #dsB.close()

    # GOALS -----------
    # comparer échelles caractéristiques SAR vs LES.
    # > SAR
    #       utiliser code PE Brilouet pour avoir une estimation de l'ellipse comme dans son papier.
    #       utiliser fonction de structure en polaire pour récupérer un R.
    #         
    # > LES
    #       récupérer distribution de la taille des updrafts en fonction de l'altitude( 2 méthodes, ellipse ou rayon eq.)
    # ajuster colorbar des fonctions de structure

    # - Compute L Obukhov for both LES and SAR
    # LES : use L_Obukhov(tht0,u_star,surf_wtht_flx) from module_tool
    # SAR : from O'Driscoll paper (10.1029/2023GL104228)

    # - Compute geometric statistics for both LES and SAR
    # L and K statistics
    # -----------------


    # Old but to be reminded ?
    # Régler problème de résolution du SAR -> done
    #   check conversion deg vs km -> ok 
    #   récupérer dataset depuis copernicus puis ouverture avec le package S1 xarray -> non trop long
    #   récupérer données Alexis Mouche ou PE Brilouet ? -> j'ai, à voir les diff avec les données OVL: oui c'est ok