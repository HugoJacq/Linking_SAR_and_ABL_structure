import xarray as xr
import numpy as np
from module_tools import *
import os
import scipy as sp
"""
Objectif : 
Compare spectra/mean quantities for several simulations

Based on: https://gitlab.com/cerfacs/notebooks/-/blob/master/spectral_analysis/1-FFT_correction.ipynb
but also with "Evaluating Mesoscale NWP Models Using Kinetic Energy Spectra " Skamarock 2004, https://doi.org/10.1175/MWR2830.1

Description:

- PLOT_FLUC and PLOT_SPECTRA will plot the desired quantities along SPECTRA_DIR, at atX3 or atY1 position (depending on the direction chosen).
- PLOT_TURB_VALIDATION will plot several spectra of the quantity CHOICE_VAR_TURB_VALID along Y, at several X positions (Validation_atX).
- PLOT_TURB_VALIDATION_GRAD is the same as PLOT_TURB_VALIDATION but for the simu with a SST gradient
- PLOT_PROFILES will plots profiles from the 'nest' and 'homo' simulations (U,V,W,THT,RV,THTV,ET,Er)
    please note that the average operator is not the same:
        -> for 'homo', <> = average over X and Y.
        -> for 'nest', <> = average over full Y but only between atX1 and atX2 in the X direction.
        
"""


# INPUT ZONE ==============================================================================
png_out_path = 'PNGs_verifTurb/'
data_path = {'nest':'/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/test_HR_nesting_NoSSTgradient/',
             'homo_c':'/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/simu_SST_298K_convection_libre/',
             'homo_w':'/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/simu_SST_299.5K_convection_libre/',
             'nest_withgrad':'/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/test_HR_nesting/06_mesonh_2models/test_same_dt_over_dx/',
             'nest_withgrad_openDad':'/mnt/60df13e6-86a2-4cac-b207-817433cddfe5/WORKDIR2/57MNH/test_nesting_dadOpen/06_mesonh_2models/'}
nhalo = 1
atZ = 200       # m, height of spectra
atX1 = 6000     # m, first X index for averages of nested simulations with SST gradient 
atX2 = 10000    # m, second index for averages of nested simulations with SST gradient 
atX = 6000     # m, X position where to plot along Y spectra (for all sim)
atY = 6000     # m, Y position where to plot along Y spectra (for all sim)
atY1 = 0        # m, first Y index for averages of nested simulations with SST gradient 
atY2 = 11000    # m, second Y index for averages of nested simulations with SST gradient 

#------------------------------------------------------------------------------------------
SPECTRA_DIR             = 'X'       # X or Y
PLOT_FLUC               = False     # for 'nest' only, with SPECTRA_DIR
PLOT_FLUC_GRAD          = False      # for 'nest_withgrad' only
PLOT_SPECTRA            = False     # for 'nest' only, with SPECTRA_DIR
#------------------------------------------------------------------------------------------
PLOT_TURB_VALIDATION    = False     # uses 'nest' and 'homo_c'
PLOT_TURB_VALIDATION_GRAD = False    # for 'nest_withgrad' only
CHOICE_VAR_TURB_VALID   = 'W'       # Et or U or W
#   to be used with PLOT_TURB_VALIDATION or PLOT_TURB_VALIDATION
Validation_atX = [1,2,3,4,10] # km, [k for k in range(1,11,1)]
#   to be used with PLOT_TURB_VALIDATION or PLOT_TURB_VALIDATION
#   last X position is the reference of comparison (good if =10km)
#------------------------------------------------------------------------------------------
PLOT_PROFILES           = False     # A nest configuration vs homo_c/homo_w
PLOT_PROFILE_VAR = 'nestgrad' # nest or nestgrad
#------------------------------------------------------------------------------------------
WHICH_GRAD = 'nest_withgrad' # 'nest_withgrad' or 'nest_withgrad_openDad',
#   to be used with PLOT_TURB_VALIDATION_GRAD 
#   or PLOT_PROFILES with PLOT_PROFILE_VAR='nestgrad'
#------------------------------------------------------------------------------------------
VAR_OVERVIEW = False            # plots a top down view of CHOICE_VAR_TURB_VALID at a height of atZ
liste_sim = ['nest','homo_c']   # paper is ['nest','homo_c']
#------------------------------------------------------------------------------------------
WIND_DIRECTION = False   # plot 10m mean wind norm and direction, and ER5 wind and SAR wind
#------------------------------------------------------------------------------------------
CHECK_LOSSY_COMPRESSION = False # determine the number of significan digits to keep 
TEST_LOSSY_COMPRESSION = False
#------------------------------------------------------------------------------------------
CHECK_SPINUP = False # plotting 4h of dad surface wind
#------------------------------------------------------------------------------------------


A=0.001 # y=log(x)**(-coeff)+log(A) for Kolmogorov Law in inertial subrange
coeff = -5/3
B_KPSD = True # plots spectra with k*F(k) if True, F(k) else
dad_over_son_res = 4 # = dx_d/dx_s

dpi = 200 

PAPER = True # plots figures for paper.
# END INPUT ZONE ==========================================================================






if SPECTRA_DIR in ['X']:
    direction='1'
    string_save = 'atY'+str(atY1)+'km'
elif SPECTRA_DIR in  ['Y']:
    direction='2'
    string_save = 'atX'+str(atX)+'km'
else:
    raise Exception('Your choice of SPECTRA_DIR='+SPECTRA_DIR+' is not recognized')

os.system('mkdir -p '+png_out_path)

# Opening files
print('* Opening and initializing fields')
dsmaster = {'nest':xr.open_dataset(data_path['nest'] +'06_mesonh_2models/TEST2.2.001.001.nc'),
            'homo_c':xr.open_dataset(data_path['homo_w'] +'HOMOG.1.001.002.nc'),
            'homo_w':xr.open_dataset(data_path['homo_c'] +'HOMOG.1.001.002.nc'),
            'nest_withgrad':xr.open_dataset(data_path['nest_withgrad']+'TEST2.2.002.001.nc'),
            'nest_withgrad_openDad':xr.open_dataset(data_path['nest_withgrad_openDad']+'TEST2.2.001.001.nc')}

SST = dsmaster['nest_withgrad'].SST[nhalo:-nhalo,nhalo:-nhalo]
ds = {'nest':xr.open_mfdataset([data_path['nest']+'06_mesonh_2models'+'/FICHIERS_OUT/TEST2.2.001.OUT.00'+str(k)+'.nc' for k in range(4,8)]), 
                        # here not the first files bc not yet turbulent
            'homo_c':xr.open_mfdataset(([data_path['homo_c']+'FICHIERS_OUT/'+'HOMOG.1.001.OUT.00'+str(k)+'.nc' for k in range(1,10)]+
                    [data_path['homo_c']+'FICHIERS_OUT/'+'HOMOG.1.001.OUT.0'+str(k)+'.nc' for k in range(10,14)])),
            'homo_w':xr.open_mfdataset(([data_path['homo_w']+'FICHIERS_OUT/'+'HOMOG.1.001.OUT.00'+str(k)+'.nc' for k in range(1,10)]+
                    [data_path['homo_w']+'FICHIERS_OUT/'+'HOMOG.1.001.OUT.0'+str(k)+'.nc' for k in range(10,14)])),
            'nest_withgrad':xr.open_mfdataset([data_path['nest_withgrad']+'/FICHIERS_OUT/TEST2.2.002.OUT.00'+str(k)+'.nc' for k in range(3,8)]),
            'nest_withgrad_openDad':xr.open_mfdataset([data_path['nest_withgrad_openDad']+'/FICHIERS_OUT/TEST2.2.001.OUT.00'+str(k)+'.nc' for k in range(3,8)])}
# Note : nest_withgrad is from +2h to +2h30, but if not compared with other sim its ok

# for simulations with a SST gradient, we need to distinguish cold and warm side
ds['nest_withgrad_c'] = ds['nest_withgrad']
ds['nest_withgrad_w'] = ds['nest_withgrad']
ds['nest_withgrad_openDad_c'] = ds['nest_withgrad_openDad']
ds['nest_withgrad_openDad_w'] = ds['nest_withgrad_openDad']
dsmaster['nest_withgrad_c'] = dsmaster['nest_withgrad']
dsmaster['nest_withgrad_w'] = dsmaster['nest_withgrad']
dsmaster['nest_withgrad_openDad_c'] = dsmaster['nest_withgrad_openDad']
dsmaster['nest_withgrad_openDad_w'] = dsmaster['nest_withgrad_openDad']

name_conversion={'nest':'VALson',
                 'homo_c':'RefVal',}

idnt = {}
ABLH = {}
U,V,W,TKE,THT,RV,THTV = {},{},{},{},{},{},{}
X,X_u,Y,Y_v,Z,Z_w,time = {},{},{},{},{},{},{}
UW_HFLX,UW_VFLX = {},{}
VW_HFLX,VW_VFLX = {},{}
THW_FLX,THVW_FLX,RCONSW_FLX = {},{},{}
UW_FLX,VW_FLX, = {},{}

meanU,meanV,meanW,meanTKE = {},{},{},{}
meanTHT,meanRV,meanTHTV = {},{},{}

u_fluc,v_fluc,w_fluc = {},{},{}
Er,Et,thtv_fluc = {},{},{}
uw,vw,wthtv = {},{},{}
mean_uw,mean_vw,mean_wthtv = {},{},{}
mean_Er,mean_Et = {},{}

# First we find the right time to compare the simulations
for name in ds.keys():
    #ds[name] = xr.open_mfdataset(ds[name])
    X[name],X_u[name] = ds[name].ni[nhalo:-nhalo],    ds[name].ni_u[nhalo:-nhalo]
    Y[name],Y_v[name] = ds[name].nj[nhalo:-nhalo],    ds[name].nj_v[nhalo:-nhalo]
    Z[name],Z_w[name] = ds[name].level[nhalo:-nhalo], ds[name].level_w[nhalo:-nhalo]
    time[name] = ds[name].time

    THT[name] = ds[name].THT[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    RV[name] = ds[name].RVT[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    THTV[name] = Compute_THTV(THT[name],RV[name])

indz = nearest(Z['nest'].values,atZ)
indx1 = nearest(X['nest'].values,atX1)
indx2 = nearest(X['nest'].values,atX2)
indx = nearest(X['nest'].values,atX)
indy = nearest(X['nest'].values,atY)
indy1 = nearest(X['nest'].values,atY1)
indy2 = nearest(X['nest'].values,atY2) # numbering like isnt clean, to be looked at

Xslice = {'nest':slice(indx1,indx2),
        'homo_c':slice(None,None),
        'homo_w':slice(None,None),
        'nest_withgrad':slice(indx1,indx2),
        'nest_withgrad_c':slice(indx1,indx2),
        'nest_withgrad_w':slice(indx1,indx2),
        'nest_withgrad_openDad':slice(indx1,indx2),
        'nest_withgrad_openDad_c':slice(indx1,indx2),
        'nest_withgrad_openDad_w':slice(indx1,indx2)}

Yslice = {'nest':slice(None,None),
        'homo_c':slice(None,None),
        'homo_w':slice(None,None),
        'nest_withgrad':slice(None,None),
        'nest_withgrad_c':slice(indy1,indy),
        'nest_withgrad_w':slice(indy,indy2),
        'nest_withgrad_openDad':slice(None,None),
        'nest_withgrad_openDad_c':slice(indy1,indy),
        'nest_withgrad_openDad_w':slice(indy,indy2)}

indt = {'nest': -1,
        'homo_c': -1, # np.argmin( np.abs(ABLH1-ABLH2) )
        'homo_w':-1,
        'nest_withgrad':-1,
        'nest_withgrad_c':-1,
        'nest_withgrad_w':-1,
        'nest_withgrad_openDad':-1,
        'nest_withgrad_openDad_c':-1,
        'nest_withgrad_openDad_w':-1} 

# Note pour Hugo (13/09/24 16h08):
# Il faut réduire le range sur lequel on cherche l'ABLH pour les cas avec gradient de SST 
# par exemple indz  = [ 10,120] semble ok (donne ABLH5=710m)
shiftB = 50 # m, search ABLH after first level
shiftT = 1100 # m, search ABLH berore top

    
# Computing ABLH
#   - compute <THTV> over ni,nj,2x time steps
#   - compute dTHT/dz
#   - find argmax
#   - ABLH = Z[<argmax>]
for name in ds.keys():
    DTHTVDZ = THTV[name].isel(time=slice(-2,-1), # indt[name]
                            ni=Xslice[name],
                            nj=Yslice[name]
                            ).mean(['ni','nj','time']
                            ).differentiate('level'
                            ).sel(level=slice(shiftB,shiftT))
    z = Z[name]
    itemp = nearest(z.values,shiftB)
    ABLH[name] = z[itemp+int(DTHTVDZ.sel(level=slice(shiftB,shiftT)).argmax('level').mean())].values
    
    # max ABLH detection doesnt work on thoses, i need more points in the real experiment
    if name=='nest_withgrad_c':
        ABLH[name] = 650.
    elif name=='nest_withgrad_w':
        ABLH[name] = 800.
     
     
    # Verification : THTV at ABLH for each domains.
    if False:
        print(name,ABLH[name])
        x,y,z = X[name],Y[name],Z[name]
        indablh = nearest(z.values,ABLH[name])
        fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=200)
        sc = ax.pcolormesh(x.isel(ni=Xslice[name])/1000,y.isel(nj=Yslice[name])/1000,
            THTV[name].isel(time=-1,level=indablh,ni=Xslice[name],nj=Yslice[name]),cmap='jet',vmin=298,vmax=299)
        plt.colorbar(sc,ax=ax,label='THTV K')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Z='+str(ABLH[name])+' '+name)
    
print('* Computing fluctuations,means for simus')
for name in ds.keys():
    # I NEED TO INTERP...
    U[name] = ds[name].UT.interp({'ni_u':ds[name].ni})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    U[name] = U[name].rename(new_name_or_name_dict={'nj_u':'nj'})
    V[name] = ds[name].VT.interp({'nj_v':ds[name].nj})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    V[name] = V[name].rename(new_name_or_name_dict={'ni_v':'ni'})
    W[name] = ds[name].WT.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    TKE[name] = ds[name].TKET[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]

    UW_HFLX[name] = ds[name].UW_HFLX.interp({'level_w':ds[name].level,'ni_u':ds[name].ni})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo] 
    UW_VFLX[name] = ds[name].UW_VFLX.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]					
    VW_HFLX[name] = ds[name].VW_HFLX.interp({'level_w':ds[name].level,'nj_v':ds[name].nj})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo] 	
    VW_VFLX[name] = ds[name].VW_VFLX.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]					
    THW_FLX[name] = ds[name].THW_FLX.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]					
    RCONSW_FLX[name] = ds[name].RCONSW_FLX.interp({'level_w':ds[name].level})[:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    UW_HFLX[name] = UW_HFLX[name].rename(new_name_or_name_dict={'nj_u':'nj'})
    VW_HFLX[name] = VW_HFLX[name].rename(new_name_or_name_dict={'ni_v':'ni'})
    UW_FLX[name] = UW_HFLX[name] + UW_VFLX[name]
    VW_FLX[name] = VW_HFLX[name] + VW_VFLX[name]
    THVW_FLX[name] = THW_FLX[name]*THTV[name]/THT[name] + 0.61*THT[name]*RCONSW_FLX[name] # subgrid  
    
    meanU[name] = Mean2(U[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanV[name] = Mean2(V[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanW[name] = Mean2(W[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanTKE[name] = Mean2(TKE[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanTHT[name] = Mean2(THT[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanRV[name] = Mean2(RV[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    meanTHTV[name] = Mean2(THTV[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    
    u_fluc[name] = U[name] - meanU[name]
    v_fluc[name] = V[name] - meanV[name]
    w_fluc[name] = W[name] - meanW[name]
    thtv_fluc[name] = THTV[name] - meanTHTV[name]
    Er[name] = 0.5*(u_fluc[name]**2+v_fluc[name]**2+w_fluc[name]**2)
    Et[name] = Er[name] + TKE[name]
    uw[name] = u_fluc[name]*w_fluc[name] + UW_FLX[name]
    vw[name] = v_fluc[name]*w_fluc[name] + VW_FLX[name]
    wthtv[name] = w_fluc[name]*thtv_fluc[name] + THVW_FLX[name]
    
    mean_uw[name] = Mean2(uw[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    mean_vw[name] = Mean2(vw[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    mean_wthtv[name] = Mean2(wthtv[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    mean_Er[name] = Mean2(Er[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    mean_Et[name] = Mean2(Et[name].isel(time=indt[name],ni=Xslice[name],nj=Yslice[name]))
    
# pratical for plots
dico_unit = {'U':r'$m.s^{-1}$','W':r'$m.s^{-1}$','Et':r'$m^2.s^{-2}$',
             'THTV':'K','RV':r'$g.kg^{-1}$',
             'uw':r'$m^2.s^{-2}$','wthtv':r'$m.s^{-1}.K^{-1}$',}
dico_spectra_unit = {'U':r'$m^2.s^{-2}$','W':r'$m^2.s^{-2}$)','Et':r'$m^2.s^{-2}$'}
dico_var = {'U':U,'W':W,'Et':Et}
dico_bornes = {'U':[5,7.5],'W':[-2,2],'Et':[0.,0.9],
               'THTV':[297.5,298.5],'RV':[9.5,11.5],
               'uw':[-0.07,0.01],'wthtv':[-0.01,0.03]}
dico_cmap = {'U':'jet','W':'seismic','Et':'plasma'}
    
if PLOT_SPECTRA or PLOT_FLUC or PLOT_TURB_VALIDATION:
    # this is looking at the nested configuration
    if SPECTRA_DIR=='X':
        print('* I need to compute spectra along X for nest at Y =',atY/1000,'km')
        # atZ and atY1
        u_fluc1D = u_fluc['nest'].isel({'level':indz,'nj':indy,'ni':Xslice['nest']})
        v_fluc1D = v_fluc['nest'].isel({'level':indz,'nj':indy,'ni':Xslice['nest']})
        w_fluc1D = w_fluc['nest'].isel({'level':indz,'nj':indy,'ni':Xslice['nest']})
        Et1D = Et['nest'].isel({'level':indz,'nj':indy,'ni':Xslice['nest']})
        # Conditioning
        #   Detrending, this is different from computing fluctuations !!
        Et1D = Et1D - Et['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('nj')
        u_fluc1D = u_fluc1D - u_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('nj')
        v_fluc1D = v_fluc1D - v_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('nj')
        w_fluc1D = w_fluc1D - w_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('nj')
        ABCISSES = X['nest'].isel(ni=Xslice['nest'])
        NAME_ABCISSES = 'X'
            
    if SPECTRA_DIR=='Y':
        print('* I need to compute spectra along Y for nest at X =',atX/1000,'km')
        # atZ and atX3
        u_fluc1D = u_fluc['nest'].isel({'level':indz,'ni':indx})
        v_fluc1D = v_fluc['nest'].isel({'level':indz,'ni':indx})
        w_fluc1D = w_fluc['nest'].isel({'level':indz,'ni':indx})
        Et1D = Et['nest'].isel({'level':indz,'ni':indx})
        # Conditioning
        #   Detrending, this is different from computing fluctuations !!
        Et1D = Et1D - Et['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('ni')
        u_fluc1D = u_fluc1D - u_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('ni')
        v_fluc1D = v_fluc1D - v_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('ni')
        w_fluc1D = w_fluc1D - w_fluc['nest'].isel({'level':indz,'ni':Xslice['nest']}).mean('ni')
        ABCISSES = Y['nest']
        NAME_ABCISSES = 'Y'
    
    # power_spectral_density
    PSD_U = np.zeros((len(time['nest']),int(len(ABCISSES) / 2)-1))
    PSD_V = np.zeros((len(time['nest']),int(len(ABCISSES) / 2)-1))
    PSD_W = np.zeros((len(time['nest']),int(len(ABCISSES) / 2)-1))
    PSD_Et = np.zeros((len(time['nest']),int(len(ABCISSES) / 2)-1))
    for t in range(len(time['nest'])):
        f_u,PSD_U[t,:] = PSD(ABCISSES.values, u_fluc1D[t,:])
        f_v,PSD_V[t,:] = PSD(ABCISSES.values, v_fluc1D[t,:])
        f_w,PSD_W[t,:] = PSD(ABCISSES.values, w_fluc1D[t,:])
        f_w,PSD_Et[t,:] = PSD(ABCISSES.values, Et1D[t,:])    
    
    print('* I need to compute spectra for the homogeneous case')
    print('     if along Y, at X =',atX/1000,'km')
    print('     if along X, at Y =',atY/1000,'km')
    print('     then averaged over several time steps')
    # PSD from homogeneous simulation
    casename = 'homo_c'
    dico_fluc = {'U':u_fluc[casename],'V':v_fluc[casename],'W':w_fluc[casename],'Et':Et[casename]}
    PSD_ref = {}
    for var in dico_fluc.keys():
        fluc_ref = dico_fluc[var]
        PSD_ref[var]= {'X':np.zeros(int(len(X['homo_c'])/2)-1),
                        'Y':np.zeros(int(len(X['homo_c'])/2)-1)}
        for t in range(len(time['homo_c'])): # Adding the spectra from ref (at X = atX3)
            f_ref,temp = PSD(Y['homo_c'].values, fluc_ref[t,indz,:,indx]) 
            PSD_ref[var]['Y']+= temp
        for t in range(len(time['homo_c'])): # Adding the spectra from ref (at Y = atY1)
            f_ref,temp = PSD(X['homo_c'].values, fluc_ref[t,indz,indy,:]) 
            PSD_ref[var]['X']+= temp    
        PSD_ref[var]['Y'] = PSD_ref[var]['Y'] / len(time['homo_c']) # avering over time to smooth spectra
        PSD_ref[var]['X'] = PSD_ref[var]['X'] / len(time['homo_c'])
        
    if PLOT_FLUC:
        print('- Plotting raw data: fluctuations (conditioned) ...')
        fig, ax = plt.subplots(4,1,figsize = (5,10),constrained_layout=True,dpi=200)
        L_VAR = [u_fluc1D,v_fluc1D,w_fluc1D,Et1D]
        L_name = ["u'","v'","w'","Et"]
        for i in range(len(L_VAR)):
            for t in range(len(time['nest'])):
                ax[i].plot(ABCISSES/1000,L_VAR[i][t,:],c='b',alpha=0.1)
                ax[i].set_ylabel(L_name[i])   
                ax[i].set_xlim([ABCISSES[0]/1000,ABCISSES[-1]/1000])
        ax[0].set_title('nest at '+NAME_ABCISSES+'='+str(atX/1000)+' km')
        fig.savefig(png_out_path+'fluc'+SPECTRA_DIR+'.png')   

    if PLOT_SPECTRA:
        print('- Plotting spectra ...')
        figsize = (5,5)
        L_var = [PSD_U,PSD_V,PSD_W,PSD_Et]
        L_name = [r'F$^'+direction+r'_{11}$',r'F$^'+direction+r'_{22}$',r'F$^'+direction+r'_{33}$',r'E$^'+direction+r'$']
        L_namesave = ['U','V','W','Et']
        L_color = ['b','g','r','k']
        # 1 plot for each spectra
        for i in range(len(L_var)):
            fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=200)
            ax.set_xlabel(r'$\frac{1}{\lambda_'+direction+r'}$ (m$^{-1}$)')
            if B_KPSD:
                for t in range(len(time['nest'])):
                    ax.plot(f_u,f_u*L_var[i][t,:],c='k',alpha=0.1)
                ax.plot(f_u,f_u*L_var[i].mean(axis=0),c='b',label='nest') # mean of several spectra
                ax.plot(f_ref,f_ref*PSD_ref[L_namesave[i]][SPECTRA_DIR],c='k',label='homo_c') # PSD from homogeneous sim
                ax.set_ylabel(r'$\frac{1}{\lambda_'+direction+r'}$.'+L_name[i])
                ax.plot(f_u,f_u*(f_u)**(coeff)*A,c='k',ls='--')
            else:
                for t in range(len(time['nest'])):
                    ax.plot(f_u,L_var[i][t,:],c='k',alpha=0.1)
                ax.plot(f_u,L_var[i].mean(axis=0),c='b',label='nest') # mean of several spectra
                ax.plot(f_ref,PSD_ref[L_namesave[i]][SPECTRA_DIR],c='k',label='homo_c') # PSD from homogeneous sim
                ax.set_ylabel(L_name[i])
                ax.plot(f_u,(f_u)**(coeff)*A,c='k',ls='--')
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylim([0.00001,10])
            ax.legend()
            fig.savefig(png_out_path+'PSD'+SPECTRA_DIR+'_'+L_namesave[i]+'_atZ'+str(atZ)+'_'+string_save+'.png')
        
        # spectra combined
        fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=200)
        for i in range(len(L_var)):
            ax.plot(f_u,f_u*L_var[i].mean(axis=0),c=L_color[i],label=L_name[i],alpha=0.5)
        ax.set_ylabel(r'$\frac{1}{\lambda_'+direction+r'}.F^'+direction+r'_{ii}(\frac{1}{\lambda_'+direction+r'})$')
        ax.set_xlabel(r'$\frac{1}{\lambda_'+direction+r'}$ (m$^{-1}$)')
        ax.plot(f_u,f_u*(f_u)**(coeff)*A,c='k',ls='--',label='Inertial subrange')
        ax.set_ylim([0.00001,10])
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title('nest')
        ax.legend()
        fig.savefig(png_out_path+'PSD'+SPECTRA_DIR+'_UVW_atZ'+str(atZ)+'_at'+string_save+'.png')

    if PLOT_TURB_VALIDATION:
        """
        Here we plot spectra (along Y) at some X positions and we compare them with the spectra at the West border.
        Goal : find the fetch of the transition zone from the turbulent scales from dad domain to son domain.
        
        Conclusions : 
        - spectra at X=4km is like at the West border, so the turbulence spinup transition zone is X=[0,4] km.
        """
        print('- Plotting spectra at several X positions for turbulence validation ...')
        print('     z='+str(atZ)+'m, x='+str(Validation_atX)+' km')
        ABCISSES = Y['nest']
        NAME_ABCISSES = 'Y'
        dad_over_son_res = 4 # = dx_d/dx_s
        res = (ABCISSES[1] - ABCISSES[0]).values
        figsize = (5,5)
        figsize2 = (10,3)
        INDIV_PLOTS = False
        GLOBAL_PLOTS = PAPER # True for paper
        
        # Computing PSD of several var at several X
        casename = 'nest'
        dico_fluc = {'U':u_fluc[casename],'W':w_fluc[casename],'Et':Et[casename]}
        dico_name = {'U':r'$F^2_{11}(\frac{1}{\lambda_2})$','W':r'$F^2_{33}(\frac{1}{\lambda_2})$','Et':r'$E^2(\frac{1}{\lambda_2})$'}
        PSD_VAL = np.zeros((len(time['nest']),int(len(ABCISSES) / 2)-1,len(Validation_atX)))
        var_fluc = dico_fluc[CHOICE_VAR_TURB_VALID]
        for i,atX in enumerate(Validation_atX):
            indx = nearest(X['nest'].values,atX*1000)
            var_fluc1D = var_fluc.isel({'level':indz,'ni':indx}) # selecting altitude and X position
            var_fluc1D = var_fluc1D - var_fluc.isel({'level':indz,'ni':Xslice['nest']}).mean('ni') # detrending
            for t in range(len(time['nest'])):
                f_val,PSD_VAL[t,:,i] = PSD(ABCISSES.values, var_fluc1D[t,:])
            
        lw = 1.0
        # Plot
        if B_KPSD:
            coeffX_val,coeffX_ref = f_val,f_ref
            NameY = r'$\frac{1}{\lambda_2}$.'
            KPSD_save = 'K'
            KolmoCoeff = '-2/3' 
        else:
            coeffX_val,coeffX_ref = 1,1
            NameY = ''
            KPSD_save = ''
            KolmoCoeff = '-5/3'
             
        if INDIV_PLOTS:
            # Individual Plots
            for i in range(len(Validation_atX)-1):
                fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=200)
                ax.plot(f_val,coeffX_val*PSD_VAL[:,:,i].mean(axis=0),c='blue',label=str(Validation_atX[i])+'km',alpha=1)
                ax.plot(f_val,coeffX_val*PSD_VAL[:,:,-1].mean(axis=0),c='red',label=str(Validation_atX[-1])+'km',alpha=1)
                ax.plot(f_ref,coeffX_ref*PSD_ref[CHOICE_VAR_TURB_VALID]['Y'],c='k',label='SST=cst',alpha=1)
                ax.set_ylabel(NameY+dico_name[CHOICE_VAR_TURB_VALID]+' '+dico_spectra_unit[CHOICE_VAR_TURB_VALID])
                ax.set_xlabel(r'$\frac{1}{\lambda_2}$ (m$^{-1}$)')
                ax.vlines(1/(4*res),0,100,colors='grey',ls='--',label=r'4$\Delta_x$') # effective resolution from son
                ax.plot(f_val,f_val*(f_val)**(coeff)*A,c='k',ls='--',label=KolmoCoeff)
                ax.set_ylim([0.00001,10])
                ax.set_yscale('log')
                ax.set_xscale('log')
                ax.set_title('Turbulence validation')
                ax.legend()
                fig.savefig(png_out_path+KPSD_save+'PSD_'+CHOICE_VAR_TURB_VALID+'_atZ'+str(atZ)+'_atX'+str(Validation_atX[i])+'.png')
        
        # Paper plot   
        if GLOBAL_PLOTS:
            fig, ax = plt.subplots(1,len(Validation_atX)-1,figsize = figsize2,constrained_layout=True,dpi=200)
            stringsave = ''
            for i in range(len(Validation_atX)-1):
                ax[i].plot(f_val,coeffX_val*PSD_VAL[:,:,i].mean(axis=0),c='blue',label=str(Validation_atX[i])+'km',alpha=1,lw=lw)
                ax[i].plot(f_val,coeffX_val*PSD_VAL[:,:,-1].mean(axis=0),c='red',label=str(Validation_atX[-1])+'km',alpha=1,lw=lw)
                ax[i].plot(f_ref,coeffX_ref*PSD_ref[CHOICE_VAR_TURB_VALID]['Y'],c='k',label='SST=cst',alpha=1,lw=lw)
                ax[i].vlines(1/(4*res),0,100,colors='grey',ls='--',label=r'4$\Delta_x$') # effective resolution from son
                ax[i].plot(f_val,f_val*(f_val)**(coeff)*A,c='k',ls='--',label=KolmoCoeff)
                ax[i].set_ylim([0.00001,10])
                ax[i].set_yscale('log')
                ax[i].set_xscale('log')
                ax[i].set_xlabel(r'$\frac{1}{\lambda_2}$ (m$^{-1}$)')
                stringsave = stringsave + str(Validation_atX[i])
                ax[i].set_title('X = '+str(Validation_atX[i])+' km',loc='right')
                if i!=0:
                    ax[i].tick_params(axis='y',labelleft=False)
            ax[0].set_ylabel(NameY+dico_name[CHOICE_VAR_TURB_VALID]+' '+dico_spectra_unit[CHOICE_VAR_TURB_VALID])
            ax[0].legend()
            #fig.suptitle('Turbulence validation')
            fig.savefig(png_out_path+KPSD_save+'PSD_'+CHOICE_VAR_TURB_VALID+'_atZ'+str(atZ)+'_X'+stringsave+'vs'+str(Validation_atX[-1])+'km.png')    
                            
if PLOT_TURB_VALIDATION_GRAD:
    """Description
    Here we do the same as PLOT_TURB_VALIDATION but without the reference homogeneous simulations (2 would be needed to showcase both cold and warm side of the SST front). 
    What is plotted ? spectra of [U,V,W] along  at several X, compared to the spectra at X = East border.
    """
    casename = WHICH_GRAD
    print('- Plotting spectra at several X positions for turbulence validation ('+casename+') ...')
    print('     z='+str(atZ)+'m, x='+str(Validation_atX)+' km')
    ABCISSES = Y[casename]
    NAME_ABCISSES = 'Y'
    res = (ABCISSES[1] - ABCISSES[0]).values
    figsize = (5,5)
    
    # Computing PSD of several var at several X
    
    dico_fluc = {'U':u_fluc[casename],'W':w_fluc[casename],'Et':Et[casename]}
    dico_name = {'U':r'$F^2_{11}(\frac{1}{\lambda_2})$','W':r'$F^2_{33}(\frac{1}{\lambda_2})$','Et':r'$E^2(\frac{1}{\lambda_2})$'}
    dico_unit = {'U':r'($m^2.s^{-2}$)','W':r'($m^2.s^{-2}$)','Et':r'($m^2.s^{-2}$)'}
    PSD_VAL = np.zeros((len(time[casename]),int(len(ABCISSES) / 2)-1,len(Validation_atX)))
    var_fluc = dico_fluc[CHOICE_VAR_TURB_VALID]
    for i,atX in enumerate(Validation_atX):
        indx = nearest(X[casename].values,atX*1000)
        var_fluc1D = var_fluc.isel({'level':indz,'ni':indx}) # selecting altitude and X position
        var_fluc1D = var_fluc1D - var_fluc.isel({'level':indz,'ni':Xslice[WHICH_GRAD]}).mean('ni') # detrending
        for t in range(len(time[casename])):
            f_val,PSD_VAL[t,:,i] = PSD(ABCISSES.values, var_fluc1D[t,:])
    
    
    if PLOT_FLUC_GRAD:
        #print('- Plotting fluctuations for nest_withgrad')
         # atZ and atX3
        u_fluc1D = u_fluc[casename].isel({'level':indz,'ni':indx})
        v_fluc1D = v_fluc[casename].isel({'level':indz,'ni':indx})
        w_fluc1D = w_fluc[casename].isel({'level':indz,'ni':indx})
        Et1D = Et[casename].isel({'level':indz,'ni':indx})
        # Conditioning
        #   Detrending, this is different from computing fluctuations !!
        Et1D = Et1D - Et[casename].isel({'level':indz,'ni':Xslice[WHICH_GRAD]}).mean('ni')
        u_fluc1D = u_fluc1D - u_fluc[casename].isel({'level':indz,'ni':Xslice[WHICH_GRAD]}).mean('ni')
        v_fluc1D = v_fluc1D - v_fluc[casename].isel({'level':indz,'ni':Xslice[WHICH_GRAD]}).mean('ni')
        w_fluc1D = w_fluc1D - w_fluc[casename].isel({'level':indz,'ni':Xslice[WHICH_GRAD]}).mean('ni')
        
        fig, ax = plt.subplots(5,1,figsize = (5,10),constrained_layout=True,dpi=dpi)
        L_VAR = [u_fluc1D,v_fluc1D,w_fluc1D,Et1D,SST]
        L_name = ["u'","v'","w'","Et","SST"]
        for i in range(len(L_VAR)):
            if i<len(L_VAR)-1:
                for t in range(len(time[casename])):
                    ax[i].plot(ABCISSES/1000,L_VAR[i][t,:],c='b',alpha=0.1)
            else:
                ax[i].plot(ABCISSES/1000,SST[:,indx],c='k')
                ax[i].set_xlabel(NAME_ABCISSES+' (km)')
            ax[i].set_ylabel(L_name[i])   
            ax[i].set_xlim([ABCISSES[0]/1000,ABCISSES[-1]/1000])
        ax[0].set_title(casename+' at '+NAME_ABCISSES+'='+str(atX/1000)+' km')
        fig.savefig(png_out_path+'fluc'+SPECTRA_DIR+'_'+casename+'.png') 
    
    # Plot  
    if B_KPSD:
        coeffX_val,coeffX_ref = f_val,f_ref
        NameY = r'$\frac{1}{\lambda_2}$.'
        KPSD_save = 'K'
        KolmoCoeff = '-2/3'
    else:
        coeffX_val,coeffX_ref = 1,1
        NameY = ''
        KPSD_save = ''
        KolmoCoeff = '-5/3'
            
    for i in range(len(Validation_atX)-1):
        fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
        ax.plot(f_val,coeffX_val*PSD_VAL[:,:,i].mean(axis=0),c='blue',label=str(Validation_atX[i])+'km',alpha=1)
        ax.plot(f_val,coeffX_val*PSD_VAL[:,:,-1].mean(axis=0),c='red',label=str(Validation_atX[-1])+'km',alpha=1)
        #ax.plot(f_ref,f_ref*PSD_ref[CHOICE_VAR_TURB_VALID]['Y'],c='k',label='SST=cst',alpha=1) # non sens when dSST/dy =/= 0
        ax.set_ylabel(NameY+dico_name[CHOICE_VAR_TURB_VALID]+' '+dico_spectra_unit[CHOICE_VAR_TURB_VALID])
        ax.set_xlabel(r'$\frac{1}{\lambda_2}$ (m$^{-1}$)')
        ax.vlines(1/(4*res),0,100,colors='grey',ls='--',label=r'4$\Delta_x$') # effective resolution from son
        ax.plot(f_val,coeffX_val*(f_val)**(coeff)*A,c='k',ls='--',label=KolmoCoeff)
        ax.set_ylim([0.00001,10])
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title('Turbulence validation')
        ax.legend()
        fig.savefig(png_out_path+KPSD_save+'PSD_'+casename+'_'+CHOICE_VAR_TURB_VALID+'_atZ'+str(atZ)+'_atX'+str(Validation_atX[i])+'.png')

if PLOT_PROFILES:
    print('- Plotting profiles')
    colors = {"homo_c":'b','homo_w':'r',WHICH_GRAD+'_c':'b',WHICH_GRAD+'_w':'r','nest':'g'}
    ls = {"homo_c":'--','homo_w':'--',WHICH_GRAD+'_c':'-',WHICH_GRAD+'_w':'-','nest':'-'}
    
    if PLOT_PROFILE_VAR=='nest':
        list_plot = ['nest','homo_c']
        namesave = 'nest'
    elif PLOT_PROFILE_VAR=='nestgrad':
        list_plot = ["homo_c",'homo_w',WHICH_GRAD+'_c',WHICH_GRAD+'_w']
        namesave = WHICH_GRAD
        
    # Plots
    if PAPER:
        list_plot = ['nest','homo_c']
        
        namesave = 'nest'
        figsize = (8,5)
        fig, ax = plt.subplots(2,3,figsize = figsize,constrained_layout=True,dpi=dpi)
        axe = ax.flatten()
        for name in list_plot:
            Znormed = Z[name]/ABLH[name]
            axe[0].plot(meanU[name],Znormed,c=colors[name],label=name_conversion[name],ls=ls[name])
            axe[0].set_xlabel('U ('+dico_unit['U']+')')
            axe[0].set_xlim(dico_bornes['U'])
            axe[1].plot(meanTHTV[name],Znormed,c=colors[name],ls=ls[name])
            axe[1].set_xlabel('THTV ('+dico_unit['THTV']+')')
            axe[1].set_xlim(dico_bornes['THTV'])
            axe[2].plot(meanRV[name]*1000,Znormed,c=colors[name],ls=ls[name])
            axe[2].set_xlabel('RV ('+dico_unit['RV']+')')
            axe[2].set_xlim(dico_bornes['RV'])
            axe[3].plot(mean_uw[name],Znormed,c=colors[name],ls=ls[name])
            axe[3].set_xlabel(r'$<uw>$ ('+dico_unit['uw']+')')
            axe[3].set_xlim(dico_bornes['uw'])
            axe[4].plot(mean_wthtv[name],Znormed,c=colors[name],ls=ls[name])
            axe[4].set_xlabel(r'$<w\theta_v>$ ('+dico_unit['wthtv']+')')
            axe[4].set_xlim(dico_bornes['wthtv'])
            axe[5].plot(mean_Et[name],Znormed,c=colors[name],ls=ls[name])
            axe[5].set_xlabel(r'Et ('+dico_unit['Et']+')')
            axe[5].set_xlim(dico_bornes['Et'])
        for axis in axe:
            axis.set_ylim([0,1.2])
            axis.grid()
        ax[0,0].set_ylabel(r'z/$z_i$')
        ax[1,0].set_ylabel(r'z/$z_i$')
        ax[0,0].legend()
        ax[0,1].set_xticks([297.5,297.75,298,298.25,298.5])
        ax[0,1].set_xticklabels(['297.5','','298','','298.5'])
        fig.savefig(png_out_path+'Verif_combined_nestVShomo_c.png')
    else:
        figWind, axWind = plt.subplots(1,3,figsize = (10,5),constrained_layout=True,dpi=dpi)
        figSc, axSc = plt.subplots(1,3,figsize = (10,5),constrained_layout=True,dpi=dpi)
        figVar, axVar = plt.subplots(1,3,figsize = (10,5),constrained_layout=True,dpi=dpi)
        for name in list_plot:
            Znormed = Z[name]/ABLH[name]
            print(name,ABLH[name])
            
            # Wind
            axWind[0].plot(meanU[name],Znormed,c=colors[name],label=name,ls=ls[name])
            axWind[1].plot(meanV[name],Znormed,c=colors[name],ls=ls[name])
            axWind[2].plot(meanW[name],Znormed,c=colors[name],ls=ls[name])
            axWind[0].set_xlim([5,7.5])
            axWind[0].set_xlabel('U m/s')
            axWind[1].set_xlim([-0.5,0.5])
            axWind[1].set_xlabel('V m/s')
            axWind[2].set_xlim([-0.05,0.05])
            axWind[2].set_xlabel('W m/s')
            for axe in axWind.flatten():
                axe.set_ylim([0,1.2])
                axe.set_ylabel(r'z/$z_i$')
            axWind[0].legend()
            
            # Scalars
            axSc[0].plot(meanTHT[name],Znormed,c=colors[name],label=name,ls=ls[name])
            axSc[1].plot(meanRV[name]*1000,Znormed,c=colors[name],ls=ls[name])
            axSc[2].plot(meanTHTV[name],Znormed,c=colors[name],ls=ls[name])
            axSc[0].set_xlim([295.5,296.5])
            axSc[0].set_xlabel('THT K')
            axSc[1].set_xlim([7,12.5])
            axSc[1].set_xlabel('RV g/kg')
            axSc[2].set_xlim([297.5,298.5])
            axSc[2].set_xlabel('THTV K')
            for axe in axSc.flatten():
                axe.set_ylim([0,1.2])
                axe.set_ylabel(r'z/$z_i$')
            axSc[0].legend()
            
            # variances  
            axVar[0].vlines(0,0,2,colors='grey',ls='--')
            axVar[0].plot(mean_uw[name],Znormed,c=colors[name],label=name,ls=ls[name])
            #axVar[0].plot(mean_vw[name],Znormed,c=colors[name],alpha=0.5,ls=ls[name])
            axVar[1].vlines(0,0,2,colors='grey',ls='--')
            axVar[1].plot(mean_wthtv[name],Znormed,c=colors[name],ls=ls[name])
            #axVar[2].plot(mean_Er[name],Znormed,c=colors[name],alpha=0.5,ls=ls[name])
            axVar[2].plot(mean_Et[name],Znormed,c=colors[name],ls=ls[name])
            axVar[0].set_xlim([-0.07,0.01])
            axVar[0].set_xlabel('uw(-)m2/s2') # ,vw(--) 
            axVar[1].set_xlim([-0.01,0.035])
            axVar[1].set_xlabel('wthtv K.m/s')
            axVar[2].set_xlim([-0.01,0.9])
            axVar[2].set_xlabel('Et m2/s2')
            for axe in axVar.flatten():
                axe.set_ylim([0,1.2])
                axe.set_ylabel(r'z/$z_i$')
            axVar[0].legend()
        figWind.savefig(png_out_path+'Verif_'+namesave+'_wind.png')
        figSc.savefig(png_out_path+'Verif_'+namesave+'_scalar.png')
        figVar.savefig(png_out_path+'Verif_'+namesave+'_turbulence.png')
           
if VAR_OVERVIEW:
    """Description
    Plotting a top down view of CHOICE_VAR_TURB_VALID at a height of 'atZ',
        for specific simulations
    """
    print('* Plotting overview of '+CHOICE_VAR_TURB_VALID)
    VAR = CHOICE_VAR_TURB_VALID
    figsize = (6,3)
    fig, ax = plt.subplots(1,len(liste_sim),figsize = figsize,constrained_layout=True,dpi=dpi)
    savename = ''
    for i,name in enumerate(liste_sim):
        sc = ax[i].pcolormesh(X[name]/1000,Y[name]/1000,dico_var[VAR][name].isel(time=indt[name],level=indz),
                cmap=dico_cmap[VAR],vmin=dico_bornes[VAR][0],vmax=dico_bornes[VAR][-1])
        ax[i].set_xlabel('X (km)')
        savename =  savename + '_'+ name
    if True: ax[0].vlines(4,X['nest'][0]/1000,X['nest'][-1]/1000,colors='k',linestyles='--')
    plt.colorbar(sc,ax=ax[-1],label=VAR+' '+dico_unit[VAR],aspect=50,pad=0.01)
    ax[0].set_ylabel('Y (km)')
    fig.savefig(png_out_path+'Overview_'+VAR+'_atZ'+str(atZ)+'m'+savename+'.png')

if WIND_DIRECTION:
   # plot 10m mean wind norm and direction, and ER5 wind and SAR wind
   
    fig = plt.figure(figsize = (10,10),constrained_layout=True,dpi=dpi,)
    ax = fig.add_subplot(projection='polar')
    print('name, angle, norme')
    for name in ds.keys():
        indz_10m = nearest(Z[name].values,10)
        angle = np.arctan( meanV[name].isel(level=indz_10m) / meanU[name].isel(level=indz_10m) )
        angledeg = angle * 180 / np.pi
        norme = np.sqrt( meanU[name].isel(level=indz_10m)**2 + meanV[name].isel(level=indz_10m)**2 )
        ax.scatter(angle,norme,label=name)
        
        
   
        #print(name,angledeg.values,norme.values)
     
    # ERA5 profiles at -35N 25.5E (indexes 17 and 11)
    indLAT = 17
    indLON = 11
    z_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_z_AGULLAS.nc')['z'][0,:,indLAT,indLON]/9.81
    q_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_q_AGULLAS.nc')['q'][0,:,indLAT,indLON]
    p_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_q_AGULLAS.nc')['level']*100
    T_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_t_AGULLAS.nc')['t'][0,:,indLAT,indLON]
    u_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_u_AGULLAS.nc')['u'][0,:,indLAT,indLON]
    v_era5 = xr.open_dataset('/home/jacqhugo/scripts/simu_alex/Contexte_grande_echelle/ERA5/ERA5_v_AGULLAS.nc')['v'][0,:,indLAT,indLON]
	
    rv_era5 = q_era5/(1-q_era5)
    TH_era5 = T_to_Theta(T_era5,p_era5) # T_era5/(p_era5/100000)**(287.05/1004)
    THv_era5 = Compute_THTV(TH_era5,rv_era5) #TH_era5*(1+1.61*rv_era5)/(1+rv_era5) 
    
    print('first level era5 (m) = ',z_era5[-1].values)
    angle = np.arctan( v_era5[-1] / u_era5[-1] )
    angledeg = angle * 180 / np.pi
    norme = np.sqrt( u_era5[-1]**2 + v_era5[-1]**2 )
    
    print(norme.values,angledeg.values)
    ax.scatter(angle,norme,label='ERA5',marker='x')
    
    ax.set_ylim([0,6.5])
    ax.legend(loc='upper left')
    fig.savefig(png_out_path+'Wind_direction.png')
   
if CHECK_LOSSY_COMPRESSION:
    print('* Checking lossy compression ...')
    
    """
    Goal : determine the number of "significant digit" for each prognostic variables.
    Then : using nco command line, compress OUTPUT files to reduce file size on disk.
    
    ncks -O -7 -L 1 --ppc default=3#PABST,THT=5 filename newfilename
    
    Working on one of the double periodic simulation (homo_w)
    
    roadmap:
    - pdf of fluctuations at z = 100m for each var in OUTPUT file
    - test : 10th and 90th percentil to determine cutoff.
    """
    
    path_save = 'PNGs_check_compression/'
    
    ID = 'homo_w'
    altZ = 100 # m
    atT = -1
    Z = ds[ID].level[nhalo:-nhalo]
    LM = ds[ID].LM[atT,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    P = ds[ID].PABST[atT,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    SV1 = ds[ID].SVT001[atT,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    SV3 = ds[ID].SVT003[atT,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
    
    u_f = U[ID][atT] - U[ID][atT].mean(['ni','nj'])
    v_f = V[ID][atT] - V[ID][atT].mean(['ni','nj'])
    w_f = W[ID][atT] - W[ID][atT].mean(['ni','nj'])
    tht_f = THT[ID][atT] - THT[ID][atT].mean(['ni','nj'])
    rv_f = RV[ID][atT] - RV[ID][atT].mean(['ni','nj'])
    sv1_f = SV1 - SV1.mean(['ni','nj'])
    sv3_f = SV3 - SV3.mean(['ni','nj'])
    lm_f = LM - LM.mean(['ni','nj'])
    P_f = P - P.mean(['ni','nj'])
    #Lvar = [u_f,v_f,w_f,tht_f,rv_f,sv1_f,sv2_f,sv3_f,lm_f]
    Lvar = ["w"]
    d_var = {"u":{'f':u_f,'mean':U[ID][atT].mean(['ni','nj'])},
             "v":{'f':v_f,'mean':V[ID][atT].mean(['ni','nj'])},
             "w":{'f':w_f,'mean':W[ID][atT].mean(['ni','nj'])},
             "tht":{'f':tht_f,'mean':THT[ID][atT].mean(['ni','nj'])},
             "rv":{'f':rv_f,'mean':RV[ID][atT].mean(['ni','nj'])},
             "sv1":{'f':sv1_f,'mean':SV1.mean(['ni','nj'])},
             "sv3":{'f':sv3_f,'mean':SV3.mean(['ni','nj'])},
             "lm":{'f':lm_f,'mean':LM.mean(['ni','nj'])},
             "P":{'f':P_f,'mean':P.mean(['ni','nj'])}
            }
    
    valPDFmax = {"u":2,
             "v":2,
             "w":2,
             "tht":0.2,
             "rv":0.0015,
             "sv1":8,
             "sv3":10,
             "lm":30,
             "P":1
            }
    units = {"u":'m/s',
             "v":'m/s',
             "w":'m/s',
             "tht":'K',
             "rv":'kg/kg',
             "sv1":'kg/kg',
             "sv3":'kg/kg',
             "lm":'m',
             "P":'Pa'
            }
    print('At z='+str(altZ)+'m, more than 90% of the population has a |fluc| > than :')
    indz = nearest(Z.values,altZ)
        
    L_10perc = {} # absolute amplitude of fluctuations as a function of Z for each variable
    for name in d_var.keys(): # d_var.keys()
        var = d_var[name]
        X3D = var['f']
        X = np.abs(X3D.sel(level=altZ,method='nearest').values.flatten()) # ampltitude des fluctuations
        kde = sp.stats.gaussian_kde(X)
        N = 400
        X_eval = np.linspace(0,valPDFmax[name],N)
        maxvalue = valPDFmax[name] # np.max([np.amin(X3D),np.amax(X3D)])
        L_10perc[name] = np.percentile(np.abs(X3D),10,axis=[1,2])
        perc10 = L_10perc[name][indz]
        print(' ',name,perc10,units[name])
        
        fig,ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=dpi)
        ax.plot(X_eval, kde.pdf(X_eval),c='k',label='density')
        ax.set_xlabel('|'+name+"'|")
        ax.set_title('PDF at z='+str(altZ)+'m\n'+name+"'(10th)="+str(perc10))
        ax.vlines(perc10,0,kde.pdf(X_eval).max(),
                colors='grey',linestyles='--',label='10th percentile')
        ax.legend()
        ax.set_xlim([0,maxvalue])
        fig.savefig(path_save+'PDF_'+name+'_fluc_atZ'+str(altZ)+'.png')
        
        perc10 = np.percentile(X3D,10,axis=(1,2))
        perc90 = np.percentile(X3D,90,axis=(1,2))
        
        fig,ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=dpi)
        ax.plot(var['mean']+perc10,Z,c='b',label='10th')
        ax.plot(var['mean']+perc90,Z,c='r',label='90th')
        ax.plot(var['mean'],Z,c='k',label='mean')
        ax.set_ylim([0,1000])
        ax.set_xlabel(name)
        ax.set_ylabel('Z')
        fig.savefig(path_save+'profile_'+name+'.png')
    
        fig,ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=dpi)
        ax.plot(L_10perc[name],Z,c='k')
        ax.set_xlabel('|'+name+"'|")
        ax.set_ylabel('Z')
        ax.set_title('10th percentile of absolute fluctuation amplitude')
        fig.savefig(path_save+'profile_10th_abs_fluc_'+name+'.png')
      
if TEST_LOSSY_COMPRESSION:
    """
    test compression/packing with xarray
    usefull doc : 
    [1] https://www.unidata.ucar.edu/blogs/developer/entry/compression_by_scaling_and_offfset
    [2] ]https://www.unidata.ucar.edu/software/netcdf/workshops/2010/bestpractices/Packing.html
    [3] : https://nco.sourceforge.net/nco.html#Precision_002dPreserving-Compression
    base size is 4.5Go
    /home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas/06_mesonh_2models/FICHIERS_OUT/DAD_600s_OUT.nc
    
    1. encoding all var with (=packing) : {'dtype': 'int16', 'scale_factor': 0.1, '_FillValue': -9999}
    result : 2.3Go (x1.95), needs to be tuned to keep the right precision.
    
    2. encoding all var with : {'zlib': True} (default complevel = 4)
    result : 2.6Go (x1.7)
    
    3. encoding (=packing) as P = scale*X + offset with int16
    this give a scale (= resolution) that is different for each variables
    result : 2.3Go (x1.95)
    needs to be tweaked more to get nicer results...
    
    4. zlib but with complevel 7
    very long to compress... and not better than complevel 4
    results 2.6Go (x1.7)
    
    5. save as Zarr (complevel = 5)
    > need to use engine='zarr' when opening again the file 
    2.8Go (x1.6)
    
    6. Using NCO (with Nathan's number of NSD)
    quit long to do but memory is ok
    1.3Go (x3.46) but errors in the approximations !!
    
    7. Using NCO with tuned parameters (see CHECK_LOSSY_COMPRESSION)
    should preserve at least 90% of fluctuations of each variables.
    quit long to do but memory is ok
    1.8Go (x2.64)
    This seems a good choice !
    
    8. NCO and zlib ? -> Nathan did that

    -----------------------------------------------
    CONCLUSION:

    * The recommanded command is : 
        ncks -O -7 -L 1 --ppc default=4#PABST=7#RVT=5#THT=6

        which gives NSD = 4 by default
                            7 for pressure
                            5 for RVT
                            6 for THT
    
         writing speed is ~ 50Mo/s  


    * The compression with nco is good. You need to figure out the Number of Significant Digit (NSD). This can be done as in CHECK_LOSSY_COMPRESSION.
        For ex:
            UT coded in float32 : 5.658712 is NSD = 7, but if you change units it has stil NSD=7 : 0.005658712 km/s
                then, lets say you want at least 0.01 m/s fluctuations to be well represented, you need at minimum NSD = 3,
                which gives 0.005 as the smallest number well represented.

    * NSD is different from DSD (Decimal Significant Digit):
            When using DSD, units matter and so 5.658712 is DSD=+6 but 0.005658712 has DSD=+9 and 5.0 has DSD=-1            
            Differences between NSD and DSD are well explained in [3].

    * How to know which NSD is necessary for your case ?

            - get order of magnitude of variable
            - get the amplitude of minimum fluctuation 


            Example : atmospheric boundary layer, zonal wind speed (U)

                > U is ranging from 3 to 8 m/s, so lets say 1 to 10 m/s
                > u' = U - <U>, plot pdf of absolute amplitude of u' at each altitude.
                    set a threshold : i used 10th percentile
                        this threshold has no physical source, it could be improved i think.
                        (idea : represente at least X amount of energy ?)

                    u'(10th) is ranging from 0.01 -> 0.07 m/s

                > i need to represent at minima 0.01 m/s so i chose NSD = 3.
                    (just to be sure i set NSD = 4, it doesnt take much more space !)

                    

    
    -----------------------------------------------
    """
    
    CHOICE = '7'
    
    basepath = '/home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas/06_mesonh_2models/FICHIERS_OUT/'
    namein = 'DAD_600s_OUT.nc'
    nameout= 'DAD_600s_OUT_scale_offset_'+CHOICE+'.nc'
    
    
    dic_encode = {'1':{'dtype': 'int16', 'scale_factor': 0.1, '_FillValue': -9999},
                  '2':{'zlib': True},
                  '3':{'dtype': 'int16', 'scale_factor':0, 'add_offset':0,'_FillValue': -9999},
                  '4':{'zlib': True,'complevel':7},
                  '5':None,
                  '6':'ncks -O -7 -L 1 --ppc default=3#PABST,THT=5 '+basepath+namein+' '+basepath+nameout,
                  '7':'ncks -O -7 -L 1 --ppc default=4#PABST=7#RVT=5#THT=6 '+basepath+namein+' '+basepath+nameout
                  }
    ds = xr.open_dataset(basepath+namein)
    L_encode = dic_encode[CHOICE]
    L_notvar = ['MNHVERSION','MASDEV','BUGFIX', 'BIBUSER', 'PROGRAM', 'STORAGE_TYPE', 'MY_NAME', 'DAD_NAME','ni','nj','ni_u','nj_u','ni_v','nj_v','level','level_w']
    L_namevar = [name for name in list(ds.keys()) if name not in L_notvar]
    
    #
    if CHOICE=='3': 
        # following the usefull doc
        encode = {}
        for name in L_namevar: # L_namevar
            encode[name] = {} 
            VAR = ds[name][:,nhalo:-nhalo,nhalo:-nhalo,nhalo:-nhalo]
            dataMax = VAR.max().values
            dataMin = VAR.min().values
            scale = (dataMax - dataMin) / (2**16 - 1)
            encode[name]['dtype'] = L_encode['dtype']
            encode[name]['_FillValue'] = L_encode['_FillValue']
            encode[name]['scale_factor'] = scale
            encode[name]['add_offset'] = dataMin
            print(type(scale))
    else:
        encode = {key:L_encode for key in L_namevar}
    
    
    #for name in encode.keys():
        #print(encode[name])
    
    # saving
    if nameout in os.listdir(basepath):
        print(' ->File already exist !\n        Exiting ...')
    else:
        print('saving')
        if CHOICE=='5':
            ds.to_zarr(basepath+nameout)
        elif CHOICE=='6' or CHOICE=='7':
            os.system(dic_encode[CHOICE])
        else:
            ds.to_netcdf(basepath+nameout,encoding=encode)
     
    # checking that it works
    if CHOICE=='5':
        ds2 = xr.open_dataset(basepath+nameout,engine='zarr')
    else:
        ds2 = xr.open_dataset(basepath+nameout)
    print(ds2.UT.encoding)
    fig,ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
    s = ax.pcolormesh(ds2.ni/1000,ds2.nj/1000,ds2.UT[0,nhalo+2,:,:],cmap='Greys_r',vmin=3,vmax=7)
    plt.colorbar(s,ax=ax)
    ax.set_xlabel('ni')
    ax.set_ylabel('nj')
    ax.set_title('UT after decoding')
    
    fig,ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
    s = ax.pcolormesh(ds2.ni/1000,ds2.nj/1000,ds2.UT[0,nhalo+2,:,:] - ds.UT[0,nhalo+2,:,:],cmap='bwr',vmin=-0.1,vmax=0.1)
    plt.colorbar(s,ax=ax)
    ax.set_xlabel('ni')
    ax.set_ylabel('nj')
    ax.set_title('anomaly UT (after - before) ')
    for name in L_namevar:
        print(name,' DIFF max =', (np.abs(ds2[name][0,nhalo+2,:,:] - ds[name][0,nhalo+2,:,:])).max().values)
    
if CHECK_SPINUP:

    basepath = '/home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas_1/02_run_dad_spinup/FICHIERS_OUT/'
    basepath2 = '/home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas_1/02bis_run_dad_spinup/FICHIERS_OUT/'
    pathspin = ([basepath + k for k in list(os.listdir( '/home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas_1/02_run_dad_spinup/FICHIERS_OUT/'))] +
                [basepath2 + k for k in list(os.listdir( '/home/jacqhugo/WORKDIR/MNH570/simu_nest_Agulhas_1/02bis_run_dad_spinup/FICHIERS_OUT/'))] )

    ds = xr.open_mfdataset(pathspin)
    os.system('mkdir -p VERIF_SPINUP_DAD')
    for k in range(len(pathspin)-1):
        fig,ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
        s = ax.pcolormesh(ds.ni/1000,ds.nj/1000,ds.UT[k,1],cmap='Greys_r',vmin=3,vmax=7)
        plt.colorbar(s,ax=ax)
        ax.set_aspect(1)
        ax.set_xlabel('X km')
        ax.set_ylabel('Y km')
        ax.set_title('UT at '+str(ds.time[k].values)[:-10])
        fig.savefig('VERIF_SPINUP_DAD/'+str(k)+'png')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#plt.show()
