import numpy as np
from module_cst import *
from module_tools import *
import pathlib

"""
To be used with analyse.py
"""

def Reynolds_avg(array):
    """
	Reynolds averaged operator (spatial, X and Y)
	
    IN : 
        - xr.DataArray
    OUT:
        - meaned array
    """
    if 'coord_x' in array.dims:
        DIMS = ['time','coord_x','coord_y']
    else:
        DIMS =['time','ni','nj']
    return array.mean(dim=DIMS)

def mass_interpolator(dsO):
    """
    This operator is interpolating variables at mass point from a C-grid used in MesoNH

    INPUT : 
        - ds is a xarray dataset, typically a OUTPUT or BACKUP file from MesoNH simulation
    OUTPUT
        - a xarray dataset with interpolated variables
    
    Note :
        See end of MesoNH user doc for a schematic of the code grid
        Grid X -> Grid 1

    """
    ds = dsO
    ds['UT'] = dsO.UT.interp({'ni_u':dsO.ni})				            # Grid : 2
    ds['VT'] = dsO.VT.interp({'nj_v':dsO.nj})				            # Grid : 3
    ds['WT'] = dsO.WT.interp({'level_w':dsO.level})				        # Grid : 4
    ds['UW_HFLX'] = dsO.UW_HFLX.interp({'level_w':dsO.level,'ni_u':dsO.ni}) 	# Grid : 6
    ds['UW_VFLX'] = dsO.UW_VFLX.interp({'level_w':dsO.level})				    # Grid : 4
    ds['VW_HFLX'] = dsO.VW_HFLX.interp({'level_w':dsO.level,'nj_v':dsO.nj}) 	# Grid : 7
    ds['VW_VFLX'] = dsO.VW_VFLX.interp({'level_w':dsO.level})				    # Grid : 4
    ds['THW_FLX'] = dsO.THW_FLX.interp({'level_w':dsO.level})					# Grid : 4
    ds['RCONSW_FLX'] = dsO.RCONSW_FLX.interp({'level_w':dsO.level})			    # Grid : 4
    ds['UV_FLX'] = dsO.UV_FLX.interp({'ni_u':dsO.ni,'nj_v':dsO.nj})             # Grid : 5
    # no need to interp on the following dimensions, 
    #  they have a different name but are the same as mass point values
    ds['UT'] = ds.UT.rename(new_name_or_name_dict={'nj_u':'nj'})
    ds['VT'] = ds.VT.rename(new_name_or_name_dict={'ni_v':'ni'})
    ds['UW_HFLX'] = ds.UW_HFLX.rename(new_name_or_name_dict={'nj_u':'nj'})
    ds['VW_HFLX'] = ds.VW_HFLX.rename(new_name_or_name_dict={'ni_v':'ni'})

    # dropping dimensions not needed anymore
    ds = ds.drop_dims(['level_w','ni_u','nj_u','ni_v','nj_v'])

    return ds.unify_chunks() # this is needed otherwise chunked are messed up

def build_mean(dsO,function_avg=lambda a : a.mean(dim=['ni','nj'])):
    """
    this function builds mean profiles. By default, average is over the full domain in time,X and Y.

    INPUT:
        - dsO : xarray dataset = OUTPUT files from MNH simulation
        - function_avg : user made function for Reynolds averaging.
            defaut is X and Y average
    OUTPUT:
        - mean values for :
        Um,Vm,Wm,THTm,RVTm,Pm,SV1m,SV3m,SV4m,THTvm,ETm,TKEm,Tsfx,rho_sfx,Q,Qv,E0,TAU,u_star
    """

    X = dsO.ni
    Y = dsO.nj
    Z = dsO.level
    Time = dsO.time

    # removing halo
    U = dsO.UT				
    V = dsO.VT				
    W = dsO.WT 	
    P = dsO.PABST
    THT = dsO.THT
    RVT = dsO.RVT		
    TKET = dsO.TKET 			
    UW_HFLX = dsO.UW_HFLX 
    UW_VFLX = dsO.UW_VFLX				
    VW_HFLX = dsO.VW_HFLX 	
    VW_VFLX = dsO.VW_VFLX				
    THW_FLX = dsO.THW_FLX					
    RCONSW_FLX = dsO.RCONSW_FLX			

    UW_FLX = UW_HFLX + UW_VFLX
    VW_FLX = VW_HFLX + VW_VFLX

    # averaging
    Um = function_avg(U)
    Vm = function_avg(V)
    Wm = function_avg(W)
    THTm = function_avg(THT)
    RVTm = function_avg(RVT)
    THTvm = Compute_THTV(THTm,RVTm)
    SV1m,SV3m,SV4m = np.zeros(RVTm.shape),np.zeros(RVTm.shape),np.zeros(RVTm.shape)
    Pm = function_avg(P)
    if 'SVT001' in dsO.keys():
        SV1m = function_avg( dsO.SVT001 )
    #SV2m = function_avg( dsO.SVT002 )
    if 'SVT003' in dsO.keys():
        SV3m = function_avg( dsO.SVT003 )
    if 'SVT004' in dsO.keys():
        SV4m = function_avg( dsO.SVT004 )

    # avering subgrid value for surface fluxes    
    UW_FLXm = function_avg( UW_FLX ) 		# subgrid
    VW_FLXm = function_avg( VW_FLX ) 		# subgrid
    THW_FLXm = function_avg( THW_FLX ) 		# subgrid
    RCONSW_FLXm = function_avg( RCONSW_FLX )# subgrid
    THvW_FLXm = THW_FLXm*THTvm/THTm + 0.61*THTm*RCONSW_FLXm # subgrid
    
    # Turbulent kinetic nrj
    u_fluc = (U - Um)
    v_fluc = (V - Vm)
    w_fluc = (W - Wm)
    Em = function_avg( 0.5*( u_fluc**2 + v_fluc**2 + w_fluc**2 ) )
    TKEm = function_avg( TKET )
    ETm = Em + TKEm

    # surface values
    Tsfx = Theta_to_T(THTm[0],Pm[0])
    rho_sfx = Pm[0]/(Rd*Tsfx)
    Q = THW_FLXm[0] 		# at ground level
    Qv = THvW_FLXm[0] 	# at ground level
    E0 = RCONSW_FLXm[0]	# at ground level
    TAU = rho_sfx*np.sqrt( UW_FLXm[0]**2 + VW_FLXm[0]**2 )
    u_star = np.sqrt( TAU / rho_sfx )

    return Um,Vm,Wm,THTm,RVTm,Pm,SV1m,SV3m,SV4m,THTvm,ETm,TKEm,Tsfx,rho_sfx,Q,Qv,E0,TAU,u_star
     
def build_mean_flux(dsO,function_avg=lambda a : a.mean(dim=['ni','nj'])):
    """
    this function builds mean profiles for fluxes. By default, average is over the full domain in time,X and Y.

    INPUT:
        - dsO : xarray dataset = OUTPUT files from MNH simulation
        - function_avg : user made function for Reynolds averaging.
            defaut is X and Y average
    OUTPUT:
        - mean values for (total and subgrid):
        total : uv_r,uw_r,vw_r,uu_r,vv_r,ww_r,wtht_r,wrv_r,wthtv_r,ETm
        subgrid : uv_s,uw_s,vw_s,uu_s,vv_s,ww_s,wtht_s,wrv_s,wthtv_s,TKEm
    """
    X = dsO.ni
    Y = dsO.nj
    Z = dsO.level
    Time = dsO.time

    # removing halo
    U = dsO.UT				
    V = dsO.VT				
    W = dsO.WT 
    P = dsO.PABST
    THT = dsO.THT
    RVT = dsO.RVT	
    THTv = Compute_THTV(THT,RVT)
    TKET = dsO.TKET 			
    UW_HFLX = dsO.UW_HFLX 
    UW_VFLX = dsO.UW_VFLX				
    VW_HFLX = dsO.VW_HFLX 	
    VW_VFLX = dsO.VW_VFLX				
    THW_FLX = dsO.THW_FLX					
    RCONSW_FLX = dsO.RCONSW_FLX	
    UV_FLX = dsO.UV_FLX
    U_VAR = dsO.U_VAR
    V_VAR = dsO.V_VAR
    W_VAR = dsO.W_VAR
    UW_FLX = UW_HFLX + UW_VFLX
    VW_FLX = VW_HFLX + VW_VFLX
    THvW_FLX = THW_FLX*THTv/THT + 0.61*THT*RCONSW_FLX 

    # averaging
    Um = function_avg(U)
    Vm = function_avg(V)
    Wm = function_avg(W)
    RVm = function_avg(RVT)
    THTm = function_avg(THT)
    RVTm = function_avg(RVT)
    THTvm = function_avg(THTv)
    Pm = np.zeros(RVTm.shape)
    SV1m,SV3m,SV4m = np.zeros(RVTm.shape),np.zeros(RVTm.shape),np.zeros(RVTm.shape)
    Pm = function_avg(P)
    if 'SVT001' in dsO.keys():
        SV1m = function_avg( dsO.SVT001 )
    #SV2m = function_avg( dsO.SVT002 )
    if 'SVT003' in dsO.keys():
        SV3m = function_avg( dsO.SVT003 )
    if 'SVT004' in dsO.keys():
        SV4m = function_avg( dsO.SVT004 )

    # resolved fluctuations
    u_fluc = (U - Um).compute()
    v_fluc = (V - Vm).compute()
    w_fluc = (W - Wm).compute()
    tht_fluc = (THT - THTm).compute()
    rv_fluc = (RVT - RVm).compute()
    thtv_fluc = (THTv - THTvm).compute()

    uv_r,uv_s = function_avg( u_fluc*v_fluc ),function_avg( UV_FLX )
    uw_r,uw_s = function_avg( u_fluc*w_fluc ),function_avg( UW_FLX )
    vw_r,vw_s = function_avg( v_fluc*w_fluc ),function_avg( VW_FLX )
    uu_r,uu_s = function_avg( u_fluc*u_fluc ),function_avg( U_VAR )
    vv_r,vv_s = function_avg( v_fluc*v_fluc ),function_avg( V_VAR )
    ww_r,ww_s = function_avg( w_fluc*w_fluc ),function_avg( W_VAR )
    wtht_r,wtht_s = function_avg( tht_fluc*w_fluc ),function_avg( THW_FLX )
    wthtv_r,wthtv_s = function_avg( thtv_fluc*w_fluc ),function_avg( THvW_FLX )
    wrv_r,wrv_s = function_avg( w_fluc*rv_fluc ),function_avg( RCONSW_FLX)

    return uv_r,uw_r,vw_r,uu_r,vv_r,ww_r,wtht_r,wrv_r,wthtv_r,uv_s,uw_s,vw_s,uu_s,vv_s,ww_s,wtht_s,wrv_s,wthtv_s

def build_mean_file(path_in,path_out,path_out_flx,nhalo,dsO_i,d_boxes):
    """
    GOAL:
        This function computes mean profiles on user defined boxes.
    INPUT:
        - path_in : string path of the OUTPUT files (saved in the output netcdf)
        - path_out : where to write the mean file
        - path_out_flx : where to write the mean flux file
        - nhalo : MNH halo
        - dsO_i : xarray dataset that contains OUTPUT files data from MNH sim, 
                    interpolated at mass points
        - d_boxes : dict with position of boxes (origin,Lx,Ly)
    OUTPUT:
        - a netcdf file with 1D profiles (along Z) with mean values inside boxes

    Note : to change the mean operator, you need to change
        Reynolds_avg,coordinates of written variables
    """
    print('* Computing mean profiles of prognostic variables and fluxes')
    Z = dsO_i.level
    N = len(d_boxes['boxes']) # number of boxes
    Lx = d_boxes['Lx']
    Ly = d_boxes['Ly']

    Bmeanfile = pathlib.Path(path_out).is_file()
    Bmeanfluxfile = pathlib.Path(path_out_flx).is_file()

    if Bmeanfile:
        print('     - File is already here : '+path_out) 
    if Bmeanfluxfile:
        print('     - File is already here : '+path_out_flx) 

    if not Bmeanfile or not Bmeanfluxfile: 
        # one of the asked file is not here 
        # We interp at mass point before any compute

        ##################
        ##   PROG VAR   ##
        ##################

        if not Bmeanfile: # we need to build the mean profile 
            print('build file for prognostic var')
            Um,Vm,Wm = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            THTm,RVTm,Pm = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            THTvm,ETm,TKEm = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            SV1m,SV3m,SV4m = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            Tsfx,rho_sfx,TAU,u_star = np.zeros(N),np.zeros(N),np.zeros(N),np.zeros(N)
            H0,HL0,L0 = np.zeros(N),np.zeros(N),np.zeros(N)
            

            for k,boxe in enumerate(d_boxes['boxes'].keys()):
                O = d_boxes['boxes'][boxe]['O']
                
                dsINPUT = dsO_i.sel(ni=slice(O[0],O[0]+Lx),
                                nj=slice(O[1],O[1]+Ly),)
                # lazy
                U,V,W,THT,RVT,P,SV1,SV3,SV4,THTv,ET,TKE,Ts,rho_s,Q0,Qv0,E0,tau,u_st = build_mean(dsINPUT,Reynolds_avg)
                # compute dask array
                print('     -> computing boxe n° '+boxe+' ...')
                Um[k,:],Vm[k,:],Wm[k,:] = U.values,V.values,W.values
                THTm[k,:],RVTm[k,:],Pm[k,:] = THT.values,RVT.values,P.values
                THTvm[k,:],ETm[k,:],TKEm[k,:] = THTv.values,ET.values,TKE.values
                SV1m[k,:],SV3m[k,:],SV4m[k,:] = SV1.values,SV3.values,SV4.values
                Tsfx[k],rho_sfx[k],TAU[k],u_star[k] = Ts.values,rho_s.values,tau.values,u_st.values
                H0[k],HL0[k],L0[k] = Q0.values,Qv0.values,E0.values

            L_dim = ['nboxe','level']
            L_dim0D = ['nboxe']
            coords={'nboxe':np.arange(1,N+1),'level': Z}
            data_vars = {'Um':(L_dim,Um,{'long_name':'Mean zonal wind',
                                'units':'m s-1',
                                'grid location':'mass_center'}),
                    'Vm':(L_dim,Vm,{'long_name':'Mean meridional wind',
                                'units':'m s-1',
                                'grid location':'mass_center'}),
                    'Wm':(L_dim,Wm,{'long_name':'Mean vertical wind',
                                'units':'m s-1',
                                'grid location':'mass_center'}),
                    'THTm':(L_dim,THTm,{'long_name':'Mean potential temperature',
                                'units':'K',
                                'grid location':'mass_center'}),
                    'RVTm':(L_dim,RVTm,{'long_name':'Mean mixing ratio',
                                'units':'kg/kg',
                                'grid location':'mass_center'}),
                    'THTvm':(L_dim,THTvm,{'long_name':'Mean virtual potential temperature',
                                'units':'K',
                                'grid location':'mass_center'}),
                    'Pm':(L_dim,Pm,{'long_name':'Mean pressure',
                                'units':'Pa',
                                'grid location':'mass_center'}),			
                    'ETm':(L_dim,ETm,{'long_name':'Total mean turbulent kinetic energy',
                                'units':'m2 s-1',
                                'grid location':'mass_center'}),
                    'TKEm':(L_dim,TKEm,{'long_name':'Subgrid mean turbulent kinetic energy',
                                'units':'m2 s-1',
                                'grid location':'mass_center'}),
                    'SV1m':(L_dim,SV1m,{'long_name':'mean tracer1 concentration',
                                'units':'kg kg-1',
                                'grid location':'mass_center'}),
                    'SV3m':(L_dim,SV3m,{'long_name':'mean tracer3 concentration',
                                'units':'kg kg-1',
                                'grid location':'mass_center'}),
                    'SV4m':(L_dim,SV4m,{'long_name':'mean tracer4 concentration',
                                'units':'kg kg-1',
                                'grid location':'mass_center'}),
                    'Q_star':(L_dim0D,H0,{'long_name':'Sensible heat flux',
                                'units':'m K s-1',
                                'grid location':'mass_center'}),
                    'Qv_star':(L_dim0D,HL0,{'long_name':'Buoyancy flux',
                                'units':'m K s-1',
                                'grid location':'mass_center'}),
                    'E0':(L_dim0D,L0,{'long_name':'Latent heat flux',
                                'units':'kg m kg-1 s-1',
                                'grid location':'mass_center'}),
                    'Tau':(L_dim0D,TAU,{'long_name':'Tau',
                                'units':'m-1 s-2',
                                'grid location':'mass_center'}),
                    'u_star':(L_dim0D,u_star,{'long_name':'Friction velocity',
                                'units':'m s-1',
                                'grid location':'mass_center'}),
                    'Tsfx':(L_dim0D,Tsfx,{'long_name':'Surface temperature',
                                'units':'K',
                                'grid location':'mass_center'}),
                    'rho_sfx':(L_dim0D,rho_sfx,{'long_name':'Surface density',
                                'units':'K',
                                'grid location':'mass_center'}),
                    }

            
            ds_mean = xr.Dataset(data_vars=data_vars,coords=coords,
                        attrs={'input_files':path_in,'file_created_by':'build_mean_file'})
            
            print(ds_mean)
            ds_mean.to_netcdf(path=path_out,mode='w')
            ds_mean.close()
            
        ################
        ##   FLUXES   ##
        ################

        if not Bmeanfluxfile: # we need to build the mean flux profile 
            print(' building file for fluxes')
            UV_r,UW_r,VW_r = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            UU_r,VV_r,WW_r = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            WTHT_r,WRV_r,WTHTV_r = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            UV_s,UW_s,VW_s = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            UU_s,VV_s,WW_s = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))
            WTHT_s,WRV_s,WTHTV_s = np.zeros((N,len(Z))),np.zeros((N,len(Z))),np.zeros((N,len(Z)))

            for k,boxe in enumerate(d_boxes['boxes'].keys()):
                O = d_boxes['boxes'][boxe]['O']
                dsINPUT = dsO_i.sel(ni=slice(O[0],O[0]+Lx),
                                nj=slice(O[1],O[1]+Ly),)
                # lazy
                uv_r,uw_r,vw_r,uu_r,vv_r,ww_r,wtht_r,wrv_r,wthtv_r,uv_s,uw_s,vw_s,uu_s,vv_s,ww_s,wtht_s,wrv_s,wthtv_s = build_mean_flux(dsINPUT,Reynolds_avg)
                # compute dask array
                print('     -> computing boxe n° '+boxe+' ...')
                UV_r[k,:],UW_r[k,:],VW_r[k,:] = uv_r.values,uw_r.values,vw_r.values
                UU_r[k,:],VV_r[k,:],WW_r[k,:] = uu_r.values,vv_r.values,ww_r.values
                WTHT_r[k,:],WRV_r[k,:],WTHTV_r[k,:] = wtht_r.values,wrv_r.values,wthtv_r.values
                UV_s[k,:],UW_s[k,:],VW_s[k,:] = uv_s.values,uw_s.values,vw_s.values
                UU_s[k,:],VV_s[k,:],WW_s[k,:] = uu_s.values,vv_s.values,ww_s.values
                WTHT_s[k,:],WRV_s[k,:],WTHTV_s[k,:] = wtht_s.values,wrv_s.values,wthtv_s.values

            L_dim = ['nboxe','level']
            coords={'nboxe':np.arange(1,N+1),'level': Z}
            data_vars = {'FLX_UV':(L_dim,UV_r + UV_s,{'long_name':'Horizontal wind turbulent flux (total)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_UW':(L_dim,UW_r + UW_s,{'long_name':'Turbulent vertical flux in x direction (total)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_VW':(L_dim,VW_r + VW_s,{'long_name':'Turbulent vertical flux in y direction (total)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_UU':(L_dim,UU_r + UU_s,{'long_name':"u'2 covariance (total)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_VV':(L_dim,VV_r + VV_s,{'long_name':"v'2 covariance (total)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_WW':(L_dim,WW_r + WW_s,{'long_name':"w'2 covariance (total)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_THW':(L_dim,WTHT_r + WTHT_s,{'long_name':'Turbulent vertical flux of heat (dry) (total)',
                                    'units':'K m s-1',
                                    'grid location':'mass_center'}),
                        'FLX_THvW':(L_dim,WTHTV_r + WTHTV_s,{'long_name':'Turbulent vertical flux of heat (moist) (total)',
                                    'units':'K m s-1',
                                    'grid location':'mass_center'}),
                        'FLX_RvW':(L_dim,WRV_r + WRV_s,{'long_name':'Turbulent vertical flux of vapor mixing ratio (total)',
                                    'units':'kg kg-1 m s-1',
                                    'grid location':'mass_center'}),
                        'FLX_UV_s':(L_dim,UV_s,{'long_name':'Horizontal wind turbulent flux (sgs)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_UW_s':(L_dim,UW_s,{'long_name':'Turbulent vertical flux in x direction (sgs)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_VW_s':(L_dim,VW_s,{'long_name':'Turbulent vertical flux in y direction (sgs)',
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_UU_s':(L_dim,UU_s,{'long_name':"u'2 covariance (sgs)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_VV_s':(L_dim,VV_s,{'long_name':"v'2 covariance (sgs)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_WW_s':(L_dim,WW_s,{'long_name':"w'2 covariance (sgs)",
                                    'units':'m2 s-2',
                                    'grid location':'mass_center'}),
                        'FLX_THW_s':(L_dim,WTHT_s,{'long_name':'Turbulent vertical flux of heat (dry) (sgs)',
                                    'units':'K m s-1',
                                    'grid location':'mass_center'}),
                        'FLX_THvW_s':(L_dim,WTHTV_s,{'long_name':'Turbulent vertical flux of heat (moist) (sgs)',
                                    'units':'K m s-1',
                                    'grid location':'mass_center'}),
                        'FLX_RvW_s':(L_dim,WRV_s,{'long_name':'Turbulent vertical flux of vapor mixing ratio (sgs)',
                                    'units':'kg kg-1 m s-1',
                                    'grid location':'mass_center'}),	
                        }

            
            ds_meanflx = xr.Dataset(data_vars=data_vars,coords=coords,
                        attrs={'input_files':path_in,'file_created_by':'build_mean_file'})
            
            print(ds_meanflx)
            ds_meanflx.to_netcdf(path=path_out_flx,mode='w')
            ds_meanflx.close() 
        else:
            print('File is already here : '+path_out_flx) 
        
        dsO_i.close()      

        




        
