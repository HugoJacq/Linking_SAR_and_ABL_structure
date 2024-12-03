import numpy as np
import pathlib
import xarray as xr
from mod_build_mean import mass_interpolator
from module_tools import *
"""
To be used with analyse.py
"""

def Integ_min(name,dsIN,dsmean):
    """
    This function integrate the standard deviation of a tracer following eq (3) of 
        Brient, F., Couvreux, F., Rio, C., & Honnert, R. (2024). 
        Coherent subsiding structures in large eddy simulations of atmospheric boundary layers. 
        Quarterly Journal of the Royal Meteorological Society, 1–23. https://doi.org/10.1002/qj.4625

    METHOD : 
        if the tracer is emitted at surface :
            we integrate from surface to model top
        else if the tracer is emitted at a specific layer : 
            we integrate from +2 levels above the max of the mean concentration 
                up to the surface.

    INPUT:
        - name : str, name of the tracer (SV1,SV3,SV4)
        - dsIN : xarray dataset from MesoNH output (interpolated at mass point)
        - dsmean : xarray dataset with mean profiles
    OUTPUT:
        - 1D (along Z) integrated standard deviation for each passive tracer


    Note : 
        In the case where name='SV3', everything above SV3.max('level) + Top_ABL_offset is set to 999 to filter out non turbulent motion.
            This is done because else in the next step (condition 1) too much air is kept, altough it is very weakly turbulent.
            As a consquence, the value of stdmin above the boundary layer is gamma*999.
            For the case where name='SV1' or 'SV4', this is not necessary as the value of stdmin is enough to prevent this behavior.
            
    """
    L_correct_s = ['SV1','SV3','SV4']
    Z = dsIN.level
    Top_ABL_offset = 10
    dim_std = ['time','ni','nj']

    indZbot = 0
    indZtop = len(Z)
    if name == 'SV1':
        s = dsIN.SVT001
    elif name =='SV3':
        s = dsIN.SVT003
        indZtop = dsmean[name+'m'].argmax('level').values + Top_ABL_offset
    elif name =='SV4':
        s = dsIN.SVT004
    else:
        raise Exception('The chosen tracer '+name+' is not recognized, should be one of '+str(L_correct_s))

    # zone de test
    # le but ici est de voir si on mettant stdmin pareil pour SV1 et SV4 on réussi à
    # s'affranchir des problèmes qui apparaissent dans les domaine où sv1 et sv4 sont tous les 2 présents.
    if name=='SV1' or name=='SV4':
        s = xr.where( dsIN.SVT001 > dsIN.SVT004, dsIN.SVT001,dsIN.SVT004)
    #

    std = s.std(dim_std)
    integ = xr.ones_like(std) * 999
    for indz in range(indZbot,indZtop):
        if name =='SV3':
            values =  - std.isel(level=slice(indz,indZtop))
            norme = Z[indz] - Z[indZtop]
        else:
            values = std.isel(level=slice(0,indz))
            norme = Z[indz]
        integ[dict(level=indz)] = values.integrate('level') / norme

    return integ

def condition1(name,dsIN,dsmean,gamma,mCS):
    """
    This function returns a mask where, for each passive tracer, the following condition is met:

                s' > m.max(std,std_min)

    METHOD : see eq (1),(2) and (3) from 
        Brient, F., Couvreux, F., Rio, C., & Honnert, R. (2024). 
        Coherent subsiding structures in large eddy simulations of atmospheric boundary layers. 
        Quarterly Journal of the Royal Meteorological Society, 1–23. https://doi.org/10.1002/qj.4625
    
    INPUT:
        - name : str, name of the tracer (SV1,SV3,SV4)
        - dsIN : xarray dataset from MesoNH output (interpolated at mass point)
        - dsmean : xarray dataset with mean profiles
        - gamma : tunable parameter 
        - mCS : strength of the turbulent filtering
    OUTPUT:
        - turbulent mask (=1) that filter non turbulent motion (=0)

    """
    L_correct_s = ['SV1','SV3','SV4']
    dim_std = ['time','ni','nj']

    if name == 'SV1':
        name2 = 'SVT001'
    elif name =='SV3':
        name2 = 'SVT003'
    elif name =='SV4':
        name2 = 'SVT004'
    else:
        raise Exception('The chosen tracer '+name+' is not recognized, should be one of '+str(L_correct_s))
    
    s = dsIN[name2]
    mean = dsmean[name+'m']
    s_fluc = s-mean
    std_min = gamma * Integ_min(name,dsIN,dsmean)
    std = s.std(dim_std)

    # figure that plots gamma*integ
    fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
    ax.plot(std_min,dsIN.level,c='k',ls='--',label='stdmin')
    ax.plot(std,dsIN.level,c='k',label='std')
    ax.set_xlabel(name)
    ax.set_ylabel('Z')
    ax.legend()
    
   
    max_cond = mCS * xr.where(std > std_min,std,std_min) # eq (2) of ref
    cond1 = xr.where( s_fluc > max_cond, 1,0)                        # eq (1) of ref

    # print('std.shape,std_min.shape',std.shape,std_min.shape)
    # print('s_fluc.shape',s_fluc.shape)
    # print('max_cond.shape,cond1.shape',max_cond.shape,cond1.shape)
    return cond1

def condition2(name,dsIN,dsmean,condition1_mask):
    """
    This function returns a numbered mask that identifies coherent structures for a specified tracer (name).
    Two categories are considered : upward moving tracer and downward moving tracers.

    METHOD :  see table (2) from 
        Jacquet, H., Ayet, A., & Couvreux, F. (2025). 
        Atmosphere Response to an Oceanic Sub-mesoscale SST Front: A Coherent Structure Analysis [preprint]. 
        Journal of Geophysical Research: Atmospheres. https://hal.science/hal-04729066

    INPUT:
        - name : str, name of the tracer (SV1,SV3,SV4)
        - dsIN : xarray dataset from MesoNH output (interpolated at mass point) 
        - dsmean : xarray dataset with mean profiles
        - mask_condition_1 : 3D mask that filter out non turbulent motions
    OUTPUT:
        - masks for each categories (same dimensions as UT, 4D in most cases) 
            1 for the upward moving structure
            2 for the downward moving structure 
            0 else

    """
    L_correct_s = ['SV1','SV3','SV4']
    if name == 'SV1':
        name2 = 'SVT001'
    elif name =='SV3':
        name2 = 'SVT003'
    elif name =='SV4':
        name2 = 'SVT004'
    else:
        raise Exception('The chosen tracer '+name+' is not recognized, should be one of '+str(L_correct_s))
    
    sv_fluc = dsIN[name2] - dsmean[name+'m']
    w_fluc =  dsIN['WT'] - dsmean['Wm']

    mask_up = np.logical_and(np.logical_and(condition1_mask,sv_fluc>0),w_fluc > 0)
    mask_down = np.logical_and(np.logical_and(condition1_mask,sv_fluc>0),w_fluc < 0)

    Z,Y,X = dsIN.level,dsIN.nj,dsIN.nj
    mask_obj = xr.where( mask_up,1,0 )              # upward moving structure
    mask_obj = xr.where( mask_down,2,mask_obj )    # downward moving structure
    return mask_obj

def CS_for_1_tracer(name,dsIN,dsmean,gamma,mCS):
    """
    This function gather condition 1 and condition 2 of the conditional sampling operation.
    This has to be run for each tracers.
    This function is just a wrapper !

    INPUT:
        - name : str, name of the tracer (SV1,SV3,SV4)
        - dsIN : xarray dataset from MesoNH output (interpolated at mass point) 
        - dsmean : xarray dataset with mean profiles
        - gamma : tunable parameter 
        - mCS : strength of the turbulent filtering
    OUTPUT:
        - a xarray DataArray with a 4D (time+space) mask
    """
    C1 = condition1(name,dsIN,dsmean,gamma,mCS)
    C2 = condition2(name,dsIN,dsmean,C1)
    print('C1.shape,C2.shape',C1.shape,C2.shape)
    return C2

def get_global_mask(dsIN,mask_condition_2):
    """
    This function returns a mask with a number for each coherent structures.
    
    METHOD :
        This step is needed as sometimes two categories overlaps.

    INPUT:
        - dsIN : xarray dataset from MesoNH output (interpolated at mass point) 
        - mask_condition_2 : 3D mask that filter coherent structures
    OUTPUT:
        - a mask with numbered categories
            1 = updrafts from SV1
            2 = sub. shells from SV1
            3 = updrafts from SV4
            4 = sub. shells from SV4
            5 = downdrafts


    Note: 
        - At some cell, both surface tracer can be present and so we keep the most concentrated one
        - In another case, it can happen that both surface emitted tracer AND top-ABL emitted tracer are present.
            In such case, we chose to keep the top-ABL tracer information.
        - Keep in mind that surface emitted tracer concentration CANNOT be compared to top-ABL tracer concentration as they are not emitted with the same method.

    """
    SV1 = dsIN.SVT001.data
    SV3 = dsIN.SVT003.data
    SV4 = dsIN.SVT004.data

    #bothSVup = xr.where( np.logical_and(mask_condition_2[0]==1,mask_condition_2[1]==1) )
    #bothSVss = xr.where( np.logical_and(mask_condition_2[0]==2,mask_condition_2[1]==2) )

    # global_mask = xr.where( mask_condition_2[0]==1 ,1,0)            # up from SV1
    # global_mask = xr.where( mask_condition_2[0]==2 ,2,global_mask)  # ss from SV1
    # global_mask = xr.where( mask_condition_2[1]==1 ,3,global_mask)  # up from SV4
    # global_mask = xr.where( mask_condition_2[1]==2 ,4,global_mask)  # ss from SV4
    global_mask = xr.where( mask_condition_2[0]==1 ,1,0)            # up from SV1
    global_mask = xr.where( mask_condition_2[0]==2 ,2,global_mask)  # ss from SV1
    global_mask = xr.where( np.logical_and(mask_condition_2[1]==1,SV1<SV4) ,3,global_mask)  # up from SV4
    global_mask = xr.where( np.logical_and(mask_condition_2[1]==2,SV1<SV4) ,4,global_mask)  # ss from SV4
    global_mask = xr.where( mask_condition_2[2]==2,5,global_mask)  # down from SV3

    return global_mask

def build_CS(path_in,path_out,dsO_i,dsmean,d_boxes,gamma=0.005,mCS=1.0):
    """
    This function builds a netcdf file with coherent structures identified.

    METHOD : following 
        Jacquet, H., Ayet, A., & Couvreux, F. (2025). 
        Atmosphere Response to an Oceanic Sub-mesoscale SST Front: A Coherent Structure Analysis [preprint]. 
        Journal of Geophysical Research: Atmospheres. https://hal.science/hal-04729066

    INPUT:
        - path_in : string path of the OUTPUT files (saved in the output netcdf)
        - path_out : where to write the mean file
        - nhalo : MNH halo
        - dsO_i : xarray dataset that contains OUTPUT files data from MNH sim, interpolated at mass point
        - dsmean : xarray dataset with mean profiles
        - d_boxes : dict with position of boxes (origin,Lx,Ly)
        - gamma : parameter that can be tuned to filter non turbulent cells
        - mCS : scaling factor that quantifies the strenght of the conditional sampling
    OUTPUT:
        - a netcdf with 3D volume, for each boxes.
            Variables are U,V,W,RV,THT,THTV,SV1,SV3,SV4

    Note : beware that the variable SV4 is not in MNH files if you have not modified the code.
            please look at https://github.com/HugoJacq/Linking_SAR_and_ABL_structure for more information.
    """
    print('* Computing conditional sampling file')
    Z = dsO_i.level
    Y = dsO_i.nj
    X = dsO_i.ni
    N = len(d_boxes['boxes']) # number of boxes
    Isfile = pathlib.Path(path_out).is_file()
    Nsv = 3 # number of different passive tracers, SV1 SV3 SV4

    if Isfile:
       print('     - File is already here : '+path_out) 
    else:
        print('i need to work ...')
        # lets get the chunks from dsO_i
        chunks_i = {}
        for dim in dsO_i.chunks.keys():
            chunks_i[dim] = dsO_i.chunks[dim][0]
        
        Um,Vm,Wm = dsmean.Um,dsmean.Vm,dsmean.Wm
        THTm,RVTm,Pm = dsmean.THTm,dsmean.RVTm,dsmean.Pm
        THTvm,ETm,TKEm = dsmean.THTvm,dsmean.ETm,dsmean.TKEm
        
        # initialise arrays
        Lx = d_boxes['Lx']
        Ly = d_boxes['Ly']
        res = int((dsO_i.ni[2] -dsO_i.ni[1]).values)

        coords1={'nboxe':np.arange(1,N+1),
                    'time':dsO_i.time,
                    'level': Z,
                    'coord_x':np.arange(0,Lx,res),
                    'coord_y':np.arange(0,Ly,res)}
        coords2={'nboxe':np.arange(1,N+1),
                    'Nsv':np.arange(1,Nsv+1),
                    'time':dsO_i.time,
                    'level': Z,
                    'coord_x':np.arange(0,Lx,res),
                    'coord_y':np.arange(0,Ly,res)}
        chunks = {'level':chunks_i['level'],'coord_x':chunks_i['ni'],'coord_y':chunks_i['nj']}


        object_mask = xr.DataArray( np.zeros((N,Nsv,len(dsO_i.time),len(Z),Ly//res,Lx//res),dtype=np.int16),
                                 dims=('nboxe','Nsv','time','level','coord_y','coord_x'),
                                 coords=coords2).chunk(chunks)
        global_mask = xr.DataArray( np.zeros((N,len(dsO_i.time),len(Z),Ly//res,Lx//res),dtype=np.int16),
                                 dims=('nboxe','time','level','coord_y','coord_x'),
                                 coords=coords1).chunk(chunks)
        X = xr.DataArray( np.zeros((N,Lx//res),dtype=np.float32),
                         dims=('nboxe','coord_x'),
                         coords={'coord_x':np.arange(0,Lx,res)})
        Y = xr.DataArray( np.zeros((N,Lx//res),dtype=np.float32),
                         dims=('nboxe','coord_y'),
                         coords={'coord_y':np.arange(0,Lx,res)})
        
        THT = xr.DataArray( np.zeros((N,len(dsO_i.time),len(Z),Ly//res,Lx//res),dtype=np.float32),
                                 dims=('nboxe','time','level','coord_y','coord_x'),
                                 coords=coords1).chunk(chunks)
        THTV,RVT,TKE = xr.zeros_like(THT),xr.zeros_like(THT),xr.zeros_like(THT)
        U,V,W = xr.zeros_like(THT),xr.zeros_like(THT),xr.zeros_like(THT)
        SV1,SV3,SV4 =  xr.zeros_like(THT),xr.zeros_like(THT),xr.zeros_like(THT)
        ET =  xr.zeros_like(THT)

        for k,boxe in enumerate(d_boxes['boxes'].keys()): # 
            print('     -> computing boxe n° '+boxe+' ...')
            O = d_boxes['boxes'][boxe]['O']
            # selecting boxe area
            dsINPUT = dsO_i.sel(ni=slice(O[0],O[0]+Lx),
                            nj=slice(O[1],O[1]+Ly),)
            
            # compute fluctuations ?
            THT[k,:,:,:] = dsINPUT.THT.data
            RVT[k,:,:,:] = dsINPUT.RVT.data
            U[k,:,:,:] = dsINPUT.UT.data		
            V[k,:,:,:] = dsINPUT.VT.data			
            W[k,:,:,:] = dsINPUT.WT.data	
            TKE[k,:,:,:] = dsINPUT.TKET.data 	
            SV1[k,:,:,:,:] = dsINPUT.SVT001.data
            SV3[k,:,:,:,:] = dsINPUT.SVT003.data
            SV4[k,:,:,:,:] = dsINPUT.SVT004.data
            THTV[k,:,:,:] = Compute_THTV(THT[k,:,:,:],RVT[k,:,:,:])

            X[k,:] = dsINPUT.ni.data
            Y[k,:] = dsINPUT.nj.data

            for j,name in enumerate(['SV1','SV4','SV3']):
                print('working on '+name)
                object_mask[k,j,:,:,:,:] = CS_for_1_tracer(name,dsINPUT,dsmean.isel(nboxe=k),gamma,mCS).data
            global_mask[k,:,:,:,:] = get_global_mask(dsINPUT,object_mask[k,:,:,:,:,:])
            
            # TEST FOR GLOBAL_MASK
            # object_mask[k,0,:,:,:,0:10] = 1
            # object_mask[k,0,:,:,:,12:20] = 2
            # object_mask[k,1,:,:,:,5:15] = 1
            # object_mask[k,1,:,:,:,32:40] = 2
            # object_mask[k,2,:,:,:,5:10] = 1
            # object_mask[k,2,:,:,:,12:15] = 2
            # global_mask[k,:,:,:,:] = get_global_mask(dsINPUT,object_mask[k,:,:,:,:,:])
            # fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
            # x = np.arange(0,len(X[k,:]))
            # ax.scatter(x,xr.where( SV1[k,0,1,1,:]>SV4[k,0,1,1,:],1,0),c='c',label='SV1>SV4',marker='x')
            # ax.plot(x,xr.where(object_mask[k,0,1,1,1,:]==1,1,0),c='r',label='up1')
            # ax.plot(x,xr.where(object_mask[k,0,1,1,1,:]==2,1,0),c='purple',label='ss1')
            # ax.plot(x,xr.where(object_mask[k,1,1,1,1,:]==1,1,0),c='orange',label='up2')
            # ax.plot(x,xr.where(object_mask[k,1,1,1,1,:]==2,1,0),c='pink',label='ss2')
            # ax.plot(x,xr.where(object_mask[k,2,1,1,1,:]==2,1,0),c='green',label='down')
            # ax.plot(x,global_mask[k,1,1,1,:],c='k',ls='--')
            # ax.set_xlabel('X')
            # ax.set_ylabel('num of obj')
            # ax.legend()
            # plt.show()
            # raise Exception

        num_to_ob = np.array(['updraft_1','sub._shells_1','updrafts_2','sub._shells_2','downdrafts'])
        
        L_dim = ['nboxe','time','level','coord_x','coord_y']
        L_dim_2 = ['nboxe','Nsv','time','level','coord_x','coord_y']

        

        #print(object_mask.shape,num_to_ob.shape,X.shape)
        data_vars = {
                'U':(L_dim,U.data,{'long_name':'Mean zonal wind',
                            'units':'m s-1',
                            'grid location':'mass_center'}),
                'V':(L_dim,V.data,{'long_name':'Mean meridional wind',
                            'units':'m s-1',
                            'grid location':'mass_center'}),
                'W':(L_dim,W.data,{'long_name':'Mean vertical wind',
                            'units':'m s-1',
                            'grid location':'mass_center'}),
                'THT':(L_dim,THT.data,{'long_name':'Mean potential temperature',
                            'units':'K',
                            'grid location':'mass_center'}),
                'RVT':(L_dim,RVT.data,{'long_name':'Mean mixing ratio',
                            'units':'kg/kg',
                            'grid location':'mass_center'}),
                'THTV':(L_dim,THTV.data,{'long_name':'Mean virtual potential temperature',
                            'units':'K',
                            'grid location':'mass_center'}),			
                'SV1':(L_dim,SV1.data,{'long_name':'mean tracer1 concentration',
                            'units':'kg kg-1',
                            'grid location':'mass_center'}),
                'SV3':(L_dim,SV3.data,{'long_name':'mean tracer3 concentration',
                            'units':'kg kg-1',
                            'grid location':'mass_center'}),
                'SV4':(L_dim,SV4.data,{'long_name':'mean tracer4 concentration',
                            'units':'kg kg-1',
                            'grid location':'mass_center'}),
                'object_mask':(L_dim_2,object_mask.data,{'long_name':'Condition 2 : individual object mask',
                            'units':'',
                            'grid location':'mass_center'}),
                'global_mask':(L_dim,global_mask.data,{'long_name':'Global mask (no overlap)',
                            'units':'',
                            'grid location':'mass_center'}),
                'obj_name':(['nobj'],num_to_ob,{'long_name':'Name of the coherent structures',
                            'units':'',
                            'grid location':'mass_center'}),
                'X':(['nboxe','coord_x'],X.data,{'long_name':'ni dimension in the boxe',
                            'units':'m',
                            'grid location':'mass_center'}),
                'Y':(['nboxe','coord_y'],Y.data,{'long_name':'nj dimension in the boxe',
                            'units':'m',
                            'grid location':'mass_center'}),
                }
        # ajouter une variable qui fait le lien entre numéro objet et nom de l'objet
        

        encoding = {'U':{'dtype':'float32'},
                    'V':{'dtype':'float32'},
                    'W':{'dtype':'float32'},
                    'THT':{'dtype':'float32'},
                    'THTV':{'dtype':'float32'},
                    'RVT':{'dtype':'float32'},
                    'SV1':{'dtype':'float32'},
                    'SV3':{'dtype':'float32'},
                    'SV4':{'dtype':'float32'},
                    'object_mask':{'dtype':'int16'},
                    'global_mask':{'dtype':'int16'},
                    'obj_name':{'dtype':'str'},
                    'X':{'dtype':'float32'},
                    'Y':{'dtype':'float32'},
                    }

        ds_CS = xr.Dataset(data_vars=data_vars,coords=coords2,
                    attrs={'input_files':path_in,'file_created_by':'build_CS'})
        
        print(ds_CS)
        ds_CS.to_netcdf(path=path_out,mode='w',encoding = encoding)
        ds_CS.close()
        dsO_i.close()