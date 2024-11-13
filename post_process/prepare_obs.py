import xarray as xr
import numpy as np
from module_tools import *
import os
import matplotlib as mpl
import warnings


"""
This scripts explore SST data to chose a nice SST for the Agulhas current

-> BUILD_SAR_netcdf
    Transform the .npy from OVL to a clean netcdf
-> SAR_FIRST_LOOK
    Plot the swath from the clean netcdf of SAR data
-> SST_FIRST_LOOK
    Plot the L4 Regional SST Reanalysis from Ifremer, with user specified bounds.
    This is useful to check the position of a simulation in the LAT/LON space.
-> SST_SAR_FIRST_LOOK
    Combine SST/SAR data on one plot. LES domains are added with user specified dimensions
-> CHANGE_ENV
    Section used to get the input value to give a MNH given a idealized initial profile.



"""
# INPUT ZONE ================================================================
githubrepo = 'https://github.com/HugoJacq/Linking_SAR_and_ABL_structure'
pwd_path = '/home/jacqhugo/scripts/simu_nest_aggv3/'
png_save = 'PNGs_observations/'

BUILD_SAR_netcdf = False # Custom function to get netcdf from SAR .npy (OVL)
SAR_FIRST_LOOK = False # SAR (OVL)
SAR_VIGNETTES = False # SAR (OVL)
SST_FIRST_LOOK = False # SST from Ifremer first look
SST_SAR_FIRST_LOOK = False # combined SAR (OVL) and SST
CHANGE_ENV = False # small tool to setup different initial conditions
ERA5_FIRST_LOOK = False # a first look at ERA5 from Copernicus
COMPARE_SAR_SOURCES = True # compare OVL and Ifremer SAR images
dpi = 200
# END INPUT ZONE ============================================================

# we convert lat-lon in km
# https://www.lexilogos.com/calcul_distances.htm
# OVL extract bounds
#LEFT = 22.5     # °E
#RIGHT = 26      # °E
#BOTT = -37      # °N
#TOP = -34.5     # °N
DegLat = 111.2 # = 277.987 / (TOP-BOTT)# 1° lat in km
DegLon = 91.6 # = 320.719 / (RIGHT-LEFT)# 1° lon in km at ~ -35°N

os.system('mkdir -p '+png_save)

if BUILD_SAR_netcdf or SAR_FIRST_LOOK or SAR_VIGNETTES:
	"""
	SAR DATA
	See 'DataOrigin'
	"""
	base_path = pwd_path+ 'SAR_from_OVL/'
	path_file = {"SAR1":('3857_SAR_roughness-s1a-iw-grd-vv-20151210t170856-20151210t170921-008982-00cdef-001-C641-' +
		                    'lat_-37.0000_-34.5005_lon_22.5000_26.0005.npy'),
		         "SAR2":('3857_SAR_roughness-s1a-iw-grd-vv-20151210t170827-20151210t170856-008982-00cdef-001-9AE1-' +
		                    'lat_-37.0000_-34.5005_lon_22.5000_26.0005.npy'),
		         "SST":('3857_ODYSSEA_SAF_SST-20151211-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0-' +
		                    'lat_-37.0000_-34.5005_lon_22.5000_26.0005.npy')}

	# note : 33th to 64th character
	start_time = {'SAR1':'20151210t170856',
		          'SAR2':'20151210t170827',}
	start_time['SAR'] = start_time['SAR2']
	end_time = {'SAR1':'20151210t170921',
		          'SAR2':'20151210t170856',}
	end_time['SAR'] = end_time['SAR1']
	name_netcdf = 'SAR_roughness_' + start_time['SAR'] + '-' + end_time['SAR']+'.nc'

	if name_netcdf in os.listdir(base_path):
		FILE_is_HERE = True
	else:
		FILE_is_HERE = False
	if (SAR_FIRST_LOOK) and (not BUILD_SAR_netcdf) and (not FILE_is_HERE):
		BUILD_SAR_netcdf = True

	#print("SAR_FIRST_LOOK",SAR_FIRST_LOOK)
	#print("BUILD_SAR_netcdf",BUILD_SAR_netcdf)
	#print("FILE_is_HERE",FILE_is_HERE)

	if BUILD_SAR_netcdf:
		print('* Building SAR netcdf...')
		 # Box coordinates
		LEFT = 22.5     # °E
		RIGHT = 26      # °E
		BOTT = -37      # °N
		TOP = -34.5     # °N

		dataSAR1 = np.load(base_path + path_file['SAR1']) / 255
		dataSAR2 = np.load(base_path + path_file['SAR2']) / 255
		dataSST = np.load(base_path + path_file['SST']) 

		SAR1 = np.where(dataSAR1 < 1,dataSAR1,0)
		SAR2 = np.where(dataSAR2 < 1,dataSAR2,0)
		SAR = SAR1 + SAR2
		#SAR = np.

		MASK1 = np.where( dataSAR1>0,np.where( dataSAR1 < 1, 0,1),1)
		MASK2 = np.where( dataSAR2>0,np.where( dataSAR2 < 1, 0,1),1)
		MASK = MASK1*MASK2



		print(DegLat,DegLon)

		NX = dataSAR1.shape[1]
		NY = dataSAR1.shape[0]

		dx = (RIGHT-LEFT) / NX * DegLon
		dy = (TOP-BOTT) / NY * DegLat

		print('resolution of data is (km) :')
		print('dx=',dx,' dy=',dy)
		print('shape of data is')
		print(dataSAR1.shape)

		LON = np.linspace(LEFT,RIGHT,NX)
		LAT = np.linspace(TOP,BOTT,NY)

		data_vars = {'SAR1':(['lat','lon'],SAR1,
				            {'long_name':'Normalized SAR backscatter',
				            'units':'',
				            'timespan':start_time['SAR1'] + '-' + end_time['SAR1']}),
				    'SAR2':(['lat','lon'],SAR2,
				            {'long_name':'Normalized SAR backscatter',
				            'units':'',
				            'timespan':start_time['SAR2'] + '-' + end_time['SAR2']}),
				    'SAR':(['lat','lon'],SAR,
				        {'long_name':'Normalized SAR backscatter',
				            'units':'',
				            'timespan':start_time['SAR'] + '-' + end_time['SAR']}),
				    'MASK1':(['lat','lon'],MASK1,
				            {"lon_name":'mask of invalid data for SAR1 variable',
				            'units':'',
				            'timespan':start_time['SAR1'] + '-' + end_time['SAR1']}),
				    'MASK2':(['lat','lon'],MASK2,
				            {"lon_name":'mask of invalid data for SAR2 variable',
				            'units':'',
				            'timespan':start_time['SAR2'] + '-' + end_time['SAR2']}),
				    'MASK':(['lat','lon'],MASK,
				            {"lon_name":'mask of invalid data for SAR variable',
				            'units':'',
				            'timespan':start_time['SAR'] + '-' + end_time['SAR']})
				    }
				            
		coords={'lat': LAT,'lon':LON}
		attrs={'origin':'Extracted from OceanDataLab, with the .npy extract feature',
			'process':'LON = np.linspace(LEFT,RIGHT,NX)',
			'x resolution (m)':str(dx*1000),
			'y resolution (m)':str(dy*1000),
			'file1':path_file['SAR1'],
			'file2':path_file['SAR2'],
			'website':'https://ovl.oceandatalab.com/?date=2015-12-10T12:00:13&timespan=1d&extent=1918521.2293534_-4773856.3112472_3483951.5684158_-3989918.1492636&center=2701236.3988846_-4381887.2302554&zoom=8&products=3857_SAR_roughness!3857_ODYSSEA_REG_SST!3857_ODYSSEA_SST!3857_GlobCurrent_CMEMS_geostrophic_streamline!3857_AMSR_sea_ice_concentration!3857_GIBS_MODIS_Terra_CorrectedReflectance_TrueColor!3857_GIBS_MODIS_Aqua_CorrectedReflectance_TrueColor!900913_User_Shapes&opacity=80_100_100_60_100_49.498_100_100&stackLevel=100.02_50.03_30.02_120.07_40.03_50.25_50.22_10000&filter=+A,+B,+IW,+EW,+SM,+WV,+VV,+HH&selection=01000000',
			'Left border':LEFT,
			'Right border':RIGHT,
			'Bottom border':BOTT,
			'Top borde':TOP,
			'github repo':githubrepo}
		ds = xr.Dataset(data_vars=data_vars,coords=coords,attrs=attrs)
		ds.to_netcdf(path=base_path+name_netcdf,mode='w') 

	if SAR_FIRST_LOOK:
		print('* Plotting SAR first look')
		ds = xr.open_dataset(base_path+name_netcdf)
		LON = ds.lon
		LAT = ds.lat
		SAR1 = ds.SAR1
		SAR2 = ds.SAR2
		SAR = ds.SAR
		MASK1 = ds.MASK1
		MASK2 = ds.MASK2
		MASK = ds.MASK

		figsize=(5,5)

		#        fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
		#        ax.pcolormesh(LON,LAT,np.ma.masked_where(MASK1,SAR1),shading='nearest',cmap='Greys')
		#        ax.set_aspect(DegLat/DegLon)
		#        ax.set_ylabel('LAT')
		#        ax.set_xlabel('LON')

		#        fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
		#        ax.pcolormesh(LON,LAT,np.ma.masked_where(MASK2,SAR2),shading='nearest',cmap='Greys')
		#        ax.set_ylabel('LAT')
		#        ax.set_aspect(DegLat/DegLon)
		#        ax.set_xlabel('LON')

		fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
		s = ax.pcolormesh(LON,LAT,np.ma.masked_where(MASK1 * MASK2,SAR),shading='nearest',cmap='Greys')
		plt.colorbar(s,ax=ax,pad=0.05, aspect=50)
		ax.set_ylabel('LAT')
		ax.set_aspect(DegLat/DegLon)
		ax.set_xlabel('LON')
		fig.savefig(png_save+'SAR.png')
	   
	if SAR_VIGNETTES:
		print('* Plotting SAR vignettes')
		dico_loc = {'box1':{'O':(24.75,-35.5), # [LEFT,BOTT]
						'Lx':20,
						'Ly':20
						},
					'box2':{'O':(24.5,-35.8), # [LEFT,BOTT]
						'Lx':20,
						'Ly':20
						},
					'box3':{'O':(24.75,-36), # [LEFT,BOTT]
						'Lx':20,
						'Ly':20
						},
				}
		vmin,vmax = 0.3,0.6
		ds = xr.open_dataset(base_path+name_netcdf)
		LON = ds.lon
		LAT = ds.lat
		SAR = ds.SAR
		MASK = ds.MASK

		figGLO, axGLO = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
		s = axGLO.pcolormesh(LON,LAT,np.ma.masked_where(MASK,SAR),
				shading='nearest',cmap='Greys_r',vmin=vmin,vmax=vmax)
		plt.colorbar(s,ax=axGLO,pad=0.05, aspect=50)
		axGLO.set_ylabel('LAT')
		axGLO.set_xlabel('LON')
		axGLO.set_aspect(DegLat/DegLon)

		for name in dico_loc.keys():
			O = dico_loc[name]['O']
			Lx = dico_loc[name]['Lx']
			Ly = dico_loc[name]['Ly']

			# add rectangle on global figure
			axGLO.add_patch(mpl.patches.Rectangle(O, Lx/DegLon, Ly/DegLat,
				 edgecolor='k',
				 fill=False,
				 lw=2))
			# solo plot of the box
			fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=100)
			s = ax.pcolormesh(LON,LAT,np.ma.masked_where(MASK,SAR),
					shading='nearest',cmap='Greys_r',vmin=vmin,vmax=vmax)
			plt.colorbar(s,ax=ax,pad=0.05, aspect=50)
			ax.set_ylabel('LAT')
			ax.set_aspect(DegLat/DegLon)
			ax.set_xlabel('LON')
			ax.set_xlim([O[0],O[0]+Lx/DegLon])
			ax.set_ylim([O[1],O[1]+Ly/DegLat])
			ax.set_title(name+', O='+str(O)+' Lx='+str(Lx)+'km, Ly='+str(Ly)+'km')
			fig.savefig(png_save+'SAR_'+name+'.png')
		figGLO.savefig(png_save+'SAR_all_boxes.png')
		
if SST_FIRST_LOOK:
    """
    SST DATA L4 Ifremer
    See 'DataOrigin'
    """
    print('* Plotting SST first look')
    path_sst = pwd_path + 'SST_from_Ifremer/'
    filename = '20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc'
    ds = xr.open_dataset(path_sst + filename)
    
    figsize = (8,8)
    #LEFT = 22.5     # °E
    #RIGHT = 26      # °E
    #BOTT = -37      # °N
    #TOP = -34.5     # °N
    LEFT = 24.5     # °E
    RIGHT = 25.2      # °E
    BOTT = -36      # °N
    TOP = -35.19     # °N
    vmin,vmax = 295,298 # K
    cmap = 'plasma'
    lat = ds.lat
    lon = ds.lon
    sst = ds.analysed_sst.isel(time=0)
    levels = [295.65] # 295.65 K = 22.5°C
    SST = sst.sel(lat=slice(BOTT,TOP),lon=slice(LEFT,RIGHT))
    maxSST = SST.max()
    minSST = SST.min()
    meanSST = SST.mean()
    mediSST = (maxSST+minSST)/2
    print('     Contours at ',levels,'K')
    print('minSST,maxSST,meanSST,mediSST',minSST.values,maxSST.values,meanSST.values,mediSST.values)
    
    
    gradXsst = SST.differentiate('lon')/DegLon
    gradYsst = SST.differentiate('lat')/DegLat
    
    gradSST = np.sqrt( gradXsst**2 + gradYsst**2)
    
    # Plot
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=100)
    s = ax.pcolormesh(lon,lat,sst,shading='nearest',cmap=cmap,vmin=vmin,vmax=vmax)
    ax.contour(lon,lat,sst,levels=levels,colors='k')
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50,label='SST (K)')
    ax.set_xlim([LEFT,RIGHT])
    ax.set_ylim([BOTT,TOP])
    ax.set_aspect(1)
    ax.set_ylabel('LAT')
    ax.set_xlabel('LON')
    ax.set_title('SST reanalysis from L4 Ifremer')
    fig.savefig(png_save+'SST.png')
    
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=100)
    s = ax.pcolormesh(lon.sel(lon=slice(LEFT,RIGHT)),lat.sel(lat=slice(BOTT,TOP)),
                      gradSST,shading='nearest',cmap=cmap) # ,vmin=vmin,vmax=vmax
    ax.contour(lon,lat,sst,levels=levels,colors='k')
    plt.colorbar(s,ax=ax,pad=0.01,aspect=50,label='gradSST (°K/km)')
    ax.set_xlim([LEFT,RIGHT])
    ax.set_ylim([BOTT,TOP])
    ax.set_aspect(1)
    ax.set_ylabel('LAT')
    ax.set_xlabel('LON')
    ax.set_title('gradSST (°K/km) from reanalysis L4 Ifremer')
    fig.savefig(png_save+'gradSST.png')
     
if SST_SAR_FIRST_LOOK:
    print('* Plotting SST+SAR first look, with LES domains definition')
    
    # LES boundaries
    O_s = (24.5,-36) # 24.5,-36
    eps_x = 200 # in dad number of points, away from left
    eps_y = 100 # in dad number of points, away from bottom
    dx_d = 200 # km
    dx_s = 50 # km
    Lx_s = 64 # km
    Ly_s = 90 # km
    ratio_dx = 4
    
    O_d = ( O_s[0] - eps_x*dx_d/1000/DegLon, O_s[1] - eps_y*dx_d/1000/DegLat)
    Ly_d = Ly_s + 2*eps_y*dx_d/1000
    Lx_d = Lx_s + eps_x*dx_d/1000
    Nx_d = int(Lx_d/dx_d * 1000)
    Ny_d = int(Ly_d/dx_d * 1000)
    
    Nx_s = int(Lx_s/dx_s*1000)
    Ny_s = int(Ly_s/dx_s*1000)
    
    
    print('')
    print('=> Input dimensions for simulations:')
    print('DAD: Nx, Ny =',Nx_d,Ny_d) # 540,640
    print('SON: Nx, Ny =',Nx_s,Ny_s)
    print('SON in DAD numbers: Nx, Ny =',Nx_s/ratio_dx,Ny_s/ratio_dx)
    print('YOU NEED TO MAKE SURE THAT Nx Ny CAN BE WRITTENT \n AS 2^n.3^m.5^l with m,n,l INTEGERS')
    print('IXOR =',eps_x)
    print('IYOR =',eps_y)
    print('dx dad / dx son =',int(dx_d/dx_s))
    print('')
    print('=> Informations about the simulations')
    print('DAD')
    print('- dx',dx_d)
    print('- Lx, Ly = ',Ly_d,Lx_d)
    print('- Nx, Ny =',Nx_d,Ny_d)
    print('- origin at',O_d[0],'°E,',O_d[1],'°N')
    print('')
    Plot_CLI_sim_dim(O_d,Lx_d,Ly_d,DegLat,DegLon)
    
    print('SON')
    print('- dx',dx_s)
    print('- Lx, Ly = ',Ly_s,Lx_s)
    print('- Nx_s, Ny_s =',Nx_s,Ny_s)
    print('- origin at',O_s[0],'°E,',O_s[1],'°N')
    print('')
    Plot_CLI_sim_dim(O_s,Lx_s,Ly_s,DegLat,DegLon)
    
    # Opening files
    sstname = pwd_path + 'SST_from_Ifremer/' + '20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc'
    sarname = pwd_path+ 'SAR_from_OVL/' + 'SAR_roughness_20151210t170827-20151210t170921.nc'
    
    dsSST = xr.open_dataset(sstname)
    dsSAR = xr.open_dataset(sarname)
    
    # figure boundaries
    LEFT = 22.5     # °E
    RIGHT = 26      # °E
    BOTT = -37      # °N
    TOP = -34.5     # °N
    vmin,vmax = 295,298 # K
    Nlevel = 5
    Levels = np.linspace(vmin,vmax,Nlevel)
    dT = Levels[1] - Levels[0]
    
    brighter = 1
    
    figsize = (8,6)
    cmap = plt.cm.plasma  # define the colormap
    bounds = np.arange(Levels[0] - dT/2,Levels[-1] + dT/2, dT) # define the bins
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N) # normalize
    
    # Plot
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=100)
    s = ax.contour(dsSST.lon,dsSST.lat,dsSST.analysed_sst.isel(time=0),cmap=cmap,levels=Levels)
    ax.pcolormesh(dsSAR.lon,dsSAR.lat,np.ma.masked_where(dsSAR.MASK1 * dsSAR.MASK2,dsSAR.SAR)*brighter,shading='nearest',
                cmap='Greys',vmin=0,vmax=1)
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),ax=ax, cmap=cmap, norm=norm,ticks=Levels,
                 aspect=50,pad=0.01,label='SST (K)')#
    # dad
    ax.add_patch(mpl.patches.Rectangle(O_d, Lx_d/DegLon, Ly_d/DegLat,
                 edgecolor='k',
                 fill=False,
                 lw=2))
    #son
    ax.add_patch(mpl.patches.Rectangle(O_s, Lx_s/DegLon, Ly_s/DegLat,
                 edgecolor='green',
                 fill=False,
                 lw=2))
    ax.set_xlim([LEFT,RIGHT])
    ax.set_ylim([BOTT,TOP])
    ax.set_aspect(DegLat/DegLon)
    ax.set_ylabel('LAT')
    ax.set_xlabel('LON')
    fig.savefig(png_save+'SST_and_SAR.png')

if CHANGE_ENV:
    # this can be wrapped in a function that gives z1,z2,tht1,tht2 etc ..
    print('* Input values to check different environmental conditions')
    """
    This small tool gives the values to enter in the PRE_IDEAL.nam for MNH simulation, 
        from an idealized profile shape and initial value for the mixed layer as well 
        as the stratification intensity for thermodynamic variables.

    Shapes of the profiles:

     Z  ^                               ^                                   ^
        |           |                   |            /                      |       \ 
        |           |                   |           /     dtht/dz           |        \     drv/dz
        |           |                   |          /                        |         \ 
        |           |                  _|_        /                        _|_         \ 
        |           |                zi |        |                       zi |           |
        |           |                   |        |                          |           |  
        |           |                   |        |                          |           |                               
        X--------------------> U        X--------------------> tht          X--------------------> rv
                    U0                          thtML                                   rvML

    Note : the potential temperature computed here use P0 = 10^5 Pa as in MNH.
           To compare directly to the SST (and to have tht(z=0)=T(z=0)), you need to multiply theta 
               by a factor that takes into account this change of reference.
    
            THT1 = THT2 * a
     
            with    THT1: the potential temperature with Pref1 = P(z=0)
                    THT2: the potential temperature used in MNH with Pref2 = 10^5 Pa
                    a: the ratio of reference pressure
        
                            a = (Pref1/Pref2)**(Rd/Cpd)

            Another possibility is to transforme the SST value into a tht(SST) to compare with THT2

    """

    # INPUT ZONE ------------------------------------------------
    Psurface = 101250 # Pa, used to compute T
    END_Z = 2000 # end of LES domain  
    # thermo initial conditions
    d_case = {'ERA5':{'strat':{"tht":8.2,    # K/km
                            'rv':-8.36   # g/kg/km
                        },
                'iniML' : {"tht":294.25, # K/km
                            'rv':10.5     # g/kg/km
                        },
                'ABLHini' : 750 # m
                        },   
            '0':{'strat':{"tht":3.0,    # K/km
                            'rv':-4.6   # g/kg/km
                        },
                'iniML' : {"tht":295.5, # K/km
                            'rv':10     # g/kg/km
                        },
                'ABLHini' : 250 # m
                        }
                }
                
    units = {'tht':'K',
            'rv':'g/kg',
            'U':'m/s'}

    # wind initial conditions
    U = 7.5 # m/s
    V = 0 # m/s
    sub = 0 # m/s, subsidence
    # END INPUT ZONE ------------------------------------------------

    rounding = 5
    b = {}

    def f(x,a,b):
        return a*x+b
    def give_b(y,a,x):
        return y-a*x

    values = {'tht':np.zeros(2),
                'rv':np.zeros(2)}
    for case in d_case.keys():
        print('')
        print('=> CASE '+case,'ABLHini =',d_case[case]['ABLHini'],'m -----------------------')
        print('')
        print('wind is :')
        print('     U=',U,units['U'])
        print('     V=',V,units['U'])
        for var in d_case[case]['strat'].keys():
            
            print(var,'MLvalue = '+str(d_case[case]['iniML'][var])+' '+units[var]+', d'+var+'/dz = '+str(d_case[case]['strat'][var])+' '+units[var]+'/km')
            a = d_case[case]['strat'][var]/1000
            b[var] = give_b(d_case[case]['iniML'][var],a,d_case[case]['ABLHini'])

            values[var][0] = f(d_case[case]['ABLHini'],a,b[var])
            print('     ',np.round(values[var][0],rounding),units[var],'at z =',d_case[case]['ABLHini'],'m')
            values[var][1] = f(END_Z,a,b[var])
            print('     ',np.round(values[var][1],rounding),units[var],'at z =',END_Z,'m')
            if var=='tht':
                print(' T at surface = '+str(np.round(Theta_to_T(values['tht'][0],Psurface),rounding))+' K')

            values['thtv']=Compute_THTV(values['tht'],values['rv']/1000)
            dthtvdz = (values['thtv'][-1] - values['thtv'][0]) / (END_Z-d_case[case]['ABLHini'])
            print('thtv MLvalue = '+str(np.round(values['thtv'][0],rounding))+' '+units['tht']+
                    ', dtht/dz = '+str(np.round(dthtvdz*1000,rounding))+' '+units['tht']+'/km')
            print('     ',np.round(values['thtv'][0],rounding),'K at z =',d_case[case]['ABLHini'],'m',)
            print('     ',np.round(values['thtv'][1],rounding),'K at z =',END_Z,'m',)

if ERA5_FIRST_LOOK:
    # lets have a look at ERA5 from ECMWF
    print(' * ERA5 first look')

    # Note : the potential temperature computed here use P0 = 10^5 Pa as in MNH.
    #       To compare directly to the SST (and to have tht(z=0)=T(z=0)), you need to multiply theta 
    #           by a factor that takes into account this change of reference.
    #
    #   THT1 = THT2 * a
    # 
    #   with THT1: the potential temperature with Pref1 = P(z=0)
    #        THT2: the potential temperature used in MNH with Pref2 = 10^5 Pa
    #        a: the ratio of reference pressure
    #
    #           a = (Pref1/Pref2)**(Rd/Cpd)

    surface_path = 'ERA5/7f3b1d80e6888f8b9b31e0c5d8cee1a9.nc' # retrieved with API script single_level_retrieve.py
    atm_path = 'ERA5/dddc7db19be813e3ed36c9b2b41262f7.nc' # retrieved with API script pressure_level_retrieve.py

    # LES boundaries (dad)
    LEFT = 24     # °E
    RIGHT = 25.2      # °E
    BOTT = -36.5      # °N
    TOP = -35     # °N

    POINT_LOC = [(24.5,-36.25), # this is on cold side of SST
                (24.5,-35.25) 	# this is on warm side of SST
                ]
    ZONE_LOC = [(24,26),(-35,-37)]
    colors_loc = ['b','r']

    indt = 12 # 12H00

    # Gathering data
    dsSurf = xr.open_dataset(surface_path).isel(valid_time=indt)
    ds3D = xr.open_dataset(atm_path)
    seltime = ds3D.valid_time.sel(valid_time=dsSurf.valid_time.values,method='nearest')
    ds3D = ds3D.sel(valid_time=seltime)
    
    SST = dsSurf.sst
    landmaskSurf = np.ma.masked_invalid(SST).mask
    U10 = np.ma.masked_where(landmaskSurf,dsSurf.u10)
    V10 = np.ma.masked_where(landmaskSurf,dsSurf.v10)
    with warnings.catch_warnings(action="ignore"):
        M10 = np.sqrt(U10**2+V10**2)
    

    ### Surface values
    print(' > at surface')
    figsize = (6,5)
    dpi=100
    pathsave = 'PNGs_observations/ERA5_'
    # U10
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    s = ax.pcolormesh(dsSurf.longitude,dsSurf.latitude,U10,
                    cmap='jet',vmin=2,vmax=15)
    plt.colorbar(s,ax=ax,label='U10 m/s',aspect=50,pad=0.01)
    ax.set_xlabel('LON')
    ax.set_ylabel('LAT')
    ax.set_aspect(1)
    ax.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    fig.savefig(pathsave + 'U10.png')
    # V10
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    s = ax.pcolormesh(dsSurf.longitude,dsSurf.latitude,V10,
                    cmap='jet',vmin=-5,vmax=5)
    plt.colorbar(s,ax=ax,label='V10 m/s',aspect=50,pad=0.01)
    ax.set_xlabel('LON')
    ax.set_ylabel('LAT')
    ax.set_aspect(1)
    ax.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    fig.savefig(pathsave + 'V10.png')
    # wind speed norm
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    s = ax.pcolormesh(dsSurf.longitude,dsSurf.latitude,M10,
                    cmap='jet',vmin=0,vmax=15)
    plt.colorbar(s,ax=ax,label='M10 m/s',aspect=50,pad=0.01)
    ax.barbs(dsSurf.longitude.values[::2],dsSurf.latitude.values[::2], 
                U10[::2,::2], V10[::2,::2],length=3)
    ax.add_patch(mpl.patches.Rectangle((LEFT,BOTT), RIGHT-LEFT, TOP-BOTT,
                    edgecolor='k',
                    fill=False,
                    lw=2))
    ax.set_xlabel('LON')
    ax.set_ylabel('LAT')
    ax.set_aspect(1)
    ax.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    fig.savefig(pathsave + 'M10.png')
    # wind direction (streamlines)
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    s = ax.pcolormesh(dsSurf.longitude,dsSurf.latitude,SST,
                    cmap='jet',vmin=290,vmax=300)
    plt.colorbar(s,ax=ax,label='SST K',aspect=50,pad=0.01)
    for k in range(len(POINT_LOC)):
        ax.scatter(POINT_LOC[k][0],POINT_LOC[k][1],marker='x',c='grey')
    ax.barbs(dsSurf.longitude.values[::2],dsSurf.latitude.values[::2], 
                U10[::2,::2], V10[::2,::2],length=3)
    ax.add_patch(mpl.patches.Rectangle((LEFT,BOTT), RIGHT-LEFT, TOP-BOTT,
                    edgecolor='k',
                    fill=False,
                    lw=2))
    ax.set_xlabel('LON')
    ax.set_ylabel('LAT')
    ax.set_aspect(1)
    ax.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    fig.savefig(pathsave + 'SST.png')
    # SST zoomed with profile location
    fig, ax = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    s = ax.pcolormesh(dsSurf.longitude,dsSurf.latitude,SST,
                    cmap='jet',vmin=290,vmax=300)
    plt.colorbar(s,ax=ax,label='SST K',aspect=50,pad=0.01)
    for k in range(len(POINT_LOC)):
        ax.scatter(POINT_LOC[k][0],POINT_LOC[k][1],marker='x',c='grey')
    ax.barbs(dsSurf.longitude.values[::2],dsSurf.latitude.values[::2], 
                U10[::2,::2], V10[::2,::2],length=5)
    ax.add_patch(mpl.patches.Rectangle((LEFT,BOTT), RIGHT-LEFT, TOP-BOTT,
                    edgecolor='k',
                    fill=False,
                    lw=2))
    ax.set_xlabel('LON')
    ax.set_ylabel('LAT')
    ax.set_aspect(1)
    ax.set_xlim([23,26])
    ax.set_ylim([-37,-34])
    ax.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    fig.savefig(pathsave + 'SST_zoomed.png')


    ### PROFILES ##
    ZMAX = 2000 #m

    print(' > profiles')
    figU, axU = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    figq, axq = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    figt, axt = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    figtht, axtht = plt.subplots(1,1,figsize = figsize,constrained_layout=True,dpi=dpi)
    
     ### > Mean profiles in the ZONE_LOC box1
    ds = ds3D.sel(longitude=slice(ZONE_LOC[0][0],ZONE_LOC[0][1]),
                latitude=slice(ZONE_LOC[1][0],ZONE_LOC[1][1])
                ).mean(dim=['longitude','latitude'])
    dss = dsSurf.sel(longitude=slice(ZONE_LOC[0][0],ZONE_LOC[0][1]),
                latitude=slice(ZONE_LOC[1][0],ZONE_LOC[1][1])
                ).mean(dim=['longitude','latitude'])
    Z = ds.z/9.81
    # U,V
    axU.plot(ds.u,Z,c='k',ls='-',label=str(ZONE_LOC[0][0])+'/'+str(ZONE_LOC[0][1])+'E_'+str(ZONE_LOC[1][0])+'/'+str(ZONE_LOC[1][1])+'N')
    axU.plot(ds.v,Z,c='k',ls='--')
    axU.scatter(dss.u10,10,c='k',marker='x')
    axU.scatter(dss.v10,10,c='k',marker='+')
    # q
    axq.plot(ds.q*1000,Z,c='k',label=str(ZONE_LOC[0][0])+'/'+str(ZONE_LOC[0][1])+'E_'+str(ZONE_LOC[1][0])+'/'+str(ZONE_LOC[1][1])+'N')
    rv = rv_from_Td(dss.d2m,dss.t2m,dss.msl)
    axq.scatter(rv*1000,2,c='k',marker='x')
    # t 
    axt.plot(ds.t,Z,c='k',label=str(ZONE_LOC[0][0])+'/'+str(ZONE_LOC[0][1])+'E_'+str(ZONE_LOC[1][0])+'/'+str(ZONE_LOC[1][1])+'N')
    axt.scatter(dss.t2m,2,c='k',marker='x')
    # tht
    tht = T_to_Theta(ds.t.values,ds.pressure_level.values*100)
    tht2m = T_to_Theta(dss.t2m,dss.msl.values)
    axtht.plot(tht,Z,c='k',label=str(ZONE_LOC[0][0])+'/'+str(ZONE_LOC[0][1])+'E_'+str(ZONE_LOC[1][0])+'/'+str(ZONE_LOC[1][1])+'N')
    axtht.plot(tht2m,2,c='k',marker='x')
    ### > Profiles at some locations
    #       with a cross at these locations on the streamline plot
    for k in range(len(POINT_LOC)):
        ds = ds3D.sel(longitude=POINT_LOC[k][0],latitude=POINT_LOC[k][1])
        dss = dsSurf.sel(longitude=POINT_LOC[k][0],latitude=POINT_LOC[k][1])
        Z = ds.z/9.81
        # U,V
        axU.plot(ds.u,Z,c=colors_loc[k],ls='-',
                    label=str(POINT_LOC[k][0])+'E_'+str(POINT_LOC[k][1])+'N')
        axU.plot(ds.v,Z,c=colors_loc[k],ls='--')
        axU.scatter(dss.u10,10,c=colors_loc[k],marker='x')
        axU.scatter(dss.v10,10,c=colors_loc[k],marker='+')
        # q
        axq.plot(ds.q*1000,Z,c=colors_loc[k],
                    label=str(POINT_LOC[k][0])+'E_'+str(POINT_LOC[k][1])+'N')
        rv = rv_from_Td(dss.d2m,dss.t2m,dss.msl)
        axq.scatter(rv*1000,2,c=colors_loc[k],marker='x')
        # t
        axt.plot(ds.t,Z,c=colors_loc[k],
                    label=str(POINT_LOC[k][0])+'E_'+str(POINT_LOC[k][1])+'N')
        axt.scatter(dss.t2m,2,c=colors_loc[k],marker='x')
        # tht
        tht = T_to_Theta(ds.t.values,ds.pressure_level.values*100)
        axtht.plot(tht,Z,c=colors_loc[k],
                    label=str(POINT_LOC[k][0])+'E_'+str(POINT_LOC[k][1])+'N')
        tht2m = T_to_Theta(dss.t2m,dss.msl.values)
        axtht.plot(tht2m,2,c=colors_loc[k],marker='x')
        print(dss.msl.values)
    axU.set_xlabel('U(-), V(--) m/s')
    axU.legend()
    axU.set_ylabel('Z m')
    axU.set_ylim([0,ZMAX])
    axU.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    axU.legend()
    	
    axq.set_xlabel('q g/kg')
    axq.set_ylabel('Z m')
    axq.set_ylim([0,ZMAX])
    axq.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    axq.legend()
    
    axt.set_xlabel('T K')
    axt.set_ylabel('Z m')
    axt.set_ylim([0,ZMAX])
    axt.set_xlim([280,297])
    axt.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    axt.legend()
    
    axtht.set_xlabel(r'$\theta$ K')
    axtht.set_ylabel('Z m')
    axtht.set_ylim([0,ZMAX])
    axtht.set_xlim([292,305])
    axtht.set_title('ERA5 '+str(dsSurf.valid_time.values)[:-10])
    axtht.legend()
        
    # saving
    figU.savefig(pathsave +'u.png')
    figq.savefig(pathsave +'q.png')
    figt.savefig(pathsave +'t.png')
    figtht.savefig(pathsave +'tht.png')

if COMPARE_SAR_SOURCES:
    dsOVL =  xr.open_dataset('SAR_from_OVL/SAR_roughness_20151210t170827-20151210t170921.nc')
    dsIFRE = xr.open_dataset('SAR_from_IFREMER/S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_9AE1.nc')

    print(dsOVL) # res is 125m
    print(dsIFRE) # res is 100m

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=200)
    ax.pcolormesh(dsOVL.lon,dsOVL.lat,dsOVL.SAR,cmap='Greys_r',vmin=0.3,vmax=0.6)
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    ax.set_title('SAR from OVL')
    ax.set_xlim([24,24.5])
    ax.set_ylim([-36,-35.5])

    fig, ax = plt.subplots(1,1,figsize = (7,5),constrained_layout=True,dpi=200)
    ax.pcolormesh(dsIFRE.longitude,dsIFRE.latitude,dsIFRE.sigma0_detrend.isel(pol=0),cmap='Greys_r',vmin=0.,vmax=0.1)
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    ax.set_title('SAR from Ifremer, pol=VV')
    ax.set_xlim([24,24.5])
    ax.set_ylim([-36,-35.5])




plt.show()













