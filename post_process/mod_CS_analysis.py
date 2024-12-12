from module_tools import *
from scipy.ndimage import label as SPlabel

def mean_vertical_contrib(flx_i,flx_mean,indzi):
	"""This function is computing the mean contribution of flx_i to the total
		flx_mean from surface to Z[indzi]
		
		eq 8 of Brient 2018 "Object-Oriented Identification of Coherent Structures in
					Large Eddy Simulations: Importance of Downdrafts
					in Stratocumulus"
	INPUT:			
		- flx_i    : field of the contribution from structure i to flux flx (xr.DataArray), Z is first dim
		- flx_mean : same size as flx_i, domain average flx (xr.DataArray)
		- indzi    : Z index to stop the average (usually 1.1zi)
    OUTPUT:
        - scalar, ratio of the contribution from ith object
	"""
	return np.abs(flx_i.isel(level=slice(2,indzi))).integrate('level') / np.abs(flx_mean.isel(level=slice(2,indzi))).integrate('level')

def compute_flx_contrib(flx,L_mask,function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
	"""This function is computing the contribution of the flux a'b' (='flx')
		relative to 'mask'. This is the sum of the "top-hat" and "inter/intra variability" 
		
		F_i = alpha_i * flx
		
			where alpha_i is the area cover of structure i
				alpha_i = 1/N * sum_over_ipoints(mask)
	INPUT:	
		- flx : flux to be partionned
		- L_mask : list of mask for the current structure, 1=belong to structure
		- meanDim : dimensions to be averaged
    OUTPUT:
        - 1D array of contribution from ith object to flx
	"""	
	L_out = []
	for i,mask in enumerate(L_mask):
		flx_i_m = function_avg(flx.where(mask,other=0)) # 
		L_out.append( flx_i_m )
	if len(L_mask)==1:
		return L_out[0]
	else:
		return tuple(flx_i for flx_i in L_out)
      
def var_conditional_object(var,L_mask,function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
    """
    Computes the values of 'var' conditioned by objects masks.

    INPUT:
        - var : xarray DataArray
        - L_mask: liste of mask for each object (1=object, 0=other)
        - function_avg : reynolds operator

    OUTPUT:
        - tuple of profiles (length = length L_mask)
    """

    L_out = []
    for i,mask in enumerate(L_mask):
        var_i = function_avg(var.where(mask))
        L_out.append( var_i )
    if len(L_mask)==1:
        return L_out[0]
    else:
        return tuple(var_i for var_i in L_out)

def Plot_mean_profiles(dsCS,dsmean,                    
                        ax,NORM,
                        function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
    """
    This function plots the profiles of prognostic variables with object decomposition

    INPUT:
        - dsCS : xarray dataset with prognostic variables at mass point and object masks
        - dsmean : xarray dataset with mean var profiles
        - dsflx : xarray dataset with mean flux profiles
        - ax : 2D axis, the matplotlib axis where to plot
        - NORM : boolean to normalise fluxes or not
    OUTPUT:
        - 1 figure for each variables with domain averaged variables as well as contribution from each object
        (variables are : U,V,W,THT,THTV,RVT,E)
    """
    K = 0.01 # plot only if coverage > K 
    Z = dsCS.level
    zi = get_ABLH(Z,dsmean.THTvm) 
    indzmax = nearest(Z.values,1.1*zi)

    U,V,W = dsCS.U,dsCS.V,dsCS.W
    Um,Vm,Wm = dsmean.Um,dsmean.Vm,dsmean.Wm
    THTV = dsCS.THTV
    THTVm = dsmean.THTvm

    # object mask
    is_up1 = xr.where( dsCS.global_mask==1,1,0 )
    is_ss1 = xr.where( dsCS.global_mask==2,1,0 )
    is_up2 = xr.where( dsCS.global_mask==3,1,0 )
    is_ss2 = xr.where( dsCS.global_mask==4,1,0 )
    is_down = xr.where( dsCS.global_mask==5,1,0 )
    A_up1 = function_avg(is_up1)
    A_ss1 = function_avg(is_ss1)
    A_up2 = function_avg(is_up2)
    A_ss2 = function_avg(is_ss2)
    A_down = function_avg(is_down)

    Lcolors = {'up1':'r',
                'ss1':'purple',
                'up2':'orange',
                'ss2':'pink',
                'down':'green'}
    if NORM:
        norm_U = 1
        norm_wthtv = 1
    else:
        norm_U = 1
        norm_wthtv = 1

    
    ### PLOTS
    # > area
    ax[0].plot(A_up1,Z/zi,c=Lcolors['up1'],label='up1')
    ax[0].plot(A_ss1,Z/zi,c=Lcolors['ss1'],label='ss1')
    ax[0].plot(A_up2,Z/zi,c=Lcolors['up2'],label='up2')
    ax[0].plot(A_ss2,Z/zi,c=Lcolors['ss2'],label='ss2')
    ax[0].plot(A_down,Z/zi,c=Lcolors['down'],label='down')



    # > U
    ax[1].plot(Um/norm_U,Z/zi,c='k',label='<U>',ls='-')
    # obj contrib
    U_up1,U_ss1,U_up2,U_ss2,U_down = var_conditional_object(U,[is_up1,is_ss1,is_up2,is_ss2,is_down],function_avg)
    ax[1].plot(U_up1.where(A_up1>K)/norm_U,Z/zi,c=Lcolors['up1'],label='up1')
    ax[1].plot(U_ss1.where(A_ss1>K)/norm_U,Z/zi,c=Lcolors['ss1'],label='ss1')
    ax[1].plot(U_up2.where(A_up2>K)/norm_U,Z/zi,c=Lcolors['up2'],label='up2')
    ax[1].plot(U_ss2.where(A_ss2>K)/norm_U,Z/zi,c=Lcolors['ss2'],label='ss2')
    ax[1].plot(U_down.where(A_down>K)/norm_U,Z/zi,c=Lcolors['down'],label='down')

    # > W
    ax[2].plot(Wm/norm_U,Z/zi,c='k',label='<W>',ls='-')
    # obj contrib
    W_up1,W_ss1,W_up2,W_ss2,W_down = var_conditional_object(W,[is_up1,is_ss1,is_up2,is_ss2,is_down],function_avg)
    ax[2].plot(W_up1.where(A_up1>K)/norm_U,Z/zi,c=Lcolors['up1'],label='up1')
    ax[2].plot(W_ss1.where(A_ss1>K)/norm_U,Z/zi,c=Lcolors['ss1'],label='ss1')
    ax[2].plot(W_up2.where(A_up2>K)/norm_U,Z/zi,c=Lcolors['up2'],label='up2')
    ax[2].plot(W_ss2.where(A_ss2>K)/norm_U,Z/zi,c=Lcolors['ss2'],label='ss2')
    ax[2].plot(W_down.where(A_down>K)/norm_U,Z/zi,c=Lcolors['down'],label='down')

    # > THTV
    ax[3].plot(THTVm,Z/zi,c='k',label=r'<$\theta_v$>',ls='-')
    # obj contrib
    thtv_up1,thtv_ss1,thtv_up2,thtv_ss2,thtv_down = var_conditional_object(THTV,[is_up1,is_ss1,is_up2,is_ss2,is_down],function_avg)
    ax[3].plot( (thtv_up1-THTVm).where(A_up1>K)/norm_U,Z/zi,c=Lcolors['up1'],label='up1')
    ax[3].plot( (thtv_ss1-THTVm).where(A_ss1>K)/norm_U,Z/zi,c=Lcolors['ss1'],label='ss1')
    ax[3].plot( (thtv_up2-THTVm).where(A_up2>K)/norm_U,Z/zi,c=Lcolors['up2'],label='up2')
    ax[3].plot( (thtv_ss2-THTVm).where(A_ss2>K)/norm_U,Z/zi,c=Lcolors['ss2'],label='ss2')
    ax[3].plot( (thtv_down-THTVm).where(A_down>K)/norm_U,Z/zi,c=Lcolors['down'],label='down')


def Plot_mean_flux(dsCS,dsmean,dsflx,
                    ax,NORM,
                    function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
    """
    This function plots the profiles of fluxes (total) with object decomposition (only resolved part)

    INPUT:
        - dsCS : xarray dataset with prognostic variables at mass point and object masks
        - dsmean : xarray dataset with mean var profiles
        - dsflx : xarray dataset with mean flux profiles
        - ax : 2D axis, the matplotlib axis where to plot
        - NORM : boolean to normalise fluxes or not
    OUTPUT:
        - 1 figure for each flux with domain averaged variables as well as contribution from each object
        (variables are : uw,wthtv)
    """
    Z = dsCS.level
    zi = get_ABLH(Z,dsmean.THTvm) # Boundary layer height
    indzmax = nearest(Z.values,1.1*zi) # for integrated vertical contribution

    U,V,W = dsCS.U,dsCS.V,dsCS.W
    THTV = dsCS.THTV

    u_fluc = dsCS.U - dsmean.Um
    w_fluc = dsCS.W - dsmean.Wm
    thtv_fluc = dsCS.THTV - dsmean.THTvm

    uw_t = dsflx.FLX_UW
    uw_s = dsflx.FLX_UW_s
    uw_r = uw_t - uw_s
    uw = u_fluc*w_fluc
    
    wthtv_t = dsflx.FLX_THvW
    wthtv_s = dsflx.FLX_THvW_s
    wthtv_r = wthtv_t - wthtv_s
    wthtv = w_fluc*thtv_fluc

    # object mask
    is_up1 = xr.where( dsCS.global_mask==1,1,0 )
    is_ss1 = xr.where( dsCS.global_mask==2,1,0 )
    is_up2 = xr.where( dsCS.global_mask==3,1,0 )
    is_ss2 = xr.where( dsCS.global_mask==4,1,0 )
    is_down = xr.where( dsCS.global_mask==5,1,0 )
    A_up1 = function_avg(is_up1)
    A_ss1 = function_avg(is_ss1)
    A_up2 = function_avg(is_up2)
    A_ss2 = function_avg(is_ss2)
    A_down = function_avg(is_down)
    Lcolors = {'up1':'r',
                'ss1':'purple',
                'up2':'orange',
                'ss2':'pink',
                'down':'green'}

    if NORM:
        norm_uw = dsmean.u_star.values**2
        norm_wthtv = dsmean.Qv_star
    else:
        norm_uw = 1
        norm_wthtv = 1
    
    ### PLOTS
    # > uw
    ax[0].plot(uw_r/norm_uw,Z/zi,c='k',label='uw(res)',ls='--')
    ax[0].plot(uw_t/norm_uw,Z/zi,c='k',label='uw(total)')
    # obj contrib
    uw_up1, uw_ss1, uw_up2, uw_ss2, uw_down = compute_flx_contrib(uw,[is_up1,is_ss1,is_up2,is_ss2,is_down],function_avg)
    somme = uw_up1+uw_ss1+uw_up2+uw_ss2+uw_down
    part_up1 = mean_vertical_contrib(uw_up1,uw_r,indzmax).values
    part_up2 = mean_vertical_contrib(uw_up2,uw_r,indzmax).values
    part_sub1 = mean_vertical_contrib(uw_ss1,uw_r,indzmax).values
    part_sub2 = mean_vertical_contrib(uw_ss2,uw_r,indzmax).values
    part_down = mean_vertical_contrib(uw_down,uw_r,indzmax).values
    obj_over_all = mean_vertical_contrib(somme,uw_r,indzmax).values
    ax[0].plot(uw_up1/norm_uw,Z/zi,c=Lcolors['up1'],label='up1 ('+str(np.round(part_up1*100,1))+'%)')
    ax[0].plot(uw_ss1/norm_uw,Z/zi,c=Lcolors['ss1'],label='ss1 ('+str(np.round(part_sub1*100,1))+'%)')
    ax[0].plot(uw_up2/norm_uw,Z/zi,c=Lcolors['up2'],label='up2 ('+str(np.round(part_up2*100,1))+'%)')
    ax[0].plot(uw_ss2/norm_uw,Z/zi,c=Lcolors['ss2'],label='ss2 ('+str(np.round(part_sub2*100,1))+'%)')
    ax[0].plot(uw_down/norm_uw,Z/zi,c=Lcolors['down'],label='down ('+str(np.round(part_down*100,1))+'%)')
    ax[0].plot(somme/norm_uw,Z/zi,c='gray',label='all ('+str(np.round(obj_over_all*100,1))+'%)')

    # > wthtv
    ax[1].plot(wthtv_r/norm_wthtv,Z/zi,c='k',label='wthtv(res)',ls='--')
    ax[1].plot(wthtv_t/norm_wthtv,Z/zi,c='k',label='wthtv(total)')
    # obj contrib
    wthtv_up1, wthtv_ss1, wthtv_up2, wthtv_ss2, wthtv_down = compute_flx_contrib(wthtv,[is_up1,is_ss1,is_up2,is_ss2,is_down],function_avg)
    somme = wthtv_up1+wthtv_ss1+wthtv_up2+wthtv_ss2+wthtv_down
    part_up1 = mean_vertical_contrib(wthtv_up1,wthtv_r,indzmax).values
    part_up2 = mean_vertical_contrib(wthtv_up2,wthtv_r,indzmax).values
    part_sub1 = mean_vertical_contrib(wthtv_ss1,wthtv_r,indzmax).values
    part_sub2 = mean_vertical_contrib(wthtv_ss2,wthtv_r,indzmax).values
    part_down = mean_vertical_contrib(wthtv_down,wthtv_r,indzmax).values
    obj_over_all = mean_vertical_contrib(somme,wthtv_r,indzmax).values
    ax[1].plot(wthtv_up1/norm_wthtv,Z/zi,c=Lcolors['up1'],label='up1 ('+str(np.round(part_up1*100,1))+'%)')
    ax[1].plot(wthtv_ss1/norm_wthtv,Z/zi,c=Lcolors['ss1'],label='ss1 ('+str(np.round(part_sub1*100,1))+'%)')
    ax[1].plot(wthtv_up2/norm_wthtv,Z/zi,c=Lcolors['up2'],label='up2 ('+str(np.round(part_up2*100,1))+'%)')
    ax[1].plot(wthtv_ss2/norm_wthtv,Z/zi,c=Lcolors['ss2'],label='ss2 ('+str(np.round(part_sub2*100,1))+'%)')
    ax[1].plot(wthtv_down/norm_wthtv,Z/zi,c=Lcolors['down'],label='down ('+str(np.round(part_down*100,1))+'%)')
    ax[1].plot(somme/norm_wthtv,Z/zi,c='gray',label='all ('+str(np.round(obj_over_all*100,1))+'%)')

def Plot_mean_progvar_allboxes(dsCS,dsmean,path_save,function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
    """
    This function plots the profiles of mean profiles with object decomposition for all boxes defined by dsCS.nboxe
    
    INPUT:
        - dsCS : xarray dataset with prognostic variables at mass point and object masks
        - dsmean : xarray dataset with mean var profiles
        - path_save : where to save figures
        - function_avg : reynolds average operator
    OUTPUT:
        - 1 figure for each flux with domain averaged variables as well with contribution from each object
            and for each boxes
        (variables are : Area,U,W,THTV)
    """

    Nboxe = len(dsCS.nboxe)
    NORM = True
    Z = dsCS.level
    figsize = (10,5)
    dpi = 200
    YLIM = [0,1.1]
    fig, ax = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # Area
    fig2, ax2 = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # U
    fig3, ax3 = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # W
    fig4, ax4 = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # THTV

    for k in range(Nboxe):
        Plot_mean_profiles(dsCS.isel(nboxe=k),dsmean.isel(nboxe=k),
                    [ax[k],ax2[k],ax3[k],ax4[k]],NORM,function_avg)
        # Area
        ax[k].set_xlim([-0.1,0.3])
        ax[k].set_xlabel(r"Area")
        ax[k].grid()
        ax[k].set_ylim(YLIM)
        ax[k].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
        ax[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        # U
        ax2[k].set_xlabel(r"U (m/s)")
        ax2[k].grid()
        ax2[k].set_ylim(YLIM)
        ax2[k].set_xlim([5,7.5])
        ax2[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        ax2[k].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax2[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        # W
        ax3[k].set_xlabel(r"W (m/s)")
        ax3[k].grid()
        ax3[k].set_ylim(YLIM)
        ax3[k].set_xlim([-1,1])
        ax3[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        ax3[k].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax3[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        # THTV
        ax4[k].set_xlabel(r"$\theta_v$-<$\theta_v$> (K)")
        ax4[k].grid()
        ax4[k].set_ylim(YLIM)
        ax4[k].set_xlim([-0.5,0.5])
        ax4[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        ax4[k].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax4[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[0].set_ylabel(r'z/$z_i$')
    ax2[0].set_ylabel(r'z/$z_i$')
    ax3[0].set_ylabel(r'z/$z_i$')
    ax4[0].set_ylabel(r'z/$z_i$')

    fig.savefig(path_save + 'Area.png')
    fig2.savefig(path_save + 'U_split.png')
    fig3.savefig(path_save + 'W_split.png')
    fig4.savefig(path_save + 'THTV_split.png')

def Plot_mean_flux_allboxes(dsCS,dsmean,dsflx,path_save,function_avg=lambda a : a.mean(dim=['time','ni','nj'])):
    """
    This function plots the profiles of fluxes (total) with object decomposition (only resolved part)
     for all boxes defined by dsCS.nboxe
    
    INPUT:
        - dsCS : xarray dataset with prognostic variables at mass point and object masks
        - dsmean : xarray dataset with mean var profiles
        - dsflx : xarray dataset with mean flux profiles
        - path_save : where to save figures
    OUTPUT:
        - 1 figure for each flux with domain averaged variables as well with contribution from each object
            and for each boxes
        (variables are : uw,wthtv)
    """

    Nboxe = len(dsCS.nboxe)
    NORM = True
    Z = dsCS.level
    figsize = (10,5)
    dpi = 200
    YLIM = [0,1.1]
    fig, ax = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # uw
    fig2, ax2 = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi) # wthtv

    for k in range(Nboxe):
        Plot_mean_flux(dsCS.isel(nboxe=k),dsmean.isel(nboxe=k),dsflx.isel(nboxe=k),
                    [ax[k],ax2[k]],NORM,function_avg)
        
        ax[k].set_xlim([-1.1,0.1])
        ax[k].set_xlabel(r"<u'w'>/u$^{*2}$")
        ax[k].grid()
        ax[k].set_ylim(YLIM)
        ax2[k].set_xlim([-0.6,1.1])
        ax[k].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        ax[k].legend()
        ax2[k].set_xlabel(r"<w'$\theta_v$'>/$Q_{v,0}$")
        ax2[k].grid()
        ax2[k].set_ylim(YLIM)
        ax2[k].set_title('Boxe '+str(dsCS.nboxe[k].values))
        ax2[k].xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        ax2[k].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax2[k].legend()
    ax[0].set_ylabel(r'z/$z_i$')
    ax2[0].set_ylabel(r'z/$z_i$')

    fig.savefig(path_save + 'uw_split.png')
    fig2.savefig(path_save + 'wthtv_split.png')

def Plot_top_view_var(dsCS,VAR,atZ,path_save):
    """
    This function plots a top view at 'atZ' altitude of 'var' with coherent structures in contours
    
    INPUT:
    - dsCS : xarray DataArray 
    - VAR : string with name of background variable
    - atZ : Altitude of the plot (in meters)
    - path_save : where to save figures
    OUTPUT:
    - 1 plot with Nboxe columns, background = var, contours = CS
    """
    figsize = (10,5)
    dpi = 200
    Nboxe = len(dsCS.nboxe)
    
    if VAR in ['U','V','W']:
        var = dsCS[VAR+'T']
    
    
    
    #fig, ax = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi)


# METHOD 1:
# fit une ellipse pour chaque structure (updrafts ?)

# à une altitude donnée
# pour chaque structure coherente 
#       possible que downdrafts ca marche pas super
#       attention près de la surface

def find_objects():
     """
     """
     # scipy.ndimage.label(input, structure=None, output=None)


def find_center_of_mass():
     """
     """
    # scipy.ndimage.center_of_mass(input, labels=None, index=None)
    # chercher le barycentre pour tous les index.

def reduce_to_labeled_object():
    """
    """
    # remove the rest of the array, keep only the labeled zone, with a 1 halo of zeros

def step_fit_ellipse():
    """
    """
    # fit ellipse on the smaller array



# METHOD 2:
# Req = A / Npixels

def equivalent_radius_object(mask,resolutionH):
    """
    This function computes the area of the mask for one object, then assume a circular shape to get an equivalent radius
    
    INPUTS:
        - mask: xarray DataArray 3D (level,coord_y,coord_x) field with 1 where a structure is present, 0 else.
        - resolutionH: horizontal resolution (m)

    OUTPUTS:
        - a scalar, the equivalent radius for the 'mask'
    
    """
    # here i will sum on dim that are not 'level'
    dims = mask.dims
    dimsum = [dim for dim in dims if dim!='level']
    A = mask.sum(dimsum)*resolutionH**2
    R = np.sqrt(A/np.pi)
    return R

def R_equivalent_for_all_objects(dsCS,dsmean,resolutionH,path_save):
    """
    This function computes a equivalent radius for each coherent structure, at every altitude.

    INPUTS:
        -

    OUTPUTS:
        - a plot of R(z) for each structures

    """
    dsCS = dsCS.isel(time=-1) # this can be removed when time is taken into account
    global_mask = dsCS.global_mask
    Z = dsCS.level
    Nboxe = len(dsCS.nboxe)

    # is_up1 = xr.where( dsCS.global_mask==1,1,0 )
    # is_ss1 = xr.where( dsCS.global_mask==2,1,0 )
    # is_up2 = xr.where( dsCS.global_mask==3,1,0 )
    # is_ss2 = xr.where( dsCS.global_mask==4,1,0 )
    # is_down = xr.where( dsCS.global_mask==5,1,0 )
    is_up1 = global_mask.where(global_mask==1,0)
    is_ss1 = global_mask.where(global_mask==2,0)/2
    is_up2 = global_mask.where(global_mask==3,0)/3
    is_ss2 = global_mask.where(global_mask==4,0)/4
    is_down = global_mask.where(global_mask==5,0)/5
    list_mask = [is_up1,is_ss1,is_up2,is_ss2,is_down]
    colors_obj = ['r','purple','orange','pink','g']
    labels_obj = ['up1','ss1','up2','ss2','down']
    Req = np.zeros((Nboxe,len(list_mask),len(Z)))
    labels = xr.zeros_like(is_up1.isel(nboxe=0,level=0))

    # unique object identification with scipy
    for k,boxe in enumerate(dsCS.nboxe[:1]):
        print('I work on boxe '+str(boxe.values))
        for imask in range(len(list_mask)): #  1 len(list_mask)
            print(' '+labels_obj[imask])
            for iz in range(len(Z)): # len(Z) 5
                labels.data, Nlabel = SPlabel(list_mask[imask].isel(nboxe=k,level=iz))
                if iz==50:
                    print('Z= ',Z[iz].values)
                    fig, ax = plt.subplots(1,1,figsize = (3,3),constrained_layout=True,dpi=200)
                    ax.pcolormesh(labels)
                    plt.show()
                for i in range(Nlabel):
                    Req[k,imask,iz] += equivalent_radius_object(xr.where( labels==i,1,0),resolutionH)
                if Nlabel>0:
                    Req[k,imask,iz] = Req[k,imask,iz] / Nlabel    
                else:
                     Req[k,imask,iz] = 0
                print('         Z = ',Z[iz].values,' at this level, '+str(Nlabel)+' obj, Req = ',Req[k,imask,iz])
                

    # for k,boxe in enumerate(dsCS.nboxe):
    #     for imask in range(len(list_mask)):
    #         Req[k,imask,:] = equivalent_radius_object(list_mask[imask][k],resolutionH)

    fig, ax = plt.subplots(1,Nboxe,figsize = (3*Nboxe,5),constrained_layout=True,dpi=200)
    for k,boxe in enumerate(dsCS.nboxe):
        zi = get_ABLH(Z,dsmean.THTvm.isel(nboxe=k)) # Boundary layer height
        for i,obj in enumerate(labels_obj):
            ax[k].plot(Req[k,i,:],Z/zi,c=colors_obj[i],label=obj)
            ax[k].set_xlabel('Req (m)')
            ax[k].set_title('boxe '+str(boxe.values),loc='right')
            ax[k].set_ylim([0,1.2])
            #ax[k].set_xlim([0,100])
    ax[0].set_ylabel('Z (m)') 
    fig.savefig(path_save + 'Requivalent_Z_allobjects.png')

# A faire: fonction qui regarde l'aire de chaque structure, plot une distribution de Req en fonction de Z.


