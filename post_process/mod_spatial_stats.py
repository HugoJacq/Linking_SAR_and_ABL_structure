from module_tools import *
import dask.array as da
import dask
import pathlib
import cv2
from scipy.linalg import inv, eigh, solve
from scipy import signal

"""
To be used with analyse.py

thank you Pierre Etienne Brilouet (pierre-etienne.brilouet@ird.fr) for :
_find_max_eigval, _find_max_eigvec, polar_axis, field_2D_polar, fit_ellipse, get_parameters
"""


# TOOL FUNCTIONS
def indicator(d,t):
    if d<t:
        return 1
    else:
        return 0
    
def Ripley_K_L(A,data,t):
    """
    This function return Ripley's K and L values for a given dataset

    METHOD: 
        following https://en.wikipedia.org/wiki/Spatial_descriptive_statistics

    INPUT:
        - A: Area of the region (m2)
        - data: dataset of size N
        - t: search radius
    OUTPUT:
        - K(t)
        - L(t)
    """

    """
    Pas encore clair : il faut spécifier des points particuliers !
            -> gros gradient de sigma0 / vent ?

    """  

def polar_axis(data, n_theta=-1):
    '''    
    This function returns 2 numpy arrays characterizing the polar coordinates 
    of the output field of the function field_2D_polar:
 
    axis_t, axis_r, size = polar_axis(data, n_theta=-1)

    Input arguments :
        
        data  : image to analyze
        n_theta : nb points desired in the theta direction if =-1: auto-detect (as R_max*pi) 
        
    Output arguments :
        
        axis_t : vector of angles
        axis_r : vector of radius
        size : shape of field in polar coordinates
        
    Usage example :
        
        image=np.random.randn((1000,1000))
        axis_t, axis_r, size = polar_axis(image, n_theta=-1)

   
    ##
    N.B Garnier (ENS Lyon) and C. Granero-Belinchon (IMT- Atlantique), September 2021
 
    '''  

    Nx,Ny = data.shape
    R_max = min(Nx//2,Ny//2)
    N_theta = n_theta if n_theta>0 else (int)(R_max*np.pi)
    size  = (R_max+1, N_theta)      # N_theta pts in [0,2*pi[ and R (tau) in [O, R_max]
    axis_t= np.arange(N_theta)/N_theta*2*np.pi - np.pi/2
    axis_r= np.arange(R_max+1)
    return axis_t, axis_r, size
  
def field_2D_polar(H, n_theta=-1):    
    '''    
    This function returns the cartesian field H in polar coordinates. The center of the field is used as origin.
    First dimension is the angle in [0, 2*pi[, second dimension is the radius in [0, R_max]
 
    H_polar = field_2D_polar(H, n_theta=-1)

    Input arguments :
        
        H  : field to analyze
        n_theta : nb points desired in the theta direction if =-1: auto-detect (as R_max*pi) 
        
    Output arguments :
        
        H_polar : H field in polar coordinates
        
    Usage example :
        
        image=np.random.randn((1000,1000))
        image_polar = field_2D_polar(image, n_theta=-1)

   
    ##
    N.B Garnier (ENS Lyon) and C. Granero-Belinchon (IMT- Atlantique), September 2021
 
    '''  

    Nx,Ny = H.shape
    H_t= np.transpose(H) # cv2 works on images, which are transposed of matrices
    origin= (Nx//2, Ny//2) 
    axis_t, axis_r, size = polar_axis(H, n_theta) # overkill here...
    R_max = size[0]-1
    return cv2.warpPolar(H_t, size, origin, R_max, cv2.WARP_POLAR_LINEAR+cv2.INTER_CUBIC) # possible option: cv2.WARP_FILL_OUTLIERS

def _find_max_eigval(S):
    """
    Finds the biggest generalized eigen value of the system
        
    S * x = l * C * x
        
    where
    
    ::
        
        C = | 0  0 2 |
            | 0 -1 0 |
            | 2  0 1 |
 
    Parameters:
    -----------    
    S : 3x3 matrix
    
    Returns:
    --------
    the highest eigen value
    """
 
    a = S[0,0]
    b = S[0,1]
    c = S[0,2]
    d = S[1,1]
    e = S[1,2]
    f = S[2,2]
 
    # computes the coefficients of the caracteristique polynomial
    # det(S - x * C) = 0
    # Since the matrix is 3x3 we have a 3rd degree polynomial
    # _a * x**3 + _b * x**2 + _c * x + _d
    _a = -4
    _b = 4 * (c - d)
    _c = a * f - 4 * b * e + 4 * c * d - c * c
    _d = a * d * f - b * b * f - a * e * e + 2 * b * c * e - c  * c * d
 
    # computes the roots of the polynomial
    # there must be 2 negative roots and one
    # positive, i.e. the biggest one.
    x2, x1, x0 = sorted(np.roots([_a, _b, _c, _d] ))
    return x0

def _find_max_eigvec(S):
    """
    Computes the positive eigen value and the corresponding
    eigen vector of the system:
        
        S * x = l * C * x
    
        where
        ::
        
            C = | 0  0 2 |
                | 0 -1 0 |
                | 2  0 1 |
                
    Parameters:
    -----------    
    S : 3x3 matrix
    
    Returns:
    --------
        (l, u)
    
    l : float
        the positive eigen value
    
    u : the corresponding eigen vector
"""
 
    l = _find_max_eigval(S)
 
    a11 = S[0,0]
    a12 = S[0,1]
    a13 = S[0,2]
    a22 = S[1,1]
    a23 = S[1,2]
 
    u = np.array([
        a12 * a23 - (a13  - 2*l) * (a22 + l),
        a12 * (a13  - 2*l) - a23 * a11,
        a11 * (a22 + l) - a12 * a12
    ])
 
    c = 4 * u[0] * u[2] - u[1] * u[1]
 
    return l, u/np.sqrt(c)

def fit_ellipse(X):
  """ Fit an ellipse.
    Computes the best least squares parameters of an ellipse  expressed as:
        a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0
    Parameters
    ----------
    X : N x 2 array
        an array of N 2d points.
    
    Returns:
    --------
    an array containing the parameters:
        [ a , b, c, d, e, f]
  """
  x = X[:,0]
  y = X[:,1]
 
  # building the design matrix
  D = np.vstack([ x*x, x*y, y*y, x, y, np.ones(X.shape[0])]).T
  S = np.dot(D.T, D)
 
  S11 = S[:3][:,:3]
  S12 = S[:3][:,3:]
  S22 = S[3:][:,3:]
 
  S22_inv = inv(S22)
  S22_inv_S21 = np.dot(inv(S22), S12.T)
 
  Sc =  S11 - np.dot(S12, S22_inv_S21)
  l, a = _find_max_eigvec(Sc)

  b = - np.dot(S22_inv_S21, a)
 
  return np.hstack([a,b])

def get_parameters(x):
  """
    Computes 'natural' parameters of an ellipse given the parameters
    of the canonical equation:
        
        a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0
    
    Parameters:
    -----------
    x : array_like
        An array of 6 elements corresponding to the coefficients of the
        canonical equation (see above)
    
    Returns:
    --------
        tuple (rx, ry), (xc, yc), alpha
        
    (rx, ry) : tuple
        Radii of the major and minor axes
    
    (xc, yc) : tuple
        coordinates of the center
    
    alpha : float
        angle between the x axis and the major axis
    
    :Note:
 
        Computed the parameters of the ellipse when it is expressed as:
            
            x'^2/rx^2 + y'/ry^2 = 1
            
        where x' and y' correpsond to the rotated coordinates:    
            
            x' =  cos(alpha)(x-xc) + sin(alpha)(y-yc)
            y' = -sin(alpha)(x-xc) + cos(alpha)(y-yc)
        
        Which can be put in matrix form as
        
            (X-Xc)' R D R' (X-Xc) = 1
        
        where
        ::
            
              X = [x y] and Xc = [xc yc]
              
              R = [ cos(alpha) -sin(alpha)]
                  [ sin(alpha) cos(alpha) ]
                  
              D = [ 1/rx^2           0    ]
                  [    0          1/ry^2  ]
                    
        Parameters are given as the parameter of the conic:
            
            a * x^2 + b * x * y + c * y^2 + d * x + e * y + f = 0
        
        In matrix form we have:
            
            X' A X + B' X + f = 0
            
            where
            ::
                  X = [ x  y ]'
            
                  A = [ a b/2]
                      [b/2 c ]
                      
                  B = [ d  e ]'
              
        Any ellipse can be written as:
            
                (X - Xc)' A (X - Xc)  = r^2
                
        which develops in:
            
            X'A X - 2 * Xc' A X + Xc' A Xc - r^2 = 0
            
            So we have:
                
                B = - 2 * A Xc
                
            and
            
                f = Xc' A Xc - r^2
                
            and thus:
                
                Xc = -1/2 * A^(-1) B
                r^2 = Xc' A Xc - f
            
            We also see that
            
            1/r^2 * (X - Xc)' A (X - Xc) = (X-Xc)' R D R' (X-Xc) = 1
            
            By performing eigen decomposition on A = U L U', we obtain
            
                R = U
                
                and 
                
                lx / r^2 = 1/rx^2
                ly / r^2 = 1/ry^2
                
            hence
            
                rx^2 = r^2 / lx
                ry^2 = r^2 / ly
                
            the angle alpha is finally determined using
            ::
                
                U = | u11,  u12 | = | cos(alpha)  sin(alpha)|
                    |-u12,  u22 |   | -sin(alpha) cos(alpha)|
                alpha = sign(u12) * arccos(u11)
  """
  a,b,c,d,e,f = x

  A = np.array([
      [ a, b/2 ],
      [b/2, c  ]
  ])    
 
  B = np.array([d,e])    

  w,u = eigh(A)

  Xc = solve(-2*A,B)
  r2 = -0.5 * np.inner(Xc,B) - f
 
  rr2 = r2 / w
 
  alpha = np.arccos(u[0,0])
  if alpha > np.pi/2:
    alpha = alpha - np.pi
 
  alpha *= np.sign(u[0,1])
 
  return tuple(np.sqrt(rr2)), tuple(Xc), alpha

def Nth_structure_function_2D(Array,resolution,N=2):
    """
    This function computes the Nth order structure function of the user-provided Array of the form A(x,y)

        S_n(lx,ly) = < ( A(x+lx,y+ly) - A(x,y) )**N >

    METHOD:
        following 
        * Granero Belinchon, C., Roux, S. G., Garnier, N. B., Tandeo, P., Chapron, B., & Mouche, A. (2022). 
          "Two-dimensional structure functions for characterizing convective rolls in the marine atmospheric boundary layer from Sentinel-1 SAR images"
          Remote Sensing Letters, https://doi.org/10.1080/2150704X.2022.2112107

        see jupyter notebook : https://github.com/cgranerob/2D-Structure-Functions/blob/main/Script.ipynb  

    INPUT:
        - Array : xarray DataArray, 2D
        - N : order of the structure function
    OUTPUT:
        - lx,ly,nth_variance

    Note:
        - the average operator <> is taken as the mean over the full Array
        - dim1,dim2 are local dimensions (m)
    """
    SERIAL = False

    # getting local dimensions
    dims = Array.dims
    dim1 = Array[dims[0]]
    dim2 = Array[dims[1]]
    name1 = dim1.name
    name2 = dim2.name

    # building lx,ly
    #resx = dim1[1]-dim1[0]
    Nx = len(dim1)
    #lx = np.arange(-(Nx//2-1)*resx,(Nx//2)*resx,resx)
    lx = np.arange(-(Nx//2-1)*resolution,(Nx//2)*resolution,resolution)
    #resy = dim2[1]-dim2[0]
    Ny = len(dim2)
    #ly = np.arange(-(Ny//2-1)*resy,(Ny//2)*resy,resy)
    ly = np.arange(-(Ny//2-1)*resolution,(Ny//2)*resolution,resolution)

    # building DataArray of S_n
    coords={'lx':lx,'ly':ly}
    S_n = xr.DataArray( np.zeros((len(lx),len(ly))),
                                 dims=('lx','ly'),
                                 coords=coords)# .chunk(chunks)
    if not SERIAL:
        # parallel version
        # method : we create an array with only 4 values filled and 0 elswhere.
        #   then all of these array are summed to get S_n
        Array = Array.compute() # going with numpy, chunk({'coord_x':200,'coord_y':200})
        def shift_and_mean(i,j,Array):
            # this is a function to build a // version
            shifted = Array.roll(shifts={name1:i,name2:j})
            results = xr.DataArray( np.zeros((len(lx),len(ly))),
                                        dims=('lx','ly'),
                                        coords=coords).chunk({'lx':200,'ly':200})
            results.loc[dict(lx=i*resolution,ly=j*resolution)] = ((shifted - Array)**N ).mean()

            shifted = Array.roll(shifts={name1:i,name2:-j})
            results.loc[dict(lx=i*resolution,ly=-j*resolution)] = ((shifted - Array)**N ).mean()

            shifted = Array.roll(shifts={name1:-i,name2:-j})
            results.loc[dict(lx=-i*resolution,ly=-j*resolution)] = ((shifted - Array)**N ).mean()
            
            shifted = Array.roll(shifts={name1:-i,name2:j})
            results.loc[dict(lx=-i*resolution,ly=j*resolution)] = ((shifted - Array)**N ).mean()
            return results

        arrays = [da.from_delayed(
                        dask.delayed(shift_and_mean)(i,j,Array),
                        shape=(len(lx),len(ly)),
                        dtype=float
                                )
                            for i in np.arange((len(dim1)//2)) 
                            for j in np.arange((len(dim2)//2))     ]

        computed = dask.compute(*arrays) # this is where computations are triggered
        S_n = da.stack(computed).sum(axis=0).compute() # this .compute is for the sum

    # serial version
    if SERIAL:
        print(' starting loop')
        for i in range(len(dim1)//2):
            for j in range(len(dim2)//2):
                # we make the vignette periodic
                shifted = Array.roll(shifts={name1:i,name2:j})
                S_n.loc[dict(lx=i*resolution,ly=j*resolution)] = ((shifted - Array)**N ).mean()

                shifted = Array.roll(shifts={name1:i,name2:-j})
                S_n.loc[dict(lx=i*resolution,ly=-j*resolution)] = ((shifted - Array)**N ).mean()

                shifted = Array.roll(shifts={name1:-i,name2:-j})
                S_n.loc[dict(lx=-i*resolution,ly=-j*resolution)] = ((shifted - Array)**N ).mean()
                
                shifted = Array.roll(shifts={name1:-i,name2:j})
                S_n.loc[dict(lx=-i*resolution,ly=j*resolution)] = ((shifted - Array)**N ).mean()
            if True:
                if i == len(dim1)//8:
                    print('     1/4 is done ...')
                elif i == len(dim1)//4:
                    print('     half is done ...')
                elif i == 3*(len(dim1)//8):
                    print('     3/4 is done ...')
    print('     done')
    return lx,ly,S_n

# FUNCTIONS FOR MY ANALYSIS 
def save_S_n_SAR(dsSAR,N,d_boxes,path_save):
    """
    This function saves the values of the Nth stucture function for SAR data.
    this is done for each boxe from 'd_boxes'

    INPUT:
        - dsSAR : xarray dataset with SAR data
        - N : order of structure function
        - d_boxes : where are the boxes
        - path_save : where to save file
    OUTPUT:
        - a netcdf file with S_n,lx,ly
    """
    print(' * Building S_'+str(N)+' for SAR')
    
    nameA = 'sig0'
    name_out = path_save+'S_'+str(N)+'_'+nameA+'.nc'
    IS_HERE = pathlib.Path(name_out).is_file()
    #IS_HERE = False # debug
    if IS_HERE:
        tmp = xr.open_dataset(name_out)
        if len(tmp.nboxe) != len(d_boxes['boxes']):
            print('     - File is here but with the wrong number of boxes ('+str(len(tmp.nboxe))+'instead of '+str(len(d_boxes['boxes']))+')')
            IS_HERE = False
        else:
            print('     - File is already here : '+name_out)
    else:
        lon = dsSAR.longitude
        lat = dsSAR.latitude

        res = dsSAR.sampleSpacing.values
        Nx = d_boxes['Lx']*1000//res
        Ny = d_boxes['Ly']*1000//res
        Nmax =  int(min(Nx,Ny)) # reduce this for debug
        Nboxe = len(d_boxes['boxes'])


        
        S_n = np.zeros((Nboxe,Nmax-1,Nmax-1))
        variance = np.zeros(Nboxe)
        origin = np.zeros((Nboxe,2))
        origin_sat = np.zeros((Nboxe,2))
        for k,boxe in enumerate(d_boxes['boxes']):
            print('Boxe='+boxe)
            O = d_boxes['boxes'][boxe]['O']
            indlineO,indsampleO = find_indx_indy_from_2D_LatLon(lat,lon,O)
            if indlineO==None or indsampleO==None:
                continue
            indlineRight,indsampleRight =  indlineO + int(Nx), indsampleO + int(Ny)
            
            sig0 = dsSAR.sigma0_detrend.isel(line=slice(indlineO,indlineRight),sample=slice(indsampleO,indsampleRight))
            lx,ly,S_n[k,:,:] = Nth_structure_function_2D(sig0,res,N)
            variance[k] = sig0.var()
            origin[k,0],origin[k,1] = O[0],O[1]
            origin_sat[k,0],origin_sat[k,1] = indsampleO,indlineO
        coords={'nboxe':np.arange(1,Nboxe+1),'lx': lx,'ly':ly,'dim_origin':['lon','lat'],'dim_origin_sat':['sample','line']}
        data_vars = {'S_'+str(N):(['nboxe','lx','ly'],S_n.data,{'long_name':str(N)+'th order structure function of sigma0'}),
                     'var':(['nboxe'],variance,{'long_name':'variance of sigma0_detrend'}),
                     'origin':(['nboxe','dim_origin'],origin,{'long_name':'origine of the boxe in 0:lon 1:lat'}),
                     'origin_sat':(['nboxe','dim_origin_sat'],origin_sat,{'long_name':'origin of the boxe in satellite 0:sample 1:line'}),
                }
        ds = xr.Dataset(data_vars=data_vars,coords=coords,
                    attrs={'Nmax':Nmax})
        ds.to_netcdf(path=name_out,mode='w')
        ds.close()

def save_S_n_SAR_OVL(dsSAR,N,d_boxes,path_save):
    """
    This function saves the values of the Nth stucture function for SAR data.
    this is done for each boxe from 'd_boxes'

    INPUT:
        - dsSAR : xarray dataset with SAR data
        - N : order of structure function
        - d_boxes : where are the boxes
        - path_save : where to save file
    OUTPUT:
        - a netcdf file with S_n,lx,ly
    """
    print(' * Building S_'+str(N)+' for SAR')
    
    nameA = 'sig0'
    name_out = path_save+'S_'+str(N)+'_'+nameA+'.nc'
    IS_HERE = pathlib.Path(name_out).is_file()
    #IS_HERE = False # debug
    if IS_HERE:
        print('     - File is already here : '+name_out)
    else:
        lon = dsSAR.lon
        lat = dsSAR.lat

        xres=(lon[1]-lon[0]).values*DegLon*1000
        yres=(lat[0]-lat[1]).values*DegLat*1000
        Nx = d_boxes['Lx']*1000//xres
        Ny = d_boxes['Ly']*1000//yres
        Nmax = 20 # int(min(Nx,Ny)) # reduce this for debug
        Nboxe = len(d_boxes['boxes'])


        
        S_n = np.zeros((Nboxe,Nmax-1,Nmax-1))
        

        for k,boxe in enumerate(d_boxes['boxes']):
            print('Boxe='+boxe)
            

            left = nearest(lon.values,d_boxes['boxes'][boxe]['O'][0])
            bott = nearest(lat.values,d_boxes['boxes'][boxe]['O'][1])
            
            sig0 = dsSAR.SAR.isel(lon=slice(left,left+Nmax),lat=slice(bott-Nmax,bott)) # to get square array
        #     fig, ax = plt.subplots(1,1,figsize = (5,5),constrained_layout=True,dpi=200)
        #     print(sig0[::-1,:])
        #     ax.set_title('boxe '+boxe)
        #     ax.pcolormesh(sig0[::-1,:],vmin=0.3,vmax=0.6,cmap='Greys_r')
        # plt.show()
        # raise Exception

            lx,ly,S_n[k,:,:] = Nth_structure_function_2D(sig0,N)

        coords={'nboxe':np.arange(1,Nboxe+1),'lx': lx,'ly':ly}
        data_vars = {'S_'+str(N):(['nboxe','lx','ly'],S_n.data,{'long_name':str(N)+'th order structure function of sigma0'}),
                }
        ds = xr.Dataset(data_vars=data_vars,coords=coords,
                    attrs={'Nmax':Nmax})
        ds.to_netcdf(path=name_out,mode='w')
        ds.close()

def save_S_n_LES(dsCS,VAR,atZ,N,path_save):
    """
    This function saves the values of the Nth stucture function for SAR data.
    this is done for each boxe from 'd_boxes'

    INPUT:
        - dsCS : xarray dataset with LES data
        - VAR : str, Choice of var
        - atZ : altitude (m)
        - N : order of structure function
        - path_save : where to save file
    OUTPUT:
        - a netcdf file with S_n,lx,ly
    """
    print(' * Building S_'+str(N)+' for LES')
    print('     var is : '+VAR+' at z='+str(atZ)+'m')
    
    
    indt = -1 
    nameA = VAR+str(atZ)
    name_out = path_save+'S_'+str(N)+'_'+nameA+'.nc'
    atZ = int(atZ) 

    IS_HERE = pathlib.Path(name_out).is_file()
    if IS_HERE:
        print('     - File is already here : '+name_out)
    else:
        Nboxe = len(dsCS.nboxe)
        Z = dsCS.level
        dsIN = dsCS.isel(time=indt)
        indz10 = nearest(Z.values,atZ)
        
        # variable selector
        if VAR=='M':
            A = np.sqrt( dsIN.U.isel(level=indz10)**2 + dsIN.V.isel(level=indz10)**2 )



        if False: # reduce the size for debug
            Nmax=10
        else:
            Nmax = len(A.coord_x)
        A = A.isel(coord_x=slice(0,Nmax),coord_y=slice(0,Nmax)) 
        S_n = np.zeros((Nboxe,Nmax-1,Nmax-1))
        variance = np.zeros(Nboxe)
        origin = np.zeros((Nboxe,2))
        for k in range(len(dsCS.nboxe)):
            print('Boxe=',k+1)
            var = A.isel(nboxe=k)
            lx,ly,S_n[k,:,:] = Nth_structure_function_2D(var,LES_res,N)
            variance[k] = var.var()
            origin[k,0],origin[k,1] = dsIN.X.isel(nboxe=k,coord_x=0),dsIN.Y.isel(nboxe=k,coord_y=0)
        coords={'nboxe':np.arange(1,Nboxe+1),'lx': lx,'ly':ly,'dim_origin':['X','Y']}
        data_vars = {'S_'+str(N):(['nboxe','lx','ly'],S_n.data,{'long_name':str(N)+'th order structure function of '+VAR+str(atZ)}),
                     'variance':(['nboxe'],variance.data,{'long_name':'variance of '+VAR+str(atZ)}),
                     'origin':(['nboxe','dim_origin'],origin,{'long_name':'origine of the boxe in 0:X 1:Y'}),
                }
        ds = xr.Dataset(data_vars=data_vars,coords=coords,
                    attrs={'Nmax':Nmax})
        ds.to_netcdf(path=name_out,mode='w')
        ds.close()

def Plot_S_n(dsS_n,nameA,N,path_save):
    """
    This function plots the 2D spatial structure function of Array, with a overview of the array

    INPUT:
        - dsS_n : xarray Dataset with S_n,lx,ly
        - nameA : name of the array
        - N : order of the structure function
        - path_save : where to save figures
    OUTPUT:
        - 1 saved figure at 'path_save' with : Nth order structure function
    """
    S_n = dsS_n['S_'+str(N)]
    lx = S_n.lx
    ly = S_n.ly
    Nboxe = len(S_n.nboxe)

#     bornesNth_str = {'2':{'M10':[[3,7],[0.2,0.8]],
# #                     'sig0':[[0.3,0.6],[0,0.003]],},
# #                 '3':{'M10':[[3,7],[-0.6,0.6]],
# #                     'sig0':[[0.3,0.6],[-0.0001,0.0001]],}  }

    # plot settings
    aspect = 1
    figsize = (20,5)
    figsize2 = (15,5)
    dpi = 200
    cmap = 'plasma'
    if nameA[:3] =='M10':
        coeffx = 1/1000
        coeffy = 1/1000
        vmin,vmax = 0.2,0.8
        resolution = LES_res
    elif nameA[:3]  =='sig':
        coeffx = 1/1000
        coeffy = 1/1000
        vmin,vmax = 0.0001,0.0002
        resolution = SAR_res
    else:
        coeff = 1.
        nameX,nameY = '',''
        vmin,vmax = 0,1

    # cartesian plot
    fig, ax = plt.subplots(1,Nboxe,figsize = figsize,constrained_layout=True,dpi=dpi)
    for k in range(Nboxe):
        s = ax[k].pcolormesh(lx*coeffx,ly*coeffy,S_n[k,:,:],cmap=cmap,vmin=vmin,vmax=vmax,shading='nearest')
        plt.colorbar(s,ax=ax[k],aspect=50,pad=0.01)
        ax[k].set_ylabel('ly (km)')
        ax[k].set_xlabel('lx (km)')
        ax[k].set_title(r'$S_'+str(N)+r'^{lx,ly}$'+' of '+nameA+' boxe'+str(k+1))
        XMAX = min(ax[k].get_ylim()[1],ax[k].get_xlim()[1])
        ax[k].set_xlim([-XMAX,XMAX])
        ax[k].set_ylim([-XMAX,XMAX])
        ax[k].set_aspect(aspect)
        #ax[k].xaxis.set_major_locator(ticker.MultipleLocator(2))
        #ax[k].yaxis.set_major_locator(ticker.MultipleLocator(2))
    fig.savefig(path_save+'S'+str(N)+'_lxly_'+nameA)


    # polar plot
    fig, ax = plt.subplots(1,Nboxe,figsize = figsize2,constrained_layout=True,dpi=dpi)
    for k in range(Nboxe):
        axis_t, axis_r, size = polar_axis(S_n[k,:,:].values, 180)
        angle_from_north = 90-(-axis_t*180/np.pi)+azimuth_sat
        S2new_p = field_2D_polar(S_n[k,:,:].values, 180)
        s = ax[k].pcolormesh(axis_r*resolution/1000,angle_from_north,S2new_p,cmap=cmap,shading='nearest') # ,vmin=vmin,vmax=vmax
        plt.colorbar(s,ax=ax[k],aspect=50,pad=0.01)
        ax[k].set_ylabel(r'$\theta$')
        ax[k].set_xlabel('r (km)')
        ax[k].set_title(r'$S_'+str(N)+r'^{r,\theta}$'+' of '+nameA+' boxe'+str(k+1))
        #XMAX = min(ax[k].get_ylim()[1],ax[k].get_xlim()[1])
        #ax[k].set_xlim([-XMAX,XMAX])
        #ax[k].set_ylim([-XMAX,XMAX])
        #ax[k].set_aspect(aspect)
        #ax[k].xaxis.set_major_locator(ticker.MultipleLocator(2))
        #ax[k].yaxis.set_major_locator(ticker.MultipleLocator(2))
    fig.savefig(path_save+'S'+str(N)+'_lxly_'+nameA+'_polar')

def compute_integral_scale_at_tht(var,S2_polar,r,atTheta,resolution):
    """
    This function computes the integral length scale of 'S2_polar' for a specific 'atTheta'
    It also provides the data that is need to fit an ellipse.


    INPUT:
        - var : variance of the signal
        - S2_polar : 2nd order structure function in polar coordinate (theta,r)
        - atTheta : value of theta (in rad)
        - resolution  : along r dimension, in meters
    OUTPUT:
        - float : integral length scale
        - DATA_TO_FIT : array of size 2 with data to be fitted with ellipse

    Note :  j'ai vérifié que ca marche
    """
    DATA_TO_FIT = np.zeros(2)
    autocorr = 1-1/(2*var)*S2_polar
    autocorr = autocorr/autocorr[0]
    iddd = np.where(np.sign(autocorr)!=1) # look for first negative autocorr
    if len(iddd[0])==0:
        integral_scale = np.nan
        DATA_TO_FIT[0] = 999.
        DATA_TO_FIT[1] = 999.
    else:
        id_zero1 = iddd[0][0]
        integral_scale = np.trapz(autocorr[0:id_zero1],r[0:id_zero1]*resolution)
        DATA_TO_FIT[0] = integral_scale*np.cos(-atTheta)
        DATA_TO_FIT[1] = integral_scale*np.sin(-atTheta)
    return integral_scale,DATA_TO_FIT

def S2_analysis(source,dsS2,d_boxes,path_out):
    """
    This function uses diagnostics from [1] to analyse the 2nd order structure function.

    METHOD:

        [1] : Brilouet, P. ‐E., Bouniol, D., Couvreux, F., Ayet, A., Granero‐Belinchon, C., Lothon, M., & Mouche, A. (2023). 
                "Trade Wind Boundary Layer Turbulence and Shallow Precipitating Convection: New Insights Combining SAR Images, 
                Satellite Brightness Temperature, and Airborne In Situ Measurements"
                https://doi.org/10.1029/2022GL102180

    INPUT:
        - source : string that specifie the source ('LES' or 'SAR')
        - S2_polar : xarray dataset of 2nd order structure function in cartesian coordinates (lx,ly)
        - d_boxes : dict with boxe location and width and height
        - path_out : where to save the texte file
    OUTPUT:
        A text file with :
        - L_E_MIN,L_E_MAX : min and max integral scales
        - SMALL_RAD_ELLIPSE : small radius of fitted ellipse
        - BIG_RAD_ELLIPSE : big radius of fitted ellipse
        - FLATNESS_ELLIPSE : flatness parameter of the fitted ellipse
        - DIR_ELLIPSE_FROM_NORTH : direction of the ellipse with respect to North
        - DIR_ROLL_FROM_NORTH : direction of rolls with respect to North
        - L_OS :
        - DELTA_AUTOCORR :
    """
    # in module_cst : azimuth_sat = -17 # degree
    
    d_boxes = d_boxes[source]
    Nboxe = len(dsS2.nboxe)
    if source=='LES':
        resolution = LES_res
    elif source=='SAR':
        resolution = SAR_res
    else:
        raise Exception('You want to use data from '+source+' but it doesnt exist')
    
    # PLAN
    # boucle sur les boxes
    #   calcul longueur intégrales
    #       j'ai Le(theta)
    #   calcul direction roll
    #   fit ellipse
    #       reduction des Le aux structures


    # initialisation de ce que je retourne : parametres ellipse, ...
    DIR_ROLL_FROM_NORTH = np.zeros(Nboxe)
    L_E_MIN,L_E_MAX = np.zeros(Nboxe),np.zeros(Nboxe)
    SMALL_RAD_ELLIPSE,BIG_RAD_ELLIPSE = np.zeros(Nboxe),np.zeros(Nboxe)
    FLATNESS_ELLIPSE = np.zeros(Nboxe)
    DIR_ELLIPSE_FROM_NORTH = np.zeros(Nboxe)
    L_OS = np.zeros(Nboxe)
    DELTA_AUTOCORR = np.zeros(Nboxe)

    for k in range(Nboxe):

        # get S2 in polar coordinates
        ds = dsS2.isel(nboxe=k)
        S2 = ds.S_2
        var = ds.var # this i need to add to the files
        theta, r, size = polar_axis(S2.values, 180)
        S2_polar = field_2D_polar(S2.values, 180)
        
        # direction is given by minimum of S2
        CUMUL_S2 = np.sum(S2_polar[:,1:],axis=1) # sum over all r, after first
        indice_para = np.argmin(CUMUL_S2)
        angle_para = theta[indice_para]*180/(np.pi)
        angle_90 = theta[indice_para]*180/(np.pi)+90
        indice_90 = np.argmin(np.abs(theta[:]*180/(np.pi)-angle_90))

        # Computing integral length scale
        #      this can be vectorized !
        Le = np.zeros(len(theta))
        DATA_TO_FIT = np.zeros((len(theta),2))
        for i,angle in enumerate(theta):
            Le[i],DATA_TO_FIT[i] = compute_integral_scale_at_tht(var,S2_polar,r,resolution)
            
        # Ellipse fitting
        #       first : remove data where no integral length scale has been computed
        DATA_TO_FIT = np.delete(DATA_TO_FIT,np.unique(np.where(DATA_TO_FIT==999.)[0]),axis=0) 
        #       then : fit the data(x,y) with custom function
        a = fit_ellipse(DATA_TO_FIT)
        rf, xcf, alpha_f = get_parameters(a)

        # Organised structure lengthscale
        if (rf[0]-rf[1])/rf[0]<0.7:
            roll_size = np.nan
            diff_tempo = np.nan
        else:  
            temp = np.nanargmin(Le)
            AUTOCORR_perp = 1-1/(2*var)*S2_polar[temp,:]
            AUTOCORR_perp = AUTOCORR_perp/AUTOCORR_perp[0]
            autocorr_smoothed = AUTOCORR_perp # maybe i need to smooth, to be looked at
            peaks, _ = signal.find_peaks(autocorr_smoothed, height=None,width=3)
            mins, _ = signal.find_peaks(-autocorr_smoothed, height=0)
            if len(peaks)==0:
                peaks = np.array([0])
            if len(mins)==0:
                mins = np.array([0])
            # critere sur diff min - max autocorr
            diff_tempo = np.nan
            roll_size = np.nan
            for ppp in peaks:
                diff_tempo = AUTOCORR_perp[ppp]-AUTOCORR_perp[mins][0]
                #if diff_tempo>=0.115:
                if (diff_tempo>=0.14) & (r[ppp]*SAR_res<1000*np.sqrt(2)*5):
                    roll_size = r[ppp]*resolution
                    break
                else:
                    diff_tempo = np.nan
                    #continue
        # -- 
        dir_roll_from_north = 90-(-theta[indice_para])*180/np.pi+azimuth_sat
        if dir_roll_from_north>180:
            dir_roll_from_north = dir_roll_from_north-180
        dir_ellipse_from_north = 90-(alpha_f)*180/np.pi+azimuth_sat
        
        # gathering results
        DIR_ROLL_FROM_NORTH[k] = dir_roll_from_north
        if np.isnan(Le).all():
            L_E_MIN[k],L_E_MAX[k] = np.nan,np.nan
        else:
            L_E_MIN[k],L_E_MAX[k] = 2*Le[np.nanargmin(Le)],2*Le[np.nanargmax(Le)]
        SMALL_RAD_ELLIPSE[k] = rf[1]
        BIG_RAD_ELLIPSE[k] = rf[0]
        FLATNESS_ELLIPSE[k] = (rf[0]-rf[1])/rf[0]
        DIR_ELLIPSE_FROM_NORTH[k] = dir_ellipse_from_north
        L_OS[k] = roll_size
        DELTA_AUTOCORR[k] = diff_tempo


    print('DIR_ROLL_FROM_NORTH')
    print(DIR_ROLL_FROM_NORTH)
    print('L_E_MIN')
    print(L_E_MIN)
    print('L_E_MAX')
    print(L_E_MAX)
    print('SMALL_RAD_ELLIPSE')
    print(SMALL_RAD_ELLIPSE)
    print('BIG_RAD_ELLIPSE')
    print(BIG_RAD_ELLIPSE)
    print('FLATNESS_ELLIPSE')
    print(FLATNESS_ELLIPSE)
    print('DIR_ELLIPSE_FROM_NORTH')
    print(DIR_ELLIPSE_FROM_NORTH)
    print('L_OS')
    print(L_OS)
    print('DELTA_AUTOCORR')
    print(DELTA_AUTOCORR)

    file_a_moi = 'SAR_clear_sky_parameters.txt' 
    fullname = path_out + file_a_moi
    file_final = open (fullname, 'w') 
    file_final.write(
    '#-------------------------------------------------------------- \n'\
    '# C00: num of boxe \n' \
    '# C01: longitude of left,bottom corner \n' \
    '# C02: lattitude of left,bottom corner \n' \
    '# C05: roll direction (degree from North) as min cumul S2  \n'\
    '# C06: 2x min of integral length scale (m)  \n'\
    '# C07: 2x max of integral length scale (m)  \n'\
    '# C08: radius of the major axe (m) \n'\
    '# C09: radius of the minor axe (m) \n'\
    '# C10: ellipse flatness (major-minor)/major   \n'\
    '# C11: Angle between the North and the major axis of the ellipse (degree from North)  \n'\
    '# C12: diff 2nd max R(r) - min R(r) for L_OS calculation  \n'\
    '# C13: Length organized structure L_OS (m)   \n'\
    '#-------------------------------------------------------------- \n')
    for index in range(Nboxe):
        file_final.write('%.0f\t %.3f\t %.3f\t %.0f\t %.0f\t %.0f\t %.0f\t %.0f\t %.2f\t %.1f\t %.2f\t %.0f\t\n' %
                        (index+1,d_boxes['boxes'][str(index)]['O'][1],d_boxes['boxes'][str(index)]['O'][1],DIR_ROLL_FROM_NORTH[index],L_E_MIN[index],L_E_MAX[index],
                        BIG_RAD_ELLIPSE[index],SMALL_RAD_ELLIPSE[index],FLATNESS_ELLIPSE[index],DIR_ELLIPSE_FROM_NORTH[index],DELTA_AUTOCORR[index],
                        L_OS[index]))

    #--- 
    file_final.close() 
    print('writing ok')

        # -> attention à l'orientation de l'array : grille lat-lon ou grille pixelX/pixelY ?

            
            