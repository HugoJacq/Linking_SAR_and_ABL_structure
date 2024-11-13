import matplotlib as mpl

from module_tools import *

"""
To be used with analyse.py
"""

def plot_energy_spectrum(dsO,altZ_convergence,liste_X,B_KPSD,A,Kcoeff,path_save):
    """
    This function plots the energy spectrum at the East side of the domain of dsO.
    This is to know  the distance at which the turbulence is at the scale of the Son domain.

    METHOD:
        At every X position, the Y dimension is split into 'Split' parts.
        A spectrum is computed on this smaller domain, for each X position and time step.
        At the end, spectrum are averaged over theses parts and in time.

    INPUT:
        - dsO : xarray dataset interpolated at mass point
        - altZ_convergence : altitude where to plot the spectrum (m)
        - liste_X : a liste of X position away from the East border (km)
        - B_KPSD : switch to plot k*PSD(k) or just PSD(k)
        - A : y=log(x)**(-coeff)+log(A) for Kolmogorov Law in inertial subrange
        - Kcoeff : -5/3 in inertial range
        - path_save : where to save figures
    OUTPUT:
        - A figure with spectrum at different X positions
        - A figure with 
    """
    Split = 4 # how many Y split ?
    
    X,Y,Z = dsO.ni.values,dsO.nj.values,dsO.level.values
    res = Y[1]-Y[0]
    indz_convergence = nearest(Z,altZ_convergence)

    Ny = len(Y)
    quarter = Ny//Split
    cmap_str = 'jet'
    cmap = mpl.colormaps[cmap_str]
    colors = cmap(np.linspace(0,1,len(liste_X)))
    
    U = dsO.UT
    V = dsO.VT
    W = dsO.WT
    E = 0.5*(U**2+V**2+W**2).isel(level=indz_convergence)
    PSD_E = np.zeros((Split,len(liste_X),len(dsO.time),quarter//2-1))

    # Computing PSD
    for s in range(Split):
        for k,atX in enumerate(liste_X):
            indx = nearest(X,X[0]+atX*1000)
            for t in range(len(dsO.time)):
                borne1 = s*quarter
                borne2 = (s+1)*quarter
                f,PSD_E[s,k,t,:] = detrended_PSD(Y[borne1:borne2], E.isel(time=t,ni=indx,nj=slice(borne1,borne2)).values )

    if B_KPSD:
        coeff_f = f
        nameY = r'$\frac{1}{\lambda_2}$.$F^{22}_{E}$ (m$^2$.s$^{-2}$)'
        addsave = 'KPSD'
        ylim = [1e-5,1e2]
    else:
        coeff_f = 1
        nameY = r'$F^{22}_{E}$ (m.s$^{-2}$)' 
        addsave = 'PSD'
        ylim = [1e-4,1e7]

    # 1 plot at each X position,
    #   to show how the average of 'Split' number of spectrum
    for k,atX in enumerate(liste_X):
        fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=200)
        ax.vlines(1/(4*res),0,10000,colors='Grey',linestyles='--')
        for s in range(Split):
            for t in range(len(dsO.time)):
                ax.plot(f,coeff_f*PSD_E[s,k,t,:],c='b',alpha=0.1)
        ax.plot(f,coeff_f*PSD_E.mean(axis=(0,2))[k,:],c=colors[k],alpha=1)
        ax.set_xlabel(r'$\frac{1}{\lambda_3'+r'}$ (m$^{-1}$)')
        ax.plot(f,coeff_f*(f)**(Kcoeff)*A,c='k',ls='--')
        ax.set_ylim(ylim)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_ylabel(nameY)
       
        ax.set_title('E spectrum at X='+str(atX)+'km')
        ax.grid()
        fig.savefig(path_save + 'E_'+addsave+'_'+str(Split)+'avg_atX'+str(atX)+'km_atZ'+str(altZ_convergence)+'m.png')

    # 1 plot with only averaged spectrum (time and Y split)
    fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=200)
    ax.vlines(1/(4*res),0,10000,colors='Grey',linestyles='--')
    for k in range(len(liste_X)):
        ax.plot(f,coeff_f*PSD_E.mean(axis=(0,2))[k,:],c=colors[k],label='atX='+str(liste_X[k])+'km')
    ax.set_xlabel(r'$\frac{1}{\lambda_3'+r'}$ (m$^{-1}$)')
    ax.plot(f,coeff_f*(f)**(Kcoeff)*A,c='k',ls='--')
    ax.set_ylim(ylim)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylabel(nameY)
    ax.set_title('E spectrum')
    ax.grid()
    ax.legend()
    fig.savefig(path_save + 'E_'+addsave+'_'+str(Split)+'avg_atZ'+str(altZ_convergence)+'m.png')
    return None