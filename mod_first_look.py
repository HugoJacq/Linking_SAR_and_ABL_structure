from module_tools import *

def Wind_10m_in_boxes_vs_SAR(dsCS,dsSAR,indt,d_boxes,path_save):
    """
    This function show the modulus of z=10m, wind in the boxes
    
    INPUT:
        - dsCS: xarray dataset with variables in boxes
        - dsSAR: xarray dataset with SAR data (from prepare_obs.py)
        - indt : chosen time step
        - d_boxes: boxes localisation for LES and SAR
        - path_save : where to save the figures
    OUTPUT:
        - 1 figure with 10m wind in the boxes
    """
    figsize = (5,7)
    dpi = 200
    cmap = 'Greys_r'
    vmin,vmax = 3,7 # wind
    sigmin,sigmax = 0.3,0.6 # sigma0
    indz = nearest(dsCS.level.values,10)
    ds = dsCS.isel(level=indz,time=indt)
    Nboxe = len(ds.nboxe)
    M = np.sqrt(ds.U**2+ds.V**2+ds.W**2)



    fig, ax = plt.subplots(Nboxe,2,figsize = figsize,constrained_layout=True,dpi=dpi)
    # LES
    for kboxe in range(Nboxe):
        s = ax[kboxe,0].pcolormesh(ds.X[kboxe]/1000,
                                   ds.Y[kboxe]/1000,
                                   M.isel(nboxe=kboxe),cmap=cmap,vmin=vmin,vmax=vmax)
        ax[kboxe,0].set_ylabel('Y (km)')
        ax[kboxe,0].set_aspect(1)
        ax[kboxe,0].xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax[kboxe,0].yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax[0,0].set_title('LES')
    ax[Nboxe-1,0].set_xlabel('X (km)')
    plt.colorbar(s,ax = ax[Nboxe-1,0],label=r'M$_{10}$ m/s',
                  orientation='horizontal',location='bottom',aspect=50)
    # SAR
    ax[0,1].set_title('SAR')
    Lx = d_boxes['SAR']['Lx'] # width of boxe in km
    Ly = d_boxes['SAR']['Ly'] # height of boxe in km
    if 'SAR' in dsSAR.keys(): # SAR from OVL
        Lx = Lx/DegLon # width of boxe in °Lon
        Ly = Ly/DegLat # height of boxe in °Lat
        for kboxe,boxe in enumerate(d_boxes['SAR']['boxes']):
            O = d_boxes['SAR']['boxes'][boxe]['O']
            ds2 = dsSAR.sel(lon=slice(O[0],O[0]+Lx),lat=slice(O[1]+Ly,O[1]))
            s = ax[kboxe,1].pcolormesh(ds2.lon,ds2.lat,np.ma.masked_where(ds2.MASK,ds2.SAR),cmap=cmap,vmin=sigmin,vmax=sigmax)
            ax[kboxe,1].set_aspect(DegLat/DegLon)
            ax[kboxe,1].set_ylabel('LAT')
            ax[kboxe,1].xaxis.set_major_locator(ticker.MultipleLocator(0.1))
            ax[kboxe,1].yaxis.set_major_locator(ticker.MultipleLocator(0.1))
            ax[kboxe,1].yaxis.tick_right()
            ax[kboxe,1].yaxis.set_label_position("right")

    elif 'sigma0_detrend' in dsSAR.keys(): # SAR from Ifremer
        res = dsSAR.sampleSpacing.values # and this is = to lineSpacing
        for kboxe,boxe in enumerate(d_boxes['SAR']['boxes']):
            O = d_boxes['SAR']['boxes'][boxe]['O']
            sigmin,sigmax = 0,0.1
            sig0 = dsSAR.sigma0_detrend
            lon = dsSAR.longitude
            lat = dsSAR.latitude
            # lets find line and sample index at origin of boxe (bottom left corner)
            indlineO,indsample0 = find_indx_indy_from_2D_LatLon(lat,lon,O)
            if indlineO==None or indsample0==None:
                continue
            indlineRight,indsampleRight =  int(indlineO+Lx*1000//res), int(indsample0+Lx*1000//res)  

            sig0 = sig0.isel(line=slice(indlineO,indlineRight),sample=slice(indsample0,indsampleRight))
            x,y = np.arange(0,len(sig0.line)*res,res),np.arange(0,len(sig0.sample)*res,res)
            s = ax[kboxe,1].pcolormesh(x/1000,y/1000,sig0,cmap=cmap,vmin=sigmin,vmax=sigmax) # sig0.sample,sig0.line,
            ax[kboxe,1].set_aspect(1)
            ax[kboxe,1].xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax[kboxe,1].yaxis.set_major_locator(ticker.MultipleLocator(5))
            ax[kboxe,1].yaxis.tick_right()
            ax[kboxe,1].yaxis.set_label_position("right")
    
    
    plt.colorbar(s,ax = ax[Nboxe-1,1],label=r'$\sigma_0$',
                 orientation='horizontal',location='bottom',aspect=50)
    ax[Nboxe-1,1].set_xlabel('X (km)')



    fig.savefig(path_save + '10m_wind_vs_SAR.png')

def Where_boxes(dsO,dsSAR,dsSST,indt,d_boxes,SAR_SIZE,path_save):
    """
    This function plots the location of boxes defined in d_boxes

    INPUT:
        - dsO: xarray dataset with variables for full domain
        - dsSAR: xarray dataset with SAR data (from prepare_obs.py)
        - dsSST : xarray dataset with ODYSEA L4 product
        - indt : chosen time step
        - d_boxes: boxes localisation for LES and SAR
        - path_save : where to save the figures
    OUTPUT:
        - 1 figure, M10 and sigma0 in the backdrop with rectangles where the boxes are
    """
    figsize = (7,5)
    dpi = 200
    cmap = 'Greys_r'
    vmin,vmax = 3,7 # wind
    sigmin,sigmax = 0.3,0.6 # sigma0

    # to plot SST with LES :
    Lx_s = 64 # km
    Ly_s = 90 # km
    O_s = (24.5,-36)
    SSTvmin,SSTvmax = 295,298 # K
    Nlevel = 5
    Levels = np.linspace(SSTvmin,SSTvmax,Nlevel)
    dT = Levels[1] - Levels[0]
    cmap2 = plt.cm.plasma  # define the colormap
    bounds = np.arange(Levels[0] - dT/2,Levels[-1] + dT/2, dT) # define the bins
    norm = mpl.colors.BoundaryNorm(bounds, cmap2.N) # normalize
    XSST = dsSST.lon.sel(lon=slice(O_s[0],O_s [0] + Lx_s/DegLon+0.02))
    XSST = (XSST - XSST[0])*DegLon + dsO.ni[0]/1000                         
    YSST = dsSST.lat.sel(lat=slice(O_s[1],O_s [1] + Ly_s/DegLat+0.02) )
    YSST = (YSST - YSST[0])*DegLat + dsO.nj[0]/1000 
    SST = dsSST.analysed_sst.isel(time=0).sel(lon=slice(O_s[0],O_s [0] + Lx_s/DegLon + 0.02),
                                              lat=slice(O_s[1],O_s [1] + Ly_s/DegLat + 0.02))
    # to plot SAR
    if 'SAR' in dsSAR.keys():
        azimuth = 0 # dont rotate boxes
        SARmin,SARmax=0.3,0.6
        XSAR = dsSAR.lon
        YSAR = dsSAR.lat
        MASK = dsSAR.MASK
        sig0 = dsSAR.SAR
        SAR = np.ma.masked_where(MASK,sig0)
        cbarTicks = 1
    elif 'sigma0_detrend' in dsSAR.keys():
        azimuth = 90 - azimuth_sat # degree, approximately
        SARmin,SARmax=0.,0.07
        SAR = dsSAR.sigma0_detrend
        XSAR = dsSAR.longitude
        YSAR = dsSAR.latitude
        cbarTicks = 0.02
    SARxlim = [24.2,25.2]
    SARylim = [-36.2,-35.1]
    # to plot LES
    indz = nearest(dsO.level.values,10)
    UT = dsO.UT.isel(time=indt,level=indz)
    VT = dsO.VT.isel(time=indt,level=indz)
    M = np.sqrt(UT**2+VT**2)
    # PLOTTING
    fig, ax = plt.subplots(1,2,figsize = figsize,constrained_layout=True,dpi=dpi)
    # LES
    s = ax[0].pcolormesh(dsO.ni/1000,dsO.nj/1000,M,cmap=cmap,vmin=vmin,vmax=vmax)
    b = plt.colorbar(s,ax=ax[0],label=r'M$_{10}$ (m/s)',orientation='horizontal',aspect=50)
    b.locator = ticker.MultipleLocator(1)
    b.update_ticks()
    CS = ax[0].contour(XSST,YSST,SST,cmap=cmap2,levels=Levels)
    ax[0].clabel(CS, inline=True, fontsize=10)
    ax[0].set_xlabel('X (km)')
    ax[0].set_ylabel('Y (km)')
    ax[0].set_title('LES')
    ax[0].set_aspect(1)
    # SAR
    s = ax[1].pcolormesh(XSAR,YSAR,SAR,cmap=cmap,vmin=SARmin,vmax=SARmax)
    b = plt.colorbar(s,ax=ax[1],label=r'$\sigma_0$',orientation='horizontal',aspect=50)
    CS = ax[1].contour(dsSST.lon,dsSST.lat,dsSST.analysed_sst.isel(time=0),cmap=cmap2,levels=Levels)
    ax[1].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax[1].yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    b.locator = ticker.MultipleLocator(cbarTicks)
    b.update_ticks()
    ax[1].set_aspect(1) # DegLat/DegLon
    ax[1].set_xlim(SARxlim)
    ax[1].set_ylim(SARylim)
    #ax[1].clabel(CS, inline=True, fontsize=10) #, fmt=contour_fmt
    ax[1].set_xlabel('lon')
    ax[1].set_ylabel('lat')
    ax[1].set_title('SAR')
    # draw boxes
    for boxe in d_boxes['LES']['boxes'].keys():
        left = d_boxes['LES']['boxes'][boxe]['O'][0]/1000
        width = d_boxes['LES']['Lx']/1000
        right = left + width
        bott = d_boxes['LES']['boxes'][boxe]['O'][1]/1000
        height = d_boxes['LES']['Ly']/1000
        top = bott + height
        rectangleLES = mpl.patches.Rectangle(
                (left,bott),width,height,
                 edgecolor='k',fill=False,lw=2)

        ax[0].add_patch(rectangleLES)
        ax[0].annotate(boxe,(0.5*(left+right), 0.5*(bott+top)),
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=10, color='black',) # transform=fig.transFigure

    for boxe in d_boxes['SAR']['boxes'].keys():
        if boxe in ['1','2','3'] and SAR_SIZE not in ['small_only']:
            lw = 2
        else:
            lw = 1
        left = d_boxes['SAR']['boxes'][boxe]['O'][0]
        width = d_boxes['SAR']['Lx']/DegLat
        right = left + width
        bott = d_boxes['SAR']['boxes'][boxe]['O'][1]
        height = d_boxes['SAR']['Ly']/DegLat
        top = bott + height
        xpos = 0.5*(left+right)
        ypos = 0.5*(bott+top)

        if 'sigma0_detrend' in dsSAR.keys():
            top = bott
            bott = top - height
            xpos = xpos -0.2*width
            ypos = ypos +0.18*height
        
        rectangleSAR = mpl.patches.Rectangle(
            (left,bott),width,height,angle=azimuth, rotation_point=(left,top),
                edgecolor='k',fill=False,lw=lw)
        ax[1].add_patch(rectangleSAR)
    

        ax[1].annotate(boxe,(xpos, ypos),
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=10, color='black',)




    fig.savefig(path_save + 'Boxe_location_SAR_LES.png')

def SAR_boxes_generator(dsSAR,d_boxe,size,Nx=-1,Ny=-1):
    """
    this function generates a dictionnary that specifies location and dimension of analysis boxes for SAR data

    d_boxe['SAR'] = {'boxes':'k':{'O':(Ox,Oy)} }

    The width and heigth of each boxe are given by d_boxe['SAR']['Lx'] and d_boxe['SAR']['Ly']

    METHOD:
        We start from bottom left corner, checking if there is data there.
        Then, depending on 'size', we look for boxes of width Lx and height Ly

    INPUT:
        - dsSAR : xarray dataset with SAR data (from Ifremer)
        - d_boxe : dict with boxe location and dimensions
        - Nx : number of boxes in along track direction (-1 for the maximum)
        - Ny : number of boxes in across track direction (-1 for the maximum)
    OUTPUT:
        - d_boxe with added boxes dimensions/position
    """
    sig0 = dsSAR.sigma0_detrend
    Lx = d_boxe['Lx'] # km, along 'sample'
    Ly = d_boxe['Ly'] # km, along 'line'
    nameBoxe = len(d_boxe['boxes']) + 1

    if size=='all': # add Lx*Ly boxes all over SAR data
        indx = 0
        indy = 10
    elif size=='small': # reducing dataset
        indx = 1000 
        indy = 500 
    elif size=='minimal': # do not add automaticaly detected boxes
        return d_boxe
    elif size=='small_only': # remove user given boxes
        nameBoxe = 1
        indx = 1000 
        indy = 500 
    else:
        raise Exception('The size '+size+' of the batch of boxe is not recognized')

    if Nx==-1:
        Nx = int( len(dsSAR.sample[indx:]) // (Lx*1000/SAR_res) )
    if Ny==-1:
        Ny = int( len(dsSAR.line[indy:]) // (Ly*1000/SAR_res) )
    
    masked = sig0.where( sig0 > 0 )

    print('SAR: number of boxe in sample,line directions:',Nx,',',Ny,'. Total=',Nx*Ny)
    offsetx = 10
    offsety = 10
    stepx = Lx*1000 // SAR_res
    stepy = Ly*1000 // SAR_res

    

    # looking for first non nan value along 'sample'
    while np.isnan(masked.isel(line=indy,sample=indx)):
        x = 0
        while np.isnan(masked.isel(line=indy,sample=x)) and x<len(dsSAR.sample)-1:
            x = x+1
        indx = x
        if x == len(dsSAR.sample)-1:
            indx += 1

    indx,indy = int(indx+offsetx),int(indy+offsety) # let's get far from boundaries
    ix,iy = 0,0
    # loop over Nx,Ny boxes
    for kx in range(indx,len(dsSAR.sample),stepx):
        if (len(dsSAR.sample)-kx) >= stepx:
            for ky in range(indy,len(dsSAR.line),stepy):
                if  (len(dsSAR.line)-ky) >= stepy:
                    origin = ( float(dsSAR.longitude[ky,kx].values),float(dsSAR.latitude[ky,kx].values) )
                    # add the new origin
                    d_boxe['boxes'][str(nameBoxe)] = {'O':origin}
                    nameBoxe += 1
    return d_boxe




