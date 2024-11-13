#import xarray_sentinel

if False:
    """
    Searching if VIIRS has better resolution than Ifremer reanalysis
    """
    basename= 'new_sst_data/VIIRS/'

    #path_file = {'VIIRS':{
                        #'2015-12-10T11:50:00.00':basename + '20151210115000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.80-v02.0-fv01.0.nc',
                        #'2015-12-10T23:10:01.00':basename + '20151210231000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.80-v02.0-fv01.0.nc',
                        #'2015-12-11T11:30:01.00':basename + '20151211113000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.80-v02.0-fv01.0.nc',
                        #'2015-12-11T13:10:00.00':basename + '20151211131000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.80-v02.0-fv01.0.nc',
                        #'2015-12-11T22:50:01.00':basename + '20151211225000-STAR-L2P_GHRSST-SSTsubskin-VIIRS_NPP-ACSPO_V2.80-v02.0-fv01.0.nc',
                        #}
                #'SAR':{},
                #}

    LATlim = [-38,-33]
    LONlim = [22,27]
    minSST = 295
    maxSST = 299
    selection = {'VIIRS':['2015-12-10T11:50:00.000000000']}


    """
    Looking at SST from VIIRS to see if we can get higher resolution.
    CCL : no, there is clouds...
    """
    source = 'VIIRS'
    for date in path_file[source].keys():
        if date in selection[source]:
            ds = xr.open_dataset(path_file['VIIRS'][date])
            lat,lon =  ds.lat,ds.lon
            sst = ds.sea_surface_temperature[0,:,:]
            valmask = ds.quality_level[0,:,:]
            
            good_quality_sst = sst.where(valmask==5)
            
            fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
            s = ax.pcolormesh(lon,lat,good_quality_sst,cmap = 'jet',vmin=minSST,vmax=maxSST)
            plt.colorbar(s,ax=ax,label='SST K',pad=0.01,aspect=50)
            ax.set_ylabel('lat')
            ax.set_xlabel('lon')
            ax.set_title('VIIRS S-NPP at '+date)
            ax.set_xlim(LONlim)
            ax.set_ylim(LATlim)
            
            
if False:
    """
    Trying to open SAR-S1 data (retrieved from copernicus)
    """
    # A essayer : utiliser le L2 OCN qui est en .nc, avec des coordonnées lat/lon + ENW + ENW direction.
    # Téléchargé depuis ici : https://browser.dataspace.copernicus.eu/?zoom=7&lat=-36.23541&lng=27.68555&themeId=OCEAN&visualizationUrl=U2FsdGVkX1%2F7%2B%2FxHRwotQv6G8VFZ7K4mBoF8xx%2FBC06JW1%2BTfxZB7dJYjEnOIUPuplD%2FJXnvSgye9YNiG55eYLNEE9iSxTV98r8A%2F33XKF617lIkElofoLwl%2FzEWvVAG&datasetId=S2_L2A_CDAS&fromTime=2024-09-24T00%3A00%3A00.000Z&toTime=2024-09-24T23%3A59%3A59.999Z&layerId=1_TRUE_COLOR

    basepath = 'new_sst_data/SAR-Sentinel1/L1/S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_4145_COG.SAFE/measurement/'
    tiffname = 's1a-iw-grd-vv-20151210t170827-20151210t170856-008982-00cdef-001-cog.tiff'
    
    
    data_sar = {'L1':('new_sst_data/SAR-Sentinel1/L1/' +
                      'S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_4145_COG.SAFE/measurement/' + 
                        's1a-iw-grd-vv-20151210t170827-20151210t170856-008982-00cdef-001-cog.tiff'),
                'L2':('new_sst_data/SAR-Sentinel1/L2/' +
                      'S1A_IW_OCN__2SDV_20151210T170827_20151210T170856_008982_00CDEF_0043.SAFE/measurement/' + 
                      's1a-iw-ocn-vv-20151210t170827-20151210t170856-008982-00CDEF-001.nc') }
    # Good signal L1: s1a-iw-grd-vv-20151210t170827-20151210t170856-008982-00cdef-001-cog.tiff 
    # Good signal L2: 
    #   /home/jacqhugo/scripts/simu_test_aggv3/new_sst_data/SAR-Sentinel1/L2/      
    #    S1A_IW_OCN__2SDV_20151210T170827_20151210T170856_008982_00CDEF_0043.SAFE/measurement/
    #    s1a-iw-ocn-vv-20151210t170827-20151210t170856-008982-00CDEF-001.nc
    
    #choice = 'L2'
    
    #if choice=='L1':
        #engine="rasterio" 
    #else:
        #engine='netcdf4'
    #ds = xr.open_dataset(data_sar[choice],engine=engine)

    product_path = 'new_sst_data/SAR-Sentinel1/L1/S1A_IW_GRDH_1SDV_20151210T170827_20151210T170856_008982_00CDEF_4145_COG.SAFE/'
    swath_group = "IW"
    swath_polarisation_group = "IW/VV"
    measurement_group = "IW/VV"
    measurement_block_slices = (slice(7000, 9000), slice(20000, 23000))
    digital_number_max = 600
    #print(xr.open_dataset(product_path, engine='sentinel-1'))
    print(xr.open_dataset(product_path, engine="sentinel-1", group=swath_group))
    print(xr.open_dataset(product_path, engine="sentinel-1", group=swath_group).attrs['subgroups'])
    #print(xr.open_dataset(product_path, engine="sentinel-1", group=swath_polarisation_group))
    
    measurement = xr.open_dataset(product_path, engine="sentinel-1", group=measurement_group, chunks=2048)
    measurement_block = measurement.measurement[measurement_block_slices]
    #print(measurement)
    #print(measurement_block)
    #abs(measurement_block).plot(y="azimuth_time", vmax=digital_number_max)
    
    gcp = xr.open_dataset(product_path, engine="sentinel-1", group=f"{swath_polarisation_group}/gcp")
    #print(gcp)
    #gcp.height.plot(y="azimuth_time")
    #gcp.plot.scatter(x="longitude", y="latitude", hue="height")
    
    orbit = xr.open_dataset(product_path, engine="sentinel-1", group=f"{swath_polarisation_group}/orbit")
    #print(orbit)
    #orbit.plot.scatter(y="azimuth_time", x="position", hue="velocity")
    
    calibration = xr.open_dataset(product_path, engine="sentinel-1", group=f"{swath_polarisation_group}/calibration")
    #print(calibration)
    #print(calibration.betaNought.mean().item(), "+-", calibration.betaNought.std().item())
    #print(calibration.dn.mean().item(), "+-", calibration.dn.std().item())
    #calibration.sigmaNought.plot(x="pixel")
    #calibration.gamma.plot(x="pixel")
    
    betaNought_block = xarray_sentinel.calibrate_amplitude(measurement_block, calibration.betaNought)
    #print(betaNought_block)
    #abs(betaNought_block).plot(y="azimuth_time", vmax=1)
    
    betaNought_block_db = xarray_sentinel.calibrate_intensity(measurement_block, calibration.betaNought, as_db=True)
    #print(betaNought_block_db)
    betaNought_block_db.plot(y="azimuth_time", vmin=-20, vmax=5)

    #noise_range = xr.open_dataset(product_path, engine="sentinel-1", group=f"{swath_polarisation_group}/noise_range")
    #print(noise_range)
    #noise_range.noiseRangeLut.plot(x="pixel")

    coordinate_conversion = xr.open_dataset(
    product_path,
    engine="sentinel-1",
    group=f"{swath_polarisation_group}/coordinate_conversion",)
    #print(coordinate_conversion)




    #fig, ax = plt.subplots(1,1,figsize = (10,10),constrained_layout=True,dpi=100)
    # s = ax.imshow(ds.band_data.isel(band=0)/255,vmin=0,vmax=1,cmap='Greys_r',origin='lower')
    # plt.colorbar(s,ax=ax)
    
    #tiff = rasterio.open(basepath+tiffname)
    #rasterio.plot.show(tiff, title = "SAR")
            
            
            
            
            
            
            
            
            
            
            
            
plt.show()
