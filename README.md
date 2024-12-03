# Linking_SAR_and_ABL_structure

<p align="center">
  <img width="600" height="450" src="https://github.com/HugoJacq/Linking_SAR_and_ABL_structure/blob/main/SST_and_SAR.png">
</p>

This repo gather the material needed for the study 'Linking_SAR_and_ABL_structure'.
In this README, you will find a description of the files (and their origin), how to reproduce the simulations and how to use the post-process scripts

This README is organised as follow:

[0. File Tree](#0-file-tree)

[1. How to run the simulations](#1-how-to-run-the-simulations)

[2. How to use post-process scripts](#2-how-to-use-post-process-scripts)

[3. Where to find SAR, SST, ERA5 data](#3-where-to-find-sar-sst-era5-data)

[4. What remains to do ?](#4-what-remains-to-do-)

[5. Contact](#5-contact)

## 0. File tree

The folder `user_modification` contains the necessary piece of code the modify MesoNH v5.7.0. See [section 1](#1-how-to-run-the-simulations) for a description on how to use them.
```
user_modification/
│    condsamp.f90
│    modd_blankn.f90
│    mode_io_write_nc4.f90
│    mode_turb_hor_uw.f90
│    mode_turb_ver_thermo_flux.f90
│    default_desfmn.f90
│    modd_condsamp.f90
│    mode_turb_hor_dyn_corr.f90
│    mode_turb_hor_vw.f90
│    modn_condsamp.f90
│    ini_modeln.f90
│    mode_io_field_write.f90
│    mode_turb_hor_uv.f90
│    mode_turb_ver_dyn_flux.f90
│    turb.f90
```

The folder `Namlists` is what you need to reproduce the simulations. The main file is `setup_big.py`, this is a python script that builds the MesoNH namlists. It is somewhat modular and you can for example easily change the initial conditions, numerical schemes and dimensions.
```
Namlists/
└───Namlist_injector/
│   |   PRE_IDEA1.nam
│   |   txt_00.py
│   |   ...
│   |   txt_09.py
│   |   setup_big.py
│   |   function_replaceSST.py
│   |   clean
│   |   compressor.py
│   |   run_compressor
```

The folder `post_process` is where are all python script for post-processing the simulation. The main file is `analyse.py`, where you can ask for specific plots with boolean switches variables.
```
post_process/
|   verif_turb.py
|   prepare_obs.py
|   mod_turb_stats.py
|   mod_build_mean.py
|   module_tools.py
|   module_cst.py
|   mod_build_CS.py
|   mod_CS_analysis.py
|   mod_spatial_stats.py
|   mod_first_look.py
|   analyse.py
```

## 1. How to run the simulations
![alt text](http://mesonh.aero.obs-mip.fr/mesonh57/Welcome?action=AttachFile&do=get&target=LogoMesoNH.jpg)
### 1.1 Install user modification of MNH
We used MesoNH version 5.7.0 to produce the simulations ([Lafore et al. 2018](https://doi.org/10.5194/gmd-11-1929-2018)) and it can be downloaded from [here](http://mesonh.aero.obs-mip.fr/mesonh56/Download).
You will need to modify the code. The modification allows the emission of two surface tracer and also the output of some subgrid fluxes.
Please follow the procedure described in the section 2.1 of [this](https://github.com/HugoJacq/ABL_response_to_SST_front/edit/main/), with the files from the folder `user_modification` and by replacing the version of MesoNH with the version 5.7.0

### 1.2 Overview of the simulation

The goal is to produce a simulation that can be compared to SAR data. We use a SST as forcing (no coupling with the ocean). The case studied is the Agulhas current on the 10th of December 2015. We used a semi-realistic configuration: we use a grid-nesting setup with the father having cyclic conditions in the mean (zonal) wind direction  and open boundary conditions at North and South boundaries. The horizontal resolution is 200m with a timestep of 4s. The father domain is represented by the dark rectangle on the first image. After a 4h spinup, we spawn a smaller domain (the 'son') at 50m resolution and with a timestep of 1s (green rectangle on the first image).

Numerical schemes and parametrisations are the same for both domains, the only differences are the resolution (spatial and temporal). We use the 3D version of the turbulence scheme, with Deardorff mixing length.

### 1.3 Performing the simulation

On your preferred super-computer, run the the command `python setup_big.py`. In the upper level of where `setup_big.py` is located, this will create folders with MesoNH namlists. They are order by number (00,01,..) and should be run in this order.
If you need to modify something about the simulation (initial conditions, number of CPUs, ...) it is better to do it on `setup_big` and then to run the script again. A saved version of the previous namelist will be generated in case you messed up.


`00_prep_ideal_dad`: this is the initialisation step for the dad domain. This is where the atmospheric state is set up. Dimensions and initial conditions are given in `setup_big.py`

`01_replace_dad_sst`: replacing the SST in the initial file for the dad domain. This step interpolate the 0.02° SST from LAT-LON grid onto the LES grid. The link between the two grid is made with the `dico_origin` entry in `setup_big.py`.

`02_run_dad_spinup`: Spinup run for the dad domain. This is 4h of simulated time. This is done in about 3h of real time with 1024 Cpus.

`03_spawning`: Spawning of the son domain. Informations about dimensions and where to spawn the domain are given in `setup_big.py`. This step writes vertical slices of the domain, this is necessary due to memory constrains. As a consequence, a lot more RAM is needed for this step to work (i used 512Go of RAM splited for 16 Cpus).

`04_prep_real_son`: Initialisation step for the son domain. This is the vertical interpolation step. As we work above the ocean, this should not be necessary but if this is not done, the code crashes... This step re-combine the many files from last step into one. It also needs a lot of memory (i used 1024Go of RAM splited for 16 Cpus).

`05_replace_son_sst`: This step is only useful if you use a idealized SST. If you set `SST_TYPE = 'REAL'` in `setup_big.py`, then you can skip this step as the change of SST has been made at the dad level (step 01).

`06_mesonh_2models`: This is the run of the two domains together. One hour of simulated time is roughly 10h of Cpu time with 2048 Cpus. Ouput files for the son domain are huge: backup files are 276Go, and frequent output are 16Go each. This is not ok for storage and this is why you will need to compress the data is you need to transfer it. You will find in the `Namlist_injector` folder two compressor files (`compressor.py` and `run_compressor`). `run_compressor` is really just the command `python compressor.py` but for super-computer batch run. In `compressor.py` the important line is the command `ncks -O -7 -L 1 --ppc default=4#RVT=5#THT=6#PABST=7 oldname newname`: this tell the NCO tool to compress the data with different significant digit (see [here](https://www.unidata.ucar.edu/blogs/developer/entry/compression_by_scaling_and_offfset) for more information about what are significant digit. Have a look, it is much more complicated than what it seems ...)

## 2. How to use post-process scripts

### 2.1 Pre-process

`where_is_my_energy.py`, `verif_turb.py` and `prepare_obs.py` are pre-process scripts. 

`prepare_obs.py` checks some caracteristics of the fields from ERA5, SAR and SST. `where_is_my_energy.py` and `verif_turb.py` are the validation scripts for the simulation. More specifically, `verif_turb.py` is checking that the grid nesting setup is working ok by computing spectra at different distance from the inflow boundary. In `verif_turb.py`, you will also find a tool to know how many significant digit number you need for each prognostic variable, and how to compress file to reduce diskspace. `where_is_my_energy.py` measure the loss of kinetic energy at the change of resolution between the two domains.

### 2.2 Post-process

All `mod_*`  and `analyse.py` are post-process scripts. The main file is `analyse.py` and it calls functions that are in the `mod_*` files. You can ask for different plots via boolean switches: for example, setting `PLOT_10m_WIND` to `True` will plot the 10m modulus of wind in the LES and the sigma0_detrend from the SAR. Each function in the `mod_*` files are documented with a description of the function, a liste of inputs and what is the ouput.

Here i provide informations about the custom files needed to plot figures. I use a Reynolds average that is an average in time, X and Y directions.

#### Before plotting, building
- `build_mean_file` in `mod_build_mean.py`: This function builds mean profiles. By default, average is over the full domain in time,X and Y. Available variables are MNH prognostic variables, and surface flux. 
- `build_CS` in `mod_build_CS.py`: This function builds a netcdf file with coherent structures identified.
- `save_S_n_SAR` in `mod_spatial_stats.py`:
- `save_S_n_LES` in `mod_spatial_stats.py`:

Here i provide more details for each switches:

#### First look
- `PLOT_10m_WIND`: Plot the 10 m modulus of the wind from the LES simulations and the SAR, with rectangle for boxes (boxes are defined in `d_boxes`). LES coordinates are in km, SAR is in lat-lon. For the SAR, boxes are defined in satellite coordinate (along track and across track pixels are called line and sample respectively).
- `WHERE_ARE_THE_BOXES`: For each boxes, plot 10 m modulus of wind in the local coordinates. For both LES and SAR. For the SAR, only the first few boxes are plotted.

#### Statistics on the textures
- `PLOT_2D_COVARIANCE`: 
- `PLOT_2D_SKEWNESS`:
- `S2_ANALYSIS`:

#### Turbulence convergence
- `VERIFY_TURB_CONVERGENCE`:
- `B_KPSD`:

#### Coherent structure analysis
- `PLOT_MEAN_PROGVAR`:
- `PLOT_MEAN_FLX`:
- `PLOT_TopView_with_CS`:

### 2.3 Known bugs
(open issue if there is a new one)

TBD

## 3. Where to find SAR, SST, ERA5 data

#### SAR
The SAR data has been provided by Alexis Mouche from Ifremer, but you can extract it from the [OVL](https://ovl.oceandatalab.com/?date=2015-12-10T12:00:00&timespan=1d&extent=1608235.0302722_-5098655.4894003_4739095.7083972_-3530779.1654331&center=3173665.3693347_-4314717.3274167&zoom=7&products=3857_SAR_roughness!3857_ODYSSEA_REG_SST!3857_ODYSSEA_SST!3857_GlobCurrent_CMEMS_geostrophic_streamline!3857_AMSR_sea_ice_concentration!3857_GIBS_MODIS_Terra_CorrectedReflectance_TrueColor!3857_GIBS_MODIS_Aqua_CorrectedReflectance_TrueColor&opacity=80_100_100_60_100_49.498_100&stackLevel=100.02_50.03_30.02_120.07_40.03_50.25_50.22&filter=+A,+B,+IW,+EW,+SM,+WV,+VV,+HH&selection=1000000) tool. For this you will need to switch your browser to "phone" view (via F11) and then increase as much as you can the resolution. This is needed because the extraction feature of OVL is screen-resolution dependent. To extract data from OVL, click to the small gear icon on the top bar. This will bring the extraction tool. Then enter the desired LAT-LON coordinates to extract and process. This will give you a numpy array file and you will need to get back the X and Y dimensions with the number of points of the file and the coordinates you entered in the last step. After this, you will need to remove yourself the non valid values from the data. The file is 1.1Go.

#### SST
The SST data is also from Ifremer but it is also available on Copernicus website (SST L4 Ifremer GHRSST Odyssea): 20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc. The file is about 2.6Mo

#### ERA5
ERA5 data can be retrieved with `python ERA5_pressure_level_retrieve.py`: this will dowload a file with ERA5 variables on pressure levels. The file is about 26Mo.

## 4. What remains to do ?
This section has been written in December of 2024 by Hugo Jacquet.

### 4.1 What are the scientific questions ?

Please have a look at the latex file.

### 4.2 What I have done

TBD

### 4.3 What I think needs to be done

I split this section in two: technical needs are some improvment that I see that would be nice to have in the current state of the scripts; scientific needs are functionnalities that are not yet present but that would provide some part of the answer of the [scientific questions](#41-what-are-the-scientific-questions-).
#### Technical needs

- modify `Where_boxes` to allow selection of which SAR boxes are plotted

#### Scientific needs

## 5. Contact

This work has been done with Fleur Couvreux (CNRM) and Alex Ayet (Gipsa-lab). You can contact me at hugo.jacquet.work[at]gmail.com
