# Linking_SAR_and_ABL_structure

This repo gather the material needed for the study 'Linking_SAR_and_ABL_structure'.
In this README, you will find a description of the files (and their origin), how to reproduce the simulations and how to use the post-process scripts

## 0. File tree
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
I used MesoNH version 5.7.0 to produce the simulations ([Lafore et al. 2018](https://doi.org/10.5194/gmd-11-1929-2018)) and it can be downloaded from [here](http://mesonh.aero.obs-mip.fr/mesonh56/Download).
Please follow the procedure described in the section 2.1 of [this](https://github.com/HugoJacq/ABL_response_to_SST_front/edit/main/), with the files from the folder `user_modification` and by replacing the version of MesoNH with the version 5.7.0

## 1.2 Overview of the simulation


## 2. How to use post-process scripts

