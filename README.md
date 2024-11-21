# Linking_SAR_and_ABL_structure

This repo gather the material needed for the study 'Linking_SAR_and_ABL_structure'.
In this README, you will find a description of the files (and their origin), how to reproduce the simulations and how to use the post-process scripts

## 0. File tree
This repo contains the following files:

Namlists/
│   PRE_IDEA1.nam
│   EXSEG1.nam.spinup
│   EXSEG1.nam.run
|   replaceSST.py
|   run_mesonh
|   run_prep_ideal
|   run_spinup

post_process/
|


## 1. How to run the simulations
![alt text](http://mesonh.aero.obs-mip.fr/mesonh57/Welcome?action=AttachFile&do=get&target=LogoMesoNH.jpg)
### 1.1 Install user modification of MNH
I used MesoNH version 5.7.0 to produce the simulations ([Lafore et al. 2018](https://doi.org/10.5194/gmd-11-1929-2018)) and it can be downloaded from [here](http://mesonh.aero.obs-mip.fr/mesonh56/Download).
Please follow the procedure described in the section 2.1 of [this](https://github.com/HugoJacq/ABL_response_to_SST_front/edit/main/), with the files from the folder `user_modification` and by replacing the version of MesoNH with the version 5.7.0

## 2. How to use post-process scripts
