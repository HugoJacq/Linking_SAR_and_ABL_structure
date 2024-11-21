# python setup.py
import os
from txt_00 import return_txt_00
from txt_01 import return_txt_01
from txt_02 import return_txt_02
from txt_03 import return_txt_03
from txt_04 import return_txt_04
from txt_05 import return_txt_05
from txt_06 import return_txt_06
"""
This script is building the namlist and executables for a MesoNH simulation.

The case studied here is a grid nesting case.

Some values of namlist can be modified within the 'INPUT ZONE',
	but keep in mind that if you intend to use this script for something else then
	you'll have to adapt some things ! 
	
	
Workflow:

00_prep_ideal_dad   : initialize dad domain.
01_replace_dad_sst  : replace dad SST
02_run_dad_spinup   : dad spinup
03_spawning         : horizontal interpolation of son
04_prep_real_son    : vertical interpolation of son
05_replace_son_sst  : replace son SST (if needed)
06_mesonh_2models   : main run, both domains run 
"""

# ============================================================================
# INPUT ZONE
# ============================================================================

ON_BELENOS = True # True if on HPC

# common to all steps --------------------------------------------------------
profile_mnh_local = "/home/jacqhugo/MNH-V5-7-0/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-7-0-MODIF_condsamp_and_OUTfluxes_v2-MPIAUTO-DEBUG" 
profile_mnh_belenos = "/home/cnrm_other/ge/mrmc/jacqueth/MNH-V5-7-0/conf/profile_mesonh-LXifort-R8I4-MNH-V5-7-0-MODIF_condsamp_and_OUTfluxes_v2-MPIAUTO-O2"

# Numerical configuration ----------------------------------------------------
if ON_BELENOS:
    ADV_SCHEME = 'CEN4TH'           # WENO5 or CEN4TH, temporal scheme is changed accordingly
    LENGTH_SPINUP = '14400.' 	   # float in seconds, dad spinup segment length
    OUT_SPINUP = '600.' 		   # float, in seconds. OUT files frequency
    BACK_SPINUP = LENGTH_SPINUP	    # float, in seconds BACKUP files frequency
    LENGTH_RUN2 = '3600.' 		   # float in seconds, run of dad and son models
    OUT_RUN2 = '600.' 			   # float, in seconds.	
    BACK_RUN2 = '3600.'		       # float, in seconds

    # MAKE SURE THAT Nx Ny CAN BE WRITTEN
    #   AS 2^n.3^m.5^l with m,n,l INTEGERS !!
    # DAD
    dt_d = '4.' # float in seconds
    NI_d = 540 	# Nx
    NJ_d = 640 	# Ny
    NZ = 160 	# Nz, common with son
    dx = 200
    XT4DIFU_d = '800.'     # float in seconds
    # SON
    dt_s = '1.' # float in seconds
    IXOR = 200 		# first X point (in DAD coordinates) for spawning Son domain
    IYOR = 100 		# first Y point (in DAD coordinates) for spawning Son domain
    NI_s = 320 		# X Number of dad cells to define son domain
    NJ_s = 450 		# Y Number of dad cells to define son domain
    dx_ratio = 4 	# = dx_dad / dx_son
    XT4DIFU_s = '100.' # Float in seconds

else: # local machine, smaller domain
    ADV_SCHEME = 'CEN4TH'           # WENO5 or CEN4TH, temporal scheme is changed accordingly
    LENGTH_SPINUP = '12.' 	   # float in seconds, dad spinup segment length
    OUT_SPINUP = '12.' 		   # float, in seconds. OUT files frequency
    BACK_SPINUP = LENGTH_SPINUP	    # float, in seconds BACKUP files frequency
    LENGTH_RUN2 = '12.' 		   # float in seconds, run of dad and son models
    OUT_RUN2 = '12.' 			   # float, in seconds.	
    BACK_RUN2 = '12.'		       # float, in seconds

    # MAKE SURE THAT Nx Ny CAN BE WRITTEN
    #   AS 2^n.3^m.5^l with m,n,l INTEGERS !!
    # DAD
    dt_d = '4.' # float in seconds
    NI_d = 64 	# Nx
    NJ_d = 64 	# Ny
    NZ = 160 	# Nz, common with son
    dx = 200
    XT4DIFU_d = '800.'     # float in seconds
    # SON
    dt_s = '1.' # float in seconds
    IXOR = 2 		# first X point (in DAD coordinates) for spawning Son domain
    IYOR = 2 		# first Y point (in DAD coordinates) for spawning Son domain
    NI_s = 32 		# X Number of dad cells to define son domain
    NJ_s = 32 		# Y Number of dad cells to define son domain
    dx_ratio = 4 	# = dx_dad / dx_son
    XT4DIFU_s = '100.' # Float in seconds


# SST parameters -------------------------------------------------------------
SST_TYPE = 'REAL' # UNIFORM,IDEALIZED,REAL
#               IDEALIZED   <=> 1 front, N/S or E/W
#               UNIFORM_SST <=> SST = T0, override IDEALIZED
#               REAL        <=> SST = Agulhas front
# if IDEALIZED:
FUNC = 'tanh' 	# shape of SST tanh or linear
DIR = 'Y' 		# direction of change of SST
T0=298 			# K, 
deltaT=1.5		# K, 
atpos=6000 		# m : position of the front
Lx=1000 		# m : width of the front
# if REAL
#  -> resolution of reanalysis is 0.02 degree.
SST_DIR = {'Belenos':'/home/cnrm_other/ge/mrmc/jacqueth/SST_L4_IFREMER_Agulhas/20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc',
            'MCP':'/home/jacqhugo/scripts/simu_nest_aggv3/SST_from_Ifremer/20151210-IFR-L4_GHRSST-SSTfnd-ODYSSEA-SAF_002-v2.0-fv1.0.nc'}
DegLat = 111.2  #  =1° lat in km
DegLon = 91.6   #  =1° lon in km at ~ -35°N
dico_origin = {"DAD":(24.06,-36.18),"SON":(24.5,-36)}

# NAM_CONSAMP (common for the 02 and 06) --------------------------------------
LBISURF = '.TRUE.'
CBISURF = 'USER'  # or 'MEAN' or 'USER' or 'MEDI'
XSST_CRIT = 295.65 # 295.65K=22.5°C, if IDEALIZED:T0+deltaT/2

# Environment : initial conditions -------------------------------------------
#	to be used with CHANGE_ENV from prepare_obs.py
CHOICE = 'ERA5'
d_case = {'ERA5':{
				'Psurf': 101250, # Pa, pressure at ground level
     			'Zlevel':[750,2000], # m
				'tht'	:[294.25,304.5], # K
                'rv' 	:[0.0105,0.00005], # kg/kg
                'U' 	:'7.5', 	# m/s
				'V'		:'0.'
                    },   
            '0':{
                'Psurf': 100000, # Pa, pressure at ground level
                'Zlevel':[250,2000], # m
				'tht'	:[295.5,300.75], # K
                'rv' 	:[0.01,0.00195], # kg/kg
                'U' 	:'7.5', 	# m/s
				'V'		:'0.'
                    }, 
            }
# Note : refaire le cas 0 avec Psurf=101250 à la surface (et pas P00=10^5)
# 			mais que si on utilise cette simu !

""" 
HOW TO CHOOSE CPU NUMBER

-> reminder of the steps:
00_prep_ideal_dad : initialize dad domain.
01_replace_dad_sst : replace dad SST
02_run_dad_spinup : dad spinup
03_spawning : horizontal interpolation of son
04_prep_real_son : vertical interpolation of son
05_replace_son_sst : replace son SST
06_mesonh_2models : main run, both domains run 
"""
# READ and WRITE number of Z split -------------------------------------------
dic_PROCIO_R = {'00':0, # number of levels to read
		'01':0,
		'02':0,
		'03':0,
		'04':64,
		'05':0,
		'06':0
		}

dic_PROCIO_W = {'00':0,  # number of levels to write
		'01':0,
		'02':0,
		'03':64,
		'04':0,
		'05':0,
		'06':0
		}

# CPUs : local ---------------------------------------------------------------
#	32 Cpu, 256Go RAM
dic_cpu_local = {'00':1,
		'01':1,
		'02':16,
		'03':1,
		'04':1,
		'05':1,
		'06':16
		}
# CPUs : Belenos (and if NI,NJ are big enough !)
#	1 node = 128 CPU with 256Go of RAM for 1 node shared to all CPUs
dic_cpu_belenos = {'00':1,
		'01':1,
		'02':1024,
		'03':16,
		'04':16,
		'05':1,
		'06':2048
		}
dic_nodes = {'00':1,
		'01':1,
		'02':8,
		'03':2, # 2 because needs a lot of memory
		'04':4, # 4 because needs a lot of memory
		'05':1,
		'06':16
		}
dic_timelimit = {
        '00':'00:10:00',
		'01':'00:10:00',
		'02':'04:00:00',
		'03':'01:00:00',
		'04':'01:00:00',
		'05':'00:10:00',
		'06':'04:00:00'
		}
MPI_BUFFER = {	# in MB
		'00':200,
		'01':-99,	# not used
		'02':100,
		'03':200, 	
		'04':200,
		'05':-99,	# not used
		'06':2000
		}		
LMNH_MPI_BSEND = {
		'00':'.TRUE.',
		'01':-99, 		# not used
		'02':'.FALSE.',
		'03':'.FALSE.',	
		'04':'.FALSE.',
		'05':-99,		# not used
		'06':'.TRUE.'
		}
			
# ============================================================================
# END OF INPUT ZONE
# ============================================================================

os.system('cp setup_big.py last_setup_big.py') # we make a save for next change
		
if ON_BELENOS:
	profile_mnh = profile_mnh_belenos
	dic_cpu = dic_cpu_belenos
	time_cmd = 'time ' 
	path_sst = SST_DIR['Belenos']
	mpicmd = 'Mpirun'
else:
	profile_mnh = profile_mnh_local
	dic_cpu = dic_cpu_local
	time_cmd = '${TIMEOUT} time '
	path_sst = SST_DIR['MCP']
	mpicmd = 'mpirun'

if ADV_SCHEME=='WENO5':
    CUVW_ADV_SCHEME = 'WENO_K'
    CTEMP_SCHEME = 'RK53'
    nhalo = 3
    JPHEXT = 3
    LNUMDIFU = 'F'
elif ADV_SCHEME=='CEN4TH':
    CUVW_ADV_SCHEME = 'CEN4TH'
    CTEMP_SCHEME = 'RKC4'
    nhalo = 1
    JPHEXT = 1
    LNUMDIFU = 'T'
else:
    raise Exception('the advective scheme '+ADV_SCHEME+' is not recognized')



if SST_TYPE=="UNIFORM":
    IDEALIZED = False # overrided
    UNIFORM_SST = True
elif SST_TYPE=="IDEALIZED":
    IDEALIZED = True
    UNIFORM_SST = False
elif SST_TYPE=="REAL":
    IDEALIZED = False
    UNIFORM_SST = False
else:
    raise Exception('the SST '+SST_TYPE+' is not recognized')


name_ini_dad = 'INIT_DAD_SST' # by default, SST is changed
name_ini_son = 'INIT_SON' # by default, SST change is done by MNH
if UNIFORM_SST:
    name_ini_dad = 'INIT_DAD'    
if SST_TYPE=="IDEALIZED":
    name_ini_son = 'INIT_SON_SST'

# 00 =================================================================


TXT_00_PREP_IDEAL,TXT_00_run = return_txt_00(T0,d_case,CHOICE,NI_d,NJ_d,NZ,
											MPI_BUFFER,nhalo,JPHEXT,dic_nodes,dic_cpu,dic_timelimit,
											profile_mnh,mpicmd,time_cmd)

# 01 =================================================================

TXT_01_replaceSST,TXT_01_run = return_txt_01(UNIFORM_SST,IDEALIZED,FUNC,DIR,T0,deltaT,atpos,Lx,
                  XSST_CRIT,dico_origin,path_sst,
                  IXOR,IYOR,NI_s,NJ_s,
                  nhalo,dic_nodes,dic_cpu,dic_timelimit)

# 02 =================================================================

TXT_02_EXSEG,TXT_02_run = return_txt_02(name_ini_dad,dt_d,XT4DIFU_d,LENGTH_SPINUP,BACK_SPINUP,OUT_SPINUP,
                  LBISURF,CBISURF,XSST_CRIT,
                    CUVW_ADV_SCHEME,CTEMP_SCHEME,LNUMDIFU,
                    profile_mnh,nhalo,JPHEXT,LMNH_MPI_BSEND,MPI_BUFFER,
                    dic_nodes,dic_cpu,dic_timelimit,mpicmd,time_cmd)


# 03 =================================================================

TXT_03_SPAWN,TXT_03_RUN = return_txt_03(IXOR,IYOR,NI_s,NJ_s,dx_ratio,
                                    LMNH_MPI_BSEND,MPI_BUFFER,dic_PROCIO_R,dic_PROCIO_W,
                                    profile_mnh,dic_nodes,dic_cpu,dic_timelimit,mpicmd,time_cmd)


# 04 =================================================================

TXT_04_PRE_REAL,TXT_04_RUN = return_txt_04(LMNH_MPI_BSEND,MPI_BUFFER,dic_nodes,dic_cpu,dic_timelimit,
                                        profile_mnh,mpicmd,time_cmd,
                                        dic_PROCIO_R,dic_PROCIO_W)

# 05 =================================================================

TXT_05_replaceSST,TXT_05_run = return_txt_05(UNIFORM_SST,IDEALIZED,FUNC,DIR,T0,deltaT,atpos,Lx,
                                        XSST_CRIT,dico_origin,path_sst,
                                        IXOR,IYOR,NI_s,NJ_s,
                                        nhalo,dic_nodes,dic_cpu,dic_timelimit)

# 06 =================================================================

TXT_06_EXSEG1,TXT_06_EXSEG2,TXT_06_RUN  = return_txt_06(XSST_CRIT,LBISURF,CBISURF,
                                                    dt_d,LENGTH_RUN2,BACK_RUN2,OUT_RUN2,
                                                    dx_ratio,dt_s,name_ini_son,
                                                    LNUMDIFU,XT4DIFU_d,XT4DIFU_s,CUVW_ADV_SCHEME,CTEMP_SCHEME,
                                                    nhalo,JPHEXT,LMNH_MPI_BSEND,MPI_BUFFER,
                                                    profile_mnh,dic_nodes,dic_cpu,dic_timelimit,mpicmd,time_cmd
                                                    )


# INJECTION =================================================================

def Check_file(liste):
	"""Check if file is here,
		if yes : rename it by adding OLD at its end
		if no : create it
	"""
	for name in liste:
		if os.path.exists(name):
			os.system('mv '+name+' '+name+'OLD')
		else:
			os.system('touch '+name)
			
			
print(' * Creating namlists...')
os.system('mkdir -p ../00_prep_ideal_dad')  
os.system('mkdir -p ../01_replace_dad_sst')       
os.system('mkdir -p ../02_run_dad_spinup') 
os.system('mkdir -p ../03_spawning')
os.system('mkdir -p ../04_prep_real_son')
os.system('mkdir -p ../05_replace_son_sst') 
os.system('mkdir -p ../06_mesonh_2models')
os.system('cp clean ../00_prep_ideal_dad ')
os.system('cp function_replaceSST.py ../01_replace_dad_sst') 
os.system('cp clean ../02_run_dad_spinup')
os.system('cp clean ../03_spawning')
os.system('cp clean ../04_prep_real_son')
os.system('cp function_replaceSST.py ../05_replace_son_sst') 
os.system('cp clean ../06_mesonh_2models')

# List of file to check if exist and if not create it
#	and if yes then rename the older one as OLD
L_namlist = ['../00_prep_ideal_dad/PRE_IDEA1.nam',
			'../01_replace_dad_sst/replaceSST.py',
			'../02_run_dad_spinup/EXSEG1.nam',
			'../03_spawning/SPAWN1.nam',
			'../04_prep_real_son/PRE_REAL1.nam',
			'../05_replace_son_sst/replaceSST.py',
			'../06_mesonh_2models/EXSEG1.nam',
			'../06_mesonh_2models/EXSEG2.nam',
			] 

L_run = ['../00_prep_ideal_dad/run_prep_phase',
		'../01_replace_dad_sst/run_replaceSST',
		'../02_run_dad_spinup/run_mesonh',
		'../03_spawning/run_spawn',
		'../04_prep_real_son/run_prep_phase',
		'../05_replace_son_sst/run_replaceSST',
		'../06_mesonh_2models/run_mesonh'] 

Check_file(L_run)
# 00	
run_00 = open('../00_prep_ideal_dad/run_prep_phase', 'w')
run_00.write(TXT_00_run)
run_00.close()
os.system('chmod +x ../00_prep_ideal_dad/run_prep_phase')
# 01
run_01 = open('../01_replace_dad_sst/run_replaceSST', 'w')
run_01.write(TXT_01_run)
run_01.close()
os.system('chmod +x ../01_replace_dad_sst/run_replaceSST')
# 02
run_02 = open('../02_run_dad_spinup/run_mesonh', 'w')
run_02.write(TXT_02_run)
run_02.close()
os.system('chmod +x ../02_run_dad_spinup/run_mesonh')
# 03
run_03 = open('../03_spawning/run_spawn', 'w')
run_03.write(TXT_03_RUN)
run_03.close()
os.system('chmod +x ../03_spawning/run_spawn')
# 04
run_04 = open('../04_prep_real_son/run_prep_phase', 'w')
run_04.write(TXT_04_RUN)
run_04.close()
os.system('chmod +x ../04_prep_real_son/run_prep_phase')
# 05
run_05 = open('../05_replace_son_sst/run_replaceSST', 'w')
run_05.write(TXT_05_run)
run_05.close()
os.system('chmod +x ../05_replace_son_sst/run_replaceSST')
# 06
run_06 = open('../06_mesonh_2models/run_mesonh', 'w')
run_06.write(TXT_06_RUN)
run_06.close()
os.system('chmod +x ../06_mesonh_2models/run_mesonh')
print('	I replaced all run files')
	

Check_file(L_namlist)
	
pre_idea_00 = open('../00_prep_ideal_dad/PRE_IDEA1.nam', 'w')
pre_idea_00.write(TXT_00_PREP_IDEAL)
pre_idea_00.close()

spawn_01 = open('../01_replace_dad_sst/replaceSST.py', 'w')
spawn_01.write(TXT_01_replaceSST)
spawn_01.close()
	
exseg1_02 = open('../02_run_dad_spinup/EXSEG1.nam', 'w')
exseg1_02.write(TXT_02_EXSEG)
exseg1_02.close()	

spawn_03 = open('../03_spawning/SPAWN1.nam', 'w')
spawn_03.write(TXT_03_SPAWN)
spawn_03.close()
	
pre_real_04 = open('../04_prep_real_son/PRE_REAL1.nam', 'w')
pre_real_04.write(TXT_04_PRE_REAL)
pre_real_04.close()
	
pre_real_05 = open('../05_replace_son_sst/replaceSST.py', 'w')
pre_real_05.write(TXT_05_replaceSST)
pre_real_05.close()
	
exseg1_06 = open('../06_mesonh_2models/EXSEG1.nam', 'w')
exseg1_06.write(TXT_06_EXSEG1)
exseg1_06.close()
exseg2_06 = open('../06_mesonh_2models/EXSEG2.nam', 'w')
exseg2_06.write(TXT_06_EXSEG2)
exseg2_06.close()
print('	I replaced all namlist files')
	
print('')
os.system('ls -lh ../00_prep_ideal_dad') 
os.system('ls -lh ../01_replace_dad_sst')       
os.system('ls -lh ../02_run_dad_spinup') 
os.system('ls -lh ../03_spawning')
os.system('ls -lh ../04_prep_real_son')
os.system('ls -lh ../05_replace_son_sst')
os.system('ls -lh ../06_mesonh_2models')


print('	done!')	
	
	
	
	
	
	
	
	
	
	
	
	
