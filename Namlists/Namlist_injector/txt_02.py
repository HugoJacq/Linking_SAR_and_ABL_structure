def return_txt_02(name_ini_dad,dt_d,XT4DIFU_d,LENGTH_SPINUP,BACK_SPINUP,OUT_SPINUP,
                  LBISURF,CBISURF,XSST_CRIT,
                    CUVW_ADV_SCHEME,CTEMP_SCHEME,LNUMDIFU,
                    profile_mnh,nhalo,JPHEXT,LMNH_MPI_BSEND,MPI_BUFFER,
                    dic_nodes,dic_cpu,dic_timelimit,mpicmd,time_cmd):

    TXT_02_run = """#!/bin/bash
#SBATCH -J MNH02
#SBATCH -N """+str(dic_nodes['02'])+"""           # nodes number (=NBP)   
#SBATCH -n """+str(dic_cpu['02'])+"""       # CPUs number (on all nodes) (=NBP*TPN) 
#SBATCH -o MNH02.eo%j   #
#SBATCH -e MNH02.eo%j   #
#SBATCH -t """+ dic_timelimit['02']+"""    # time limit

# Echo des commandes
set -x
ulimit -c 0
ulimit -s unlimited
# Arret du job des la premiere erreur

# Nom de la machine
hostname

# links
ln -sf ../00_prep_ideal_dad/DAD_PGD* .
ln -sf ../00_prep_ideal_dad/INIT_DAD* .


mkdir -p FICHIERS_OUT

#. /home/cnrm_other/ge/mrmc/jacqueth/MNH-V5-6-0/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-6-0-MYMODIF1-MPIAUTO-DEBUG
. """+profile_mnh+"""

set -e
export MPIRUN='"""+mpicmd+""" -np """+str(dic_cpu['02'])+"""'
export TIMEOUT="timelimit -t 18000" #259200=3j 86400=1j

rm -f file_for_xtransfer pipe_name

ls -lrt
"""+time_cmd+"""${MPIRUN} MESONH${XYZ}
#time $MPIRUN ~escobar/MPI/iping_ipong_bi_10it_check_network_intel2019_ompi4051 # pour voir si un noeud a un probleme
#${TIMEOUT} time -v ${MPIRUN} MESONH${XYZ} # local cmd
#time ${MPIRUN} MESONH${XYZ} #belenos cmd
#${TIMEOUT} time -v ${MPIRUN} idrmem MESONH${XYZ} PAS EN PRODUCTION !!
mv -f OUTPUT_LISTING1  OUTPUT_LISTING1_mnh
mv -f OUTPUT_LISTING0  OUTPUT_LISTING0_mnh

rm -f REMAP*
ls -lrt
"""



    TXT_02_EXSEG = """&NAM_CONFZ LMNH_MPI_BSEND = """+str(LMNH_MPI_BSEND['02'])+""",MPI_BUFFER_SIZE="""+str(MPI_BUFFER['02'])+""", NZ_VERB=5 /
&NAM_CONFIO LCDF4=.TRUE., LLFIOUT=.FALSE., LLFIREAD=.FALSE. / ! Type of file

&NAM_LUNITn  CINIFILE = '"""+name_ini_dad+"""',	! Name of file from which to start the segment
	     CINIFILEPGD="DAD_PGD" /      ! Name of PGD
	     
&NAM_DYNn    XTSTEP = """+dt_d+""", 		! Timestep
             CPRESOPT = "ZRESI",	! Pressure solver
             XT4DIFU="""+XT4DIFU_d+""",		! Diffusion time for momentum
             LHORELAX_UVWTH=.FALSE.,	! Relaxation in case LBC='OPEN'
	     NRIMX=1,			! Number of relaxation points in X
             NRIMY=1,			! Number of relaxation points in Y
             XRIMKMAX=0.01/		! K relaxation coefficient
             
&NAM_ADVn CUVW_ADV_SCHEME = '"""+CUVW_ADV_SCHEME+"""', 	! CEN4TH = Centered order 4 for U,V,W
	  CTEMP_SCHEME='"""+CTEMP_SCHEME+"""', 		! RungeKuta 4 for time if CEN4TH
                                    ! 4*dx of effective resolution if CEN4TH and RKC4
          NWENO_ORDER=5,
          CMET_ADV_SCHEME = "PPM_01",	! Piecewise Parabolic Method for THT,TKE,Scalar 
          CSV_ADV_SCHEME="PPM_01"/	! Piecewise Parabolic Method for tracer 
          
&NAM_PARAMn  CTURB = "TKEL",		! 1.5 order closure (TKE and LM)
             CRAD = "NONE", 		! Radiation
             CCLOUD= "ICE3", 		! Cloud 
             CSCONV= "NONE",		! Param. shallow convection
             CDCONV= "NONE" /		! Param. deep convection
             
&NAM_PARAM_ICEn LWARM=.TRUE./		! ICE3 namelist 
             
&NAM_SEAFLUXn CSEA_FLUX="COARE3",CSEA_ALB="UNIF"/  ! Surface scheme (see Surfex)

&NAM_LBCn    CLBCX = 2*"CYCL",		! LBC X direction
             CLBCY = 2*"OPEN"/		! LBC Y direction
             
&NAM_TURBn   XIMPL = 1.,		! 1=full implicit, 0=full explicit
             CTURBLEN = "DEAR",		! Turbulent mixing length
             CTURBDIM = "3DIM",		! 3D or 1D
             LRMC01 = T,		! Separate mixing and dissipative mixing lengths
             LTURB_FLX = T,		! Turbulent flux stored in BACKUPs (THW_FLX, ...)
             LTURB_DIAG = T,		! Turbulent diag (TKE_DISS, ...)
             LSIG_CONV = F,		! Compute Sigmas due to subgrid condensation
             /		! Flag for subgrid condensation
             
&NAM_CONF    CCONF="START",		! RESTA or START
             LFLAT = T,			! flat terrain
             CEQNSYS = "DUR",		! system of equation, Durran
             LFORCING = T,		! use the forcing defined in ZFRC
             NMODEL = 1,		! Number of nested models
             NVERB = 1,			! verbose, 10 is maximum
             CEXP = "DAD_S",		! name of experiment
             CSEG = "001",		! name of segment
             NHALO="""+str(nhalo)+""", 			! halo for parallel computation
             JPHEXT="""+str(JPHEXT)+"""/    ! Horizontal exterior points

&NAM_CONDSAMP 	LCONDSAMP = .TRUE.	! Flag to activate conditional sampling
	NCONDSAMP = 4,		! Number of conditional samplings
	!XRADIO = 3*900.,	! Period of radioactive decay (15min)
	!XSCAL =3*1.,		! Scaling factor (1)
	XHEIGHT_BASE = 50,	! Height below Cbase where the 2nd tracer is released
	XDEPTH_BASE = 0.,	! Depth on which the 2nd tracer is released
	XHEIGHT_TOP = 50,	! Height above Ctop where the 3rd tracer is released
	XDEPTH_TOP = 50,	! Depth on which the 3rd tracer is released
	NFINDTOP = 1		! see doc
	!XTHVP =0.25,		! if NFINDTOP = 2, threshold to detect ABLH base on thtv
	LTPLUS = .FALSE.,	! see doc
	LBISURF = """+LBISURF+""",	! Flag to use 2 tracers at surface
	CBISURF = '"""+CBISURF+"""',	! Type of emission when LBISURF=TRUE 
	XSST_CRIT = """+str(XSST_CRIT)+"""/ ! if CBISURF='USER'
             
&NAM_DYN     XSEGLEN = """+LENGTH_SPINUP+"""	, 	! length of segment
             XASSELIN = 0.2,		! Asselin temporal filter
             LCORIO = T,		! T=Earth's rotation is taken into account
             XALKTOP = 0.005, 		! Top sponge layer coefficient
             XALZBOT = 1800., 		! Altitude of the begining of sponge layer
             LNUMDIFU = """+LNUMDIFU+""" /		! Flag for num. diffusion (for CEN4TH)
 
&NAM_FRC  LTEND_THRV_FRC= F,		! Flag to use THT and RV tendencies
          LVERT_MOTION_FRC= F,		! Flag to use large scale vertical transport
          LRELAX_THRV_FRC=F,		! Flag to relax to ZFRC values of THT and RV
          LGEOST_UV_FRC=T,		! Flag to use ZFRC values of U,V as geo wind
          LRELAX_UV_FRC=F/ 		! Flag to relax to ZFRC values of U,V    
	  	
&NAM_BACKUP XBAK_TIME_FREQ(1) = """+BACK_SPINUP+""", 	! Frequency of output of all variables
            XBAK_TIME_FREQ_FIRST(1)= """+BACK_SPINUP+""" /	! First time 
&NAM_OUTPUT 				! Smaller files
	COUT_VAR(1,1)='RVT',
	COUT_VAR(1,2)='THT',
	COUT_VAR(1,3)='TKET',
	COUT_VAR(1,4)='UT',
	COUT_VAR(1,5)='VT',
	COUT_VAR(1,6)='WT',
	XOUT_TIME_FREQ(1)="""+OUT_SPINUP+""",	! Frequency of output of pronostics variables
					! 	and subgrid flux (custom modification)
	LOUT_BEG=.TRUE.,		! Save at begining of segment
	LOUT_END=.TRUE.,		! Save at end of segment
	COUT_DIR='FICHIERS_OUT',	! Directory to save the files
	XOUT_TIME_FREQ_FIRST(1) = """+str(OUT_SPINUP)+"""	! Time at which the saving starts
	LOUT_REDUCE_FLOAT_PRECISION(1)=.TRUE./	! Flag to use simple precision
&NAM_NEBn LSUBG_COND = F   /
"""
    return TXT_02_EXSEG,TXT_02_run