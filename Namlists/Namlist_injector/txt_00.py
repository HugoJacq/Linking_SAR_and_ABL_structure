def return_txt_00(T0,d_case,CHOICE,NI_d,NJ_d,NZ,
                  MPI_BUFFER,nhalo,JPHEXT,dic_nodes,dic_cpu,dic_timelimit,
                  profile_mnh,mpicmd,time_cmd):
    TXT_00_PREP_IDEAL = """&NAM_CONFZ MPI_BUFFER_SIZE="""+str(MPI_BUFFER['00'])+""" /
    &NAM_CONFIO LCDF4=.TRUE.,LLFIOUT=.FALSE., LLFIREAD=.FALSE. /
    &NAM_DIMn_PRE  NIMAX="""+str(NI_d)+""",	! Nx
                NJMAX="""+str(NJ_d)+""" / 	! Ny
    &NAM_VER_GRID NKMAX="""+str(NZ)+""",
                YZGRID_TYPE='FUNCTN', !option manual et prescrire une grille vert.
                ZDZGRD=2.,	! zGrid size at bottom
                ZDZTOP=20.,	! zGrid size at top
                ZZMAX_STRGRD=100.,
                ZSTRGRD=10.,
                ZSTRTOP=0.3 /
    &NAM_CONF_PRE  LCARTESIAN=T,
                CEQNSYS='DUR',	! Durran 
                NVERB=10,
                CIDEAL='RSOU',	! to prescribe custom profiles
                CZS='FLAT',
                LBOUSS=F,
                LPERTURB=T, 	! Initial perturbation for spinup
                LFORCING=T,	! Relaxation toward geo wind = pressure gradient
                NHALO="""+str(nhalo)+"""
                JPHEXT="""+str(JPHEXT)+"""/ Horizontal exterior points
    &NAM_PERT_PRE  CPERT_KIND='WH',	! Type of perturbation
                XAMPLIWH=0.1 /	! Amplitude of perturbation
    &NAM_CONFn     LUSERV=F,LUSERC=F/	! Flags for water and cloud mixing ratios
    &NAM_GRID_PRE  XLAT0=-35.5,	! Approximate position of the Agulhas
                XLON0=25.0 / 
    &NAM_GRIDH_PRE XDELTAX=200.,	! Grid resolution
                XDELTAY=200.  /
    &NAM_LUNITn    CINIFILE='INIT_DAD',CINIFILEPGD="DAD_PGD" /
    &NAM_DYNn_PRE  CPRESOPT='ZRESI'/ 
    &NAM_LBCn_PRE  CLBCX=2*'CYCL',	! Boundary conditions
                CLBCY=2*'CYCL' /

    &NAM_GRn_PRE   CSURF='EXTE' /
    &NAM_PGD_SCHEMES CSEA="SEAFLX", CNATURE='NONE' ,CWATER='NONE',CTOWN='NONE'/
    &NAM_PREP_SURF_ATM NYEAR=2000, NMONTH=01, NDAY=1, XTIME=0. /
    &NAM_COVER XUNIF_COVER(1)=1. /
    &NAM_PREP_SEAFLUX XSST_UNIF="""+str(T0)+"""/ ! Convective, Ta_ini=296.55K

    RSOU
    2000 01 01 0.
    'ZUVTHDMR'	! type of radiosounding
    0.		! height of ground, in meters
    """+str(d_case[CHOICE]['Psurf'])+""".		! pressure at ground level
    """+str(d_case[CHOICE]['tht'][0])+"""		! dry potential temperature at ground level
    0.01		! mixing ratio at ground level
    2		! number of wind levels
    0. """+d_case[CHOICE]['U']+""" """+d_case[CHOICE]['V']+"""	! altitude, u, v
    """+str(d_case[CHOICE]['Zlevel'][-1])+""". """+d_case[CHOICE]['U']+""" """+d_case[CHOICE]['V']+"""
    """+str(len(d_case[CHOICE]['Zlevel'])+1)+"""		! number of mass level, ground is 1
    """+str(d_case[CHOICE]['Zlevel'][0])+""". """+str(d_case[CHOICE]['tht'][0])+""" """+str(d_case[CHOICE]['rv'][0])+"""	! altitude, tht, rv
    """+str(d_case[CHOICE]['Zlevel'][1])+""". """+str(d_case[CHOICE]['tht'][1])+""" """+str(d_case[CHOICE]['rv'][1])+""" 	

    ZFRC !z u v theta q w dtheta/dt dq/dt du/dt dv/dt
    1   
    2000  01 01 0.
    0.0     
    """+str(d_case[CHOICE]['Psurf'])+""".  
    """+str(d_case[CHOICE]['tht'][0])+"""
    0. 
    2
    10.     """+d_case[CHOICE]['U']+""" """+d_case[CHOICE]['V']+"""  0   0   -0.000    0.    0.0000e-08 0. 0.
    1000.     """+d_case[CHOICE]['U']+""" """+d_case[CHOICE]['V']+"""  0   0   -0.0000   -0.     0.00000e-08 0. 0."""

    TXT_00_run = """#!/bin/bash
    #SBATCH -J MNH00
    #SBATCH -N """+str(dic_nodes['00'])+"""           # nodes number (=NBP)   
    #SBATCH -n """+str(dic_cpu['00'])+"""         # CPUs number (on all nodes) (=NBP*TPN) 
    #SBATCH -o MNH00.eo%j   #
    #SBATCH -e MNH00.eo%j   #
    #SBATCH -t """+dic_timelimit['00']+"""   # time limit
    # Echo des commandes
    set -x
    ulimit -c 0
    ulimit -s unlimited

    # Nom de la machine
    hostname

    . """+profile_mnh+"""
    # Arrete du job des la premiere erreur
    set -e

    export MPIRUN='"""+mpicmd+""" -np """+str(dic_cpu['00'])+"""'
    export TIMEOUT="timelimit -t 7200"

    ls -lrt

    # PREP_IDEAL of dad domain
    """+time_cmd+"""${MPIRUN} PREP_IDEAL_CASE${XYZ}

    rm -f file_for_xtransfer pipe_name

    ls -lrt
    echo "Exit status is" $? """

    return TXT_00_PREP_IDEAL,TXT_00_run