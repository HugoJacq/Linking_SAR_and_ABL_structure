def return_txt_03(IXOR,IYOR,NI_s,NJ_s,dx_ratio,
                  LMNH_MPI_BSEND,MPI_BUFFER,dic_PROCIO_R,dic_PROCIO_W,
                  profile_mnh,dic_nodes,dic_cpu,dic_timelimit,mpicmd,time_cmd):

    TXT_03_SPAWN = """&NAM_CONFIO  LCDF4=T, LLFIOUT=F, LLFIREAD=F /
&NAM_GRID2_SPA   IXOR = """+str(int(IXOR))+""",
				 IYOR = """+str(int(IYOR))+""",
                 IXSIZE = """+str(NI_s)+""",
                 IYSIZE = """+str(NJ_s)+""",
                 IDXRATIO = """+str(dx_ratio)+""",
                 IDYRATIO = """+str(dx_ratio)+"""/
                 
&NAM_CONFZ LMNH_MPI_BSEND = """+LMNH_MPI_BSEND['03']+""",
           MPI_BUFFER_SIZE="""+str(MPI_BUFFER['03'])+""",
           NZ_VERB=0,
           NB_PROCIO_R="""+str(dic_PROCIO_R['03'])+""",
           NB_PROCIO_W="""+str(dic_PROCIO_W['03'])+""" /
           
&NAM_LUNIT2_SPA  CINIFILE = "DAD_S.1.001.001",
                 CINIFILEPGD="DAD_PGD",
                 YSPANBR='00'/
"""

    TXT_03_RUN = """#!/bin/bash
#SBATCH -J MNH03
#SBATCH -N """+str(dic_nodes['03'])+"""           # nodes number (=NBP)   
#SBATCH -n """+str(dic_cpu['03'])+"""           # CPUs number (on all nodes) (=NBP*TPN) 
#SBATCH -o MNH03.eo%j   #
#SBATCH -e MNH03.eo%j   #
#SBATCH -t """+dic_timelimit['03']+"""    # time limit

# Echo des commandes
set -x
ulimit -c 0
ulimit -s unlimited

# Nom de la machine
hostname

#. /home/cnrm_other/ge/mrmc/jacqueth/MNH-V5-6-0/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-6-0-MYMODIF1-MPIAUTO-DEBUG
#. /home/jacqhugo/MNH-V5-7-0/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-7-0-MPIAUTO-DEBUG
. """+profile_mnh+"""

# Arret du job des la premiere erreur
set -e
export MPIRUN='"""+mpicmd+""" -np """+str(dic_cpu['03'])+"""' 
export TIMEOUT="timelimit -t 7200" #259200=3j 86400=1j

# SPAWNING (horizontal interpolation ) ----------------------------------------
# links
ln -sf ../00_prep_ideal_dad/DAD_PGD* .
#ln -sf ../00_prep_ideal_dad/INIT_DAD* .
ln -sf ../02_run_dad_spinup/DAD_S.1.001.001.{nc,des} .

#rm -f file_for_xtransfer pipe_name

ls -lrt
"""+time_cmd+"""${MPIRUN} SPAWNING${XYZ}
#${TIMEOUT} time ${MPIRUN} SPAWNING${XYZ} # local cmd
#time ${MPIRUN} MESONH${XYZ} #belenos cmd
#${TIMEOUT} time -v ${MPIRUN} idrmem MESONH${XYZ} PAS EN PRODUCTION !!
mv -f OUTPUT_LISTING1  OUTPUT_LISTING1_spa
mv -f OUTPUT_LISTING0  OUTPUT_LISTING0_spa"""
    return TXT_03_SPAWN,TXT_03_RUN