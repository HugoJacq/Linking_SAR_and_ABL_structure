def return_txt_04(LMNH_MPI_BSEND,MPI_BUFFER,dic_nodes,dic_cpu,dic_timelimit,
                  profile_mnh,mpicmd,time_cmd,
                  dic_PROCIO_R,dic_PROCIO_W):

    TXT_04_PRE_REAL = """&NAM_REAL_CONF NVERB = 10/
&NAM_CONFIO LCDF4=.TRUE.,LLFIOUT=.FALSE., LLFIREAD=.FALSE. /

&NAM_CONFZ LMNH_MPI_BSEND = """+LMNH_MPI_BSEND['04']+""",
           MPI_BUFFER_SIZE="""+str(MPI_BUFFER['04'])+""",
           NZ_VERB=0,
           NB_PROCIO_R="""+str(dic_PROCIO_R['04'])+""",
           NB_PROCIO_W="""+str(dic_PROCIO_W['04'])+""" /

&NAM_VER_GRID YZGRID_TYPE='SAMEGR'/

&NAM_FILE_NAMES HATMFILE ='ATM_SON',
                HATMFILETYPE='MESONH',
                HPGDFILE ='ATM_SON_AS_PGD',
                CINIFILE='INIT_SON' /
                
&NAM_PREP_SURF_ATM CFILE = 'ATM_SON',
                   CFILETYPE = 'MESONH',
                   CFILEPGD = 'ATM_SON_AS_PGD',
                  CFILEPGDTYPE = 'MESONH' /
"""

    TXT_04_RUN = """#!/bin/bash
#SBATCH -J MNH04
#SBATCH -N """+str(dic_nodes['04'])+"""           # nodes number (=NBP)   
#SBATCH -n """+str(dic_cpu['04'])+"""           # CPUs number (on all nodes) (=NBP*TPN) 
#SBATCH -o MNH04.eo%j   #
#SBATCH -e MNH04.eo%j   #
#SBATCH -t """+dic_timelimit['04']+"""    # time limit

# Echo des commandes

ulimit -c 0
ulimit -s unlimited

# Nom de la machine
hostname


#. /home/jacqhugo/MNH-V5-7-0/conf/profile_mesonh-LXgfortran-R8I4-MNH-V5-7-0-MPIAUTO-DEBUG
. """+profile_mnh+"""

# Arrete du job des la premiere erreur
set -e

export MPIRUN='"""+mpicmd+""" -np """+str(dic_cpu['04'])+"""'
export TIMEOUT="timelimit -t 7200"



#ln -sf ../03_spawning/DAD_S.1.001.001.spa00.nc ATM_SON.nc
#ln -sf ../03_spawning/DAD_S.1.001.001.spa00.des ATM_SON.des
#cp ATM_SON.nc ATM_SON_AS_PGD.nc # i need to change the name, idk why but else it crash at next step
#cp ATM_SON.des ATM_SON_AS_PGD.des 

export SPW_NAME='DAD_S.1.001.001.spa00'
export NEW_SPWN_NAME='ATM_SON'
export SPW_AS_PGD='ATM_SON_AS_PGD'
export number_zsplit="""+str(dic_PROCIO_R['04'])+"""

ln -sf ../03_spawning/${SPW_NAME}.des $NEW_SPWN_NAME.des
ln -sf ../03_spawning/${SPW_NAME}.nc $NEW_SPWN_NAME.nc
ln -sf ../03_spawning/${SPW_NAME}.des $SPW_AS_PGD.des
ln -sf ../03_spawning/${SPW_NAME}.nc $SPW_AS_PGD.nc

if  [ $number_zsplit -gt 1 ] # Z splitting is used if true.
then
	for count in $(seq $number_zsplit)
	do
		if [ $count -lt 10 ]
		then
			add=".Z00"$count".nc"
		elif [ $count -lt 100 -a $count -ge 10 ]
		then
			add=".Z0"$count".nc"	
		fi
		ln -sf ../03_spawning/${SPW_NAME}$add $NEW_SPWN_NAME$add
		#ln -sf ../03_spawning/${SPW_NAME}$add $SPW_AS_PGD$add
	done	
fi

ls -lrt

"""+time_cmd+"""${MPIRUN} PREP_REAL_CASE${XYZ}

rm -f file_for_xtransfer pipe_name

ls -lrt
echo "Exit status is" $? """

    return TXT_04_PRE_REAL,TXT_04_RUN