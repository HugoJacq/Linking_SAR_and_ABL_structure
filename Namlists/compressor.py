import os

PATH_OUT = './FICHIERS_OUT/'
liste=os.listdir(PATH_OUT)
N = len(liste)

for k in range(0,N):
	
	if k<10:
		num = '00'+str(k)
	elif k>=10 and k<100:
		num = '0'+str(k)
	else:
		num = str(k)
	oldname = liste[k]
	newname = liste[k][:-3]+'_comp.nc'
	print(k+1,'/',N)
	print('ncks -O -7 -L 1 --ppc default=4#RVT=5#THT=6#PABST=7 '+PATH_OUT+oldname+' '+PATH_OUT+newname)
	os.system('ncks -O -7 -L 1 --ppc default=4#RVT=5#THT=6#PABST=7 '+PATH_OUT+oldname+' '+PATH_OUT+newname)
