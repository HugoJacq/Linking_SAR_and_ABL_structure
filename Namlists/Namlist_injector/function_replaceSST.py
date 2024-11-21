""" Function file that work with replaceSST.py"""

import xarray as xr
import numpy as np

def linearSST(LX,THT0,x1,x2,deltaTHT,L1,L2):
	# 2 linear ramps  of length L1 and L2
	# and SST is SST[0] = SST[-1] cyclic
	# x1 and x2 are the starts of the ramps
	SST=np.zeros(len(LX))
	if x1>LX[-1] or x2>LX[-1]:
		raise Exception('x1 or x2 is outside the domain')
	for k,x in enumerate(LX):
		if x >= x1 and x < x1+L1:
			SST[k]=THT0 + deltaTHT/L1*(x-x1)
		elif x >= x1+L1 and x < x2:
			SST[k]= THT0 + deltaTHT
		elif x >= x2 and x < x2+L2:
			SST[k]= THT0 + deltaTHT - deltaTHT/L2*(x-x2)
		else:
			SST[k]=THT0
	return SST

def tanhSST(x,THT0,x1,x2,deltaTHT,L1,L2): 
	# tanh centered on x1+L1/2 and x2+L2/2
	SST= THT0 + deltaTHT/2*( np.tanh( (x- (x1+L1/2))/L1) - np.tanh( (x- (x2+L2/2))/L2) ) 
	return SST

def buildSST(Dir,SST1D,NX,NY):
	SST=np.zeros((NY,NX))
	if Dir=='X':
		for j in range(NY):
			SST[j,:] = SST1D
	elif Dir=='Y':
		for i in range(NX):
			SST[:,i] = SST1D
	else:
		raise Exception('NX or NY doesnt match SST1D shape')
	return SST
	
def One_Step_SST(dim,atX,Lx,T0,deltaT,form):
	"""Produce an array that is a step of SST at 'atX' position,
	 wide of 'Lx' and of amplitude 'deltaT'. It has the shape 'form'.
	"""
	SST=np.zeros(len(dim))
	if atX>dim[-1] or atX>dim[-1]:
		raise Exception('atX is outside the domain')
	if form=='tanh':
		SST= T0 + deltaT/2*(1+np.tanh( (dim- (atX+Lx/2))/Lx))
	elif form=='linear':
		print('tbd')
	else:
		raise Exception("you choice '"+form+"' is not yet coded")
	
	return SST
		
		
			
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
