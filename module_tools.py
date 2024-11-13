# To be used with analyse.py 
import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from module_cst import *
from math import factorial
from numpy.fft import fft, ifft
import netCDF4
import matplotlib.ticker as ticker

def nearest(array,value):
	"""
	Array is 1D
	value is 0D
	"""
	return np.argmin(np.abs(array-value))	

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return a.flat[idx]

def find_indx_indy_from_2D_LatLon(lat,lon,O):
	"""
	lat,lon are function of 2 variables : lat(x,y) for eg
	O is a tuple with coordinate to look for

	return : a tuple of index that gives lon(indx,indy),lat(indx,indy) = O
	"""
	lonmin = O[0] - 0.001
	lonmax = O[0] + 0.001
	latmin = O[1] - 0.001
	latmax = O[1] + 0.001
	lat_mask = np.ma.masked_where((lat<latmin) | (lat>latmax) | (lon<lonmin)  | (lon>lonmax)  ,lat)
	lon_mask = np.ma.masked_where((lat<latmin) | (lat>latmax) | (lon<lonmin)  | (lon>lonmax)   ,lon)
	#
	if len(lat_mask.nonzero()[0])==0:
		indx,indy = None,None 
	else:
		indx,indy = (lat_mask.nonzero()[0][0],lat_mask.nonzero()[1][0]) # indices lon(ix0,ix1) lat(ix0,ix1) du point central 

	return indx,indy
def T_to_Theta(T,P,P0=P00): 
	#Renvoie la température potentielle en K
	#IN : Température (K)
	#     Pression (Pa)
	#OUT : Theta (K)	
	theta=T*(P0/P)**(Rd/Cpd) #K
	return theta
	
def Theta_to_T(Theta,P,P0=P00): 
	#Renvoie la température en K
	#IN : Température potentielle (K)
	#     Pression (Pa)
	#OUT : Température (K)	
	temp=Theta*(P/P0)**(Rd/Cpd) #K
	return temp

def sat_vap_pressure(T):
	# compute saturation vapor pressure
	return e0*np.exp( Lv/Rv * (1/273.15 - 1/T) )

def Td(rv,P):
	# compute dew point temperature = 'point de rosé'
	# A saturation mixing ratio can be deduced from this on emagram
	#	if you follow iso rvsat and crosses the dry adiabatic, you get LCL
	# P in Pa
	return 1/( 1/273.15 - Rv/Lv * np.ln(rv*P/(e0*(rv+Rd/Rv))) )
	
def compute_rv_s(T,P):
	# Compute saturation mixing ratio at T and P
	# P in Pa
	es = sat_vap_pressure(T)
	return Rd/Rv * es / (P-es)

def RH_from_Td_T(Td,T):
	# RH from Td and T
	return 100*np.exp( - Lv/(Rv*T*Td) * (T-Td))

def rv_from_Td(Td,T,P):
	# Td : dew point temperature
	# T : dry temperature
	# P : Pressure
	# source : http://climate.envsci.rutgers.edu/pdf/LawrenceRHdewpointBAMS.pdf
	Psat = sat_vap_pressure(T)
	HR = RH_from_Td_T(Td,T)
	Pv = Psat * HR/100
	rv = 0.622/(P/Pv-1) 
	return rv
	
def qtoRH(RVT,THT,P):
	# Calcul l'humidité relative
	# Tiré du code de MNH
	# Constantes (déclarée : modd_cst.f90 , initialisée : ini_cst.f90)	
	XBOLTZ      = 1.380658E-23		# Boltzman constant
	XAVOGADRO   = 6.0221367E+23		# Avogadro number
	XMD    = 28.9644E-3			# Molar mass of dry air
	XMV    = 18.0153E-3			# Molar mass of vapor
	#XRD    = XAVOGADRO * XBOLTZ / XMD	# Gaz constant for dry air
	XRV    = XAVOGADRO * XBOLTZ / XMV	# Gaz constant for vapor
	#XCPD   = 7.* XRD /2.			# Cpd (dry air)
	XCPV   = 4.* XRV			# Cpv (vapor)
	#XRHOLW = 1000.				# Volumic mass of liquid water
	#XRHOLI = 900.				# Volumic mass of ice
	#XCONDI = 2.22				# Thermal conductivity of ice (W m-1 K-1)
	XCL    = 4.218E+3			# Cl (liquid)
	#XCI    = 2.106E+3			# Ci (ice)
	XTT    = 273.16				# Triple point temperature
	XLVTT  = 2.5008E+6			# Vaporization heat constant
	#XLSTT  = 2.8345E+6			# Sublimation heat constant
	#XLMTT  = XLSTT - XLVTT			# Melting heat constant
	XESTT  = 611.14				# Saturation vapor pressure  at triple point
	# Constants for saturation vapor
	XGAMW  = (XCL - XCPV) / XRV
	XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
	XALPW  = np.log(XESTT) + (XBETAW /XTT) + (XGAMW *np.log(XTT))
	# Conversion Theta en T
	T = Theta_to_T(THT,P)
	# Calcul Pression de saturation (compute_function_thermo.f90)
	PSAT = np.exp( XALPW - XBETAW/T - XGAMW*np.log(T))
	SATmr = (XMV/XMD)*PSAT/(P-PSAT)
	# Humidité relative
	RHU = RVT/SATmr
	return RHU

def Exner(P):
	P0 = 100000	# Pression de référence
	Cpd = 1004.71	# Capacité thermique à pression cste
	Rd = 287.05	# Constante des gaz parfaits air sec
	return (P/P0)**(Rd/Cpd)
	
def inv_Exner(Pi):
	P0 = 100000	# Pression de référence
	Cpd = 1004.71	# Capacité thermique à pression cste
	Rd = 287.05	# Constante des gaz parfaits air sec
	return P0*Pi**(Cpd/Rd)
	
def L_Obukhov(tht0,u_star,surf_wtht_flx):
	g = 9.81
	K = 0.41 # von karman
	return - tht0*u_star**3/( K*g*surf_wtht_flx)
	
def Compute_w_star(flx,zi):
	g=9.81 # gravity
	thtv=300 # K ref temperature
	return (g/thtv * flx*zi)**(1/3)
	
	
def Compute_THTV(THT,RVT):
	"""input can be DataArray
	"""
	return THT*(1+Rv/Rd*RVT)/(1+RVT)
	
def PSD(time_vect, signal_vect):
    """This function automates the computation of the Power Spectral Density of a signal.
    """
    # Same as amplitude spectrum START ===================
    total_time = time_vect[-1] - time_vect[0]
    dt = time_vect[1] - time_vect[0]
    freq_resolution = 1.0 / total_time
    freq_nyquist = 1. / dt
    
    # samples
    N = len(time_vect)
    
    # build frequency
    frequency = np.arange(0, freq_nyquist, freq_resolution, dtype=float)
    frequency = frequency[: int(N / 2) - 1]  # limit to the first half
    
    raw_fft = np.fft.fft(signal_vect, norm="backward")
    raw_fft /= N
    
    # Takes half, but double the content, excepted the first component
    amplitude_spectrum = 2*np.absolute(raw_fft)[:int(N/2) - 1] 
    amplitude_spectrum[0] /= 2
    # Same as amplitude spectrum END===================
    
    power_spectral_density = 1/2 * np.absolute(
        amplitude_spectrum
        *np.conjugate(amplitude_spectrum)
    )/freq_resolution 
    
    power_spectral_density[0] *= 2
    
    return frequency, power_spectral_density

def detrended_PSD(time_vect,signal_vect,order=1):
	"""
	This is a function that computes the PSD with before a detrending of signal_vect
	By default the detrending is linear (removing ax+b)
	"""
	poly = np.polynomial.Polynomial.fit(time_vect,signal_vect,order)

	detrended = signal_vect - poly(time_vect)
	return PSD(time_vect,detrended)



def LambdaToK(x):
    return 2*np.pi/x

def Mean1(var):
    """DataArray from xarray
    
    Mean1 = average along time,nj and between indx1 and indx2 in ni
    """
    return var.mean(['ni','time','nj'])
def Mean2(var):
    return var.mean(['ni','nj'])


def Plot_CLI_sim_dim(ORIGIN,Lx,Ly,DegLat,DegLon):
    """
    Plot a nice view of simulation dimensions, in the given unit
    """
    print(np.round(ORIGIN[0],2),'°E ------------------------',np.round(ORIGIN[0]+Lx/DegLon,2),'°E')
    print(np.round(ORIGIN[1]+Ly/DegLat,2),'°N                         ',np.round(ORIGIN[1]+Ly/DegLat,2),'°N')
    print('|                                   |')
    print('|                                   |')
    print('|                                   |')
    print('|                                   |')
    print(np.round(ORIGIN[0],2),'°E ------------------------',np.round(ORIGIN[0]+Lx/DegLon,2),'°E')
    print(np.round(ORIGIN[1],2),'°N                         ',np.round(ORIGIN[1],2),'°N')
    print('')
    
	
def give_NETCDF4_chunking(path):
	"""
	Return chunks from the file path
	This information can be used to set dask chunks in a coherent way.
	"""
	data = netCDF4.Dataset(path)
	print(data.variables['UT'].chunking())


def sec2hms(ss):
	(hh, ss)=divmod(ss, 3600)
	(mm, ss)=divmod(ss, 60)
	time = (hh, mm, np.round(ss,0))
	return str(time[0])+'h '+str(time[1])+'min '+str(time[2])+'sec'
	
	
	
def contour_fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return s
	#return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"
	
def get_ABLH(Z,THTV):
	"""
	Return the ABL height given a THTV profile, zi = altitude at max dthtv/dz

	THTV is xarray DataArray

	"""
	dimname = THTV.dims[0]

	indz = np.argmax(THTV.differentiate(dimname).values)
	return Z[indz].values
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
