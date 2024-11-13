"""
Here constants are declared, to be used in other modules
"""
ABLH_S1 = 600 # m, unified ABLH
g = 9.80665 	# gravity m s-2
C2M = 0.2 	# Constant for transport/diffusion term TKE budget
Ceps = 0.85 	# Constant for dissipation terms TKE budget
f = -8.469e-5 # 2*Omega*sin LAT, ici LAT=-35.5°
f_star = 1.18732481e-4 # 2*Omega*cos LAT, ici LAT=-35.5°
Cpd = 1004.71 # J kg-1 K-1
Md = 28.9644e-3 # Masse molaire air sec
Mv = 18.0153e-3 # Masse molaire vapeur eau
Rd = 287.0597 # gas constant for dry air J kg-1 K-1
Rv = 461.5245 # gas constant for moist air J kg-1 K-1
Cvd = Cpd - Rd
P00 = 100000 # reference pressure (for Exner function) Pa
Psfx = 101250. # pressure at ocean surface
nhalo = 1
nu_air = 1.516e-5 # Kinematic air viscosity, m2/s	
Lv = 2.5*10**6 # J/kg, taken as constant but varies with T
e0 = 611.3 # Pa reference vapor pressure

# Conversion latlon to km at around -35°N
DegLat = 111.2 # km/°lat
DegLon = 91.6 # km/°lon

azimuth_sat = -17 # degree
SAR_res = 100 # m
LES_res = 50 # m

