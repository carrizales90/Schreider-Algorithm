from datetime import datetime as dt
import matplotlib.pyplot as plt
from haversine import haversine
from datetime import timedelta
import datetime as DDTT
import warnings, csv
import numpy as np
#########################################################################################################################################################
font0 = {'family':'serif','color':'black','size':20}
date_format = "%Y-%m-%d %H:%M:%S"
warnings.filterwarnings("ignore")
#########################################################################################################################################################
def distance(Z1, lon1, lat1, Z2, lon2, lat2): 
    """
    This function computes the euclidean distance between two hypocenters (considering the earth as a sphere of radii 6,371 km)
    """
    Rt = 6371
    r1 = Rt - abs(Z1); phi1 = (90 - lat1) * (np.pi/180); theta1 = lon1 * (np.pi/180)
    r2 = Rt - abs(Z2); phi2 = (90 - lat2) * (np.pi/180); theta2 = lon2 * (np.pi/180)
    x1 = r1*np.cos(phi1)*np.sin(theta1); y1 = r1*np.sin(phi1)*np.sin(theta1); z1 = r1*np.cos(theta1)
    x2= r2*np.cos(phi2)*np.sin(theta2); y2 = r2*np.sin(phi2)*np.sin(theta2); z2 = r2*np.cos(theta2)

    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
#########################################################################################################################################################
def ff(n,s): 
    """
    This function computes the Laplace function for the smoothness procedure suggested by Schreider for his algorithm (1990)
    """
    A = 1/(s * np.sqrt(2 * np.pi))
    return A * np.exp( -(n**2/(2 * s**2)) )
#########################################################################################################################################################
'''
Inputs for Schreider algorithm 
'''
Mc = 4.4 #-----------------### Mc = Completeness Magnitude
Ti = '2000-01-01 00:00:00' ### Ti = Initial time (Date and Hour) for computation of the Schreider algorithm
Tf = '2025-01-01 00:00:00' ###Tf = Final time (Date and Hour) for computation of the Schreider algorithm
Center_Latitude = 14.761 #-### Center_Latitude = Latitude at which is centered the circular region of study
Center_Longitude = -94.103 ### Center_Longitude = Longitude at which is centered the circular region of study
Minimum_Depth = 30 #-------### Minimum_Depth = inferior threshold at which the depth of earthquakes are considered (the closest to the surface)
Maximum_Depth = 90 #-------### Maximum_Depth = superior threshold at which the depth of earthquakes are considered (the farthest from the surface)
Radius = 200 #-------------### Radius = radius of the circular region of study
#########################################################################################################################################################
'''
The data is loaded, filtered (according to the Inputs used) and saved in independent lists
'''
Date=[]; Hour=[]; Time=[]; Mag=[];Lat=[];Lon=[]; Depth=[]
with open('DataBase_SSN_1990-2024.dat') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for i in readCSV:
        lat = float(i[3]); lon = float(i[4]); t = dt.strptime(i[0]+' '+i[1], date_format); d = float(i[5]); m = float(i[2])
        if t >= dt.strptime(Ti, date_format) and t < dt.strptime(Tf, date_format) and  haversine((lat, lon),(Center_Latitude, Center_Longitude)) <= Radius and d >= abs(Minimum_Depth) and d <= abs(Maximum_Depth) and m >= Mc:
            Date.append(i[0]); Hour.append(i[1]); Mag.append(m); Lat.append(lat); Lon.append(lon); Depth.append(d); Time.append(dt.strptime(i[0]+' '+i[1], date_format))
#########################################################################################################################################################
'''
the smoothness parameter is chosen s=9 (see section 2.1 of the publication)
the limit l of the Laplace function f(n,s) is determined when f(n,s) < 10^{-2} (see Eq. (1) of the publication)
'''
for n in range(1000):
	if ff(n,9) < 10**(-2):
		l = n
		break
#########################################################################################################################################################
'''
Convolutions are computed (see section of Methods of the publication)
'''
Plot_Time=[]; T=[]; R=[]
for i in range(len(Time)-1):
    delta_t = DDTT.timedelta.total_seconds(Time[i+1]-Time[i])
    if delta_t > 0:
        Plot_Time.append(Time[i+1]); T.append(delta_t); R.append(distance(Depth[i+1], Lon[i+1], Lat[i+1] , Depth[i], Lon[i], Lat[i]))
V = np.log10(np.divide(R,T)); RT = np.multiply( np.divide(R, np.std(R)) , np.divide(T, np.std(T)) )
T_prima = T; T2_prima=V; T3_prima = RT; N = len(T); Tk = []; Tk2 = []; Tk3 = []

for i in range(N):
    preTk = []; preTk2 = []; preTk3 = []
    for n in range(l):
        if n <= i:
            preTk.append(T_prima[i-n] * ff(n,9)); preTk2.append(T2_prima[i-n] * ff(n,9)); preTk3.append(T3_prima[i-n] * ff(n,9))
    Tk.append(sum(preTk)); Tk2.append(sum(preTk2)); Tk3.append(sum(preTk3))
#########################################################################################################################################################
'''
Convolutions are plotted following the format of the results shown in sections 4 and 5 of the publicaction
'''
plt.figure(figsize=(8, 5))
TSigma=3
plt.subplot(3,1,1)
signal = np.divide(Tk, np.std(Tk))
MEAN = np.mean(signal); STD=np.std(signal) 
plt.plot(Plot_Time, signal, linewidth=1)
plt.hlines(y = MEAN, xmin=dt.strptime(Ti, date_format), xmax=dt.strptime(Tf, date_format), colors='r', linestyles='--', label=r'$\mu$')
plt.hlines(y = MEAN + TSigma*STD, xmin=dt.strptime(Ti, date_format), xmax=dt.strptime(Tf, date_format), colors='orange', linestyles=':', label=r'$'+str(TSigma)+'\\sigma$')
plt.legend(loc='center right', fontsize=9, ncol=2)
plt.grid(linestyle = '--', linewidth = 0.5)
plt.ylabel(r'$\widetilde{\mathcal{T}}(k)$', fontdict=font0)

plt.subplot(3,1,2)
II=10 
MEAN = np.mean(Tk2[II:]); STD=np.std(Tk2[II:]) 
plt.plot(Plot_Time[II:], Tk2[II:], linewidth=1)
plt.hlines(y = MEAN, xmin=dt.strptime(Ti, date_format), xmax=dt.strptime(Tf, date_format), colors='r', linestyles='--', label=r'$\mu$')     
plt.grid(linestyle = '--', linewidth = 0.5)
plt.ylabel(r'$\mathcal{V}(k)$', fontdict=font0)

plt.subplot(3,1,3)
MEAN = np.mean(Tk3); STD=np.std(Tk3) 
plt.plot(Plot_Time, Tk3, linewidth=1)
plt.hlines(y=MEAN, xmin=dt.strptime(Ti, date_format), xmax=dt.strptime(Tf, date_format), colors='r', linestyles='--', label=r'$\mu$')
plt.hlines(y= MEAN + TSigma*STD, xmin=dt.strptime(Ti, date_format), xmax=dt.strptime(Tf, date_format), colors='orange', linestyles=':', label=r'$'+str(TSigma)+'\\sigma$')
plt.grid(linestyle = '--', linewidth = 0.5)
plt.ylabel(r'$\widetilde{\mathcal{RT}}(k)$', fontdict=font0)

plt.tight_layout()
plt.subplots_adjust(hspace=0)
'''plt.savefig('Example_Quiescence_EQ8.2.png')''' #for saving figure, uncomment function plt.savefig()
plt.show()