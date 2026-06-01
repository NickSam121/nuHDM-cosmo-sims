import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
from scipy.integrate import quad
import sympy as sy
import math
import sys

zs = np.arange(0,5,0.01)

def lookback(z, OmegaM, OmegaL, H0):   #Carroll, Press & Turner : https://ui.adsabs.harvard.edu/abs/1992ARA%26A..30..499C/abstract
    #return ((1+z)*np.sqrt((OmegaM*(1+z)**3)+OmegaL+(OmegaK*(1+z)**2)) )
    return ((1+OmegaM*z)*(1+z)**2 - z*(2+z)*OmegaL)**(-1/2)/(H0*(1+z))

Lb_values_lcdm = np.zeros(len(zs))
for (i,z) in enumerate(zs):
	K = quad(lookback, 0, z, args=(0.3153, 1-0.3153, 67.64))[0]
	Lb_values_lcdm[i] = K
	
Lb_values_lcdm2 = np.zeros(len(zs))
for (i,z) in enumerate(zs):
	K2 = quad(lookback, 0, z, args=(0.3153, 1-0.3153, 72.0))[0]
	Lb_values_lcdm2[i] = K2	
	
Lb_values_nuHDM = np.zeros(len(zs))
for (i,z) in enumerate(zs):
	L = quad(lookback, 0, z, args=(0.3153, 1-0.3153, 67.64))[0]
	Lb_values_nuHDM[i] = L
	
Lb_values = np.zeros(len(zs))
for (i,z) in enumerate(zs):
	M = quad(lookback, 0, z, args=(0.495, 1-0.495, 55.65))[0]
	Lb_values[i] = M

font1 = {'family':'serif','color':'black','size':35}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]
plt.plot(zs, 1000*Lb_values_lcdm, label = r"$\Lambda$CDM - $h = 0.67$", color = "blue", linestyle = "-", linewidth= 16)
plt.plot(zs, 1000*Lb_values_lcdm2, label = r"$\Lambda$CDM - $h = 0.72$", color = "lightblue", linestyle = "-", linewidth= 6)
plt.plot(zs, 1000*Lb_values_nuHDM, label = r"$\nu$HDM - $h = 0.67$", color = "red", linestyle = "--", linewidth= 4)
plt.plot(zs, 1000*Lb_values, label = r"opt-$\nu$HDM - $h = 0.55$", color = "green", linestyle = "-", linewidth= 6)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel("Redshift z", fontdict=font1)
plt.ylabel("Lookback Time [Gyr]", fontdict=font1)
plt.legend(loc="best", fontsize = 35)
plt.grid()
plt.show()
