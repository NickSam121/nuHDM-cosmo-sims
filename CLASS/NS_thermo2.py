import numpy as np
import matplotlib.pyplot as plt

class_planck  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/CLASS_planck_201803_thermodynamics.dat", delimiter="")
class_nu  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/nu20_thermodynamics.dat", delimiter="")
class_opt_nu  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/opt_nu03_thermodynamics.dat", delimiter="")

pl1 = class_planck[:,0] 
pl2 = class_planck[:,6] 

nu_1 = class_nu[:,0] 
nu_2= class_nu[:,6] 

opt_nu_1 = class_opt_nu[:,0] 
opt_nu_2= class_opt_nu[:,6] 

font1 = {'family':'serif','color':'black','size':35}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]

#plt.plot(pl1, pl2,  label = r"$x_e$ - $\Lambda$CDM", color = "blue", linestyle = "-", linewidth= 10)
#plt.plot(nu_1, nu_2,  label = r"$x_e$ - $\nu$HDM", color = "red", linestyle = "-", linewidth= 3)
#plt.plot(nu_1, opt_nu_2,  label = r"$x_e$ - opt-$\nu$HDM", color = "green", linestyle = "dotted", linewidth= 7)

plt.plot(pl1, pl2,  label = r"$\Lambda$CDM", color = "blue", linestyle = "-", linewidth= 15)
plt.plot(nu_1, nu_2,  label = r"$\nu$HDM", color = "red", linestyle = "-", linewidth= 3.5)
plt.plot(nu_1, opt_nu_2,  label = r"opt-$\nu$HDM", color = "green", linestyle = "dotted", linewidth= 7)

#plt.xlabel(r"$\tau $[Mpc]", fontdict=font1)
plt.xlabel("redshift z", fontdict=font1)

plt.ylabel("Matter Temperature [K]", fontdict=font1)
plt.legend(loc="best", fontsize = 35)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
#plt.xlim([1,10**4])
#plt.ylim([10**-4,10**1])
plt.grid()
plt.xscale('log')
plt.yscale('log')
#plt.title("CLASS", fontsize = 35, fontdict=font1)
plt.show()
