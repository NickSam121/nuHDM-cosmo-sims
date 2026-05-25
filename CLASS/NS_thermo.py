import numpy as np
import matplotlib.pyplot as plt

class_planck  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/CLASS_planck_201803_thermodynamics.dat", delimiter="")
class_nu  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/nu20_thermodynamics.dat", delimiter="")
class_opt_nu  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/class/output/opt_nu03_thermodynamics.dat", delimiter="")

z = class_planck[:,0] 
xe= class_planck[:,2] 

nu_z = class_nu[:,0] 
nu_xe= class_nu[:,2] 

opt_nu_z = class_opt_nu[:,0] 
opt_nu_xe= class_opt_nu[:,2] 

font1 = {'family':'serif','color':'black','size':35}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]

plt.plot(z, xe,  label = r"$x_e$ - $\Lambda$CDM", color = "blue", linestyle = "-", linewidth= 10)
plt.plot(nu_z, nu_xe,  label = r"$x_e$ - $\nu$HDM", color = "red", linestyle = "-", linewidth= 4)
plt.plot(nu_z, opt_nu_xe,  label = r"$x_e$ - opt-$\nu$HDM", color = "green", linestyle = "dotted", linewidth= 7)

plt.xlabel("redshift z", fontdict=font1)
plt.ylabel(r"$x_e$", fontdict=font1)
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
