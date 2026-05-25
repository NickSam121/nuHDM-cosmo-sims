import numpy as np
import matplotlib.pyplot as plt

class_planck  = np.genfromtxt(fname="/home/nsamaras/Desktop/astro/nuHDM_AUUK_desktop_folder/class/output/CLASS_planck_201804_background.dat", delimiter="")
class_nu  = np.genfromtxt(fname="/home/nsamaras/Desktop/astro/nuHDM_AUUK_desktop_folder/class/output/nu20_background.dat", delimiter="")
class_opt_nu  = np.genfromtxt(fname="/home/nsamaras/Desktop/astro/nuHDM_AUUK_desktop_folder/class/output/opt_nu03_background.dat", delimiter="")

z = class_planck[:,0]
a = 1/(1+z)
conf_time = class_planck[:,2] 
dens_rad = class_planck[:,8] 
dens_b = class_planck[:,9] 
dens_cdm = class_planck[:,10] 
dens_crit = class_planck[:,13] 
dens_lambda = class_planck[:,11] 

nu_z = class_nu[:,0]
nu_a = 1/(1+nu_z)
nu_conf_time = class_nu[:,2] 
nu_dens_rad = class_nu[:,8] 
nu_dens_b = class_nu[:,9] 
nu_dens_hdm =class_nu[:,10] 
nu_dens_crit = class_nu[:,14] 
nu_dens_lambda = class_nu[:,12] 

opt_nu_z = class_opt_nu[:,0]
opt_nu_a = 1/(1+opt_nu_z)
opt_nu_conf_time = class_opt_nu[:,2] 
opt_nu_dens_rad = class_opt_nu[:,8] 
opt_nu_dens_b = class_opt_nu[:,9] 
opt_nu_dens_hdm = class_opt_nu[:,10] 
opt_nu_dens_crit = class_opt_nu[:,14] 
opt_nu_dens_lambda = class_opt_nu[:,12] 
###########################################################################

###########################################################################
fig, ax = plt.subplots(figsize=(11, 11))
font1 = {'family':'serif','color':'black','size':30}
#plt.rcParams["figure.figsize"] = [6*4., 4*4.]

plt.plot(a, dens_rad/dens_crit,  label = r"$\Lambda$CDM - rad", color = "blue", linestyle = "-", linewidth= 2)
plt.plot(a, dens_b/dens_crit,  label = r"$\Lambda$CDM - bar", color = "green", linestyle = "-", linewidth= 2)
plt.plot(a, dens_cdm/dens_crit,  label = r"$\Lambda$CDM - CDM", color = "red", linestyle = "-", linewidth= 2)
plt.plot(a, dens_lambda/dens_crit,  label = r"$\Lambda$CDM - $\Lambda$", color = "pink", linestyle = "-", linewidth= 2)

plt.plot(nu_a, nu_dens_rad/nu_dens_crit,  label = r"$\nu$HDM - rad", color = "blue", linestyle = "--", linewidth= 4.5)
plt.plot(nu_a, nu_dens_b/nu_dens_crit,  label = r"$\nu$HDM - bar", color = "green", linestyle = "--", linewidth= 4.5)
plt.plot(nu_a, nu_dens_hdm/nu_dens_crit,  label = r"$\nu$HDM - HDM", color = "red", linestyle = "--", linewidth= 4.5)
plt.plot(nu_a, nu_dens_lambda/nu_dens_crit,  label = r"$\nu$HDM - $\Lambda$", color = "pink", linestyle = "--", linewidth= 4.5)

plt.plot(opt_nu_a, opt_nu_dens_rad/opt_nu_dens_crit,  label = r"opt-$\nu$HDM - rad", color = "blue", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_a, opt_nu_dens_b/opt_nu_dens_crit,  label = r"opt-$\nu$HDM - bar", color = "green", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_a, opt_nu_dens_hdm/opt_nu_dens_crit,  label = r"opt-$\nu$HDM - HDM", color = "red", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_a, opt_nu_dens_lambda/opt_nu_dens_crit,  label = r"opt-$\nu$HDM - $\Lambda$", color = "pink", linestyle = "dotted", linewidth= 7)

plt.axvline(x = 0.005, color = 'black', label = 'Hydro sims', linewidth=2)
plt.axvline(x = 0.0003, color = 'grey', label = 'rad-mat eq', linewidth=2)

plt.plot(np.zeros(1), np.zeros([1,2]), color='w', alpha=0, label=' ')

plt.xlabel(r"$\alpha/\alpha_{today}$", fontsize=30)
plt.ylabel(r"$\Omega_i$",  rotation=0, fontsize=32)
ax.yaxis.set_label_coords(-.12, .5)

#plt.legend(loc="best", fontsize = 20)
#plt.legend(bbox_to_anchor=(0.5, 1), loc='center', fontsize=15)

plt.legend(
    ncol=4,
    loc='upper center',
    bbox_to_anchor=(0.5, 1.32),
    fontsize=19,
    frameon=True,
    edgecolor='black',
    facecolor='white',
    columnspacing=1.2
)

#, 0.03,0.08

plt.xticks(fontsize=32)
plt.yticks(fontsize=32)
plt.xlim([10**-9,1])
plt.ylim([10**-4, 1])
plt.grid()
plt.xscale('log')
plt.yscale('log')
#ax.yaxis.set_label_coords(-.1, -.9)

#plt.tight_layout(h_pad=-3)
#plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.2, wspace=0.1)

plt.tight_layout()
#plt.subplots_adjust(top=0.7)
plt.show()

###########################################################################
'''
###########################################################################
plt.plot(conf_time, dens_rad,  label = r"$\Lambda$CDM - rad", color = "blue", linestyle = "-", linewidth= 2)
plt.plot(conf_time, dens_b,  label = r"$\Lambda$CDM - bar", color = "green", linestyle = "-", linewidth= 2)
plt.plot(conf_time, dens_cdm,  label = r"$\Lambda$CDM - CDM", color = "red", linestyle = "-", linewidth= 2)
plt.plot(conf_time, dens_lambda,  label =r"$\Lambda$CDM - $\Lambda$", color = "pink", linestyle = "-", linewidth= 2)

plt.plot(nu_conf_time, nu_dens_rad,  label = r"$\nu$HDM - rad", color = "blue", linestyle = "--", linewidth=4.5)
plt.plot(nu_conf_time, nu_dens_b,  label = r"$\nu$HDM - bar", color = "green", linestyle = "--", linewidth= 4.5)
plt.plot(nu_conf_time, nu_dens_hdm,  label = r"$\nu$HDM - HDM", color = "red", linestyle = "--", linewidth= 4.5)
plt.plot(nu_conf_time, nu_dens_lambda,  label = r"$\nu$HDM - $\Lambda$", color = "pink", linestyle = "--", linewidth= 4.5)

plt.plot(opt_nu_conf_time, opt_nu_dens_rad,  label = r"opt-$\nu$HDM - rad", color = "blue", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_conf_time, opt_nu_dens_b, label = r"opt-$\nu$HDM - bar", color = "green", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_conf_time, opt_nu_dens_hdm,  label = r"opt-$\nu$HDM - HDM", color = "red", linestyle = "dotted", linewidth= 7)
plt.plot(opt_nu_conf_time, opt_nu_dens_lambda,  label = r"opt-$\nu$HDM - $\Lambda$", color = "pink", linestyle = "dotted", linewidth= 7)

plt.xlabel("conformal time [Mpc]", fontdict=font1)
plt.ylabel(r"Densities $[Mpc^2]$", fontdict=font1)
plt.legend(loc="best", fontsize = 18.4)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlim([1,10**4])
plt.ylim([10**-10,10**14])
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.show()
###########################################################################
'''