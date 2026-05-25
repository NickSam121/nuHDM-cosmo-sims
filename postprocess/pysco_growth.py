import numpy as np
import matplotlib
matplotlib.use('TkAgg') # or 'Qt5Agg', 'Qt6Agg', 'WxAgg', etc.
import matplotlib.pyplot as plt

lcdm1 = np.loadtxt('/media//nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n128_planck2024/evolution_table_pysco.txt', dtype='float', delimiter=None)
lcdm_k1 = lcdm1[:,0]
lcdm_pk1 = lcdm1[:,1]
#lcdm_pk2 = lcdm1[:,5]

'''
lcdm2 = np.loadtxt('/media//nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n256_planck2024/evolution_table_pysco.txt', dtype='float', delimiter=None)
lcdm_k2 = lcdm2[:,0]
lcdm_pk2 = lcdm2[:,4]

lcdm3 = np.loadtxt('//media//nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n256_l8_planck2024/evolution_table_pysco.txt', dtype='float', delimiter=None)
lcdm_k3 = lcdm3[:,0]
lcdm_pk3 = lcdm3[:,4]

lcdm4 = np.loadtxt('/media//nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n256_l9_planck2024/evolution_table_pysco.txt', dtype='float', delimiter=None)
lcdm_k4 = lcdm4[:,0]
lcdm_pk4 = lcdm4[:,4]

'''
opt_nu199 = np.loadtxt('/media/nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n256_l9_opt_nuHDM_z0/evolution_table_pysco.txt', dtype='float', delimiter=None)
opt_nu_a  = opt_nu199[:,0]         
opt_nu_f1 = opt_nu199[:,1]
#opt_nu_f2 = opt_nu199[:,5]

nu199 = np.loadtxt('/media/nikolaos-samaras/T7/pysco+sterile_dm/pysco/examples/b200_n256_l9_nuHDM_z0/evolution_table_pysco.txt', dtype='float', delimiter=None)
nu_a = nu199[:,0]         
nu_f1 = nu199[:,1]
#nu_f2 = nu199[:,5]

def redshift_2_scale(z):
    return 1.0/(z+1)

def scale_2_redshift(a):
    return (1.0/a)-1

font1 = {'family':'serif','color':'black','size':25}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]

fig, ax1 = plt.subplots()	

dex_x = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9,1.0]
#dex_y = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
#dex_x2 = [16, 12, 10, 8, 6, 4 ,3, 2, 1, 0]
dex_y2 =  [10, 4]

ax1.set_xticks(dex_x)
ax1.set_xticklabels(dex_x, fontsize=25)
#ax1.set_yticks(dex_y)
#ax1.set_yticklabels(dex_y, fontsize=25)
secax = ax1.secondary_xaxis('top', functions=(scale_2_redshift, redshift_2_scale))
secax.set_xlabel("1+z", fontsize=25)
secax.set_xticks(dex_x)
secax.set_xticklabels(dex_x, fontsize=10.2)

ax1.set_xlabel("Scalefactor", fontdict=font1)
ax1.set_ylabel(r"$H/H_0$", fontdict=font1)

#ax1.set_ylabel(r"Growth factor 1st order", fontdict=font1)

ax1.plot(lcdm_k1, lcdm_pk1, label = r"1st order $\Lambda$CDM (h = 0.69, $\Omega_m$ = 0.3, $\omega_{CDM} = 0.12$)", color = "blue", linewidth= 10.0)
#ax1.plot(lcdm_k1, lcdm_pk2, label = r"2nd order $\Lambda$CDM (h = 0.69, $\Omega_m$ = 0.3, $\omega_{CDM} = 0.12$)", color = "blue", linewidth= 10.0, linestyle  = "dashed" )

ax1.plot(nu_a, nu_f1, label = r"1st order $\nu$HDM (h = 0.69, $\Omega_m = 0.3$, $\omega_{\nu} = 0.12$)", color = "red", linewidth = 5.0)
#ax1.plot(nu_a, nu_f2, label = r"2nd order $\nu$HDM (h = 0.69, $\Omega_m = 0.3$, $\omega_{\nu} = 0.12$)", color = "red", linewidth = 5.0, linestyle  = "dashed")

ax1.plot(opt_nu_a, opt_nu_f1, label = r"1st order opt-$\nu$HDM (h = 0.55, $\Omega_m = 0.49$, $\omega_{\nu} = 0.13$)", color = "green", linewidth= 5.0)
#ax1.plot(opt_nu_a, opt_nu_f2, label = r"2nd order opt-$\nu$HDM (h = 0.55, $\Omega_m = 0.49$, $\omega_{\nu} = 0.13$)", color = "green", linewidth= 5.0, linestyle  = "dashed")

#ax1.plot(lcdm_k2, lcdm_pk2, label = r"$\Lambda$CDM - pysco - b200l7_256", color = "black", linewidth= 3.5)
#ax1.plot(lcdm_k3, lcdm_pk3, label = r"$\Lambda$CDM - pysco - b200l8_256", color = "green", linewidth= 1)
#ax1.plot(lcdm_k4, lcdm_pk4, label = r"$\Lambda$CDM - pysco - b200l9_256", color = "red", linewidth= .5)

#ax1.plot(opt_nu_a, opt_nu_f1, label = r"opt-$\nu$HDM", color = "green", linewidth= 3.0)
#ax1.plot(nu_a, nu_f1, label = r"$\nu$HDM", color = "red", linewidth= 3.0)
ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.legend(loc="best", fontsize = 15)
ax1.grid()
plt.show()
