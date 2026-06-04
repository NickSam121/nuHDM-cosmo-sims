#########################################################################################################################################################
from matplotlib import pyplot as plt
import numpy as np
#########################################################################################################################################################

################    Observations           #############################################################################################################
#obs = np.loadtxt ('/local/home/nicksam/Desktop/NS+2025/data/2ds.txt', dtype='float', delimiter=None)
#obs_k = obs[:,0]
#obs_pk = obs[:,1]
#########################################################################################################################################################
nu199 = np.loadtxt('/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/angus_nils_matterpower.dat', dtype='float', delimiter=None)
nu_k = nu199[:,0]         
nu_pk = nu199[:,1]

k_opt_nuHDM  = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/output/opt-angus/matter_power_lin/k_h.txt", delimiter="")
#z_opt_nuHDM  = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/output/opt-angus/matter_power_lin/z.txt", delimiter="")
pk_opt_nuHDM  = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/output/opt-angus/matter_power_lin/p_k.txt", delimiter="")
pk_opt_nuHDM_z0 = pk_opt_nuHDM[1,:]
print(len(pk_opt_nuHDM_z0))

camb_opt_nu199 = np.loadtxt ('/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/angus_opt_matterpower_z199.dat', dtype='float', delimiter=None)
camb_opt_nu_k = camb_opt_nu199[:,0]         
camb_opt_nu_pk = camb_opt_nu199[:,1]
print(len(camb_opt_nu_k))

lcdm199 = np.loadtxt ('/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/planck2024_matterpower.dat', dtype='float', delimiter=None)
lcdm_k = lcdm199[:,0] 
lcdm_pk = lcdm199[:,1]

music_nuHDM = np.loadtxt('/local/home/nicksam/Desktop/NS+2025/data/input_powerspec_nuHDM.txt', dtype='float', delimiter=None)
music_nuHDM_k = music_nuHDM[:,0]
music_nuHDM_pk = 8.0*(np.pi**3.0)*music_nuHDM[:,5]#5 column for Ptotal

music_opt_nuHDM_b200_lmin8lmx12 = np.loadtxt('/local/home/nicksam/Desktop/NS+2025/data/input_powerspec_opt_nuHDM_b200lmn8lmx12.txt', dtype='float', delimiter=None)
music_opt_nuHDM_b200_lmin8lmx12_k = music_opt_nuHDM_b200_lmin8lmx12[:,0]
music_opt_nuHDM_b200_lmin8lmx12_pk = 8.0*(np.pi**3.0)*music_opt_nuHDM_b200_lmin8lmx12[:,5]#5 column for Ptotal

#test = np.loadtxt ('/local/home/nicksam/Desktop/NS+2025/codes/pk_test_ahf.txt', dtype='float', delimiter=None)
#test_k = lcdm199[:,0]
#test_pk = lcdm199[:,1]

#stacy = np.loadtxt ('/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/stacy_matterpower_z199.dat', dtype='float', delimiter=None)
#stacy_k = stacy[:,0]
#stacy_pk = stacy[:,1]
###################################################################################

###################################################################################
font1 = {'family':'serif','color':'black','size':30}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]
#fig, ax1 = plt.subplots(1, 1)

fig, ax1 = plt.subplots()	
dex_x = [0.0001, 0.001, 0.01, 0.1, 1, 10, 12]
dex_x2 = [0.001, 0.01, 0.1, 1]
dex_y = [0.0000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 10000]
dex_y2 = [0.1, 0.2, 0.5, 0.75, 1.2]

ax1.set_xticks(dex_x)
ax1.set_xticklabels(dex_x, fontsize=20)
ax1.set_yticks(dex_y)
ax1.set_yticklabels(dex_y, fontsize=20)
ax1.set_yscale("log")
ax1.set_xscale("log")
ax1.set_xlabel("Wavenumber k [h/cMpc]", fontdict=font1)
ax1.set_ylabel(r"P(k) $[(h^{-1}Mpc)^3]$", fontdict=font1)
ax1.set_xlim(0.0001,12)
ax1.set_xlim(0.001,12) #NS 7May2025 for zoom in to probe different neutrino free-streaming lengths
#ax1.set_ylim(bottom = 0.0000001, top= max(pk_lcdm_z0 + (1/2)*pk_lcdm_z0))

#ax1.scatter(obs_k, obs_pk, label= "2dFGRS 2005",color = "black", linewidth= 3.0)

ax1.plot(lcdm_k, lcdm_pk, label = r"$\Lambda$CDM - CAMB - $\it{Planck}$ 2024", color = "blue", linewidth= 3.0)

ax1.plot(nu_k, nu_pk, label = r"$\nu$HDM - CAMB - Wittenburg+ 2023", color = "red", linewidth= 1.0)
ax1.plot(music_nuHDM_k, music_nuHDM_pk, label = r"$\nu$HDM - MUSIC - this work", color = "red", linewidth= 3.0)

#ax1.plot(k_opt_nuHDM, pk_opt_nuHDM_z0, label = r"opt-$\nu$HDM - CosmoSIS - Samaras+ 2025", color = "green", linewidth= 1.0)
#ax1.plot(camb_opt_nu_k, camb_opt_nu_pk, label = r"opt-$\nu$HDM - CAMB - this work", color = "green", linewidth= 2.0, linestyle = "--")
#ax1.plot(music_opt_nuHDM_b200_lmin8lmx12_k, music_opt_nuHDM_b200_lmin8lmx12_pk, label = r"opt-$\nu$HDM - MUSIC - this work", color = "green", linewidth= 3.0)

#ax1.plot(stacy_k, stacy_pk, label = r"McGaugh 1999", color = "grey", linewidth= 1.0, linestyle = "-")
'''

ax1.plot(lcdm_k, (lcdm_k**3)*(lcdm_pk)/(2*(np.pi)**2), label = r"$\Lambda$CDM - CAMB - Planck 2024", color = "blue", linewidth= 3.0)
ax1.plot(nu_k, (nu_k**3)*(nu_pk)/(2*(np.pi)**2), label = r"$\nu$HDM - CAMB - Wittenburg+ 2023", color = "red", linewidth= 2.0)
#ax1.plot(k_opt_nuHDM, (k_opt_nuHDM**3)*(pk_opt_nuHDM_z0/(2*(np.pi)**2)), label = r"opt-$\nu$HDM - CosmoSIS - Samaras+ 2025", color = "green", linewidth= 3.0)
ax1.plot(camb_opt_nu_k, (camb_opt_nu_k**3)*(camb_opt_nu_pk)/(2*(np.pi)**2), label = r"opt-$\nu$HDM - CAMB - this work", color = "green", linewidth= 1.0, linestyle = "--")
ax1.hlines(1, xmin=0.001,xmax=10, color = "black")
'''
ax1.legend(loc="lower left", fontsize = 21)
ax1.grid()
ax1.set_title("Matter Power Spectrum", fontsize = 36)
plt.show()
