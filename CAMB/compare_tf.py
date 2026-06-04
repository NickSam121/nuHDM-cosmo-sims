from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize
#############################################################################################################################################################################################
#############################################################################################################################################################################################
##############################################################################
opt_x    = np.genfromtxt(fname="/media/nsamaras/T7/Archive/cosmosis-standard-library_prague_home_directory/output/opt-angus/matter_power_transfer_func/k_h.txt", delimiter="")
opt_y    = np.genfromtxt(fname="/media/nsamaras/T7/Archive/cosmosis-standard-library_prague_home_directory/output/opt-angus/matter_power_transfer_func/t_k.txt", delimiter="")

camb_opt_nu199 = np.loadtxt ('/home/nsamaras/Desktop/astro/AUUK_CAMB+CLASS/CAMB-1.3.5/fortran/angus_opt_transfer_out_z199.dat', dtype='float', delimiter=None)
camb_opt_nu_k = camb_opt_nu199[:,0]
#camb_opt_bar = camb_opt_nu199[:,1] #ATTENTION  here
camb_opt_nu_tk = 15*camb_opt_nu199[:,5] #ATTENTION  here scale massive neutrino by a factor of 15 (just multiply the entry *15)

camb_opt_nu0 = np.loadtxt ('/home/nsamaras/Desktop/astro/AUUK_CAMB+CLASS/CAMB-1.3.5/fortran/angus_opt_transfer_out_z0.dat', dtype='float', delimiter=None)
camb_opt_nu_k0 = camb_opt_nu0[:,0]
#camb_opt_bar = camb_opt_nu199[:,1] #ATTENTION  here
camb_opt_nu_tk0 = 15*camb_opt_nu0[:,5] #ATTENTION  here scale massive neutrino by a factor of 15 (just multiply the entry *15)

camb_nu199 = np.loadtxt ('/home/nsamaras/Desktop/astro/AUUK_CAMB+CLASS/CAMB-1.3.5/fortran/angus_nils_transfer_out.dat', dtype='float', delimiter=None)
camb_nu_k = camb_nu199[:,0]         
#camb_nu_bar = camb_nu199[:,2]
camb_nu_tk = camb_nu199[:,5]

lcdm199 = np.loadtxt ('/home/nsamaras/Desktop/astro/AUUK_CAMB+CLASS/CAMB-1.3.5/fortran/planck2024_transfer_out.dat', dtype='float', delimiter=None)
lcdm_k = lcdm199[:,0] 
#lcdm_bar = lcdm199[:,1]
lcdm_tk = lcdm199[:,1]

#stacy= np.loadtxt ('/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/stacy_transfer_out_z199.dat', dtype='float', delimiter=None)
#stacy_k = stacy[:,0] 
#stacy_tk = stacy[:,5]
###################################################################################

###################################################################################
def k2l(x):
    return (2 * np.pi / x)*(1/(1+199))

def l2k(x):
    return (2*np.pi/x)*(1/(1+199))
###################################################################################

###################################################################################
fig, ax = plt.subplots(figsize=(24, 24))
font1 = {'family':'serif','color':'black','size':30}
font2 = {'family':'serif','color':'black','size':20}

#fig, ax1 = plt.subplots()	
dex_x = [0.0001, 0.001, 0.01, 0.1, 1, 10, 12]
dex_x2 = [0.001, 0.01, 0.1, 1]
dex_y = [0.0000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 10000]
dex_y2 = [0.1, 0.2, 0.5, 0.75, 1.2]

#ax.plot(abs(lcdm_k),abs(lcdm_tk), label='CDM - $\Lambda$CDM - CAMB - Planck 2024', color = "blue", linestyle = "-",  linewidth= 5.0)
ax.plot(abs(lcdm_k),abs(lcdm_tk), label=r'$\rho_{CDM}$ - $\Lambda$CDM - CAMB', color = "blue", linestyle = "-",  linewidth= 1.5)
#ax.plot(abs(lcdm_k),abs(lcdm_bar), label=r'$v_{bar}$ - $\Lambda$CDM', color = "blue", alpha = 0.5, linestyle = "-",  linewidth= 3.0)

ax.plot(abs(camb_nu_k),  abs(camb_nu_tk), label = r"$\rho_{\nu}$ - $\nu$HDM - CAMB - Wittenburg+ 2023", color = "red", linestyle = "-", linewidth= 1.5)
#ax.plot(abs(camb_nu_k),  abs(camb_nu_bar), label = r"$v_{bar}$ - $\nu$HDM", color = "red", alpha=0.5, linestyle = "--", linewidth= 3.0)

ax.plot(abs(opt_x),    abs(opt_y), label = r"$\rho_{\nu}$ - opt-$\nu$HDM - CosmoSIS - Samaras+ 2025", color = "green", linestyle = "-", linewidth= 1.5)
ax.plot(abs(camb_opt_nu_k),  abs(camb_opt_nu_tk), label = r"$\rho_{\nu}$ - opt-$\nu$HDM - CAMB - this work", color = "green", linestyle = "--", linewidth= 4.0)
#ax.plot(abs(camb_opt_nu_k0),  abs(camb_opt_nu_tk0), label = r"$\rho_{\nu}$ - opt-$\nu$HDM - CAMB - this work - z = 0 ", color = "olive", linestyle = "--", linewidth= 4.0)
#ax.plot(abs(camb_opt_nu_k),  abs(camb_opt_bar), label = r"$v_{bar}$ - opt-$\nu$HDM", color = "green", alpha = 0.5,  linestyle = "--", linewidth= 3.0)

#ax.plot(abs(stacy_k),  abs(stacy_tk), label = r"McGaugh 1999", color = "grey", linestyle = "--", linewidth= 1.0)
#plt.plot(opt_x2, monoExp(opt_x2, m, t, b), label="fitted", color = "teal", linestyle = "--", linewidth= 2.0)

ax.set_xlabel(r"Wavenumber k [h/cMpc]", fontdict=font1)
ax.set_ylabel(r"|T(k)|", fontsize=30, rotation=0)

secax = ax.secondary_xaxis('top', functions=(k2l,l2k))
secax.set_xlabel(r'Wavelength [$h^{-1}$ cMpc]', fontdict=font2)
secax.set_xscale('log')

ax.axvline(x = 200/0.55, color = 'black', label = 'boxsize sidelength', linewidth = 1.5)

#ax.vlines(2*np.pi/100 * (1/200))
ax.set_xscale('log')
ax.set_yscale('log')
ax.yaxis.set_label_coords(-.1, .5)
ax.xaxis.set_label_coords(0.5, -0.085)
ax.set_xlim([0.01,10])
#ax.set_xticks([10**-2, 10**-1,10**0,10], labels = [0.01, 0.1,1,10], fontsize=20)
#secax.set_xticks([10**0, 10**1,10**0,10**2], labels = [0, 1,1,2], fontsize=15)
secax.set_xticks([10**-2, 10**-1,10**0,10], labels = [0.01, 0.1,1,10], fontsize=20)
ax.set_yticks([10**-5, 10**-3, 10**-1, 10**1, 10**3,10**5, 10**6], labels= ["-5","-3","-1","1", "3","5", "6"], fontsize=20)
ax.grid()
ax.legend(fontsize = 20,  loc='lower left')
#plt.title("Transfer function", fontsize = 36)
plt.show()

