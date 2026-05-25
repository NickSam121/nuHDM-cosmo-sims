########################################################################################
import numpy as np
from matplotlib import pyplot as plt
########################################################################################

#I have commented out the line 87 of AHF/src/defin.h/#define AHFptfocus  0.
#Then run this:
#./bin/ramses2gadget -i output_000XX
#sbatch slurm
#then halos will be slightly more massive 1.8x10**16 to 1.9x10**16, but the mass gas will be different in the output.

######nuHDM##############################################################################
nuHDM = np.loadtxt('/home/nsamaras/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00018.z2.974.AHF_halos', delimiter=None)
mass_tot_nuHDM = nuHDM[:,3]
npart_nuHDM   = nuHDM[:,11]
mass_gas_nuHDM = nuHDM[:,44]
fnu_nuHDM = 1-(mass_gas_nuHDM/mass_tot_nuHDM)

#print(mass_nuHDM[:10])
#print(gas_mass_nuHDM[:10],"\n")
#########################################################################################

###############opt-nuHDM###################################################################
#opt_nuHDM =  np.loadtxt('/Users/nicksam121/Desktop/astro/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00016.z3.961.AHF_halos', delimiter=None)

#gas_opt_nuHDM =  np.loadtxt('/local/home/nicksam/Desktop/NS+2025/data/opt_nuHDM_b200_lmin8_lmax9_output_00025.z0.000.AHF_gas_halos.txt', delimiter=None)
#mass_tot_opt_nuHDM = opt_nuHDM[:,3]
#mass_gas_opt_nuHDM = opt_nuHDM[:,44]

#print(min(mass_gas_opt_nuHDM[:5]))
#print(len(mass_gas_opt_nuHDM))
######################################################################################

'''
######################################################################################
dex     = np.linspace(11.6, 16.8, num=27, endpoint=True) #2dex
dex2     = np.linspace(13.6, 15.4, num=10, endpoint=True) #2dex

n_tot_neu_opt_nuHDM = np.zeros(len(dex))
err_opt_neu_opt_nuHDM  = np.zeros(len(dex))
mass_binned_neu_opt_nuHDM = np.zeros(len(dex))
mass_binned_gas_opt_nuHDM = np.zeros(len(dex))

n_tot_gas_opt_nuHDM = np.zeros(len(dex2))
err_opt_gas_opt_nuHDM  = np.zeros(len(dex2))

for (i,j) in enumerate(dex): 
    k = np.where( (mass_nu_opt_nuHDM>= (10**j))  & (mass_nu_opt_nuHDM <= 10**(j+0.2)))
    n = len(mass_neu_opt_nuHDM[k])
    print("Between ", j ," and ", j+0.2, ", there are", n, "Nb of galaxies = len(mass_opt_nuHDM[k]), ",len(mass_neu_opt_nuHDM[k]))
    n_tot_neu_opt_nuHDM[i] = n
    tot_neu_mass = np.log10(sum(mass_neu_opt_nuHDM[k]))
    print("And their total mass is = ", np.log10(sum(mass_neu_opt_nuHDM[k])))
    mass_binned_neu_opt_nuHDM[i] = tot_neu_mass
print(n_tot_neu_opt_nuHDM)
print("Sum of n_tot = ", sum(n_tot_neu_opt_nuHDM),"\n")
err_opt_neu_opt_nuHDM = 1/np.sqrt(n_tot_neu_opt_nuHDM)

for (i,j) in enumerate(dex2): 
    k = np.where( (mass_gas_opt_nuHDM>= (10**j))  & (mass_gas_opt_nuHDM <= 10**(j+0.2))) 
    n = len(mass_gas_opt_nuHDM[k])
    print("Between ", j ," and ", j+0.2, ", there are", n, "Nb of galaxies = len(mass_gas_nuHDM[k]), ",len(mass_gas_opt_nuHDM[k]))
    n_tot_gas_opt_nuHDM[i] = n
    tot_gas_mass = np.log10(sum(mass_gas_opt_nuHDM))
    mass_binned_gas_opt_nuHDM[i] = tot_gas_mass
print(n_tot_gas_opt_nuHDM)
print("Sum of n_tot = ", sum(n_tot_gas_opt_nuHDM),"\n")
err_opt_gas_opt_nuHDM = 1/np.sqrt(n_tot_gas_opt_nuHDM)
######################################################################################
'''
######################################################################################

######################################################################################
plt.rcParams["figure.figsize"] = [24, 24]
plt.rcParams['font.size'] = 24
plt.rcParams['font.family'] = 'serif'
fig, ax1 = plt.subplots()

#ax1.scatter(mass_tot_opt_nuHDM, 1-(mass_gas_opt_nuHDM/mass_tot_opt_nuHDM), linestyle = "solid", linewidth=4.5, color = "lightgreen",  label = r"opt-$\nu$HDM bound structures" )
#plt.hlines(0.8391,xmin=min(mass_tot_opt_nuHDM),xmax=max(mass_tot_opt_nuHDM), linestyle = "solid", linewidth=4, color = "blue",  label = r"$\Lambda$CDM - Tristram+ 2024")
#plt.hlines(0.8434,xmin=min(mass_tot_opt_nuHDM),xmax=max(mass_tot_opt_nuHDM), linestyle = "solid", linewidth=4, color = "red",  label = r"$\nu$HDM - Wittenburg+ 2023")
#plt.hlines(0.853,xmin=min(mass_tot_opt_nuHDM),xmax=max(mass_tot_opt_nuHDM), linestyle = "solid", linewidth=4, color = "green",  label = r"opt-$\nu$HDM - Samaras+ 2025")

scatter = ax1.scatter(mass_tot_nuHDM, 1-(mass_gas_nuHDM/mass_tot_nuHDM), c=npart_nuHDM, cmap='viridis', linestyle = "solid", alpha=1,  label = r"bound objects" )
cbar = plt.colorbar(scatter)
cbar.set_label(r'Radius [kpc/$h$]')
plt.hlines(0.86,xmin=min(mass_tot_nuHDM),xmax=max(mass_tot_nuHDM), linestyle = "solid", linewidth=4, color = "blue",  label = r"$\Lambda$CDM - Tristram+ 2024")
#plt.hlines(0.8434,xmin=min(mass_tot_nuHDM),xmax=max(mass_tot_nuHDM), linestyle = "solid", linewidth=4, color = "red",  label = r"$\nu$HDM - Wittenburg+ 2023")
plt.hlines(0.85,xmin=min(mass_tot_nuHDM),xmax=max(mass_tot_nuHDM), linestyle = "solid", linewidth=4, color = "green",  label = r"opt-$\nu$HDM - Samaras+ 2025")
plt.hlines(0.5,xmin=min(mass_tot_nuHDM),xmax=max(mass_tot_nuHDM), linestyle = "dashed", linewidth=4, color = "grey")
plt.text(10**13.1, 0.53, r'$f_{\nu}$ = $f_{gas}$', color = "grey",fontsize=30)


ax1.grid()
ax1.legend(loc='lower right', prop ={'size': 30})
ax1.set_ylabel(r'$f_{\nu}$', rotation = 0, fontsize = 35)
ax1.yaxis.set_label_coords(-.1, .5)

ax1.set_xlabel(r'M $[M_{\odot}/h]$', fontsize = 25)
ax1.set_xscale("log")
plt.title("Sterile neutrino mass fraction at z = 2.97", fontsize = 30)
plt.show()
