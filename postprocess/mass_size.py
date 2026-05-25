######################################################################################################################################################
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
######################################################################################################################################################

######################################################################################################################################################
#nuHDM = np.loadtxt('/home/nikolaos-samaras/Desktop/astro/investigating/nuHDM_b200_lmin8_lmax9_sfr_z0_R_vir500.txt', delimiter=None)
nuHDM = np.loadtxt('/home/nsamaras/Desktop/astro/investigating/nuHDM_b200_lmin8_lmax9_sfr_z0_R_vir500.txt', delimiter=None)
nu_radius = nuHDM[:,11]
nu_mass_tot = nuHDM[:,3]
nu_mass_gas = nuHDM[:,44]
a = np.where( nu_radius>402)#>402
nu_radius = nu_radius[a]
nu_mass_tot = nu_mass_tot[a]
nu_mass_gas = nu_mass_gas[a]
######################################################################################################################################################

######################################################################################################################################################
nuHDM_f = np.loadtxt('/home/nsamaras/Desktop/astro/investigating/nuHDM_b200_lmin8_lmax9_sfr_focus_gas_MWP.z0.001.AHF_halos.txt', delimiter=None)
nu_radius_f = nuHDM_f[:,11]
nu_mass_tot_f = nuHDM_f[:,3]
nu_mass_gas_f = nuHDM_f[:,44]
######################################################################################################################################################

######################################################################################################################################################
opt_nuHDM = np.loadtxt('/home/nsamaras/Desktop/astro/investigating/opt_nuHDM_b200_lmin8_lmax9_sfr_z0_R_vir500.txt', delimiter=None)
opt_nu_radius = opt_nuHDM[:,11]
opt_nu_mass_tot = opt_nuHDM[:,3]
opt_nu_mass_gas = opt_nuHDM[:,44]
b = np.where(opt_nu_radius>402)
opt_nu_radius = opt_nu_radius[b]
opt_nu_mass_tot = opt_nu_mass_tot[b]
opt_nu_mass_gas = opt_nu_mass_gas[b]
######################################################################################################################################################

######################################################################################################################################################
opt_nuHDM_f = np.loadtxt('/home/nsamaras/Desktop/astro/investigating/opt_nuHDM_b200_lmin8_lmax9_sfr_focus_gas_MWP.z0.000.AHF_halos.txt', delimiter=None)
opt_nu_radius_f = opt_nuHDM_f[:,11]
opt_nu_mass_tot_f = opt_nuHDM_f[:,3]
opt_nu_mass_gas_f = opt_nuHDM_f[:,44]
######################################################################################################################################################

######################################################################################################################################################
clusters = fits.open("/home/nsamaras/Desktop/astro/erass1cl_primary_v3.2.fits")
# https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/erass1cl_primary_v3.2.html

clusters_data = clusters[1].data
z = clusters_data['BEST_Z']

#print("max eRossita z = ", max(z))
mass         = clusters_data['M500']
#print("max of M500 = ", max(mass),"\n")
mass_gas     = clusters_data['MGAS500']
radius       = clusters_data['R500']
k = np.where(( z < 0.05) & ( radius>402)) # add radius cut to make point at the left of the plot disappear
#k = np.where(z<0.1)
mass = mass[k]*(10**13)
mass_gas = mass_gas[k]*(10**11)
radius = radius[k]

radius_minus = clusters_data["R500_L"]
radius_minus = radius_minus[k]

radius_plus  = clusters_data["R500_H"]
radius_plus = radius_plus[k]

mass_L = clusters_data["M500_L"]
mass_L = mass_L[k]*(10**13)
print(min(mass_L), "min of mass_L")
print(max(mass_L), "max of mass_L")
mass_H = clusters_data["M500_H"]
print(min(mass_H), "min of mass_H")
print(max(mass_H), "max of mass_H")
mass_H = mass_H[k]*(10**13)

mass_gas_L = clusters_data["MGAS500_L"]
print(min(mass_gas_L), "min of mass_gas_L")
print(max(mass_gas_L), "max of mass_gas_L")
mass_gas_L = mass_gas_L[k]*(10**11)
print(min(mass_gas_L), "min of mass_gas_L")
print(max(mass_gas_L), "max of mass_gas_L")
mass_gas_H = clusters_data["MGAS500_H"]
mass_gas_H = mass_gas_H[k]*(10**11)
asymmetric_error_y_gas = [mass_gas_L, mass_gas_H]

asymmetric_error_x = [radius_minus, radius_plus]
asymmetric_error_y = [mass_L, mass_H]


#super_clusters = fits.open("/home/nsamaras/Desktop/astro/investigating/erass1sc.fits")#superclusters
# https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/LiuA_DR1/erass1sc.html

#super_clusters_data = super_clusters[1].data

#super_z = super_clusters_data["REDSHIFT"]
#l = np.where(super_z<0.05)
#super_clusters= super_clusters[1].data
#super_mass = super_clusters["MTOT"]*10**14
#super_mass = super_mass[l]
#super_mass_error = super_clusters["emtot"]*10**14
#super_mass_error = super_mass_error[l]
#super_size = super_clusters['lproj']*1000
#super_size = super_size[l]
######################################################################################################################################################

######################################################################################################################################################
plt.rcParams["figure.figsize"] = [24, 24]#24,24
plt.rcParams['font.size'] = 25
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 19
plt.rcParams['ytick.labelsize'] = 19

fig, ax1 = plt.subplots()
ax1.set_xscale("log")
ax1.set_yscale("log")

ax1.errorbar(radius, mass, xerr=asymmetric_error_x, yerr= asymmetric_error_y, fmt='o', color = "grey", ecolor="grey", alpha=0.2, label = "SRG/eROSITA - M500")
ax1.errorbar(radius, mass_gas, xerr=asymmetric_error_x, yerr = asymmetric_error_y_gas, color= "black", fmt = '^', ecolor = "black", alpha=0.4, label = "SRG/eROSITA - MGAS500")
#ax1.scatter(super_size, super_mass, color="orange", alpha=1, label = "SRG/eROSITA - superclusters")

#ax1.errorbar(radius, mass, xerr=asymmetric_error, fmt='.', color = 'grey')
#ax1.errorbar(radius, mass, xerr=asymmetric_xerr, fmt="o", color = "grey", ecolor="grey",linestyle=None, alpha=0.05, label = "SRG/eROSITA - eRASS1 (X-ray)")
#ax1.scatter(radius, mass,  color = "grey", alpha=0.5, label = "SRG/eROSITA - eRASS1 (X-ray)")

ax1.scatter(nu_radius*0.674, nu_mass_tot*0.674, linestyle = "solid", alpha=0.8, s=120, color = "red",  label = r"$\nu$HDM - total (gas + $\nu_s$) mass" )
ax1.scatter(nu_radius*0.674, nu_mass_gas*0.674, linestyle = "solid", alpha=1, s=30, color = "red",  label = r"$\nu$HDM - only gas mass" )

ax1.scatter(nu_radius_f*0.674, nu_mass_tot_f*0.674,  linestyle = "solid", alpha=0.8, marker="^", s=300, color = "pink",  label = r"$\nu$HDM - focused on gas ( - $\nu_s$)")
ax1.scatter(nu_radius_f*0.674, nu_mass_gas_f*0.674,  linestyle = "solid", alpha=0.9, marker="^", color = "magenta",  label = r"$\nu$HDM - focused on gas mass (only)" )

ax1.scatter(opt_nu_radius*0.556, opt_nu_mass_tot*0.556,  linestyle = "solid", s=100, alpha=0.8, color = "green",  label = r"opt-$\nu$HDM - total (gas + $\nu_s$) mass" )
ax1.scatter(opt_nu_radius*0.556, opt_nu_mass_gas*0.556,  linestyle = "solid", s=30, alpha=1,  color = "green",  label = r"opt-$\nu$HDM - only gas mass" )

ax1.scatter(opt_nu_radius_f*0.556, opt_nu_mass_tot_f*0.556,  linestyle = "solid", alpha=0.15, marker='^', s=300, color = "olive",  label = r"opt-$\nu$HDM - focused on gas ( - $\nu_s$)")
ax1.scatter(opt_nu_radius_f*0.556, opt_nu_mass_gas_f*0.556,  linestyle = "solid", alpha=0.9, marker='^', color = "olive",  label = r"opt-$\nu$HDM - focused on gas mass (only)" )

ax1.set_xlim(400, 1400) # or max = 2.0*10**3)
ax1.set_ylim(10**13, 10**16.1)

ax1.grid()
ax1.legend(loc='best', prop ={'size': 10})
ax1.set_xlabel(r'$R_{500}$ [kpc]', rotation = 0, fontsize =30)
ax1.set_ylabel(r'$M_{500} [M_{\odot}]$', rotation = 90, fontsize = 30)
#ax1.yaxis.set_label_coords(-.1, -.9)
plt.title("Mass-size relation (z < 0.05)", fontsize = 35)
#plt.savefig("mass_size.pdf")
plt.show()
######################################################################################################################################################