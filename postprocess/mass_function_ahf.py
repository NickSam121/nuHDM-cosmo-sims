##############################################################################################################################################################################
import numpy as np
from matplotlib import pyplot as plt
import hmf
from hmf import MassFunction     # The main hmf class
from hmf.halos import mass_definitions as md
print("Using hmf version v%s"%hmf.__version__)
from hmf import cosmo
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from hmf import Transfer
##############################################################################################################################################################################

###### Bahcall & Cen 1993 #################################################################################################################################################
#obs =  np.genfromtxt('/Users/nicksam121/Desktop/astro/investigating/bahcall_cen1993.txt', delimiter=",")
#mass_x = obs[:,0]*(10**15)
#phi_y  = obs[:,1]
#print(" Masses of Bahcall & Cen 1993 = ",  mass_x, "\n")
#print("Phi of Bahcall & Cen 1993 = ", phi_y,"\n")

obs2 = np.genfromtxt('/home/nsamaras/Desktop/astro/investigating/Driver2022.txt', delimiter=",")
mass_x2 = 0.67*10**obs2[:,0]
phi_y2  = 10**obs2[:,3]
error   = 10**(obs2[:,7]+obs2[:,2])
print(" Masses of Driver2022 = ",  mass_x2, "\n")
print("Phi of Driver2022 = ", phi_y2,"\n")
##############################################################################################################################################################################
'''
##############################################################################################################################################################################
clusters = fits.open("/Volumes/T7/NS+2025/data/erass1cl_primary_v3.2.fits")
# https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/erass1cl_primary_v3.2.html
clusters_data = clusters[1].data
eros_mass_x         = clusters_data['M500']
z = clusters_data['BEST_Z']
print(len(eros_mass_x), "= ALL eRosita clusters" )

k = np.where( z < 0.11 )
eros_mass_x = eros_mass_x[k]*(10**13)

print(len(eros_mass_x), "= total eRosita clusters with z<0.11")
print(min(eros_mass_x), " = min mass of eRosita")
print(max(eros_mass_x), " = max mass of eRosita", "\n" )
##############################################################################################################################################################################
'''
######nuHDM_lmin8####################################################################################################################################################################
nuHDM =  np.loadtxt('/home/nsamaras/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.001.AHF_halos', delimiter=None) # nuHDM_b200_lmin8_lmax9_all_sf_noefe_z0.txt
mass_nuHDM = nuHDM[:,3]
print(min(np.log10(mass_nuHDM))," = min of mass nuHDM structure")
print(max(np.log10(mass_nuHDM))," = max of mass nuHDM structure")
print(len(mass_nuHDM), " = tot number of nuHDM structures")
print("\n")
###############################################################################################################################################################################

###############opt-nuHDM#########################################################################################################################################################
opt_nuHDM =  np.loadtxt('/home/nsamaras/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.000.AHF_halos', delimiter=None)
mass_opt_nuHDM = opt_nuHDM[:,3]
print(min(np.log10(mass_opt_nuHDM)), " = min of mass of opt nuHDM structure")
print(max(np.log10(mass_opt_nuHDM)), " = max of mass of opt nuHDM structure")
print(len(mass_opt_nuHDM), " = tot number of opt nuHDM structures")
print("\n")
############################################################################################################################################################################

###############Bohemian#########################################################################################################################################################
#bohemaian_b100_lmin7_lmax12 =  np.loadtxt('/local/home/nicksam/Desktop/NS+2025/data/bohemian_b100_lmin7_lmax12_gas_output_00030.z0.000.AHF_halos', delimiter=None)
#mass_bohem_lmin7 = bohemaian_b100_lmin7_lmax12[:,3]
#print(min(np.log10(mass_bohem_lmin7)), " = min of mass bohemian lmin7 structure")
#print(max(np.log10(mass_bohem_lmin7)), " = max of mass Bohemianlmin7 structure")
#print(len(mass_bohem_lmin7)," = total number of bohemian lmin7 gas structures")
#print("\n")
############################################################################################################################################################################

############################################################################################################################################################################
dex     = np.linspace(min(np.log10(mass_nuHDM)), max(np.log10(mass_nuHDM)), num=10, endpoint=True) #2dex
print(dex[1]-dex[0], "= dex[1]-dex[0]","\n")
dex2    = np.linspace(min(np.log10(mass_opt_nuHDM)), max(np.log10(mass_opt_nuHDM)), num=10, endpoint=True) #2dex
print(dex2[1]-dex2[0], "= dex2[0]-dex2[0]","\n")

#dex3    = np.linspace(min(np.log10(eros_mass_x)), max(np.log10(eros_mass_x)), num=10, endpoint=True) #2dex

print(dex,"\n")
print(dex2,"\n")
#print(dex3, "\n")

n_tot_nuHDM  = np.zeros(len(dex))
#err_nuHDM     = np.zeros(len(dex))

#n_tot_eros  = np.zeros(len(dex3))

n_tot_opt_nuHDM = np.zeros(len(dex2))
#err_opt_nuHDM  = np.zeros(len(dex2))
############################################################################################################################################################################

############################################################################################################################################################################
for (i,j) in enumerate(dex):
    k = np.where( (mass_nuHDM>= (10**j)) & (mass_nuHDM <= 10**(j+abs(dex[1]-dex[0]))))
    n = len(mass_nuHDM[k])
    print("Between ", j ," and ", j+(dex[1]-dex[0]), ", there are", n, "Nb of nuHDM galaxies = ), ", len(mass_nuHDM[k]))
    n_tot_nuHDM[i] = n
    #print ("i = ", i , " and j = ", j, "\n")
#print(n_tot_nuHDM)
print("Sum of n_tot = ", sum(n_tot_nuHDM),"\n")
#err_nuHDM = 1/np.sqrt(n_tot_nuHDM)

for (i,j) in enumerate(dex):
    k = np.where( (mass_opt_nuHDM >= (10**j)) & (mass_opt_nuHDM <= 10**(j+(dex[1]-dex[0]))))
    n = len(mass_opt_nuHDM[k])
    print("Between ", j ," and ", j+(dex[1]-dex[0]), ", there are", n, "Nb of opt-nuHDM galaxies = ), ", len(mass_opt_nuHDM[k]))
    n_tot_opt_nuHDM[i] = n
#print(n_tot_opt_nuHDM_lmin8)
print("Sum of n_tot = ", sum(n_tot_opt_nuHDM),"\n")
#err_opt_nuHDM = 1/np.sqrt(n_tot_opt_nuHDM)
#print(err_opt_nuHDM)
'''
for (i,j) in enumerate(dex):
    k = np.where( (eros_mass_x >= (10**j)) & (eros_mass_x <= 10**(j+(dex[1]-dex[0]))))
    n = len(eros_mass_x[k])
    print("Between ", j ," and ", j+(dex[1]-dex[0]), ", there are", n, "Nb of eRosita galaxies = ), ", len(eros_mass_x[k]))
    n_tot_eros[i] = n
print("Sum of n_tot = ", sum(n_tot_eros),"\n")
#err_eros = 1/np.sqrt(n_tot_eros)
'''

############################################################################################################################################################################

###############PLOTS########################################################################################################################################################
plt.rcParams["figure.figsize"] = [13, 13]
plt.rcParams['font.size'] = 18
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
fig, ax1 = plt.subplots()
dex_x = [10**11.5, 10**12.0, 10**13.0, 10**14.0, 10**15.0, 10**16.0, 10**16.5, 10**17]
dex_y = [10**-8,    10**-7.0, 10**-6.0, 10**-5.0, 10**-4.0, 10**-3,   10**-2,  10**-1]

#ax1.axvline(2*10**15, color = 'black', linestyle="--", label = 'El Gordo')
#ax1.plot(10**dex, n_tot/200**3, linestyle = "solid", linewidth=4.5, color = "green", marker="o", markersize = 13, label = r"opt-$\nu$HDM_lmin7_lmax9" )
#ax1.fill_between(10**dex, (n_tot-err_opt)/(200**3), (n_tot+err_opt)/(200**3), color = "green", alpha= 0.08)

mf1 = MassFunction(Mmin=12.5, Mmax=17, z = 0.0, n = 0.966, transfer_model="EH", cosmo_params={"Om0":0.3, "H0":67.4, "Ob0":0.0486} , mdef_params = {"overdensity": 180}) #LCDM

mf2 = MassFunction(Mmin=12.5, Mmax=17, z = 0.0, n = 0.966, transfer_model="EH", cosmo_params={"Om0":0.3, "H0":67.4, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 11.2]}, mdef_params = {"overdensity": 180}) #nuHDM no TF
mf3 = MassFunction(Mmin=12.5, Mmax=17, z = 0.0, n = 0.872, transfer_model="EH",cosmo_params={"Om0":0.495, "H0":55.64, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 12]}, mdef_params = {"overdensity": 180}) #opt-nuHDM no TF

Transfer.transfer_model ="/home/nsamaras/Desktop/astro/investigating/nuHDM_CAMB_transfer_z199.dat"
mf4 = MassFunction(Mmin=12.5, Mmax=17, z = 0.0, n = 0.966, cosmo_params={"Om0":0.3, "H0":67.4, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 11.2]}, mdef_params = {"overdensity": 180}) #nuHDM

Transfer.transfer_model = "/home/nsamaras/Desktop/astro/investigating/opt_nuHDM_CAMB_transfer_z199.dat"
mf5 = MassFunction(Mmin=12.5, Mmax=17, z = 0.0, n = 0.872, cosmo_params={"Om0":0.495, "H0":55.64, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 12]}, mdef_params = {"overdensity": 180}) #opt-nuHDM

#https://hmf.readthedocs.io/en/stable/_autosummary/hmf.density_field.transfer_models.html#module-hmf.density_field.transfer_models
#mf3 = MassFunction(Mmin=11.5, Mmax=17, z = 0.96, n = 0.872, cosmo_params={"Om0":0.495, "H0":55.64, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 11.2]}, mdef_params = {"overdensity": 500})

#ax1.scatter(10**dex3, n_tot_eros/(471.35**3), linestyle = "solid", linewidth=4.1, color = "grey",  marker="^", label = r"eRosita" )

#mf3 = MassFunction(z = 0, cosmo_params={"Om0":0.3, "H0":67.4, "Ob0":0.0486}, hmf_model="PS")
#mf4 = MassFunction(z = 0, cosmo_params={"Om0":0.495, "H0":55.64, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 11.2]}, hmf_model="PS")
#mf5 = MassFunction(z = 0, cosmo_params={"Om0":0.3, "H0":67.4, "Ob0":0.0486, "Neff": 4.044, "m_nu": [0, 0, 0, 11.2]}, hmf_model="PS")

#ax1.scatter(mass_x, phi_y,  color = "black",  linewidth=10, marker="^", label = r"Bahcall & Cen 1993" )

### COMOVING quantities!!!! ##########################################################################################
ax1.plot(mass_x2, phi_y2 , linewidth=5, color = "black",  markersize = 10, marker="o", label = r"Driver+ 2022")
ax1.fill_between(mass_x2, phi_y2-error, phi_y2+error, color = "grey", alpha= 0.8)

#plt.plot(mf1.m, mf1.dndlog10m/(0.67**3), color="blue",  linestyle = "dotted", linewidth= 5.0, label= r"$\Lambda$CDM (hmf)")
plt.plot(mf1.m, mf1.dndlog10m, color="blue",  linestyle = "solid", linewidth= 5.0, label= r"$\Lambda$CDM (hmf)")
#plt.plot(mf11.m, mf11.dndlog10m, color="blue",  linestyle = "dotted", linewidth= 15.0, label= r"$\Lambda$CDM (hmf)")

plt.plot(mf2.m, mf2.dndlog10m, color="red",  linestyle = "dashed", linewidth= 5.0, label= r"$\nu$HDM (hmf - EH)")#/(0.67**3)
plt.plot(mf3.m, mf3.dndlog10m, color="green",  linestyle = "dashed", linewidth= 5.0, label= r"opt-$\nu$HDM (hmf - EH)")#/(0.556**3)

plt.plot(mf4.m, mf4.dndlog10m, color="red", linestyle = "dotted", linewidth= 5.0, label=r"$\nu$HDM (hmf - TF)")
plt.plot(mf5.m, mf5.dndlog10m, color="green", linestyle = "dotted", linewidth= 5.0, label=r"opt-$\nu$HDM (hmf - TF)")#/(0.556**3)

ax1.plot(10**dex, n_tot_nuHDM/200**3, linestyle = "solid", linewidth=4.5, color = "red", label = r"$\nu$HDM (PoR)" )
ax1.plot(10**dex2, n_tot_opt_nuHDM/200**3, linestyle = "solid", linewidth=4.5, color = "green",  label = r"opt-$\nu$HDM (PoR)" )

#ax1.plot(10**dex3, (dex3[1]-dex3[0])*n_tot_bohem_lmin7/100**3, linestyle = "solid", linewidth=4.5, color = "grey", marker="o", markersize = 13, label = r"Bohemian b100_lmin7_lmax12" )
#ax1.fill_between(10**dex, (n_tot_nuHDM-err_nuHDM)/(200**3), (n_tot_nuHDM+err_nuHDM)/(200**3), color = "red", alpha= 0.5, linewidth=4.5)
#ax1.fill_between(10**dex2, (n_tot_opt_nuHDM-err_opt_nuHDM)/(200**3), (n_tot_opt_nuHDM+err_opt_nuHDM)/(200**3), color = "green", alpha= 0.5, linewidth=4.5)

ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xticks(dex_x)
ax1.set_xlim(10**12.8, 10**16.5)
ax1.set_ylim(10**-8, 10**-1.5)
ax1.set_yticks(dex_y)
ax1.grid()
ax1.legend(loc='upper right', prop ={'size': 20})
ax1.set_xlabel(r'M $[M_{\odot}/h]$', fontsize = 25)
ax1.set_ylabel(r'number density (cMpc$^{-3}$ dex$^{-1}$)', fontsize = 25)
ax1.set_title(r"Dynamical MOND mass (gas+$\nu_s$) function at z = 0.0", fontsize=30)
#ax1.set_title(r"Mass function of bound structures (gas+$\nu_s$, with SF, no EFE) z=0.11", fontsize=30)
plt.show()