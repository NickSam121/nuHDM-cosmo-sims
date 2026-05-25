########################################################################################
import numpy as np
from matplotlib import pyplot as plt
########################################################################################

#################################################################################################################################################################################################
nuHDM_z0 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.001.AHF_halos', delimiter=None)
nuHDM_z011 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00024.z0.111.AHF_halos', delimiter=None)
nuHDM_z0513 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00023.z0.513.AHF_halos', delimiter=None)
nuHDM_z0966 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00022.z0.996.AHF_halos', delimiter=None)
nuHDM_z1492 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00021.z1.492.AHF_halos', delimiter=None)
nuHDM_z2020 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00020.z2.020.AHF_halos', delimiter=None)
nuHDM_z2547 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00019.z2.547.AHF_halos', delimiter=None)
nuHDM_z2961 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00018.z2.961.AHF_halos', delimiter=None)
nuHDM_z3544 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00017.z3.544.AHF_halos', delimiter=None)
#nuHDM_z3930 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00016.z3.930.AHF_halos', delimiter=None)
#nuHDM_z4514 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00015.z4.514.AHF_halos', delimiter=None)
#nuHDM_z5001 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00014.z5.001.AHF_halos', delimiter=None)

n_struct_z0 = len(nuHDM_z0[:,1])
n_struct_z1 = len(nuHDM_z011[:,1])
n_struct_z2 = len(nuHDM_z0513[:,1])
n_struct_z3 = len(nuHDM_z0966[:,1])
n_struct_z4 = len(nuHDM_z1492[:,1])
n_struct_z5 = len(nuHDM_z2020[:,1])
n_struct_z6 = len(nuHDM_z2547[:,1])
n_struct_z7 = len(nuHDM_z2961[:,1])
n_struct_z8 = len(nuHDM_z3544[:,1])
#n_struct_z9 = len(nuHDM_z3930[:,1])
#n_struct_z10 = len(nuHDM_z4514[:,1])
#n_struct_z11 = len(nuHDM_z5001[:,1])

z = [0.0, 0.11, 0.513, 0.966, 1.492, 2.02, 2.547, 2.961, 3.544, 3.930, 4.514, 5.001]
lens_z = [n_struct_z0, n_struct_z1, n_struct_z2, n_struct_z3, n_struct_z4, n_struct_z5, n_struct_z6, n_struct_z7, 1.0, 0.0, 0.0, 0.0] #n_struct_z8, n_struct_z9. n_struct_z10, n_struct_z11]
#################################################################################################################################################################################################

#################################################################################################################################################################################################
opt_nuHDM_z0 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.000.AHF_halos', delimiter=None)
opt_nuHDM_z011 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00024.z0.110.AHF_halos', delimiter=None)
opt_nuHDM_z0515 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00023.z0.515.AHF_halos', delimiter=None)
opt_nuHDM_z0966 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00022.z0.996.AHF_halos', delimiter=None)
opt_nuHDM_z1490 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00021.z1.490.AHF_halos', delimiter=None)
opt_nuHDM_z2022 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00020.z2.022.AHF_halos', delimiter=None)
opt_nuHDM_z2544 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00019.z2.544.AHF_halos', delimiter=None)
opt_nuHDM_z2974 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00018.z2.974.AHF_halos', delimiter=None)
opt_nuHDM_z3483 = np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00017.z3.483.AHF_halos', delimiter=None)

opt_n_struct_z0 = len(opt_nuHDM_z0[:,1])
opt_n_struct_z1 = len(opt_nuHDM_z011[:,1])
opt_n_struct_z2 = len(opt_nuHDM_z0515[:,1])
opt_n_struct_z3 = len(opt_nuHDM_z0966[:,1])
opt_n_struct_z4 = len(opt_nuHDM_z1490[:,1])
opt_n_struct_z5 = len(opt_nuHDM_z2022[:,1])
opt_n_struct_z6 = len(opt_nuHDM_z2544[:,1])
opt_n_struct_z7 = len(opt_nuHDM_z2974[:,1])
opt_n_struct_z8 = len(opt_nuHDM_z3483[:,1])
#n_struct_z9 = len(nuHDM_z3930[:,1])
#n_struct_z10 = len(nuHDM_z4514[:,1])
#n_struct_z11 = len(nuHDM_z5001[:,1])

opt_z = [0.0, 0.11, 0.513, 0.966, 1.490, 2.022, 2.544, 2.974, 3.483, 3.930, 4.514, 5.001]
opt_lens_z = [opt_n_struct_z0, opt_n_struct_z1, opt_n_struct_z2, opt_n_struct_z3, opt_n_struct_z4, opt_n_struct_z5, opt_n_struct_z6, opt_n_struct_z7, 1.0, 0.0, 0.0, 0.0] #n_struct_z8, n_struct_z9. n_struct_z10, n_struct_z11]
#################################################################################################################################################################################################


#################################################################################################################################################################################################
plt.rcParams["figure.figsize"] = [6*3, 7.5]
plt.rcParams['font.family'] = ['serif']
plt.rcParams['xtick.labelsize'] = 25
plt.rcParams['ytick.labelsize'] = 25
fig, ax1 = plt.subplots()

plt.plot(z,lens_z, linestyle = "solid", linewidth=4, color = "red",  marker = "o" , markersize = 10, label = r"$\nu$HDM model (h = 0.67, $\Omega_m$ = 0.3)")
plt.plot(opt_z, opt_lens_z, linestyle = "solid", linewidth=4, color = "green",  marker = "o", markersize = 10, label = r"opt-$\nu$HDM model (h = 0.55, $\Omega_m$ = 0.49)")

ax1.grid()
ax1.legend(loc='upper right', prop ={'size': 25})
ax1.set_ylabel('Nb of structures', fontsize =30)
ax1.set_xlabel(r'redshift z', fontsize =30)
plt.show()
########################################################################################
