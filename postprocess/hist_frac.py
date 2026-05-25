###########################################################################################################################################################################################
import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.stats import norm
from scipy.optimize import curve_fit
###########################################################################################################################################################################################

###########################################################################################################################################################################################
#nuHDM =  np.loadtxt('/home/nikolaos-samaras/Desktop/astro/investigating/investigating_nuHDM/investigating_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.001.AHF_halos', delimiter=None) # nuHDM_b200_lmin8_lmax9_all_sf_noefe_z0.txt
#careful here!
nuHDM =  np.loadtxt('/home/nikolaos-samaras/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.000.AHF_halos', delimiter=None)

mass_tot_nuHDM = nuHDM[:,3]
mass_gas_nuHDM = nuHDM[:,44]
f_nu = 1.0 - (mass_gas_nuHDM/mass_tot_nuHDM)

mean_all= np.mean(mass_tot_nuHDM)
median_all = np.median(mass_tot_nuHDM)
print("mean all = ", np.log10(mean_all))
print("median all = ", np.log10(median_all))

data = nuHDM[:,3]
a = np.where(f_nu < 0.5) # no difference if "<="
data2 = data[a]

mean_gas_dom= np.mean(data2)
median_gas_dom = np.median(data2)
print("mean gas = ", np.log10(np.mean(data2)))
print("median gas = ", np.log10(np.median(data2)))
###########################################################################################################################################################################################

##Careful here
# I have the new file with the previous name for the shake of brevity.
###########################################################################################################################################################################################
#opt_nuHDM =  np.loadtxt('/Users/nicksam121/Desktop/astro/investigating/investigating_opt/investigating_opt_nuHDM_b200_lmin8_lmax9_sfr_noefe_output_00025.z0.000.AHF_halos', delimiter=None)
#mass_tot_opt_nuHDM = opt_nuHDM[:,3]
#mass_gas_opt_nuHDM = opt_nuHDM[:,44]
#f_nu_opt = 1.0 - (mass_gas_opt_nuHDM/mass_tot_opt_nuHDM)
###########################################################################################################################################################################################

###########################################################################################################################################################################################
# Ensure data has positive min/max values for log calculation
MIN_VAL = np.min(data)
MAX_VAL = np.max(data)
if MIN_VAL <= 0:
    # Handle cases where data might be non-positive, e.g., filter or shift
    # This example assumes you want to work with the positive subset
    print("Warning: Data contains non-positive values. Log binning requires positive data.")
    data = data[data > 0]
    MIN_VAL = np.min(data)
    MAX_VAL = np.max(data)

# 2. Define the log10-spaced bin edges
# Calculate the log10 of min and max data points
log10_min = np.log10(MIN_VAL)
log10_max = np.log10(MAX_VAL)

# Generate logarithmically spaced bin edges (e.g., 20 bins)
num_bins = 20
bins = np.logspace(log10_min, log10_max, num=num_bins) #

# 3. Bin the data and count occurrences using numpy.histogram
counts, bin_edges = np.histogram(data, bins=bins)
counts2, bin_edges2 = np.histogram(data2, bins=bins)

# 4. Display the results
print("Bin Edges (log10 spaced):")
print(bin_edges)
print("\nCounts per Bin:")
print(counts)
print("\nTotal data points counted:", np.sum(counts))
print("Total data points in original array:", len(data))
print("Finished with Total mass")
print("-------------------------------------------------------------------------------------------------", "\n")
print("\nCounts per Bin:")
print(counts)
print("\nTotal data points counted:", np.sum(counts))
print("Total data points in original array:", len(data2))
print("Finished with gas mass")
######################################################################################################################################################################################

###########################################################################################################################################################################################
plt.rcParams["figure.figsize"] = [6*4, 6*4]
plt.rcParams['font.size'] = 24
plt.rcParams['font.family'] = ['serif']
fig, ax1 = plt.subplots()

plt.hist(data, bins=bins, color='red', edgecolor='black', label = r"$\nu$HDM all bound objects")
plt.hist(data2, bins=bins, color='pink', edgecolor='black', label = r"$\nu$HDM gas-dom bound objects")

#plt.hist(data, bins=bins, color='green', edgecolor='black', label = r"opt-$\nu$HDM all bound objects")
#plt.hist(data2, bins=bins, color='olive', edgecolor='black', label = r"opt-$\nu$HDM gas-dom bound objects")
plt.axvline(mean_all, color = "black", linewidth = 2, linestyle = "-", label=r"mean Mass of all bound objects")
plt.axvline(median_all, color = "black",  linewidth = 2, linestyle = "--", label = r"median Mass of all bound objects")
plt.axvline(mean_gas_dom, color = "grey", linewidth = 2, linestyle = "-", label=r"mean Mass of gas-dom objects")
plt.axvline(median_gas_dom, color = "grey",  linewidth = 2, linestyle = "--", label = r"median Mass of gas-dom objects")

ax1.legend(loc='upper right', prop ={'size': 20})
plt.grid()
plt.gca().set_xscale("log") 
#plt.title("Histogram with Log10 Bins")
ax1.set_xlabel(r'$\log_{10} M [M_{\odot}/h]$', fontsize = 25)
ax1.set_ylabel(r'N counts', fontsize = 25)
plt.show()
###########################################################################################################################################################################################

'''
######################################################################################################################################################################################
# Ensure data has positive min/max values for log calculation

def gaussian(x, amplitude, mean, sigma):
    return amplitude * np.exp(-((x - mean)**2 / (2 * sigma**2)))


print("Started with gas mass")
MIN_VAL2 = np.min(data2)
MAX_VAL2 = np.max(data2)
if MIN_VAL2 <= 0:
    # Handle cases where data might be non-positive, e.g., filter or shift
    # This example assumes you want to work with the positive subset
    print("Warning: Gas mas Data contains non-positive values. Log binning requires positive data.")
    data2 = data2[data2 > 0]
    MIN_VAL2 = np.min(data2)
    MAX_VAL2 = np.max(data2)

# 2. Define the log10-spaced bin edges
# Calculate the log10 of min and max data points
log10_min2 = np.log10(MIN_VAL2)
log10_max2 = np.log10(MAX_VAL2)

# Generate logarithmically spaced bin edges (e.g., 20 bins)
num_bins = 20
bins = np.logspace(log10_min2, log10_max2, num=num_bins) #

# 3. Bin the data and count occurrences using numpy.histogram
counts2, bin_edges2 = np.histogram(data2, bins=bins)
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Initial guesses for the fit parameters
p0 = [max(counts2), np.mean(data2), np.std(data2)]  
# Perform the curve fit
params, cov_matrix = curve_fit(gaussian, bin_centers2, counts2, p0=p0, maxfev=5000)
amplitude, mean_fit, sigma_fit = params

# 4. Display the results
print("Bin Edges (log10 spaced):")
print(bin_edges2)
print("\nCounts per Bin:")
print(counts2)
print("\nTotal data points counted:", np.sum(counts2))
print("Total data points in original array:", len(data2))

#x_fit_log = np.linspace(np.log10(MIN_VAL2), np.log10(MAX_VAL2), 1000)
#y_fit = gaussian(x_fit_log, *params)
#ax1.plot(10**x_fit_log, y_fit, 'r-', label='Gaussian Fit (log-space mean: {:.2f})'.format(np.log10(mean_fit)))
###########################################################################################################################################################################################
'''
