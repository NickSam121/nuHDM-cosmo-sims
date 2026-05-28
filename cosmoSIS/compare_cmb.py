import numpy as np
import matplotlib.pyplot as plt

#nuHDM           = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/output/angus/cmb_cl/tt.txt", delimiter="")
camb_nuHDM  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/angus_nils_scalCls.dat", delimiter="")
camb_opt_nuHDM  = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/angus_opt_scalCls.dat", delimiter="")
opt_nuHDM    = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/output/opt-angus/cmb_cl/tt.txt", delimiter="")
planck            = np.genfromtxt(fname="/local/home/nicksam/cosmosis-standard-library/NS_data/COM_PowerSpect_CMB-TT-binned_R3.01.txt", delimiter="")
#camb_planck   = np.genfromtxt(fname="/local/home/nicksam/Desktop/cream/nuHDM/CAMB-1.3.5/fortran/planck_2018_scalCls.dat", delimiter="")

l_planck   = planck[:,0]
dl_planck  = planck[:,1]
dl_minus   = planck[:,2]
dl_plus    = planck[:,3]
asymmetric_error = [dl_minus, dl_plus]
#cambplanck1 = camb_planck[:,0] 
#cambplanck2 = camb_planck[:,1] 

cambnuHDM1 = camb_nuHDM[:,0] 
cambnuHDM2 = camb_nuHDM[:,1]

camb_opt_nuHDM1 = camb_opt_nuHDM[:,0] 
camb_opt_nuHDM2 = camb_opt_nuHDM[:,1]

font1 = {'family':'serif','color':'black','size':25}
plt.rcParams["figure.figsize"] = [6*4., 4*4.]

#plt.plot(cambplanck1, cambplanck2,  label = r"planck camb", color = "gold", linestyle = "--", linewidth= 1.5)
plt.plot(cambnuHDM1, cambnuHDM2,  label = r"CAMB $\nu$HDM Wittenburg+23", color = "red", linestyle = "--", linewidth= 1.5)
plt.plot(camb_opt_nuHDM1, camb_opt_nuHDM2, label = r"CAMB opt-$\nu$HDM Samaras et al 2025", color = "teal", linestyle = "--", linewidth= 1.5)
plt.plot(opt_nuHDM, label = r"CosmoSIS opt-$\nu$HDM this work", color = "green", linestyle = "-", linewidth= 3.0)

plt.errorbar(l_planck, dl_planck, yerr = asymmetric_error, label = "PLANCK 2018 data", color = "blue", linewidth= 3.0,  fmt="o")
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlim([0,2500])
plt.xlabel(r"Multiple Moment l", fontdict=font1)
plt.ylabel(r"$\ell(\ell+1) C^{TT}_{\ell}/2\pi$ [$\mu K^2$]", fontdict=font1)
plt.legend(loc="best", fontsize = 21)
plt.grid()
plt.title("CMB Power Spectrum", fontsize = 35, fontdict=font1)
plt.show()
