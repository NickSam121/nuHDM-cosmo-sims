import numpy as np
import matplotlib.pyplot as plt
import pynbody

a = pynbody.load('/Volumes/T7/NS+2025/data/output_00025/')
a.physical_units()

font1 = {'family':'serif','size':25}
plt.rcParams["figure.figsize"] = [6*3, 4*2]

#ATTENTION configure the subplot with bottom=0.079 , top = 0.952, right = 0.91, left=0, NS Aug 2025
cax = pynbody.plot.image(a.gas, width="200 Mpc", units="Msol kpc^-2", cmap="bone", show_cbar=False)
#cax = pynbody.plot.image(a.dm, width="200 Mpc", units="Msol kpc^-2", cmap="twilight", show_cbar=False)

cbar = plt.colorbar(cax)

plt.xticks([-100000,-50000,0,50000,100000], ['100','-50','0','50','100'], fontsize=30)
plt.yticks([-100000,-50000,0,50000,100000], ['100','-50','0','50','100'], fontsize=30)
plt.xlabel("x/Mpc", fontsize=30)
plt.ylabel("y/Mpc", fontsize=30)
cbar.set_label('Density', fontsize=30)
cbar.ax.tick_params(labelsize=30)
plt.title(r"$\nu$HDM model (gas) - z = 0", fontsize=30)
#plt.title(r"$\nu$HDM model ($\nu_S$) - z = 0", fontsize=30)
plt.tight_layout()
plt.show()
