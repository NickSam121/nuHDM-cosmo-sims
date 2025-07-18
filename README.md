# nuHDM-cosmo-sims
This describes how to conduct cosmological MOND simulations with Phantom of RAMSES (PoR - which is a MONDified version of RAMSES), on the nuHDM cosmological framework. Samaras, Grandis and Kroupa 2025 (10.1093/mnras/staf1041) have optimized the traditional nuHDM cosmological model to a newer version of it, the opt-nuHDM model. They have modernized it with respect to ESA's Planck 2018 data. Both nuHDM variant models are flat (omK=0), following a FLRW metric in the background. The nuHDM has the same amount of omch2 in omnuh2, but the opt-nuHDM has very different cosmological parameters (opt_nuHDM.ini). There are 3 steps to perform these simulations: first, the generation of Initial Conditions, second the actual hydrodynamical simulations and third the identification of bound structures.

To generate the Initial Conditions with CAMB and MUSIC:
1. a) These simulations are designed to commence at z = 199.0, so one needs to specify that in CAMB (transfer_redshift(1) = 199.0). Once CAMB is installed, one needs to incorporate the extra, massive but light (11eV), sterile neutrino on the namelist:
massless_neutrinos = 2.0293
nu_mass_eigenstates = 2
massive_neutrinos = 1 1
share_delta_neff = F
nu_mass_degeneracies = 1.0147 1
nu_mass_fractions = 0.0044 0.9956

I follow the Wittenburg et al 2023(10.1093/mnras/stad1371) protocol. The sterile neutrino will have to two mass eigenstates with different mass degeneracies. The opt-nuHDM sterile neutrino will have a similar mass (\approx 13eV), so identical parameters. If the accurracy boost, if set = 5, will naturally slow down the CAMB simulation.  The first thing to check is the CMB fit (want_CMB = T). What is however needed is the transfer function at z=199 (get_transfer = T, transfer_high_precision = T accurate_massive_neutrino_transfers = T). The transfer function will be the input for MUSIC, which will sample the spectrum and it will go from the Fourier space to the physical space, providing a binary file as an input for PoR. The opt-nuHDM will NOT fit the Planck CMB. While the opt-nuHDM was a nearly perfect fit to the CMB, with CosmoSIS mcmc sampler, if one writes the output posterior parameters into CAMB, one simple does not find any correspondance and the fit is bad. I solve the problem, by scaling the resulting z = 199 transfer function, so that the CAMB transfer function is the same as the CosmoSIS one. This is not the best way to do it, but one needs to understand the different physics that the two codes are making use of. In the end, one needs not the CMB but the Transfer function for MUSIC. Obviously, one obtains one transfer function for nuHDM and one for the opt-nuHDM, running the corresponding namelists. Once they are generated, each of them individually go to MUSIC.

1. b) Wittenburg et al 2023 have modified MUSIC in the source code, in order to include the negative values of the transfer function. Negative T(f) values will mean overdense regions become underdense, and since we work in log space, these will be undesirably neglected. The MUSIC file responsible for that is the src/plugins/transfer_camb.cc (lines 133-159).

2. After the IC are generated, one needs to run PoR to actually perform the hydrodynamical simulations. There is a number of parameters one needs to pay attention.
3. The ngridmax and npartmax usually are set to be equal. If ngridmax >> 2^lmax^3, then the . This happens for instance in this movie (https://www.youtube.com/watch?v=6dDqgxzIuqg), when one can basically see the grid structure. This naturally increase the integration time and the simulation slows down. Therefore ngridmax=npartmax=(2^lmax)^3, if there is no SF.

8. AHF convert end wit herror
   no mpi + makefile.config+ ahf_halos.c varies with model

9. pynbody???


This will generate a box of 200 Mpc/h with Star Formation, radiative cooling, SN feedback, with levelmin = ((2)**7)**3 particles
