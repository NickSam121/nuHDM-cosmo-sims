# νHDM and opt-νHDM cosmological simulations 
This describes how to conduct cosmological MOND simulations with Phantom of RAMSES (PoR - which is a MONDified version of RAMSES), on the nuHDM cosmological framework (no SF, no EFE). Samaras, Grandis and Kroupa 2025 (10.1093/mnras/staf1041) have optimized the traditional nuHDM cosmological model to a newer version of it, the opt-nuHDM model. They have modernized it with respect to ESA's Planck 2018 data. Both nuHDM variant models are flat (omK=0), following a FLRW metric in the background. The nuHDM has the same amount of omch2 in omnuh2, but the opt-nuHDM has very different cosmological parameters (opt_nuHDM.ini). There are 3 steps to perform these simulations: first, the generation of Initial Conditions, second the actual hydrodynamical simulations and third the identification of bound structures.

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

2. After the IC are generated, one needs to run PoR to actually perform the hydrodynamical simulations. There is a number of parameters one needs to pay attention. The ngridmax and npartmax usually are set to be equal. If ngridmax >> (2^lmin)^3, then the grid is oversized. This happens for instance in this movie (https://www.youtube.com/watch?v=6dDqgxzIuqg), when one can basically see the grid structure. This naturally increases the integration time and the simulation slows down. Therefore, I have noticed thatit should be: 2^lmin < ngridmax=npartmax, but not (2^lmin)^3 << ngridmax! (if there is neither SF nor EFE). More particularly, if ngridmax is manually set to = 2,100,000, the:
ngridmax - (2^lmin)^3 = 2,100,000 - 2,097,152 > 0 (but not >> 0). Last, also take into account this note by RAMSES: https://ramses-organisation.readthedocs.io/en/latest/wiki/Amr.html

The best spatial resolution will be boxsize/(2**lmax) [boxsize units]

Furthermore, for NO star-formation, the parameters are:
&PHYSICS_PARAMS
isothermal=.false.
cooling=.true.
g_star=1.6666D0
n_star=0.1D0
eps_star=0.0D0
t_star=0.0d0 !this basically set the SF=0
/

IF one wants to activate SF:
&PHYSICS_PARAMS
cooling=.true.
g_star=1.6666D0
n_star=0.1D0
eps_star=0.05d0
t_star=3.0d0
T2star = 1.0d4
/

For SN feedback:
&PHYSICS_PARAMS
yield = 0.1d0
etaSN = 0.1d0
rbubble = 150.0d0
fek = 0.5d0
/

If movie=.true. , then:
the movie_vars=1,1,0,0,0,0,0,0 mean:
0: temp, 1: dens, 2: vx, 3: vy, 4: vz, 5: pres, 6: dm, 7: stars

Modern computer clusters often use schedulers like SLURM. I have been using SLURM in the  CHIMERA cluster in Prague (https://gitlab.mff.cuni.cz/mff/hpc/clusters, but you probably need access for that...). 
The PoR runs with MPI in 16 tasks in the following case:
[samarasn@hpc-head b1500]$ cat slurm
#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --mail-user=nicksam@sirrah.troja.mff.cuni.cz
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="optb1500"
#SBATCH -N 2
#SBATCH -n 16
#SBATCH --mem-per-cpu=50G
#SBATCH -p ffa

##Chimera
module load oneapi/mpi 
srun ~/bonnpor/PoR_hydro/ramses/bin/NWramses3d b1500.nml

##Karolina 
ml OpenMPI/4.1.4-GCC-11.3
srun ~/bonnpor/PoR_hydro/ramses/bin/NWramses3d b1500.nml

Remember that RAMSES (https://ramses-organisation.readthedocs.io/en/latest/wiki/Amr.html) :
1.4 * (ngridmax / 1e6) + 0.7 * (npartmax / 1e7) 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

3. Start analyzing the PoR simulation:
a) With gnuplot, one can have a nice movie of the evolution, like this one: https://www.youtube.com/watch?v=XDeCCL_ln-k
i) You need to specify in the aexp.c file how many outputs you have in your movie files (lines 19,20,24 and 38)
ii) You need to specify the boxlength in the rmviter.gp, where the plot is actually happening (line 37 the ’100’ in my file)
ii) in the rmvloop.gp you need to specify half the boxlength on line 34 and you can of course change other things there
iv) run ./aexp.c
v) ./a.out
vi) gnuplot rmvloop.gp
In case this does not work, there is a stupid solution which always works. One needs to comment out from all the info_000XXX files the last lines right after unit_t with #comment the lines specifying the ordering and the DOMAIN ind_min ind_max.

8. AHF convert end with error
   no mpi + makefile.config+ ahf_halos.c varies with model

9. pynbody???
