# Tinker: Software Tools for Molecular Design

<H2><B>Introduction</B></H2>

The Tinker molecular modeling software is a complete and general package for molecular mechanics and dynamics, with some special features for biopolymers. Tinker has the ability to use any of several common parameter sets, such as Amber (ff94, ff96, ff98, ff99, ff99SB), CHARMM (19, 22, 22/CMAP), Allinger MM (MM2-1991 and MM3-2000), OPLS (OPLS-UA, OPLS-AA), Merck Molecular Force Field (MMFF), Liam Dang's polarizable model, and the AMOEBA, AMOEBA+ and HIPPO polarizable atomic multipole force fields. Parameter sets for other widely-used force fields are under consideration for future releases.

The Tinker software contains a variety of interesting algorithms such as: flexible implementation of atomic multipole-based electrostatics with explicit dipole polarizability, various continuum solvation treatments including several generalized Born (GB/SA) models, generalized Kirkwood implicit solvation for AMOEBA, an interface to APBS for Poisson-Boltzmann calculations, efficient truncated Newton (TNCG) local optimization, surface areas and volumes with derivatives, free energy calculations via the Bennett Acceptance Ratio (BAR) method, normal mode vibrational analysis, minimization in Cartesian, torsional or rigid body space, symplectic RESPA multiple time step integration for molecular dynamics, velocity Verlet stochastic dynamics, pairwise neighbor lists and splined spherical energy cutoff methods, particle mesh Ewald (PME) summation for partial charges and polarizable multipoles, a novel reaction field treatment of long range electrostatics, fast distance geometry metrization with better sampling than standard methods, Elber's reaction path algorithm, potential smoothing and search (PSS) methods for global optimization, Monte Carlo Minimization (MCM) for efficient potential surface scanning, tools for fitting charge, multipole and polarization models to QM-based electrostatic potentials and more....

<H2><B>Current Release</B></H2>

Tinker 8 is a major update of the Ponder Lab tool set for molecular mechanics and dynamics calculations. An important change in this new version is the switch from old-style common blocks to Fortran modules. Use of modules and greatly increased use of dynamic memory allocation means Tinker can now support very large molecular systems. Tinker 8 also implements improved OpenMP parallelization throughout additional parts of the code. Additional big improvements include parallel neighbor list building and updating, and big reduction in iteration needed to converge AMOEBA polarization via an efficient PCG solver. Other changes from the previous Tinker version include new and updated force field parameter sets and numerous minor additions and bug fixes, many of them suggested by users of the package. Please note that as with prior new releases, version 8 is neither backward nor forward compatible with earlier versions of Tinker. In particular, older versions of parameter files should not be used with Tinker 8 executables and vice versa. We strongly suggest users switch to Tinker 8 with its many important new features and bug fixes.

<H2><B>Building from Source Code</B></H2>

Tinker is provided as a complete source distribution. After unpacking the release, you can build a set of Tinker executables on almost any machine with a Fortran compiler. Makefiles, a cmake configure script, as well as standalone scripts to compile, build object libraries, and link executables on a wide variety of machine-CPU-operating system combinations are provided. To build using the Tinker Makefile, just copy /make/Makefile into the /source area of the distribution. Small changes are needed near the top of the file to set directory names and activate the appropriate operating system and compiler. Then issue the "make" command while in the /source directory. As noted below, the FFTW libraries must be available before building Tinker.

Tinker requires two object libraries from the FFTW Fourier transform package, libfftw3.a and libfftw3_threads.a. The FFTW libraries must be available in the one of the locations searched by the Tinker Makefile prior to building Tinker executables. While we do not provide the required FFTW libraries with the Tinker distribution, they are easy to build from the included FFTW source. Just follow the instructions in the /fftw/0README file. Optional support for APBS Poisson-Boltzmann calculations within Tinker requires object libraries from the APBS 1.3 software package, but this is not included in the default Tinker build.

If you do not want to build Tinker youself, pre-built Tinker executables for Linux, MacOS, and Windows are available for download from https://dasher.wustl.edu/tinker/. They should run on most recent vintage machines using these operating systems, and can handle a maximum of 1 million atoms provided sufficient memory is available. The Linux executables require at least glibc-2.6 or later. Note starting with Tinker 8, we no longer provide pre-built executables for any 32-bit operating systems.

<H2><B>MNDO related keywords</B></H2>
This is a listo of the keywords added or modified in the MNDO-Tinker interface.

<B>mndomm</B>: enable QM(MNDO)/MM calculation

<B>mndoexe</B>: specify the MNDO2020 command in the shell environment (default mndo2020)

<B>mndonoiterguess</B>: disable MNDO guess based on previous calculation (default is to do it). This option disable automatic keyword ipubo=1,ktrial=11,imomap=3.

<B>mndotemplate</B>: specify the MNDO template to create the MNDO input during a QM/MM calculation (default is template.inp). Note that according to the option provided, Tinker read the option of template file, modify them where needed and than reassemble the output from scratch.

<B>mndopostexe</B>: if specified, the script/excutable is run after the excution of MNDO. It can be used to archive files, extract information from MNDO logs etc.

<B>mndoallgrd</B>: ask mndo to compute gradients for all the electronic states included in the calculation (default is to only compute gradients on current state). This is only allowed with kci=5.

<B>mndostates</B>: number of electronic states to be computed in MNDO calculation (default is to only compute ground state).

<B>mndocurrentstate</B>: state to be used for the dynamic propagation (default is ground state).

<B>inactive</B>: same behavior of standard Tinker but if an MM atom is set inactive, its gradients are not computed in QM/MM calculation, actually reducing the overall cost.

<B>qmatoms</B>: specify which are the atoms in the QM part of the system. LA are automatically generated where needed. Indexes are provided following the xyz numbering.

<B>conjatoms</B>: for calculation with pipop keyword in mndo template, specify the atoms in the conjugate system. Indexes are provided following xyz numbering.

<H3><B>Link atom policies</B></H3>
Link atom is placed on the QM-MM bond at a distance of 1.1 A (controlled from variable distqmla in mndo.f at compile time) from the QM atom. Charges on MM atoms at 1 and 2 bonded of distance from the QM atoms are turned off in the QM/MM calculation (charges on atoms at two bonds from the QM part can be preserved setting mndo\_la13 to .false. in mndo.f at compile time). All the electrostatic integrals are included in the QM/MM coupling scheme in MNDO (mmlink=2).
