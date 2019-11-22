# HPSC Lab 13
2019-11-22

Click to [make your own repo](https://classroom.github.com/a/rw4AstkF).

The goals for this lab are:
* Use GROMACS

-----

[GROMACS](http://www.gromacs.org/) is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles.

## Download and install GROMACS

```
wget "http://ftp.gromacs.org/pub/gromacs/gromacs-2019.4.tar.gz"
tar xfz gromacs-2019.4.tar.gz
cd gromacs-2019.4
mkdir build & mkdir install
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DCMAKE_INSTALL_PREFIX=../install
make
make check
make install
source ../install/bin/GMXRC
```

-----

The git repo contains a .pdb file downloaded from the [Protein Data Bank](https://www.rcsb.org/structure/1AKI).  1AKI is the orthorhombic form of a hen egg-white lysozyme at 1.5 Angstroms resolution.  We are going to use GROMACS to simulate this protein in a cube of water with a handful of ions.

These instructions are heavily inspired by Justin A. Lemkul's [GROMACS Tutorials](http://www.mdtutorials.com/gmx/index.html).

## Generate Topology

Strip out crystal water entries from the file and make sure there are no "MISSING" entries:
```
grep -v HOH 1aki.pdb > 1AKI_clean.pdb
grep MISSING 1AKI_clean.pdb  ## should return nothing
```

Convert the pdb file into GROMACS files, choosing the OPLS-AA/L force field (15) when prompted:
```
gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce
```

This shoudl create three new files: 1AKI_processed.gro, topol.top, and posre.itp.  1AKI_processed.gro is a GROMACS-formatted structure file that contains all the atoms defined within the force field (i.e., H atoms have been added to the amino acids in the protein).  The topol.top file contains the system topology.  The posre.itp file contains information used to restrain the positions of heavy atoms.

Looking at the topol.top file, one can find a list of the atoms in the hen protein (Protein_chain_A), along with bonds, pairs, angles, and dihedrals.  These parameters are discussed in more detail in the GROMACS manual chapter 5.  At the bottom of the file, position restraint file, solvent, and ion information files are included.  (The solvent model is SPC/E water.)  After those, there is a list of molecules in the system.  Right now, there is just the one Protein_chain_A; we need to add other molecules.

## Define Water Cube

Define a cube (-bt) with the protein in the center (-c), with at least (-d) 0.2 nm to the edges:
```
gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 0.2 -bt cubic
```

Using the file containing the protein configuration in the cube, fill the cube with solvent (water):
```
gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top
```
"spc216.gro is a generic equilibrated 3-point solvent model."

Note that the topol.top file now has the solvent under the molecule list (4471 molecules).  "Note that if you use any other (non-water) solvent, solvate will not make these changes to your topology!"

## Adding Ions

The protein has a sharge of +8e (as summarized by the output of pdb2gmx).  We need to balance this by adding ions to the system.

To do this, we need a file that lists the details of each atom.  To obtain that file, we need another file to tell GROMACS what we want.  This file is `ions.mdp`.  Use this to create a .tpr file that is used by the ion generating tool.  When prompted, choose group 13 "SOL" in which to substitute in ions.

```
gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

Looking again at the bottom of the topol.top file, we see the addition of 8 chlorine ions to the molecule list.

-----

The system is now assembled!  Before simulating the dynamics, though, we need a good initial state.

## Energy Minimization

Similar process as adding ions, except we need to use an energy minimization tool instead.  Also, things kick off with a slightly different .mdp file:

```
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
```

Make sure that the reported potential energy is on the order of 10^5 - 10^6, and that the maximum force is less than the 1000 set earlier, otherwise the system may not be unstable.

This generates four files:
* em.log: ASCII-text log file of the EM process
* em.edr: Binary energy file
* em.trr: Binary full-precision trajectory
* em.gro: Energy-minimized structure

Examine the potential energy over time to make sure that it has decreased.  At the prompt, type `10 0` to choose just potential energy.  Then look at the time date of the generated potential.xvg file to ensure energy decreases.
```
gmx energy -f em.edr -o potential.xvg
```

## Equilibration

The energy minimization mostly optimized the solvent with respect to itself.  The solvent still needs to be brought to the proper temperature and orientation with respect to the protein.

### Phase 1

The earlier generated posre.itp file allows us to constrain the protein while equilibrating the solvent.  The first equilibration phase is conducted under a constant number of particles, volume, and temperature to stabilize the temperature.  We run for 5000 steps of 2 femptoseconds.

```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v
```

This time let's look at temperature using the energy module (type `16 9` when prompted):
```
gmx energy -f nvt.edr -o temperature.xvg
```

From this we can see that the temperature stabilizes very quickly.

### Phase 2

Continuing, the second equilibrium phase is conducted under a constant number of particles, pressure, and temperature to stabilize the pressure.

```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -v
```

Now monitor the pressure progression (`18 0` at prompt):
```
gmx energy -f npt.edr -o pressure.xvg
```

The pressure looks to be fluctuating wildly.  This is to be expected.  Is it too wild?  Let's look at density (`24 0` at prompt) instead:

```
gmx energy -f npt.edr -o density.xvg
```

This shows more stability over time.

-----

## Simulation

We are ready to release the position restraints and run a molecular dynamics simulation!

```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1 -v
```

The simulation time has been reduced from the original by a factor of 100 to complete in a reasonable amount of time.

GROMACS has several tools for analysis.  We will look at a few basic ones.

`trjconv` is used as a post-processing tool to strip out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc).  Use it to remove the physical periodicity.  At the prompts that follow, choose to center the protein (`1`) and output the system (`0`).
```
gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
```

With this corrected trajectory, we can look at the structural stability over time.  At the prompts that follow, choose backbone (`4`) for both.
```
gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
```

The resulting .xvg file won't be very interesting in this instance because the duration was so short, but we can observe some stability begin to arise.