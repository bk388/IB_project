#!/usr/local/bin/python
import numpy as np

import ase
import ase.lattice
from ase.optimize.precon import PreconLBFGS, Exp

# Initialise the atomistic energy model
import Morse
calc = Morse.MorsePotential()

# Create the structure of the bar. The directions in the function below
# orient the crystal such that the dislocations are easer to see in this
# quasi 2-dimensional geometry. "pbc" refers to "Periodic Boundary Conditions"
# and are set to be active in the Z direction, but not in X and Y
unitcell = ase.lattice.cubic.FaceCenteredCubic(directions=[[1,-1,0], [1,1,-2], [1,1,1]],
                             size=(1,1,1), symbol='Cu', pbc=(0,0,1))

# the unit cell is replicated in the X and Y directions, but kept minimal in Z
atoms = unitcell*(20,40,1)

# assign the energy model to the atoms object
atoms.set_calculator(calc)

# obtain the coordinates of the centre of the bar
cx = atoms.cell[0,0]/2.0
cy = atoms.cell[1,1]/2.0
cz = atoms.cell[2,2]/2.0

# make two notches at the top and bottom
del atoms[[atom.index for atom in atoms if abs(atom.position[0]-cx) < 5 and abs(atom.position[1]-cy) > 45 ]]

# add vaccum
#c = atoms.cell
#c[0,0] += 10.0
#c[1,1] += 10.0
#c[2,2] += 10.0
#atoms.set_cell(c, scale_atoms=False)

# calculate the extent of the bar in the X direction
left = atoms.positions[:, 0].min()
right = atoms.positions[:, 0].max()

# We will pull on the right and left of the bar, so we need the atoms on the
# two edges not to move by themselves, as if they were "gripped" by our "tension test apparatus"
# we achive this by selecting the "gripped" atoms, and calling an ASE function on them
fixed_mask = ((abs(atoms.positions[:, 0] - right) < 3.0) |
              (abs(atoms.positions[:, 0] - left) < 3.0))
fix_atoms = ase.constraints.FixAtoms(mask=fixed_mask) # this creates the constraint
atoms.set_constraint(fix_atoms) # this stores the constraint in the atoms object


# We now define some helper functions that record the experiment

# We will save the positions of the atoms during the experiment to this file
trajectory = open("traj2d_qs.xyz", mode='w')

def write_frame(atoms):
    # the energy of each atom is not saved by default, so we explicitly copy
    # it from the results of the calculation into the atoms object, this
    # way it gets saved to our file and we can use it to help visualise
    atoms.arrays['local_energy'] = atoms.calc.results['local_energy']
    ase.io.write(trajectory, atoms, format="extxyz")
    del atoms.arrays['local_energy'] # we delete this array, otherwise ASE breaks, so this is a "workaround"
    

# write out the initial atom positions, but first call the energy model so that
# the invidiual atomic energies are also saved in addition to the position. 
atoms.get_potential_energy()
write_frame(atoms)


nu = 0.25 # assumed Poisson Ratio  - replace this by what you calculated 

# record the original width of the simulation box, so that we can compute Engineering strain later
origLx = atoms.cell[0,0]



# Now apply the starting strain. We could start from zero strain, but
# it saves time to jump to a finite strain right away. We could move
# every atom position explicitly, but it is easier to mimic this by
# scaling the "cell" variable, but we need to explicitly ask ASE to
# move all the atoms with this call, using the scale_atoms argument
s = 0.05   # starting strain. 
atoms.set_cell(atoms.cell*np.diag([1.0+s,1.0-s*nu,1.0-s*nu]), scale_atoms=True)

# simple header for the output
print("        Iteration       Time     Energy     max Force")

# loop to perform the successive small strain increments
for i in range(50):
    
    # record and print the current strain
    atoms.info["strain"] = (atoms.cell[0,0]-origLx)/origLx # Engineering strain
    print "strain: ", atoms.info["strain"]

    # set up an optimizer, this is a particularly efficient one
    opt = PreconLBFGS(atoms, precon=Exp(3.0))

    # attach the optimiser to the atoms object, asking it to call our helper function 
    # that writes out the atomic positions, after every 2nd iterations of the optimisation
    opt.attach(write_frame, 2, atoms)

    # run the optimizer, until the maximum force on any atom is less than a tolerance in eV/Ang
    opt.run(fmax=0.02)

    # update the "grip" by applying the small strain increment. 
    s = 0.01 # small strain increment
    atoms.set_cell(atoms.cell*np.diag([1.0+s,1.0-s*nu,1.0-s*nu]), scale_atoms=True)

    
    
#######################################################################
# 
# Notes
#
#######################################################################


# A much faster model is available in the "quippy" package, which in turn comes with the QUIP package
#                                             http://github.com/libatoms/QUIP
# whose installation is more complicated (needs a Fortran compiler, e.g. gfortran, and the use of the
# "make" command, for more details see the README in QUIP)
#
# import quippy
# calc = Potential('IP Morse', param_filename='ip.parms.Morse.Cu_PhysicsLetters_v23_p309_y1966.xml')
