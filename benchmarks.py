"""Stuff for atomgroups to try and do!

These should all be runnable in ipython behind a %time magic function

(%timeit causes caching often, so isn't "fair")

Changing and querying stuff is done on 3 properties
 name - representing strings
 resid - representing ints
 charge - representing floats

"""
import MDAnalysis as mda

from atomtypes import convert

def query_atom_names(atomgroup):
    return atomgroup.names()

def query_atom_resids(atomgroup):
    return atomgroup.resids()

def query_atom_charges(atomgroup):
    return atomgroup.charges()

def change_atom_names(atomgroup, newnames):
    atomgroup.set_name(newnames)

def change_atom_resids(atomgroup, resids):
    atomgroup.set_resid(resids)

def change_atom_masses(atomgroup, masses):
    atomgroup.set_mass(masses)

u = mda.Universe('big.gro')

oldag = u.atoms
newag = convert(oldag)
