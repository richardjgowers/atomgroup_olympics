"""Doing AtomGroups as numpy structured arrays 2 - Electric Boogaloo

ATOMTYPE - a numpy dtype that represents an atom

MasterGroup - where all the data lives

StrucAtomGroup - a list of indices, which refer to the master list

StrucAtom - the new Atom, with a single index

"""
import numpy as np
import MDAnalysis as mda

ATOMTYPE = np.dtype([('number', 'i4'),
                     ('name', 'S10'),
                     ('type', 'S10'),
                     ('resname', 'S10'),
                     ('resid', 'i4'),
                     ('resnum', 'i4'),
                     ('segid', 'S10'),
                     ('mass', 'f8'),
                     ('charge', 'f8')])


class MasterGroup(object):
    """The master list of Atoms, only 1 should exist

    This is where the actual record of properties lives

    getitem on this only returns indices and a pointer back
    """
    def __init__(self, array):
        self._indices = np.arange(len(array))
        self._atoms = array

    def names(self):
        return self._atoms['name']

    def charges(self):
        return self._atoms['charge']

    def resids(self):
        return self._atoms['resid']

    def resnames(self):
        return self._atoms['resname']

    def set_names(self, new):
        self._atoms['name'] = new

    def set_charges(self, new):
        self._atoms['charge'] = new

    def set_resids(self, new):
        self._atoms['resid'] = new

    def set_resnames(self, new):
        self._atoms['resname'] = new

    def __getitem__(self, item):
        if isinstance(item, int):
            return StrucAtom(item, self._atoms)
        elif isinstance(item, (np.ndarray, slice)):
            return StrucAtomGroup(self._indices[item], self._atoms)


class StrucAtomGroup(object):
    def __init__(self, array, master):
        # Receives an array of indices and a pointer to the master list
        self._atoms = array
        self._master = master

    def names(self):
        return self._master['name'][self._atoms]

    def charges(self):
        return self._master['charge'][self._atoms]

    def resids(self):
        return self._master['resid'][self._atoms]

    def resnames(self):
        return self._master['resname'][self._atoms]

    def set_names(self, new):
        self._master['name'][self._atoms] = new

    def set_charges(self, new):
        self._master['charge'][self._atoms] = new

    def set_resids(self, new):
        self._master['resid'][self._atoms] = new

    def set_resnames(self, new):
        self._master['resname'][self._atoms] = new

    def __getitem__(self, item):
        if isinstance(item, int):
            return StrucAtom(self._atoms[item], self._master)
        elif isinstance(item, (np.ndarray, slice)):
            return StrucAtomGroup(self._atoms[item], self._master)


class StrucAtom(object):
    def __init__(self, idx, master):
        # Recieves a single idx, and where to look to get info
        self.number = idx
        self.name = master['name'][idx]
        self.type = master['type'][idx]
        self.resname = master['resname'][idx]
        self.resid = master['resid'][idx]
        self.resnum = master['resnum'][idx]
        self.segid = master['segid'][idx]
        self.mass = master['mass'][idx]
        self.charge = master['charge'][idx]

    def __repr__(self):
        return("<Atom {idx}: {name} of type {t} of resname {rname}, "
                "resid {rid} and segid {sid}>".format(
                    idx=self.number + 1,
                    name=self.name,
                    t=self.type,
                    rname=self.resname, rid=self.resid, sid=self.segid))


def convert(atomgroup):
    """Take a traditional AtomGroup and spit out master atomgroup"""
    new_ag = np.zeros(len(atomgroup), dtype=ATOMTYPE)

    new_ag['number'] = atomgroup.indices()
    new_ag['name'] = atomgroup.names()
    new_ag['type'] = atomgroup.types()
    new_ag['resname'] = [a.resname for a in atomgroup]
    new_ag['resid'] = [a.resid for a in atomgroup]
    new_ag['resnum'] = [a.resnum for a in atomgroup]
    new_ag['segid'] = [a.segid for a in atomgroup]
    new_ag['mass'] = atomgroup.masses()
    new_ag['charge'] = [a.charge for a in atomgroup]

    return MasterGroup(new_ag)
