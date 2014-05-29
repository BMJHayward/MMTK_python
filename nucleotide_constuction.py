""" nucleotide_construction.py creates a nucleotide chain with a ligand from PDB file """

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.NucleicAcids import NucleotideChain
from MMTK.Visualization import view

""" Load PDB entry 110d. It contains a single DNA strand with a ligand daunomycin """
configuration = PDBConfiguration('110d.pdb')

""" Construct nucleotide chain object. This also constructs positions for missing hydrogens, using geometrical criteria. """

chain = configuration.createNucleotideChains()[0]

""" Construct the ligand. There is no definition of it in the database, so it can only be constructed as a collection of atoms. The second argument of createMolecules() is set to one in order to allow this use of an unknown residue. """
ligand = configuration.createMolecules(['DM1'], 1)

# Put everything in a universe and show it graphically
universe = InfiniteUniverse()
universe.addObject(chain)
universe.addObject(ligand)

view(universe)