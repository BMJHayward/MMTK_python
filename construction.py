""" construction.py contains examples for non-standard generation of protein objects from PDB files """

from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.Proteins import PeptideChain, Protein

"""First problem: construct an all-atom model from a structure without hydrogens. This is the standard problem when using an all-atom force field with crystallographic structures.

Note: the simple solution in this case is just
    insulin = Protein('insulin.pdb')
but the explicit form shown below is necessary when any kind of modification is required."""

# Load PDB file
configuration = PDBConfiguration('insulin.pdb')

# Construct peptide chain objects. This also constructs positions for any missing hydrogens, using geometrical criteria.
chains = configuration.createPeptideChains()

# Make the protein object
insulin = Protein(chains)

# Write out structure with hydrogens to new file for use as input later
insulin.writeToFile('insulin_with_h.pdb')

""" Second problem: read a file with hydrogens and create a structure without them. This is useful for analysis; if you don't need the hydrogens it's faster without them. May be useful to compare with and without hydrogens"""

# Load the PDB file
configuration = PDBConfiguration('./insuling_with_h.pdb')

# Delete hydrogens
configuration.deleteHydrogens()

# Construct peptide chain objects without hydrogens
chains = configuration.createPeptideChains(model = 'no_hydrogens')

# make protein object
insulin = Protein(chains)

""" Third problem: cut off three residues from start of second chain before constructing protein. Useful for comparing incomplete structures. """

# Load PDB file
configuration = PDBConfiguration('insulin.pdb')

# Cut off first three residues of third chain
configuration.peptide_chains[2].removeResidues(0, 3)

# Mark the second chain as modified
for chain in configuration.peptide_chains:
    chain.modified = 0
configuration.peptide_chains[2].modified = 1

""" Construct peptide chain objects without hydrogens. For modified chains don't use N-terminal form for first residue"""
chains = []
for chain in configuration.peptide_chains:
    if chain.modified:
        chains.append(PeptideChain(chain, model='no_hydrogens', n_terminus=0))
    else:
        chains.append(PeptideChain(chain, model='no_hydrogens'))

# Make protein object
insulin = Protein(chains)

""" Stop here since 3rd problem is more an illustration. Can't be run directly because there is no 'special residue' definition. """
if 0:
    """ Fourth problem: construct a protein with a non-standard residue. New residue name is 'XXX' and definition is stored in MMTK group database under name 'special_residue'."""
    from MMTK.Proteins import defineAminoAcidResidue
    defineAminoAcidResidue('special_residue', 'xxx')
    protein = Protein('special_protein.pdb')
