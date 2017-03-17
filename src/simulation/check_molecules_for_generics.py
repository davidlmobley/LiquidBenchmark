#!/bin/env python

"""Process the molecules which density calculations were run for (yes, I already ran them before checking this! It would be better to check this BEFORE running them...) and see if any of them used generic parameters to bugs in smirff99Frosst."""

from density_simulation_parameters import DATA_PATH
from openeye.oechem import *
from smarty.forcefield import *
from smarty.forcefield_utils import get_molecule_parameterIDs
import glob
import os
import time

# Get filenames of mol2 files to check
mol2dir = os.path.join(DATA_PATH, 'monomers')
mol2files = glob.glob( mol2dir+'/*.mol2')

flavor = OEIFlavor_Generic_Default | OEIFlavor_MOL2_Default | OEIFlavor_MOL2_Forcefield

# Read OEMols
oemols = []
for mol2file in mol2files:
    istream = oemolistream(mol2file)
    istream.SetFlavor(OEFormat_MOL2, flavor)
    mol = OEMol()
    OEReadMolecule(istream, mol)
    oemols.append(mol)
    istream.close()

# load forcefield
ff = ForceField('smirff99Frosst.ffxml')

# Label molecules
init = time.time()
labels = ff.labelMolecules(oemols)
end = time.time()
print("Took %.1f seconds to label all %i molecules." % (end-init, len(oemols)))

# Check whether any of these use generics
generic_ids = ['n1', 'b1', 'a1', 't1']
for index, mol in enumerate(oemols):
    generic_used = False
    # Get parameter IDs used
    pids = [ pid for force, lists in labels[index].items() for (indices, pid, smarts) in lists ]
    # Check for use of generics
    for g in generic_ids:
        if g in pids:
            generic_used = True
            break

    if generic_used:    
        print("Generic used for molecule %s, SMILES %s..." % ( os.path.basename( mol2files[index]), OECreateIsoSmiString(mol)) ) 

