import pandas as pd
from simtk import unit as u
from liquid_tools import AmberMixtureSystem
from density_simulation_parameters import MOLECULES_PER_BOX
import sys

data = pd.read_csv("../../tables/data_dielectric.csv")
rank = int(sys.argv[1])

for k0, k1, components, smiles, cas, temperature, pressure, density, perm in data.itertuples():
    if k0 == rank:
        print(k0, k1, components, smiles, cas, temperature, pressure, density)
        builder = AmberMixtureSystem([cas], [MOLECULES_PER_BOX], temperature * u.kelvin, 'smirff99Frosst.ffxml')
        builder.run()
