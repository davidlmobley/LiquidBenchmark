from simtk import unit as u
from liquid_tools import AmberMixtureSystem
from density_simulation_parameters import MOLECULES_PER_BOX

builder = AmberMixtureSystem(["CO"], [MOLECULES_PER_BOX], 300 * u.kelvin, 'smirff99Frosst.ffxml')
builder.run()
