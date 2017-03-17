#!/bin/env python

import os
import glob
from density_simulation_parameters import MOLECULES_PER_BOX

targetdir = '/work/cluster/dmobley/liquid_benchmark/02_06_2017'

# Get list of those we are trying to run
import pandas as pd

data = pd.read_csv("../../tables/data_dielectric.csv")

rank_by_name = {}
for k0, k1, components, smiles, cas, temperature, pressure, density, perm in data.itertuples():
    #print(k0, k1, components, smiles, cas, temperature, pressure, density)
    # Compose cryptic name of job
    name = cas + '_'+str(MOLECULES_PER_BOX)+'_'+str(temperature)
    #print(name)
    rank_by_name[name]=k0


# Get list of those which have been built so far
packmolfiles = glob.glob(os.path.join(targetdir, 'packmol_boxes', '*.pdb'))
packmolnames = [ os.path.basename(filenm).replace('.pdb','') for filenm in packmolfiles]
#print("First 10 names of packmol files:")
#print(packmolnames[0:10])

# Get list of those which have successfully completed
prodcsvs = glob.glob(os.path.join(targetdir, 'production', '*.csv'))
csvnames = [ os.path.basename(filenm).replace('_production.csv','') for filenm in prodcsvs]
#print("First 10 names of csv files:")
#print(csvnames[0:10])

# Get list of those which have NOT successfully completed
incomplete = []
for filenm in packmolnames:
    if not filenm in csvnames:
        incomplete.append(filenm)
print("There are %s runs which have not completed:" % len(incomplete))
for jobname in incomplete:
    print('   %s, rank %s' % (jobname, rank_by_name[jobname]))

print("There are %s runs which have completed, out of %s total." % (len(prodcsvs), len(rank_by_name))) 
