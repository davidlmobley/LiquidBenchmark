import numpy as np
import scipy.stats
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

sns.set_palette("bright")
sns.set_style("white")
sns.set(font_scale=1.2)
sns.set(style="ticks")

x = pd.read_csv("./tables/timestep_dependence.csv")

power = 1.0
XY = x.set_index("timestep").ix[[0.5, 1.0]]
slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(XY.index.values ** power, XY.mu.values)
f = lambda timestep: slope * timestep ** power + intercept

#reference = x.set_index("timestep").mu[0.5]
reference = intercept

x.relerr = (x.mu - reference) / reference  # OVERWRITES WHAT WAS ALREADY THERE!

xgrid = np.linspace(0.0, 1.0, 100)
ygrid = np.array([f(timestep) for timestep in xgrid])
rel_grid = np.array([y / reference for y in ygrid]) - 1.0

plt.errorbar(x.timestep, x.mu, yerr=x.stderr, fmt='o')
plt.plot(xgrid, ygrid, 'k')
ticks = [0.0, 0.5, 1.0, 1.5, 2.5]
plt.xticks(ticks)
plt.ylabel("Density [g / cm^3]")
plt.xlabel("Timestep [fs]")
plt.title("Timestep dependence: Density")
plt.savefig("./manuscript/figures/timestep_dependence_density.pdf", bbox_inches="tight")

plt.figure()
plt.errorbar(x.timestep, x.relerr, yerr=x.stderr / reference, fmt='o')
plt.plot(xgrid, rel_grid, 'k')
plt.xticks(ticks)
plt.ylabel("Relative Error [unitless]")
plt.xlabel("Timestep [fs]")
plt.title("Timestep dependence: Relative Error")
plt.savefig("./manuscript/figures/timestep_dependence_relative_error.pdf", bbox_inches="tight")
