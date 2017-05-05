import numpy as np
import matplotlib
matplotlib.use('Agg')
import sklearn.metrics, sklearn.cross_validation
import statsmodels.formula.api as sm
import simtk.unit as u
import polarizability
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.pyplot as plt

from mpltools import style # mpltools from git@github.com:tonysyu/mpltools.git for better color styles
plt.style.use('ggplot')

# Extended from Beauchamp's plot_tbv.py to plot both Beauchamp's GAFF results and my SMIRNOFF results on parallel plots w same settings.

FIGURE_SIZE = (6.5, 6.5)
DPI = 1600

expt = pd.read_csv("../tables/data_with_metadata.csv")
expt["temperature"] = expt["Temperature, K"]

# Load SMIRNOFF results
pred = pd.read_csv("../tables/predictions.csv")
# Load GAFF results
pred_gaff = pd.read_csv("../tables/beauchamp_predictions.csv")

# Compute polarization correction, dielectric for both
for dataset in [pred, pred_gaff]:
    dataset["polcorr"] = pd.Series(dict((cas, polarizability.dielectric_correction_from_formula(formula, density * u.grams / u.milliliter)) for cas, (formula, density) in dataset[["formula", "density"]].iterrows()))
    dataset["corrected_dielectric"] = dataset["polcorr"] + dataset["dielectric"]

expt = expt.set_index(["cas", "temperature"])  # Can't do this because of duplicates  # Should be fixed now, probably due to the CAS / name duplication issue found by Julie.
#expt = expt.groupby(["cas", "temperature"]).mean()  # Fix a couple of duplicates, not sure how they got there.

# Tabulate some data for both sets and do plots
dataset_labels = ['SMIRNOFF', "GAFF"]
for idx,dataset in enumerate([pred, pred_gaff]):
    dataset = dataset.set_index(["cas", "temperature"])
    dataset["expt_density"] = expt["Mass density, kg/m3"]
    dataset["expt_dielectric"] = expt["Relative permittivity at zero frequency"]
    dataset["expt_density_std"] = expt["Mass density, kg/m3_uncertainty_bestguess"]
    dataset["expt_dielectric_std"] = expt["Relative permittivity at zero frequency_uncertainty_bestguess"]

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)

    for (formula, grp) in dataset.groupby("formula"):
        x, y = grp["density"], grp["expt_density"]
        xerr = grp["density_sigma"]
        yerr = grp["expt_density_std"].replace(np.nan, 0.0)
        x = x / 1000.  # Convert kg / m3 to g / mL
        y = y / 1000.  # Convert kg / m3 to g / mL
        xerr = xerr / 1000.  # Convert kg / m3 to g / mL
        yerr = yerr / 1000.  # Convert kg / m3 to g / mL
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='.', label=formula)

    plt.plot([.600, 1.400], [.600, 1.400], 'k', linewidth=1)
    plt.xlim((.600, 1.400))
    plt.ylim((.600, 1.400))
    plt.xlabel("Predicted (%s)" % dataset_labels[idx])
    plt.ylabel("Experiment (ThermoML)")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()

    x, y = dataset["density"], dataset["expt_density"]
    plt.title(r"Density [g cm$^{-3}$]")
    plt.savefig("densities_thermoml_%s.pdf" % dataset_labels[idx], bbox_inches="tight")
    plt.savefig("densities_thermoml_%s.tif" % dataset_labels[idx], bbox_inches="tight")

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)

    for (formula, grp) in dataset.groupby("formula"):
        x, y = grp["density"], grp["expt_density"]
        xerr = grp["density_sigma"]
        yerr = grp["expt_density_std"].replace(np.nan, 0.0)
        x = x / 1000.  # Convert kg / m3 to g / mL
        y = y / 1000.  # Convert kg / m3 to g / mL
        xerr = xerr / 1000.  # Convert kg / m3 to g / mL
        yerr = yerr / 1000.  # Convert kg / m3 to g / mL
        plt.errorbar(x - y, y, xerr=xerr, yerr=yerr, fmt='.', label=formula)

    plt.xlim((-0.1, 0.25))
    plt.ylim((.600, 1.400))
    plt.xlabel("Predicted (%s) - Experiment" % dataset_labels[idx])
    plt.ylabel("Experiment (ThermoML)")
    plt.gca().set_aspect('auto', adjustable='box')
    plt.draw()

    x, y = dataset["density"], dataset["expt_density"]
    plt.title(r"Density [g cm$^{-3}$]")

    plt.savefig("densities_differences_thermoml_%s.pdf" % dataset_labels[idx], bbox_inches="tight")
    plt.savefig("densities_differences_thermoml_%s.tif" % dataset_labels[idx], bbox_inches="tight")



    yerr = dataset["expt_dielectric_std"].replace(np.nan, 0.0)
    xerr = dataset["dielectric_sigma"].replace(np.nan, 0.0)

    plt.figure(figsize=FIGURE_SIZE, dpi=DPI)

    plt.xlabel("Predicted (%s)" % dataset_labels)
    plt.ylabel("Experiment (ThermoML)")
    plt.title("Inverse Static Dielectric Constant")


    plt.plot([0.0, 1], [0.0, 1], 'k')  # Guide


    x, y = dataset["dielectric"], dataset["expt_dielectric"]
    ols_model = sm.OLS(y, x)
    ols_results = ols_model.fit()
    r2 = ols_results.rsquared
    plt.errorbar(x ** -1, y ** -1, xerr=xerr * x ** -2, yerr=yerr * y ** -2, fmt='.', label=dataset_labels[idx])  # Transform xerr and yerr for 1 / epsilon plot

    plt.xlim((0.0, 1))
    plt.ylim((0.0, 1))
    plt.legend(loc=0)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.savefig("dielectrics_thermoml_nocorr_%s.pdf" % dataset_labels[idx], bbox_inches="tight")


    x, y = dataset["corrected_dielectric"], dataset["expt_dielectric"]
    ols_model = sm.OLS(y, x)
    ols_results = ols_model.fit()
    r2 = ols_results.rsquared
    plt.errorbar(x ** -1, y ** -1, xerr=xerr * x ** -2, yerr=yerr * y ** -2, fmt='.', label="Corrected")  # Transform xerr and yerr for 1 / epsilon plot

    plt.xlim((0.0, 1.02))
    plt.ylim((0.0, 1.02))
    plt.legend(loc=0)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()
    plt.savefig("dielectrics_thermoml_%s.pdf" % dataset_labels[idx], bbox_inches="tight")
    plt.savefig("dielectrics_thermoml_%s.tif" % dataset_labels[idx], bbox_inches="tight")
