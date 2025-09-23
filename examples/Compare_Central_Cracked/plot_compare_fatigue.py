r"""
.. _ref_compare_central_cracked_fatigue:

Central cracked specimen: FATIGUE
---------------------------------

"""


###############################################################################
# Import necessary libraries
# --------------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyvista as pv
import os
import sys

sys.path.insert(0, os.path.abspath('../../'))
plt.style.use('../../graph.mplstyle') 
import plot_config as pcfg

# results_folder = "results_compare_fatigue"
# if not os.path.exists(results_folder):
#     os.makedirs(results_folder)


###############################################################################
# Parameters definition
# ---------------------
# Define material and specimen parameters
E  = 210  # Young's modulus (kN/mm^2)
nu = 0.3  # Poisson's ratio (-)
Gc = 0.0027  # Critical strain energy release rate (kN/mm)
Cparis = 10e-9  # Paris' Law constant C
m = 2.5          # Paris' Law exponent m
Ep = E / (1.0 - nu**2) # Plane strain modulus (kN/mm^2)

# %
# Specimen geometry
B = 1.0  # Specimen thickness (mm)


###############################################################################
# Load results
# ------------
# Once the simulation finishes, the results are loaded from the results folder.
# The AllResults class takes the folder path as an argument and stores all
# the results, including logs, energy, convergence, and DOF files.
# Note that it is possible to load results from other results folders to compare results.
# It is also possible to define a custom label and color to automate plot labels.

results_pff = pd.read_csv("../Phase_Field_Central_Cracked/results_2_a03_l2/results.pff", delimiter="\t", comment="#", header=0)
results_pff_bourdin = pd.read_csv("../Phase_Field_Central_Cracked/results_2_a03_l2/results_corrected_bourdin.pff", delimiter="\t", comment="#", header=0)
results_pff_geometry = pd.read_csv("../Phase_Field_Central_Cracked/results_2_a03_l2/results_corrected_geometry.pff", delimiter="\t", comment="#", header=0)
label_2 = r"$a_0=0.3$ mm, $l=l_2$"
label_2_a = r"$a_0=0.3$ mm"
label_2_l = r"$l_2=0.0025$ mm"
color_2 = pcfg.color_orangered
a02 = 0.3


# %%
# From Linear elastic fracture mechanics theory
# :ref:`ref_lefm_center_cracked`

results_lefm =  pd.read_csv("../LEFM/results_central_cracked/a0_03.lefm_problem", delimiter="\t", comment="#", header=0)
LABEL_LEFM = r"LEFM"
color_var = pcfg.color_black

# %%
# From Linear elastic fracture mechanics theory
# :ref:`ref_lefm_center_cracked`
SCHEME_3 = np.loadtxt("../LEFM/results_central_cracked/center_cracked.lefm", delimiter="\t", skiprows=1)
a_lefm = SCHEME_3[:,0]
k_lefm = 1/SCHEME_3[:,2]
c_lefm = SCHEME_3[:,2]
dcda_lefm = SCHEME_3[:,3]
color_lefm = pcfg.color_black
pc_lefm = np.sqrt(2*B*Gc/(dcda_lefm))

 
markevery_2 = max(1, len(results_pff["displacement"])//20)


# %%
# From FEM elastic force controlled simulation 
# :ref:`ref_elasticity_center_cracked_force_controlled`
force_control = pd.read_csv("../Elasticity/results_center_cracked_force_control/results.elasticity", delimiter="\t", comment="#", header=0)
LABEL_FORCE = r"Force-controlled (FEM)"
a_forc = force_control["crack_length"]
k_forc = force_control["stiffness"]
c_forc = 1/k_forc
dCda_forc = np.gradient(c_forc, a_forc)
pc_forc = np.sqrt(2*B*Gc/(0.5*dCda_forc))
markevery_forc =  max(1, len(a_forc)//10)

###############################################################################
# Plot: Gc-Corrected Crack Surface Area vs Stiffness
# --------------------------------------------------
# This plot demonstrates the improved alignment between phase-field predictions
# and LEFM theory when applying the Gc-based crack surface correction method.
# The gamma_corrected_Gc value adjusts the raw crack surface area to account
# for the diffuse nature of the phase-field representation.
#
# Key observations:
# 1. Corrected results align much better with the LEFM reference curve
# 2. Both length scales (l1 and l2) produce similar results after correction
# 3. This confirms that phase-field models can accurately predict structural
#    stiffness for varying crack lengths when properly calibrated
fig, ax_reaction = plt.subplots()

# LEFM reference
ax_reaction.plot(results_lefm["a"], results_lefm["P"]/results_lefm["u"], color=pcfg.color_black, linestyle='-', label=LABEL_LEFM)

ax_reaction.plot(results_pff["gamma"], results_pff["stiffness"], color=pcfg.color_blue, linestyle='-', label=r"$\int \gamma(\phi, \nabla \phi) d\Omega$", markevery=markevery_2, marker='o')
ax_reaction.plot(results_pff_bourdin["gamma"], results_pff_bourdin["stiffness"], color=pcfg.color_orangered, linestyle='-', label=r"$1+h/(2l)$", markevery=markevery_2, marker='s')
ax_reaction.plot(results_pff_geometry["gamma"], results_pff_geometry["stiffness"], color=pcfg.color_green, linestyle='--', label=r"Skeleton", markevery=markevery_2, marker='^')

ax_reaction.set_xlabel(pcfg.crack_length_label)
ax_reaction.set_ylabel(pcfg.stiffness_label)
ax_reaction.legend()

# plt.savefig(os.path.join(results_folder, "gamma_compare_vs_k"))


###############################################################################
# Crack length vs stiffness
# -------------------------
# The stiffness as function of crack length is plotted for the three methods.
fig, ax0 = plt.subplots()

ax0.plot(a_lefm, pc_lefm, color=pcfg.color_black, linestyle='-', label=LABEL_LEFM)

ax_reaction.plot(a_lefm, k_lefm, color=pcfg.color_black, linestyle='-', label=LABEL_LEFM)

ax0.plot(results_pff["gamma"], np.sqrt(2*2*B*Gc/results_pff["dCda"]),  color=pcfg.color_blue, linestyle='-', label=pcfg.gamma_ref_label, markevery=markevery_2, marker='o')
ax0.plot(results_pff_bourdin["gamma"], np.sqrt(2*2*B*Gc/results_pff_bourdin["dCda"]), color=pcfg.color_orangered, linestyle='-', label=pcfg.gamma_bourdin_label, markevery=markevery_2, marker='s')
ax0.plot(results_pff_geometry["gamma"], np.sqrt(2*2*B*Gc/results_pff_geometry["dCda"]),  color=pcfg.color_green, linestyle='--', label=pcfg.gamma_geometry_label, markevery=markevery_2, marker='^')

ax0.set_xlim(left=0.28, right=1.1)
ax0.set_ylim(bottom=-0.07, top=2.0)
# Enhance plot aesthetics
ax0.set_xlabel(pcfg.crack_length_label)
ax0.set_ylabel(pcfg.critical_force_label)
# ax0.legend()

# Save the figure
# plt.savefig(os.path.join(results_folder, "critical_force_vs_crack_length_pff"))



###############################################################################
# Fatigue
# -------
# Once the compliance curves are obtained, it is possible to calculate the fatigue lives from the compliance respect the crack area for the different methods.
# In this case the fatigue analysis is performed for a crack that goes from an initial crack length `value_1` to a final crack length `value_2`.

def slice_array_by_values(a, value_1, value_2):
    """
    Returns a slice of the array `a` between the indices of the nearest values to `value_1` and `value_2`.

    Parameters:
    a (numpy.ndarray): The input array.
    value_1 (float): The first value to find in the array.
    value_2 (float): The second value to find in the array.

    Returns:
    numpy.ndarray: A new array sliced between the indices of the nearest values to `value_1` and `value_2`.
    """
    # Find the indices of the nearest values
    index_1 = (np.abs(a - value_1)).argmin()
    index_2 = (np.abs(a - value_2)).argmin()

    # Ensure index_1 is less than index_2
    if index_1 > index_2:
        index_1, index_2 = index_2, index_1

    # Return the sliced array
    return index_1, index_2 + 1

Ni = 0.0
AP = 33.0

# %%
# Slice the arrays to obtain the fatigue region
a0 = 0.4
af = 1.2

# i_o_lefm, i_f_lefm = slice_array_by_values(results_lefm['a'],  a0, max(results_lefm['a']))
# a_fatigue_lefm, c_fatigue_lefm = results_lefm['a'][i_o_lefm:i_f_lefm], results_lefm['gamma'][i_o_lefm:i_f_lefm]
# dCda_fatigue_lefm = 2*results_lefm['dCda'][i_o_lefm:i_f_lefm]

solution_lefm = np.loadtxt("../LEFM/results_central_cracked/a_04_099.lefm_fatigue", skiprows=1)


i_o_gamma, i_f_gamma = slice_array_by_values(results_pff["gamma"],  a0, af)
a_fatigue_gamma, c_fatigue_gamma = results_pff["gamma"][i_o_gamma:i_f_gamma], results_pff["gamma"][i_o_gamma:i_f_gamma]
dCda_fatigue_gamma = 0.5*results_pff["dCda"][i_o_gamma:i_f_gamma]
# dCda_fatigue_gamma = np.gradient(simulation_2["compliance"][i_o_gamma:i_f_gamma], 2*simulation_2["gamma"][i_o_gamma:i_f_gamma])


i_o_bourdin, i_f_bourdin = slice_array_by_values(results_pff_bourdin["gamma"],  a0, af)
a_fatigue_bourdin, c_fatigue_bourdin = results_pff_bourdin["gamma"][i_o_bourdin:i_f_bourdin], results_pff_bourdin["gamma"][i_o_bourdin:i_f_bourdin]
dCda_fatigue_bourdin = 0.5*results_pff_bourdin["dCda"][i_o_bourdin:i_f_bourdin]
# dCda_fatigue_bourdin = np.gradient(simulation_2_bourdin["compliance"][i_o_bourdin:i_f_bourdin], 2*simulation_2_bourdin["gamma"][i_o_bourdin:i_f_bourdin])

i_o_geometry, i_f_geometry = slice_array_by_values(results_pff_geometry["gamma"],  a0, af)
a_fatigue_geometry, c_fatigue_geometry = results_pff_geometry["gamma"][i_o_geometry:i_f_geometry], results_pff_geometry["gamma"][i_o_geometry:i_f_geometry]
dCda_fatigue_geometry = 0.5*results_pff_geometry["dCda"][i_o_geometry:i_f_geometry]
dCda_fatigue_geometry2 = np.gradient(results_pff_geometry["compliance"][i_o_geometry:i_f_geometry], 2*results_pff_geometry["gamma"][i_o_geometry:i_f_geometry])

i_o_force, i_f_force = slice_array_by_values(force_control["crack_length"],  a0, af)
a_fatigue_force, c_fatigue_force = force_control["crack_length"][i_o_force:i_f_force], force_control["stiffness"][i_o_force:i_f_force]

dCda_fatigue_force = np.gradient(1/force_control["stiffness"][i_o_force:i_f_force], 2*force_control["crack_length"][i_o_force:i_f_force])



from scipy.integrate import cumulative_trapezoid, cumulative_simpson

Nf_dCda_gamma = Ni + 1/(Cparis * (Ep/(2*B))**(m/2) * AP**m)*cumulative_trapezoid(1/(dCda_fatigue_gamma)**(m/2), a_fatigue_gamma, initial=0)
Nf_dCda_bourdin = Ni + 1/(Cparis * (Ep/(2*B))**(m/2) * AP**m)*cumulative_trapezoid(1/(dCda_fatigue_bourdin)**(m/2), a_fatigue_bourdin, initial=0)
Nf_dCda_geometry = Ni + 1/(Cparis * (Ep/(2*B))**(m/2) * AP**m)*cumulative_trapezoid(1/(dCda_fatigue_geometry)**(m/2), a_fatigue_geometry, initial=0)

Nf_dCda_force = Ni + 1/(Cparis * (Ep/(2*B))**(m/2) * AP**m)*cumulative_trapezoid(1/(dCda_fatigue_force)**(m/2), a_fatigue_force, initial=0)

###############################################################################
# Crack length vs stiffness
# -------------------------

# fig, ax4 = plt.subplots()
# ax4.plot(a_fatigue_lefm, dCda_fatigue_lefm, color=COLOR_LEFM, linestyle='-', marker='^', label=LABEL_LEFM, markevery=markevery_fatigue_lefm)

# # ax4.plot(a_fatigue_disp, dCda_fatigue_disp, color=COLOR_DISPLACEMENT, linestyle='-', marker='o', label=LABEL_DISPLACEMENT, markevery=markevery_fatigue_disp)
# ax4.plot(a_fatigue_forc, dCda_fatigue_forc, color=COLOR_FORCE, linestyle='--', marker='s', label=LABEL_FORCE, markevery=markevery_fatigue_forc)

# ax4.plot(a_fatigue_phas, dCda_fatigue_phas, color=COLOR_PHASEFIELD, linestyle='--', marker='D', label=r"Phase-field", markevery=markevery_phas)
# ax4.plot(a_fatigue_gc, dCda_fatigue_gc, color=COLOR_PHASEFIELD, linestyle='--', marker='D', label=r"Phase-field Gc", markevery=markevery_phas)

# # Enhance plot aesthetics
# ax4.set_xlabel(pcfg.crack_length_label)
# ax4.set_ylabel(pcfg.dCda_label)
# ax4.legend()
# plt.savefig(os.path.join(results_folder,"compare_dCda_phase_field"))


##############################################################################
# Crack length vs number of cycles
# --------------------------------

fig, ax5 = plt.subplots(figsize=(11.69, 5.85))
# ax5.plot(Nf_dCda_3, a_fatigue_lefm, linestyle='-', marker='^', label=LABEL_LEFM, markevery=markevery_fatigue_lefm)


LABEL_FORCE = r"Elasticity: Force-controlled"
ax5.plot(solution_lefm[:, 0], solution_lefm[:, 1], linestyle='--', color=pcfg.color_black, label="LEFM")
ax5.plot(Nf_dCda_force, a_fatigue_force,  color=pcfg.color_purple, linestyle='-.', marker='D', label=LABEL_FORCE, markevery=markevery_forc)


ax5.plot(Nf_dCda_gamma, a_fatigue_gamma,  color=pcfg.color_blue, linestyle='-', label=r"Phase-Field: Gamma", markevery=markevery_2, marker='o')
ax5.plot(Nf_dCda_bourdin, a_fatigue_bourdin, color=pcfg.color_orangered, linestyle='-', label=r"Phase-Field: Bourdin", markevery=markevery_2, marker='s')
ax5.plot(Nf_dCda_geometry, a_fatigue_geometry,  color=pcfg.color_green, linestyle='--', label=r"Phase-Field: Skeleton", markevery=markevery_2, marker='^')



ax5.set_xlim(left=0.28, right=11000)
ax5.set_ylim(bottom= 0.37, top=1.1)
# Enhance plot aesthetics
ax5.set_ylabel(pcfg.crack_length_label)
ax5.set_xlabel(pcfg.cycles_label)
ax5.legend()
# plt.savefig(os.path.join(results_folder, "compare_paris_law_phase_field"))

plt.show()
