r"""
.. _ref_paper_Wagner_phd:

Visualization of Experimental and Simulation Data from Wagner's PhD Thesis
--------------------------------------------------------------------------

This script reproduces several key figures from the PhD thesis by Wagner :footcite:t:`example_Wagner2018_phd_thesis`, which investigates crack propagation in complex specimens.

The plots compare experimentally measured crack paths with simulation results from different numerical methods, including FRANC3D and ZFEM. The data for each figure has been digitized from the original thesis and is stored in text files.

.. footbibliography::

The script generates the following figures:
-------------------------------------------
- **Figure 5.10**: Compares the experimental crack path with simulation data for specimen 17-157.
- **Figure 5.11**: Compares experimental results with FRANC3D and ZFEM simulations for specimen 17-159.
- **Figure 5.12a**: Compares experimental results with FRANC3D and ZFEM simulations for specimen 17-155.

"""

###############################################################################
# Import required libraries
# -------------------------
# Import standard libraries for numerical operations and plotting, and set up the
# project-specific plotting style and configuration.
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add the parent directory to the system path to allow imports from the project root.
sys.path.insert(0, os.path.abspath('../../'))
plt.style.use('../../graph.mplstyle')
import plot_config as pcfg


###############################################################################
# Define plot labels and styles
# -----------------------------
# Standardize the labels and visual styles for different data sources to ensure
# consistency across all plots.
label_experimental = "Experimental"
color_shape_experimental = 'b.'

label_franc3d = "FRANC3D"
color_shape_franc3d = 'r.'

label_zfem = "ZFEM"
color_shape_zfem = 'g.'

###############################################################################
# Figure 5.10: Crack Path Comparison for Specimen 17-157
# ------------------------------------------------------
# This plot reproduces Figure 5.10 from the thesis, comparing the experimental
# crack path with simulation and reference data for specimen 17-157.

# Load data from text files
fig510_frame = np.loadtxt("Wagner_phd/fig510/reference.txt")
fig510_experimental = np.loadtxt("Wagner_phd/fig510/experimental.txt")
fig510_simulation = np.loadtxt("Wagner_phd/fig510/simulation.txt")  # FRANC3D simulation data

# Create the plot
fig, ax0 = plt.subplots()
ax0.plot(fig510_frame[:, 0], fig510_frame[:, 1], 'k.', label="Reference")
ax0.plot(fig510_experimental[:, 0], fig510_experimental[:, 1], color_shape_experimental, label=label_experimental)
ax0.plot(fig510_simulation[:, 0], fig510_simulation[:, 1], color_shape_franc3d, label=label_franc3d)

# Configure axis labels, legend, and aspect ratio
ax0.set_xlabel("X (mm)")
ax0.set_ylabel("Y (mm)")
ax0.set_aspect('equal', adjustable='box')
ax0.legend()


###############################################################################
# Figure 5.11: Crack Path Comparison for Specimen 17-159
# ------------------------------------------------------
# This plot reproduces Figure 5.11, comparing the experimental crack path
# with results from both FRANC3D and ZFEM simulations for specimen 17-159.

# Load data from text files
fig511_frame = np.loadtxt("Wagner_phd/fig511/reference.txt")
fig511_experimental = np.loadtxt("Wagner_phd/fig511/experimental.txt")
fig511_FRANC3D = np.loadtxt("Wagner_phd/fig511/FRANC3D.txt")
fig511_ZFEM = np.loadtxt("Wagner_phd/fig511/ZFEM.txt")

# Create the plot
fig, ax0 = plt.subplots()
ax0.plot(fig511_frame[:, 0], fig511_frame[:, 1], 'k.', label="Reference")
ax0.plot(fig511_experimental[:, 0], fig511_experimental[:, 1], color_shape_experimental, label=label_experimental)
ax0.plot(fig511_FRANC3D[:, 0], fig511_FRANC3D[:, 1], color_shape_franc3d, label=label_franc3d)
ax0.plot(fig511_ZFEM[:, 0], fig511_ZFEM[:, 1], color_shape_zfem, label=label_zfem)

# Configure axis labels, legend, and aspect ratio
ax0.set_xlabel("X (mm)")
ax0.set_ylabel("Y (mm)")
ax0.set_aspect('equal', adjustable='box')
ax0.legend()


###############################################################################
# Figure 5.12a: Crack Path Comparison for Specimen 17-155
# -------------------------------------------------------
# This plot reproduces Figure 5.12a, comparing the experimental crack path
# with results from both FRANC3D and ZFEM simulations for specimen 17-155.

# Load data from text files
fig512_17_155_frame = np.loadtxt("Wagner_phd/fig512a/reference.txt")
fig512_17_155_experimental = np.loadtxt("Wagner_phd/fig512a/experimental.txt")
fig512_17_155_FRANC3D = np.loadtxt("Wagner_phd/fig512a/FRANC3D.txt")
fig512_17_155_ZFEM = np.loadtxt("Wagner_phd/fig512a/ZFEM.txt")

# Create the plot
fig, ax0 = plt.subplots()
ax0.plot(fig512_17_155_frame[:, 0], fig512_17_155_frame[:, 1], 'k.', label="Reference")
ax0.plot(fig512_17_155_experimental[:, 0], fig512_17_155_experimental[:, 1], color_shape_experimental, label=label_experimental)
ax0.plot(fig512_17_155_FRANC3D[:, 0], fig512_17_155_FRANC3D[:, 1], color_shape_franc3d, label=label_franc3d)
ax0.plot(fig512_17_155_ZFEM[:, 0], fig512_17_155_ZFEM[:, 1], color_shape_zfem, label=label_zfem)

# Configure axis labels, legend, and aspect ratio
ax0.set_xlabel("X (mm)")
ax0.set_ylabel("Y (mm)")
ax0.set_aspect('equal', adjustable='box')
ax0.legend()

# Display all generated plots
plt.show()