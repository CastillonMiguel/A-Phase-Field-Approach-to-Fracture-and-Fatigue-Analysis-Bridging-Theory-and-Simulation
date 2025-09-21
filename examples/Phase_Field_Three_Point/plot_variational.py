r"""
.. _ref_phase_field_three_point_variational:

Variational Approach
--------------------

The three-point bending specimen consists of a rectangular plate with a centrally located notch, supported at both ends. To enhance computational efficiency, a symmetric half-model is employed. The boundary conditions include a fixed vertical displacement over a small area (Asurface) at the lower-left support and an applied downward vertical force over an equal area at the top center. Additionally, the symmetry plane is constrained horizontally.

.. note::
   In this case, only half part of the model is considered due to symmetry. Furthermore, a regular mesh is utilized.

.. code-block::
      
   #                             ||||      
   #                             \/\/               
   #            ------------------------------------* 
   #            |                                   | 
   #            |                                   |
   #            |                                   |
   #            |                 /\                | 
   #            *-----------------  ----------------*
   #            /_\/_\        (0,0,0)          /_\/_\       
   #    |Y     ///////                         oo  oo
   #    |
   #    ---X
   # Z /

.. code-block::
      
   #                             ||     
   #                             \/               
   #            ------------------* o|/
   #            |                 | o|/
   #            |                 | o|/
   #            |               _ | o|/
   #            |              a0 | 
   #            *-----------------*
   #           /_\/_\  
   #    |Y     oo  oo
   #    |
   #    ---X
   # Z /
   

The Young's modulus, Poisson's ratio, and the critical energy release rate are provided in the table :ref:`Properties <table_properties_label>`. Young's modulus $E$ and Poisson's ratio $\nu$ can be expressed using the Lamé parameters as: $\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}$; $\mu=\frac{E}{2(1+\nu)}$.

.. _table_properties_label:

+----+---------+--------+
|    | VALUE   | UNITS  |
+====+=========+========+
| E  | 20.8    | kN/mm2 |
+----+---------+--------+
| nu | 0.3     | [-]    |
+----+---------+--------+
| Gc | 0.0005  | kN/mm  |
+----+---------+--------+
| l  | 0.03    | mm     |
+----+---------+--------+

"""

###############################################################################
# Import necessary libraries
# --------------------------
import numpy as np
import dolfinx
import mpi4py
import petsc4py
import os

import pyvista as pv
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath('../../'))
plt.style.use('../../graph.mplstyle')
import plot_config as pcfg

###############################################################################
# Import from phasefieldx package
# -------------------------------
# from phasefieldx.Element.Phase_Field_Fracture.Input import Input
# # from phasefieldx.Element.Phase_Field_Fracture.solver.solver_ener_variational import solve
# from phasefieldx.Boundary.boundary_conditions import bc_y, bc_x, get_ds_bound_from_marker
# from phasefieldx.PostProcessing.ReferenceResult import AllResults


###############################################################################
# Parameters Definition
# ---------------------
# `Data` is an input object containing essential parameters for simulation setup
# and result storage:
#
# - `E`: Young's modulus, set to 20.8 $kN/mm^2$.
# - `nu`: Poisson's ratio, set to 0.3.
# - `Gc`: Critical energy release rate, set to 0.0005 $kN/mm$.
# - `l`: Length scale parameter, set to 0.03 $mm$.
# - `degradation`: Specifies the degradation type. Options are "isotropic" or "anisotropic".
# - `split_energy`: Controls how the energy is split; options include "no" (default), "spectral," or "deviatoric."
# - `degradation_function`: Specifies the degradation function; here, it is "quadratic."
# - `irreversibility`: Not used/implemented for this solver.
# - `save_solution_xdmf` and `save_solution_vtu`: Specify the file formats to save displacement results.
#   In this case, results are saved as `.vtu` files.
# - `results_folder_name`: Name of the folder for saving results. If it exists,
#   it will be replaced with a new empty folder.
# Data = Input(E=20.8,
#              nu=0.3,
#              Gc=0.0005,
#              l=0.03,
#              degradation="anisotropic",
#              split_energy="spectral",
#              degradation_function="quadratic",
#              irreversibility="no",
#              fatigue=False,
#              fatigue_degradation_function="no",
#              fatigue_val=0.0,
#              k=0.0,
#              save_solution_xdmf=False,
#              save_solution_vtu=True,
#              results_folder_name="results_variational")

# %%
# The variable `a0` defines the initial crack length in the mesh. This parameter
# is crucial for setting up the simulation, as it determines the starting point
# of the crack in the domain.
a0 = 0.2  # Initial crack length in the mesh


###############################################################################
# Plot: Phase-Field
# -----------------
file_vtu = pv.read(os.path.join("results_variational", "paraview-solutions_vtu", "phasefieldx_p0_000016.vtu"))
pv.start_xvfb()
file_vtu.plot(scalars='phi', cpos='xy', show_scalar_bar=True, show_edges=False)

###############################################################################
# Plot: Mesh
# ----------
file_vtu.plot(cpos='xy', color='white', show_edges=True)

plt.show()
