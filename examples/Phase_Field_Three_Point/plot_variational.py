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

###############################################################################
# Import from phasefieldx package
# -------------------------------
from phasefieldx.Element.Phase_Field_Fracture.Input import Input
# from phasefieldx.Element.Phase_Field_Fracture.solver.solver_ener_variational import solve

def solve(Data,
          msh,
          final_gamma,
          V_u,
          V_Φ,
          bc_list_u=[],
          bc_list_phi=[],
          f_list_u=None,
          T_list_u=None,
          ds_bound=None,
          dt=0.0005,
          dt_min=1e-12,
          dt_max=0.1,
          path=None,
          bcs_list_u_names=None,
          c1=1.0,
          c2=1.0,
          threshold_gamma_save=None,
          continue_simulation=False,
          step_continue=10):
    return print("Check docs")


from phasefieldx.Boundary.boundary_conditions import bc_y, bc_x, get_ds_bound_from_marker
from phasefieldx.PostProcessing.ReferenceResult import AllResults

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
Data = Input(E=20.8,
             nu=0.3,
             Gc=0.0005,
             l=0.03,
             degradation="anisotropic",
             split_energy="spectral",
             degradation_function="quadratic",
             irreversibility="no",
             fatigue=False,
             fatigue_degradation_function="no",
             fatigue_val=0.0,
             k=0.0,
             save_solution_xdmf=False,
             save_solution_vtu=True,
             results_folder_name="results_variational")

# %%
# The variable `a0` defines the initial crack length in the mesh. This parameter
# is crucial for setting up the simulation, as it determines the starting point
# of the crack in the domain.
a0 = 0.2  # Initial crack length in the mesh
    
###############################################################################
# Mesh Definition
# ---------------
# The mesh is a structured grid with quadrilateral elements:
#
# - `divx`, `divy`: Number of elements along the x and y axes.
# - `lx`, `ly`: Physical domain dimensions in x and y.
divx, divy = 400, 100
lx, ly = 4.0, 1.0
h = ly / divy

c1 = 35
c2 = 1.0
    
run_simulation = False

if run_simulation:
    msh = dolfinx.mesh.create_rectangle(mpi4py.MPI.COMM_WORLD,
                                        [np.array([-lx, 0.0]),
                                        np.array([0.0, ly])],
                                        [divx, divy],
                                        cell_type=dolfinx.mesh.CellType.quadrilateral)

    fdim = msh.topology.dim - 1 # Dimension of the mesh facets

    ###############################################################################
    # Boundary Identification Functions
    # ---------------------------------
    # These functions identify points on the specific boundaries of the domain
    # where boundary conditions will be applied. The `bottom_left` function checks 
    # if a point lies on the bottom left boundary, returning `True` for points where 
    # `y=0` and `x` is less than `-lx + surface`, and `False` otherwise. Similarly, 
    # the `center` function identifies points on the center boundary, returning `True` 
    # for points where `x=0` and `y` is greater than or equal to `a0`. The `top` function 
    # identifies points on the top boundary, returning `True` for points where `y=ly`, 
    # and `x` is greater than or equal to `-surface`, and `False` otherwise.
    #
    # This approach ensures that boundary conditions are applied to specific parts of 
    # the mesh, which helps in defining the simulation's physical constraints.

    surface = 0.075
    def bottom_left(x):
        return np.logical_and(np.isclose(x[1], 0), np.less(x[0], -lx + surface))

    def center(x):
        return np.logical_and(np.isclose(x[0], 0), np.greater_equal(x[1], a0))

    def top(x):
        return np.logical_and(np.isclose(x[1], ly), np.greater_equal(x[0], -surface))


    # %%
    # Using the `bottom_left`, `center`, and `top` functions, we locate the facets on the respective boundaries of the mesh:
    # - `bottom_left`: Identifies facets on the bottom left side of the mesh where `y = 0` and `x < -lx + surface`.
    # - `center`: Identifies facets on the center boundary where `x = 0` and `y >= a0`.
    # - `top`: Identifies facets on the top side of the mesh where `y = ly` and `x >= -surface`.
    # The `locate_entities_boundary` function returns an array of facet indices representing these identified boundary entities.
    bottom_left_facet_marker = dolfinx.mesh.locate_entities_boundary(msh, fdim, bottom_left)
    center_facet_marker = dolfinx.mesh.locate_entities_boundary(msh, fdim, center)
    top_facet_marker = dolfinx.mesh.locate_entities_boundary(msh, fdim, top)


    # %%
    # The `get_ds_bound_from_marker` function generates a measure for applying boundary conditions 
    # specifically to the surface marker where the load will be applied, identified by `top_facet_marker`. 
    # This measure is then assigned to `ds_top`.
    ds_top = get_ds_bound_from_marker(top_facet_marker, msh, fdim)

    # %%
    # `ds_list` is an array that stores boundary condition measures along with names 
    # for each boundary, simplifying result-saving processes. Each entry in `ds_list` 
    # is formatted as `[ds_, "name"]`, where `ds_` represents the boundary condition measure, 
    # and `"name"` is a label used for saving. Here, `ds_bottom` and `ds_top` are labeled 
    # as `"bottom"` and `"top"`, respectively, to ensure clarity when saving results.
    ds_list = np.array([
                    [ds_top, "top"],
                    ])


    ###############################################################################
    # Function Space Definition
    # -------------------------
    # Define function spaces for displacement and phase-field using Lagrange elements.
    V_u = dolfinx.fem.functionspace(msh, ("Lagrange", 1, (msh.geometry.dim, )))
    V_phi = dolfinx.fem.functionspace(msh, ("Lagrange", 1))


    ###############################################################################
    # Boundary Conditions
    # -------------------
    # Dirichlet boundary conditions are defined as follows:
    #
    # - `bc_bottom_left`: Constrains both x and y displacements on the bottom left boundary, 
    #   ensuring that the leftmost bottom edge remains fixed.
    # - `bc_bottom_right`: Constrains only the vertical displacement (y-displacement) on the 
    #   bottom left boundary, while allowing horizontal movement.
    #
    # These boundary conditions ensure that the relevant portions of the mesh are correctly 
    # fixed or allowed to move according to the simulation requirements.

    bc_bottom_left = bc_y(bottom_left_facet_marker, V_u, fdim)
    bc_center = bc_x(center_facet_marker, V_u, fdim)

    # %%
    # The bcs_list_u variable is a list that stores all boundary conditions for the displacement
    # field $\boldsymbol u$. This list facilitates easy management of multiple boundary
    # conditions and can be expanded if additional conditions are needed.
    bcs_list_u = [bc_bottom_left, bc_center]
    bcs_list_u_names = ["bottom_left",  "center"]

    ###############################################################################
    # External Load Definition
    # ------------------------
    # Here, we define the external load to be applied to the top boundary (`ds_top`). 
    # `T_top` represents the external force applied in the y-direction.

    T_top = dolfinx.fem.Constant(msh, petsc4py.PETSc.ScalarType((0.0, -1.0/surface)))

    # %%
    # The load is added to the list of external loads, `T_list_u`, which will be updated
    # incrementally in the `update_loading` function.
    T_list_u = [
            [T_top, ds_top]
            ]
    f = None

    ###############################################################################
    # Boundary Conditions for phase field
    bcs_list_phi = []


    ###############################################################################
    # Solver Call for a Phase-Field Fracture Problem
    # ----------------------------------------------
    final_gamma = 0.5

    # %%
    # Uncomment the following lines to run the solver with the specified parameters.


    solve(Data,
        msh,
        final_gamma,
        V_u,
        V_phi,
        bcs_list_u,
        bcs_list_phi,
        f,
        T_list_u,
        ds_list,
        dt=0.001,
        dt_min=1e-12,
        dt_max=1.0,
        path=None,
        bcs_list_u_names=bcs_list_u_names,
        c1=c1,
        c2=c2,
        threshold_gamma_save=0.025,
        continue_simulation=False,
        step_continue=0)


##############################################################################
# Load results
# ------------
# Once the simulation finishes, the results are loaded from the results folder.
# The AllResults class takes the folder path as an argument and stores all
# the results, including logs, energy, convergence, and DOF files.
# Note that it is possible to load results from other results folders to compare results.
# It is also possible to define a custom label and color to automate plot labels.

import pyvista as pv
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.abspath('../../'))
plt.style.use('../../graph.mplstyle')
import plot_config as pcfg

S = AllResults(Data.results_folder_name)
S.set_label('Simulation')
S.set_color('b')

aux = 135
lambda_ = S.dof_files["lambda.dof"]["lambda"][:aux]


###############################################################################
# Simulation quantities
# ---------------------
force_half        = abs(S.reaction_files['bottom_left.reaction']["Ry"][:aux])
displacement_half = abs(2*S.energy_files['total.energy']["E"][:aux]/(S.reaction_files['bottom_left.reaction']["Ry"][:aux]))
stiffness_half    = abs(S.reaction_files['bottom_left.reaction']["Ry"][:aux]/displacement_half)
compliance_half   = 1/stiffness_half
dCda_half         = 2*Data.Gc/S.reaction_files['bottom_left.reaction']["Ry"][:aux]**2
gamma_half        = a0/2 + S.energy_files['total.energy']["gamma"][:aux]
lambda_half       = S.dof_files["lambda.dof"]["lambda"][:aux]


###############################################################################
# Complete model without corrections
# ----------------------------------
displacement_complete  = displacement_half
force_complete         = 2*force_half
compliance_complete    = 0.5*compliance_half
stiffness_complete     = 2*stiffness_half
dCda_complete          = 1.0/4.0*dCda_half
gamma_complete         = a0 + 2.0 * S.energy_files['total.energy']["gamma"][:aux]
gamma_phi_complete     = a0 + 2.0 * S.energy_files['total.energy']["gamma_phi"][:aux]
gamma_gradphi_complete = a0 + 2.0 * S.energy_files['total.energy']["gamma_gradphi"][:aux]
lambda_complete        = lambda_half

# header = ["displacement", "force", "gamma", "compliance", "stiffness", "dCda", "lambda"]
# data_save = np.column_stack((displacement_complete, force_complete, gamma_complete, compliance_complete,stiffness_complete, dCda_complete, lambda_complete ))
# save_path = os.path.join(Data.results_folder_name, "results.pff")
# np.savetxt(save_path, data_save, fmt="%.6e", delimiter="\t", header="\t".join(header), comments="")


###############################################################################
# Complete model equivalent Gc
# -----------------------------
alpha = 1/np.sqrt(1+c1*S.dof_files["lambda.dof"]["lambda"][:aux]/Data.Gc)

displacement_complete_equivalent = displacement_complete*alpha
force_complete_equivalent        = force_complete*alpha
compliance_complete_equivalent   = compliance_complete 
stiffness_complete_equivalent    = stiffness_complete  
dCda_complete_equivalent         = dCda_complete
gamma_complete_equivalent        = a0 + 2.0 * S.energy_files['total.energy']["gamma"][:aux]
gamma_phi_complete_equivalent    = a0 + 2.0 * S.energy_files['total.energy']["gamma_phi"][:aux]
gamma_gradphi_complete_equivalent = a0 + 2.0 * S.energy_files['total.energy']["gamma_gradphi"][:aux]
lambda_complete_equivalent        = lambda_complete

header = ["displacement", "force", "gamma", "compliance", "stiffness", "dCda", "lambda"]
data_save = np.column_stack((displacement_complete_equivalent, 
                             force_complete_equivalent,
                             gamma_complete_equivalent,
                             compliance_complete_equivalent,
                             stiffness_complete_equivalent,
                             dCda_complete_equivalent,
                             lambda_complete_equivalent))

save_path = os.path.join(Data.results_folder_name, "results_equivalent.pff")
np.savetxt(save_path, data_save, fmt="%.6e", delimiter="\t", header="\t".join(header), comments="")


###############################################################################
# Complete model with Gc corrections
# ----------------------------------
gc_factor = 1 + 2*h/(2*Data.l)
displacement_complete_corrected_gc  = displacement_complete_equivalent/np.sqrt(gc_factor)
force_complete_corrected_gc         = force_complete_equivalent/np.sqrt(gc_factor)
compliance_complete_corrected_gc    = compliance_complete_equivalent
stiffness_complete_corrected_gc     = stiffness_complete_equivalent
dCda_complete_corrected_gc          = dCda_complete_equivalent*gc_factor
gamma_complete_corrected_gc         = a0 + 2.0 * S.energy_files['total.energy']["gamma"][:aux]/gc_factor
gamma_phi_complete_corrected_gc     = a0 + 2.0 * S.energy_files['total.energy']["gamma_phi"][:aux]/gc_factor
gamma_gradphi_complete_corrected_gc = a0 + 2.0 * S.energy_files['total.energy']["gamma_gradphi"][:aux]/gc_factor
lambda_complete_corrected_gc        = lambda_complete_equivalent

# header = ["displacement", "force", "gamma", "compliance", "stiffness", "dCda", "lambda"]
# data_save = np.column_stack((displacement_complete_corrected_gc, force_complete_corrected_gc, gamma_complete_corrected_gc, compliance_complete_corrected_gc,stiffness_complete_corrected_gc, dCda_complete_corrected_gc, lambda_complete_corrected_gc))
# save_path = os.path.join(Data.results_folder_name, "results_equivalent_bourdin.pff")
# np.savetxt(save_path, data_save, fmt="%.6e", delimiter="\t", header="\t".join(header), comments="")

# Load as pandas DataFrame using header from the first line of the file
# data = pd.read_csv(save_path, delimiter="\t", comment="#", header=0)
# Example access: data["gamma"]


dCda_numerical = np.gradient(compliance_complete, gamma_complete)

###############################################################################
# Plot: Phase-Field
# -----------------
file_vtu = pv.read(os.path.join(Data.results_folder_name, "paraview-solutions_vtu", "phasefieldx_p0_000016.vtu"))
pv.start_xvfb()
file_vtu.plot(scalars='phi', cpos='xy', show_scalar_bar=True, show_edges=False)

###############################################################################
# Plot: Mesh
# ----------
file_vtu.plot(cpos='xy', color='white', show_edges=True)


###############################################################################
# Plot: Displacement vs Gamma Energy
# ----------------------------------.
fig, energyg = plt.subplots()
energyg.plot(displacement_complete, gamma_complete , 'b-')

energyg.set_xlabel('displacement - u $[mm]$')
energyg.set_ylabel('Energy')
energyg.legend()


###############################################################################
# Plot: Force vs Vertical Displacement
# ------------------------------------
fig, ax_reaction = plt.subplots()

ax_reaction.plot(displacement_complete, force_complete, 'k-')

ax_reaction.set_xlabel(pcfg.displacement_label)
ax_reaction.set_ylabel(pcfg.force_label)
ax_reaction.legend()


###############################################################################
# Plot: Force vs Vertical Displacement
# ------------------------------------
fig, ax_reaction = plt.subplots()

ax_reaction.plot(displacement_complete, lambda_, 'k-')

ax_reaction.set_xlabel(pcfg.displacement_label)
ax_reaction.set_ylabel(pcfg.lambda_label)
ax_reaction.legend()


###############################################################################
# Plot: Gamma vs stiffness
# ------------------------
fig, ax_reaction = plt.subplots()

ax_reaction.plot(gamma_complete, stiffness_complete, 'k-')

ax_reaction.set_xlabel(pcfg.gamma_label)
ax_reaction.set_ylabel(pcfg.stiffness_label)
ax_reaction.legend()

plt.show()
