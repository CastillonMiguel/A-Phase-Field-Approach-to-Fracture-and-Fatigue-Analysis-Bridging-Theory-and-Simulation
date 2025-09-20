# filepath: ArticleFatigue/examples/plot_config.py
"""
Central configuration for plot labels and styles.
Import this module in your plotting scripts to ensure consistency.
"""

# General Plotting Labels
crack_length_label   = r"Crack Length (a) [mm]"
stiffness_label      = r"Stiffness (K) [kN/mm]"
compliance_label     = r"Compliance (C) [mm/kN]"
dCda_label           = r"dC/da [(mm/kN)]"
dadN_label           = r"da/dN [mm/cycles]"
dCda_1_label         = r"1/(dC/da) [1/((mm/kN)/mm)]"
critical_force_label = r"Critical Force ($P_c$) [kN]"
displacement_label   = r"Displacement (u) [mm]"
force_label          = r"Force (P) [kN]"
critical_force_label = r"Critical Force ($P_c$) [kN]"
strain_energy_label  = r"Strain Energy [kN mm]"
cycles_label         = r"Fatigue Cycles (N)"
DeltaK_label         = r"Stress Intensity Factor Range ($\Delta K$)"
gamma_label          = r"Gamma [$mm^2$]"
iterations_label     = r"Iterations [-]"
G_effective_label    = r"Effective Fracture Energy ($G_{eff}$) [kN/mm]"
lambda_label         = r"$\lambda$ [kN/mm]"

gamma_ref_label      = r"Gamma"
gamma_bourdin_label  = r"Gamma Bourdin"
gamma_geometry_label  = r"Gamma Geometry"

# Color definitions (scientific publication style)
color_orangered = '#D55E00'   # Scientific orange-red
color_blue = '#0173B2'        # Scientific blue
color_gold = '#DE8F05'        # Scientific gold
color_green = '#029E73'       # Scientific green
color_purple = '#CC78BC'      # Scientific purple
color_brown = '#CA9161'       # Scientific brown
color_pink = '#FBAFE4'        # Scientific pink
color_grey = '#949494'        # Scientific grey
color_yellow = '#ECE133'      # Scientific yellow
color_lightblue = '#56B4E9'   # Scientific light blue
color_black = '#000000'       # Scientific black


# Example usage for plotting:
# colors = [color_blue, color_orangered, color_gold, color_green, color_purple, color_brown, color_pink, color_grey, color_yellow, color_lightblue]
# ax.plot(x, y, color=colors[i], label=labels[i])