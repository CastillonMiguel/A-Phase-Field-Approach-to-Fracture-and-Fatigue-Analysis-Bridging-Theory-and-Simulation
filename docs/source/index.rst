.. raw:: html

    <style type="text/css">
        .big-font {font-size: var(--pst-font-size-h2); color: var(--pst-color-primary)}
    </style>


.. .. line-break::
   
.. rst-class:: font-weight-bold big-font

        A Phase-Field Approach to Fracture and Fatigue Analysis: Bridging Theory and Simulation


This repository provides all simulation files, scripts, and supplementary data referenced in the paper :footcite:t:`Castillon2025_arxiv`.

.. raw:: html

    <div style="background: white; padding: 24px; border-radius: 12px; display: flex; gap: 24px; justify-content: center;">
        <img src="_static/SIMULATION_2.gif" width="200px" loop="infinite" autoplay />
        <img src="_static/SIMULATION_1.gif" width="200px" loop="infinite" autoplay />
        <img src="_static/SIMULATION_4.gif" width="200px" loop="infinite" autoplay />
    </div>

This article presents a novel, robust and efficient framework for fatigue crack-propagation that combines the principles of Linear Elastic Fracture Mechanics (LEFM) with phase-field fracture (PFF). Contrary to cycle-by-cycle PFF approaches, this work relies on a single simulation and uses standard crack propagation models such as Paris' law for the material response, simplifying its parametrization.
The core of the methodology is the numerical evaluation of the derivative of a specimen's compliance with respect to the crack area. To retrieve this compliance the framework relies on a PFF-FEM simulation, controlled imposing a monotonic crack growth. This control of the loading process is done by a new crack-control scheme which allows to robustly trace the complete equilibrium path of a crack, capturing complex instabilities. The specimen's compliance obtained from the PFF simulation enables the integration of Paris' law to predict fatigue life.

The proposed methodology is first validated through a series of benchmarks with analytical solutions to demonstrate its accuracy. The framework is then applied to more complex geometries where the crack path is unknown, showing a very good agreement with experimental results of both crack paths and fatigue life.

.. highlights::

   - **Develops a phase-field framework** for fracture and fatigue analysis.
   - New Energy-Controlled Solvers: It introduces two new, efficient solvers (variational and non-variational) that can simulate the entire fracture process, including complex behaviors like snap-back, in a single run.
   - Efficient Fatigue Analysis: The framework uses the results from a single quasi-static simulation to predict fatigue life based on Paris' law, avoiding the need for computationally expensive cycle-by-cycle simulations.
   - Improved Accuracy: It includes a validated method to correct for the overestimation of crack length that is common in phase-field models, leading to results that align well with established theoretical solutions (LEFM).
   - Complex Geometries: The method is successfully applied to complex specimens with holes where the crack path is unknown beforehand, demonstrating its ability to predict intricate crack trajectories that match experimental data.
   - Open-Source and Reproducible: The entire toolchain, including the PhaseFieldX library and all simulation scripts, is open-source to ensure full reproducibility of the findings.

This repository is designed to ensure complete reproducibility of the results by providing all simulation data, parameter sets, meshes, and detailed numerical configurations.

.. code:: latex

    @misc{castillon2025,
        title={A Phase-Field Approach to Fracture and Fatigue Analysis: Bridging Theory and Simulation}, 
        author={M. Castillón and I. Romero and J. Segurado},
        year={2025},
        eprint={2509.08939},
        archivePrefix={arXiv},
        primaryClass={cond-mat.mtrl-sci},
        url={https://arxiv.org/abs/2509.08939}, 
    }

All the files are provided in the following `GitHub Repository <https://github.com/CastillonMiguel/A-Phase-Field-Approach-to-Fracture-and-Fatigue-Analysis-Bridging-Theory-and-Simulation>`_

Since the simulations were conducted using the open-source **PhaseFieldX** :footcite:t:`code_phasefieldx` library, the implementation details of the models can be found in the **PhaseFieldX** documentation and source code.

- GitHub Repository: `https://github.com/CastillonMiguel/phasefieldx <https://github.com/CastillonMiguel/phasefieldx>`_
- Documentation: `https://phasefieldx.readthedocs.io <https://phasefieldx.readthedocs.io>`_
- Paper: `https://doi.org/10.21105/joss.07307 <https://doi.org/10.21105/joss.07307>`_

The **PhaseFieldX** project is designed to simulate and analyze material behavior using phase-field models, which provide a continuous approximation of interfaces, phase boundaries, and discontinuities such as cracks. Leveraging the robust capabilities of *FEniCSx*, a renowned finite element framework for solving partial differential equations, this project facilitates efficient and precise numerical simulations. It supports a wide range of applications, including phase-field fracture, solidification, and other complex material phenomena, making it an invaluable resource for researchers and engineers in materials science.

.. image:: _static/logo_phasefieldx.png
   :width: 600px
   :align: center

.. footbibliography::

.. toctree::

   indications/index
   crack_measurements/index
   auto_examples/index
   references/index.rst
   acknowledgements/index.rst
