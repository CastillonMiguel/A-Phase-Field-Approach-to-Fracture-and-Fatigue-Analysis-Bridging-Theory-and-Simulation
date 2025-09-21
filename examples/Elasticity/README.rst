.. _ref_examples_fem_elasticity:

FEM: Elasticity Models
======================

This section presents a set of finite element simulations based on linear elastic fracture mechanics (LEFM). The examples include:

- Center-cracked tension tests (analyzed under both displacement-controlled and force-controlled loading)
- Compact tension specimen (analyzed under force-controlled loading)

For these geometries, the crack path and associated geometry factors are well-established, allowing direct comparison between analytical solutions and numerical results. By explicitly modeling cracks in the mesh and simulating various crack lengths, the computed compliance can be quantitatively compared to LEFM predictions.

Further details on the elasticity model and solver implementation are available in the `PhaseFieldX documentation <https://phasefieldx.readthedocs.io/en/latest/theory/elasticity/main.html>`_ and the PhaseFieldX library.
