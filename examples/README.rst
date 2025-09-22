Examples
========

This directory contains the examples supporting the article, demonstrating various approaches to analyzing crack behavior and material fatigue. The analyses are structured into three main types:

*   **Theoretical LEFM Solutions**: Analytical solutions derived from Linear Elastic Fracture Mechanics (LEFM) provide a theoretical baseline for comparison.
*   **Elasticity Simulations**: Numerical analyses are performed on meshes with predefined cracks to calculate the specimen's compliance and stiffness as a function of crack length.
*   **Phase-Field Simulations**: The proposed phase-field models are used to simulate crack propagation without requiring predefined crack paths.

The examples are organized into the following directories:

1.  **LEFM Analysis**: Presents the theoretical analysis based on LEFM. :ref:`ref_examples_LEFM`
2.  **Elasticity Simulations**: Contains elasticity simulations with imposed cracks. :ref:`ref_examples_fem_elasticity`
3.  **Phase-Field Simulations**: Includes simulations for three different specimen geometries:
   
    *   Three-point bending :ref:`ref_examples_phase_field_three_point`
    *   Center-cracked specimen :ref:`ref_examples_phase_field_central_crack`
    *   Compact specimen :ref:`ref_examples_phase_field_compact_specimen`
  
4.  **Papers Data**: Provides data from published papers for comparison purposes. :ref:`ref_examples_papers_data`
5.  **Comparison**: Combines results from the different methods for direct comparison.
   
    *   Central cracked specimen :ref:`ref_examples_compare_central_cracked`
    *   Compact specimen :ref:`ref_examples_compare_compact_specimen`
    *   Three-point bending :ref:`ref_examples_compare_three_point`

6.  **GmshGeoFiles**: Contains the geometry files and instructions for generating the meshes used in the simulations. :ref:`ref_examples_geo_files`

Each ``Comparison`` directory corresponds to a validation or results section in the paper, ensuring that the presented graphs align with the discussed findings.
