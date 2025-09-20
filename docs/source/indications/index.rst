Indications
===========

The analysis presented in this article is based on four main components:

* **Theoretical Analysis:** This component focuses on the theoretical framework, using Linear Elastic Fracture Mechanics (LEFM) and Paris' law. The results are saved in files with the `.lefm` extension for comparison with other methods.

* **Elastic Simulations:** This component includes simulations of linear elastic problems. The results are saved in files with the `.elastic` extension for comparison with other methods.

* **Phase-Field Simulations:** This component involves simulations using phase-field energy-controlled schemes. The results are saved in files with the `.pff` extension for comparison with other methods.

* **Reference Data:** This component consists of additional data from articles or other references. This data is saved in files with the `.databib` extension.

For all the simulations and calculations the following units system is considered:

+-------------------------------------+---------------------------------------------+
| **Quantity**                        | **Unit**                                    |
+=====================================+=============================================+
| Length                              | $mm$                                        |
+-------------------------------------+---------------------------------------------+
| Force: $P$                          | $kN$                                        |
+-------------------------------------+---------------------------------------------+
| Area                                | $mm^2$                                      |
+-------------------------------------+---------------------------------------------+
| Energy                              | $kN \cdot mm$                               |
+-------------------------------------+---------------------------------------------+
| Young's modulus: $E$                | $kN/mm^2$                                   |
+-------------------------------------+---------------------------------------------+
| Poisson's ratio  $\nu$              | dimensionless                               |
+-------------------------------------+---------------------------------------------+
| Energy release rate: $G$            | $kN/mm$                                     |
+-------------------------------------+---------------------------------------------+
| Crack length                        | $mm$                                        |
+-------------------------------------+---------------------------------------------+
| Crack area                          | $mm^2$                                      |
+-------------------------------------+---------------------------------------------+
| Stiffness  $K$                      | $kN/mm$                                     |
+-------------------------------------+---------------------------------------------+
| Compliance  $C$                     | $mm/kN$                                     |
+-------------------------------------+---------------------------------------------+
| Paris constant $n$                  | dimensionless                               |
+-------------------------------------+---------------------------------------------+
| Paris constant $C_\text{Paris}$     | $mm^{(2+3n)/2} / \text{cycle} \cdot {kN}^n$ |
+-------------------------------------+---------------------------------------------+

All 2D simulations are performed using quadrilateral elements under plane strain conditions.
