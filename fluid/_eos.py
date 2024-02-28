"""
Van der Waals Equation of State:
Proposed by Johannes Diderik van der Waals in 1873, it corrects for
the volume occupied by gas molecules and accounts for
intermolecular attractive forces.

Peng-Robinson Equation of State:
Developed by Dukhin, Glandt, and V. V. Kafarov in the late 1960s,
it is an improvement over the Van der Waals equation and incorporates
the concept of acentric factor to better model non-ideal behavior.

Soave-Redlich-Kwong Equation of State (SRK EOS):
Proposed by Raul Soave in 1972, it is another modification of the
Van der Waals equation and includes an additional parameter
to improve accuracy in the prediction of vapor-liquid equilibria.

Redlich-Kwong Equation of State (RK EOS):
Proposed by Otto Redlich and J. N. S. Kwong in 1949, it's an early
attempt to improve upon the Van der Waals equation by introducing
temperature-dependent parameters.

Cubic Equations of State (e.g., Van der Waals, Redlich-Kwong, Peng-Robinson):
This class of equations includes several EOS that assume the volume occupied by
gas molecules can be represented by a cubic polynomial.

Lee-Kesler Equation of State:
Developed by John D. Lee and Michael Kesler in 1975, it is based on a simple
mixing rule and has been widely used in the petroleum industry.

Benedict-Webb-Rubin Equation of State (BWRS EOS):
Proposed by M. Benedict, G. B. Webb, and L. C. Rubin in 1940, it's one of the
earliest attempts to improve the accuracy of EOS by incorporating
temperature-dependent parameters.

Dranchuk-Abou-Kassem Equation of State (DAK EOS):
Developed by Michael D. Dranchuk and H. C. Abou-Kassem in 1975,
it's an extension of the Peng-Robinson equation with additional
terms to improve accuracy, particularly in the near-critical region.
"""