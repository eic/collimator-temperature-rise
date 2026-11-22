/*
 * getDeltaTempShock.C
 *
 * Elastic estimate of the thermal-shock temperature-rise threshold for the
 * aluminum demonstration:
 *
 *   DeltaT_elastic = sigma_UTS * (1 - nu) / (E * alpha),
 *
 * i.e. the temperature rise at which the thermal stress E alpha DeltaT / (1 - nu)
 * reaches the ultimate tensile strength sigma_UTS.
 *
 * The macro also prints the practical threshold that drawScoring.C actually
 * uses (Al.dT_shock in ThermalModel.h) and its ratio to the elastic estimate.
 *
 * Usage (from any directory):
 *   root -l -b -q '../ana/getDeltaTempShock.C'
 */
#include "ThermalModel.h"

#include <iostream>

void getDeltaTempShock()
{
    // Young's modulus (modulus of elasticity) [GPa]
    // Aluminum: https://www.matweb.com/search/DataSheet.aspx?MatGUID=0cd1edf33ac145ee93a0aa6fc666c0e0
    const double E = 68.0;
    // Coefficient of thermal expansion [um/m/K]
    const double alpha = 24.0;
    // Poisson's ratio
    // Al 1100: https://alloysintl.com/inventory/aluminum-alloys-supplier-old/aluminum-1100/
    const double nu = 0.330;
    // Ultimate tensile strength [MPa]
    // Al 1100: https://alloysintl.com/inventory/aluminum-alloys-supplier-old/aluminum-1100/
    const double sigma_ultimate = 90.0;

    // SI units: MPa -> Pa (1e6), GPa -> Pa (1e9), um/m/K -> 1/K (1e-6).
    const double dT_elastic =
        (sigma_ultimate * 1e6) * (1.0 - nu) / ((E * 1e9) * (alpha * 1e-6));   // [K]

    std::cout << "Al : E = " << E << " [GPa] | alpha = " << alpha << " [um/m/K] | nu = "
              << nu << " | sigma_UTS = " << sigma_ultimate << " [MPa]" << std::endl;
    std::cout << "Elastic estimate sigma_UTS (1 - nu) / (E alpha) = " << dT_elastic
              << " [K]" << std::endl;
    std::cout << "Practical threshold used by drawScoring.C (Al.dT_shock) = "
              << Al.dT_shock << " [K] = " << Al.dT_shock / dT_elastic
              << " x elastic estimate" << std::endl;
}
