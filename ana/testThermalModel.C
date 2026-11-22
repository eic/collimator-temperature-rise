/*
 * testThermalModel.C
 *
 * Self-test of the thermal model in ThermalModel.h.  The numerical
 * DeltaT(E_dep) conversion is compared with the closed-form enthalpy of the
 * NIST Shomate parameterization, and basic properties (continuity,
 * monotonicity, melting plateau, threshold inversion) are checked.
 *
 * Usage (from any directory):
 *   root -l -b -q '../ana/testThermalModel.C'
 * Exits with status 1 if any check fails.
 */
#include "ThermalModel.h"

#include <TSystem.h>

#include <cmath>
#include <cstdio>

namespace
{
int gFailures = 0;

void check(bool ok, const char* what)
{
    std::printf("  [%s] %s\n", ok ? " OK " : "FAIL", what);
    if (!ok) ++gFailures;
}

/// Antiderivative of the Shomate heat capacity with respect to t = T/1000:
/// F(t) = A t + B t^2/2 + C t^3/3 + D t^4/4 - E/t; int cp dT = 1000 [F(t2) - F(t1)].
double shomateIntegral(const std::array<double, 5>& p, double t)
{
    return p[0] * t + p[1] * t * t / 2 + p[2] * t * t * t / 3 + p[3] * t * t * t * t / 4 - p[4] / t;
}

/// Exact energy density [J/cm^3] needed to heat aluminum from T0 to T.
double exactEnergyDensityAl(double T)
{
    const std::array<double, 5> solid  = {28.08920, -5.414849, 8.560423, 3.427370, -0.277375};
    const std::array<double, 5> liquid = {31.75104, 3.935826e-8, -1.786515e-8, 2.694171e-9, 5.480037e-9};
    const double perMol = Al.rho / kAlMolarMass * 1000.0;   // [mol/cm^3] x (dT/dt = 1000 K)
    const double t0 = T0 / 1000.0, tm = Al.Tm / 1000.0, t = T / 1000.0;
    if (T <= Al.Tm)
        return perMol * (shomateIntegral(solid, t) - shomateIntegral(solid, t0));
    return perMol * (shomateIntegral(solid, tm) - shomateIntegral(solid, t0)) + Al.Lf +
           perMol * (shomateIntegral(liquid, t) - shomateIntegral(liquid, tm));
}
} // namespace

void testThermalModel()
{
    std::printf("Thermal-model self-test (Al)\n");

    // 1. Comparison with the closed-form enthalpy.
    double maxDiff = 0.0;
    for (double T : {300.5, 310.0, 374.0, 500.0, 700.0, 900.0, 933.9, 934.5, 1200.0, 2000.0, 2700.0})
    {
        const double E  = exactEnergyDensityAl(T);
        const double dT = computeDeltaT(E, Al, cp_al).deltaT;
        maxDiff = std::max(maxDiff, std::abs(dT - (T - T0)));
        std::printf("    T = %7.1f K : E = %9.3f J/cm^3 -> DeltaT = %9.4f K (exact %9.4f K)\n",
                    T, E, dT, T - T0);
    }
    check(maxDiff < 0.01, Form("DeltaT agrees with the closed-form enthalpy (max |diff| = %.2g K)", maxDiff));

    // 2. Small energy densities: no staircase / minimum step.
    bool smallOk = true;
    for (double E : {1e-6, 1e-3, 0.05, 0.15, 1.0})
    {
        const double expected = E / (Al.rho * cp_al(T0));
        smallOk &= std::abs(computeDeltaT(E, Al, cp_al).deltaT / expected - 1.0) < 1e-3;
    }
    check(smallOk, "DeltaT(E) = E / (rho cp(T0)) for small E (no discretization staircase)");
    check(computeDeltaT(0.0, Al, cp_al).deltaT == 0.0, "DeltaT(0) = 0");

    // 3. Monotonic and continuous on a fine grid across the melting plateau.
    bool monotonic = true;
    double maxJump = 0.0, prev = 0.0;
    for (double E = 0.0; E <= 4000.0; E += 0.37)
    {
        const double dT = computeDeltaT(E, Al, cp_al).deltaT;
        monotonic &= (dT >= prev);
        maxJump = std::max(maxJump, dT - prev);
        prev = dT;
    }
    check(monotonic, "DeltaT(E) is monotonically non-decreasing");
    check(maxJump < 0.37 / (Al.rho * 0.8), Form("DeltaT(E) is continuous (max step %.3f K for dE = 0.37 J/cm^3)", maxJump));

    // 4. Melting plateau.
    const double eOnset = exactEnergyDensityAl(Al.Tm);
    const auto mid = computeDeltaT(eOnset + 0.5 * Al.Lf, Al, cp_al);
    check(mid.deltaT == Al.Tm - T0 && mid.melted, "partially molten cell: DeltaT = Tm - T0, melted = true");
    check(!computeDeltaT(eOnset - 1.0, Al, cp_al).melted, "solid cell below the plateau: melted = false");

    // 5. Threshold inversion used for the shock / melt volumes.
    const double eShock      = energyDensityForDeltaT(Al.dT_shock, Al, cp_al);
    const double eMeltOnset  = energyDensityForMeltingOnset(Al, cp_al);
    const double eMeltFull   = energyDensityForDeltaT(Al.Tm - T0, Al, cp_al);
    std::printf("    E(DeltaT_shock = %.0f K) = %.4f J/cm^3 (exact %.4f)\n",
                Al.dT_shock, eShock, exactEnergyDensityAl(T0 + Al.dT_shock));
    std::printf("    melting starts at %.3f J/cm^3 (exact %.3f); completely molten at %.3f J/cm^3 (exact %.3f)\n",
                eMeltOnset, eOnset, eMeltFull, eOnset + Al.Lf);
    check(computeDeltaT(eShock * (1 + 1e-9), Al, cp_al).deltaT > Al.dT_shock &&
          computeDeltaT(eShock * (1 - 1e-6), Al, cp_al).deltaT <= Al.dT_shock,
          "shock threshold inversion is consistent with computeDeltaT()");
    check(std::abs(eShock / exactEnergyDensityAl(T0 + Al.dT_shock) - 1.0) < 1e-4,
          "shock threshold agrees with the closed-form enthalpy");
    check(computeDeltaT(eMeltOnset * (1 + 1e-9), Al, cp_al).melted &&
          !computeDeltaT(eMeltOnset * (1 - 1e-6), Al, cp_al).melted,
          "melting-onset threshold is consistent with computeDeltaT()");
    check(std::abs(eMeltOnset / eOnset - 1.0) < 1e-4,
          "melting-onset threshold agrees with the closed-form enthalpy E(Tm)");
    check(std::abs((eMeltFull - eMeltOnset) / Al.Lf - 1.0) < 1e-6,
          "complete melting = melting onset + Lf");

    std::printf("%s (%d failure%s)\n", gFailures ? "FAILED" : "ALL CHECKS PASSED", gFailures,
                gFailures == 1 ? "" : "s");
    if (gFailures) gSystem->Exit(1);
}
