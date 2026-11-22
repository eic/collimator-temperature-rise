/*
 * ThermalModel.h
 *
 * Thermal model shared by the ROOT analysis macros of the collimator
 * temperature-rise example (drawScoring.C, getCp.C, getDeltaTempShock.C,
 * testThermalModel.C).  Keeping it in one place guarantees that all macros
 * use identical material data.
 *
 * Contents
 *   T0                       initial (ambient) temperature
 *   MaterialProps / Al       material constants used by the analysis
 *   cp_al(T)                 temperature-dependent specific heat of aluminum
 *   computeDeltaT()          deposited energy density -> temperature rise
 *   energyDensityForDeltaT() inverse of computeDeltaT() for thresholds
 *   energyDensityForMeltingOnset()  energy density at which melting starts
 *   getMaterialModel()       look up a material by name ("Al")
 *
 * Units used throughout
 *   temperature [K], specific heat [J/g/K], density [g/cm^3],
 *   energy density [J/cm^3].
 *
 * The macros include this file with  #include "ThermalModel.h".  ROOT
 * resolves the path relative to the including macro, so the macros can be
 * executed from any working directory.
 */
#ifndef COLLIMATOR_THERMAL_MODEL_H
#define COLLIMATOR_THERMAL_MODEL_H

#include <TString.h>

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <limits>

// ---------------------------------------------------------------------------
// Global settings
// ---------------------------------------------------------------------------

/// Initial (ambient) temperature of the collimator [K].
const double T0 = 300.0;

/// Default energy-density step used to integrate dE = rho * cp(T) dT [J/cm^3].
/// 0.1 J/cm^3 corresponds to a temperature step of ~0.04 K in aluminum.
const double kDefaultEnergyStep = 0.1;

// ---------------------------------------------------------------------------
// Material description
// ---------------------------------------------------------------------------

/// Material constants used by the thermal analysis.
struct MaterialProps
{
    double rho;       ///< density [g/cm^3]
    double Tm;        ///< melting (or sublimation) temperature [K]
    double Lf;        ///< latent heat of fusion per unit volume [J/cm^3]; 0 = no liquid phase
    double dT_shock;  ///< thermal-shock threshold on the temperature rise [K]
    double k;         ///< thermal conductivity [W/m/K]; used only for the diffusion length
};

/// Molar heat capacity from the NIST Shomate equation [J/mol/K]:
///   cp = A + B t + C t^2 + D t^3 + E / t^2,   with t = T[K] / 1000.
inline double cp_mol(const std::array<double, 5>& par, double T)
{
    const double t  = T / 1000.0;
    const double t2 = t * t;
    return par[0] + par[1] * t + par[2] * t2 + par[3] * t2 * t + par[4] / t2;
}

// ---------------------------------------------------------------------------
// Aluminum
// ---------------------------------------------------------------------------

/// Melting temperature of aluminum used by the analysis [K] (NIST: 933.47 K).
/// The same value is used for the solid/liquid switch in cp_al(), so that the
/// heat capacity and the phase change are always consistent.
const double kAlMeltingPoint = 934.0;

/// Molar mass of aluminum [g/mol].
const double kAlMolarMass = 26.98;

/// Specific heat capacity of aluminum [J/g/K] as a function of T [K].
///
/// Source: NIST Chemistry WebBook, aluminum, condensed-phase thermochemistry
/// (Shomate parameters; solid 298-933.45 K, liquid 933.45-2790.812 K):
/// https://webbook.nist.gov/cgi/cbook.cgi?ID=C7429905&Units=SI&Mask=2#Thermo-Condensed
inline double cp_al(double T)
{
    static const std::array<double, 5> solid  = {28.08920, -5.414849, 8.560423,
                                                 3.427370, -0.277375};
    static const std::array<double, 5> liquid = {31.75104, 3.935826e-8, -1.786515e-8,
                                                 2.694171e-9, 5.480037e-9};
    const auto& par = (T < kAlMeltingPoint) ? solid : liquid;
    return cp_mol(par, T) / kAlMolarMass;   // [J/mol/K] -> [J/g/K]
}

/// Aluminum constants used by the demonstration:
///   rho      = 2.70 g/cm^3
///   Tm       = kAlMeltingPoint
///   Lf       = 1045 J/cm^3 latent heat of fusion per unit volume
///   dT_shock = 74 K practical thermal-shock threshold adopted for the
///              demonstration (see getDeltaTempShock.C for its relation to
///              the elastic estimate sigma_UTS (1 - nu) / (E alpha))
///   k        = 210 W/m/K
const MaterialProps Al = {2.70, kAlMeltingPoint, 1045.0, 74.0, 210.0};

// ---------------------------------------------------------------------------
// Energy density -> temperature rise
// ---------------------------------------------------------------------------

/// Result of computeDeltaT().
struct DeltaTResult
{
    double deltaT;       ///< temperature rise T - T0 [K]
    bool   melted;       ///< true if the melting temperature was reached
    bool   shockFailed;  ///< true if deltaT exceeds the thermal-shock threshold
};

/// Heat one phase of the material: absorb the energy density E [J/cm^3]
/// starting from temperature T [K], integrating dE = rho * cp(T) dT with
/// forward-Euler steps of at most dEmax.  The integration stops when the
/// energy is exhausted or when T reaches Tstop; T and E are updated in place.
///
/// The last step absorbs exactly the energy that is left (and a step that
/// would cross Tstop absorbs only what is needed to reach it), so the result
/// is a continuous function of E without a staircase or floating-point drift.
inline void heatSinglePhase(double& T, double& E, double Tstop, double rho,
                            const std::function<double(double)>& cp,
                            double dEmax)
{
    while (E > 0.0 && T < Tstop)
    {
        const double rhoCp = rho * cp(T);        // volumetric heat capacity [J/cm^3/K]
        const double dE    = std::min(dEmax, E); // never absorb more than is left
        if (T + dE / rhoCp >= Tstop)
        {
            E -= std::min(E, (Tstop - T) * rhoCp);
            T  = Tstop;
        }
        else
        {
            T += dE / rhoCp;
            E -= dE;
        }
    }
}

/// Convert a deposited energy density Edep [J/cm^3] into the instantaneous
/// (adiabatic) temperature rise of a material initially at T0.
///
///   1. Solid heating:  dE = rho * cp(T) dT  from T0 up to Tm.
///   2. Melting:        the next Lf [J/cm^3] are absorbed at constant T = Tm
///                      (a partially molten cell has deltaT = Tm - T0).
///   3. Liquid heating: dE = rho * cp(T) dT  above Tm.
///
/// If Lf <= 0 the material has no liquid phase in this model and the
/// temperature rise saturates at Tm - T0.  Vaporization is not modelled.
inline DeltaTResult computeDeltaT(double Edep, const MaterialProps& mat,
                                  const std::function<double(double)>& cp_mat,
                                  double dEmax = kDefaultEnergyStep)
{
    double T = T0;
    double E = std::max(Edep, 0.0);   // energy density still to be absorbed

    // 1. Solid phase.
    heatSinglePhase(T, E, mat.Tm, mat.rho, cp_mat, dEmax);
    if (T < mat.Tm)
    {
        const double dT = T - T0;
        return {dT, false, dT > mat.dT_shock};
    }

    // 2. Phase change at Tm.
    const double dTm = mat.Tm - T0;
    if (mat.Lf <= 0.0 || E <= mat.Lf)
        return {dTm, mat.Lf > 0.0, dTm > mat.dT_shock};
    E -= mat.Lf;

    // 3. Liquid phase.
    heatSinglePhase(T, E, std::numeric_limits<double>::infinity(), mat.rho,
                    cp_mat, dEmax);
    const double dT = T - T0;
    return {dT, true, dT > mat.dT_shock};
}

/// Smallest deposited energy density [J/cm^3] above which the condition
/// reached(computeDeltaT(Edep)) holds.  The condition must be monotonic in
/// Edep (false below the threshold, true above); this is the case for any
/// lower bound on deltaT, because DeltaT(Edep) is non-decreasing.
/// Returns +infinity if the condition is not met below 1e9 J/cm^3.
inline double energyDensityThreshold(const std::function<bool(const DeltaTResult&)>& reached,
                                     const MaterialProps& mat,
                                     const std::function<double(double)>& cp_mat,
                                     double relTolerance = 1e-10)
{
    auto isReached = [&](double e) { return reached(computeDeltaT(e, mat, cp_mat)); };

    // Bracket the threshold: isReached(lo) == false, isReached(hi) == true.
    double lo = 0.0, hi = 1.0;
    while (!isReached(hi))
    {
        lo = hi;
        hi *= 2.0;
        if (hi > 1e9) return std::numeric_limits<double>::infinity();
    }
    // Bisection.
    while (hi - lo > relTolerance * hi)
    {
        const double mid = 0.5 * (lo + hi);
        if (isReached(mid))
            hi = mid;
        else
            lo = mid;
    }
    return hi;
}

/// Smallest deposited energy density [J/cm^3] for which computeDeltaT()
/// gives a temperature rise strictly larger than dTthreshold, i.e.
/// "deltaT > dTthreshold" is equivalent to "Edep > energyDensityForDeltaT()".
/// For dTthreshold = Tm - T0 this is the energy density needed to melt the
/// material completely (end of the melting plateau).
inline double energyDensityForDeltaT(double dTthreshold, const MaterialProps& mat,
                                     const std::function<double(double)>& cp_mat,
                                     double relTolerance = 1e-10)
{
    return energyDensityThreshold(
        [&](const DeltaTResult& r) { return r.deltaT > dTthreshold; },
        mat, cp_mat, relTolerance);
}

/// Deposited energy density [J/cm^3] needed to heat the material from T0 to
/// its melting temperature Tm (start of the melting plateau).  Melting has
/// started in a cell with Edep above this value; the cell is completely
/// molten above this value + Lf.
inline double energyDensityForMeltingOnset(const MaterialProps& mat,
                                           const std::function<double(double)>& cp_mat,
                                           double relTolerance = 1e-10)
{
    const double dTm = mat.Tm - T0;
    return energyDensityThreshold(
        [&](const DeltaTResult& r) { return r.deltaT >= dTm; },
        mat, cp_mat, relTolerance);
}

// ---------------------------------------------------------------------------
// Material lookup
// ---------------------------------------------------------------------------

/// Material constants together with the heat-capacity function.
struct MaterialModel
{
    MaterialProps props;
    std::function<double(double)> cp;   ///< specific heat [J/g/K] vs T [K]
};

/// Look up the thermal model of a material by name.  Only aluminum ("Al")
/// is implemented in this example; returns false for any other name.
///
/// To add a material, define its MaterialProps and cp function above and add
/// a branch here (see README, "Adding a material").
inline bool getMaterialModel(const TString& name, MaterialModel& model)
{
    if (name == "Al" || name == "AL")
    {
        model = {Al, cp_al};
        return true;
    }
    std::cerr << "ERROR: unsupported material \"" << name
              << "\" (only \"Al\" is implemented)." << std::endl;
    return false;
}

#endif // COLLIMATOR_THERMAL_MODEL_H
