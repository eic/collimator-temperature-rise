#include "DetectorConstruction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Trd.hh"
#include "G4VisAttributes.hh"

#include <cmath>

// This example intentionally keeps the geometry simple.  The production
// collimator model can be substituted for the G4Trd below without changing
// the downstream scoring/thermal-analysis concept.

DetectorConstruction::DetectorConstruction(RunAction* runAct)
    : _runAction(runAct)
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // World: vacuum box, much larger than the collimator.  The beam starts
    // 1 m upstream of the collimator center (PrimaryGeneratorAction).
    const G4double worldHalfSize = 2.0 * m;
    auto* world_solid = new G4Box("world", worldHalfSize, worldHalfSize, worldHalfSize);
    auto* world_log = new G4LogicalVolume(world_solid, GetMaterial("Vacuum"), "world");
    auto* world_phys = new G4PVPlacement(nullptr, G4ThreeVector(), world_log, "world",
                                         nullptr, false, 0);
    world_log->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Collimator.  A tip length L1 <= 0 builds no collimator (empty world).
    const G4double L1 = _runAction->GetCollLength();
    if (L1 <= 0)
    {
        G4cout << "[INFO] DetectorConstruction::Construct ==> tip length is 0: "
                  "no collimator is built" << G4endl;
        return world_phys;
    }

    G4cout << "[INFO] DetectorConstruction::Construct ==> Building a simple collimator" << G4endl;
    // Side view (Z-Y plane, beam along +Z) after the rotation below:
    //
    //                    |<-L1->|
    //          y = +H/2   ------
    //    Beam =====>     /      \        H  = full height of the jaw (Y)
    //                   /        \       W  = full width of the jaw (X)
    //                  /)Th       \      tan(Th) = 2 H / (L2 - L1)
    //          y = -H/2 ------------     =>  L2 = L1 + 2 H / tan(Th)
    //                  |<--- L2 --->|
    //
    const G4double H = 30 * mm;   // full height of the jaw
    const G4double W = 30 * mm;   // full width of the jaw
    const G4double L2 = L1 + 2. * H / _runAction->GetCollTanTh();
    G4cout << "L1 = " << L1 / mm << " [mm]; L2 = " << L2 / mm << " [mm]; W = " << W / mm
           << " [mm]; H = " << H / mm << " [mm]" << G4endl;

    // The jaw must not reach the primary vertex (and must fit in the world).
    const G4double sourceDistance = std::abs(PrimaryGeneratorAction::kSourcePositionZ);
    if (L2 / 2. >= sourceDistance)
    {
        G4ExceptionDescription msg;
        msg << "The collimator base length L2 = " << L2 / mm << " mm is too long: the jaw "
            << "would reach the primary vertex at z = -" << sourceDistance / mm << " mm.\n"
            << "Increase tan(theta) or reduce the tip length.";
        G4Exception("DetectorConstruction::Construct()", "CollTempRise002",
                    FatalException, msg);
    }

    // G4Trd(dx1, dx2, dy1, dy2, dz) takes half-lengths: the y half-length is
    // L2/2 at local z = -H/2 and L1/2 at local z = +H/2.  The 90-degree
    // rotation about X maps local y -> global -z and local z -> global y, so
    // the taper lies along the beam (Z) with the tip (L1) at y = +H/2 and the
    // base (L2) at y = -H/2.
    auto* coll_trd = new G4Trd("coll_trd", W / 2., W / 2., L2 / 2., L1 / 2., H / 2.);
    auto* coll_log = new G4LogicalVolume(
        coll_trd, GetMaterial(_runAction->GetCollMaterial()), "coll_trd");
    auto* coll_rot = new G4RotationMatrix();
    coll_rot->rotateX(90 * deg);

    new G4PVPlacement(coll_rot, G4ThreeVector(0., 0., 0.), coll_log, "coll_box",
                      world_log, false, 0, true /* check overlaps */);
    coll_log->SetVisAttributes(G4VisAttributes(G4Colour::Magenta()));

    return world_phys;
}

G4Material* DetectorConstruction::GetMaterial(const G4String& matName)
{
    // The public demonstration uses only "Al", mapped to Geant4's G4_Al.
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* material = nullptr;

    if (matName == "Vacuum")
    {
        material = nistManager->FindOrBuildMaterial("G4_Galactic");
    }
    else if (matName == "Al" || matName == "AL")
    {
        material = nistManager->FindOrBuildMaterial("G4_Al");
    }
    else
    {
        G4ExceptionDescription msg;
        msg << "Unknown collimator material \"" << matName
            << "\". Only \"Al\" (G4_Al) is supported.";
        G4Exception("DetectorConstruction::GetMaterial()", "CollTempRise001",
                    FatalException, msg);
        return nullptr;
    }

    G4cout << "[MATERIAL] " << *material << G4endl;
    return material;
}
