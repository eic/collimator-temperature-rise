#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "RunAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    explicit DetectorConstruction(RunAction* runAct);
    ~DetectorConstruction() override;

    G4VPhysicalVolume* Construct() override;

private:
    G4Material* GetMaterial(G4String matName);
    RunAction* _runAction;
};

#endif
