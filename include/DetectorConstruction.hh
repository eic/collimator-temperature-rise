#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4VPhysicalVolume;
class RunAction;

/**
 * Vacuum world (4 x 4 x 4 m^3) containing a simple tapered collimator jaw
 * (G4Trd) centered at the origin.  The tip length L1, the taper tan(theta)
 * and the material are taken from the command line via RunAction.
 * See DetectorConstruction.cc for the geometry definition.
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    explicit DetectorConstruction(RunAction* runAct);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

private:
    G4Material* GetMaterial(const G4String& matName);
    RunAction* _runAction;
};

#endif
