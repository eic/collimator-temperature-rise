#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class G4ParticleDefinition;
class G4ParticleGun;
class RunAction;

/**
 * Generates one electron per event travelling along +z.
 *
 * The vertex is at z = kSourcePositionZ (1 m upstream of the collimator
 * center).  The transverse position is sampled from independent Gaussian
 * distributions with the user-specified sigma_x and sigma_y around
 * (kBeamMeanX, kBeamMeanY) = (0, -2 mm); the scoring mesh in run.mac is
 * centered at the same y.  The energy given on the command line is the total
 * energy; the corresponding kinetic energy is passed to the particle gun.
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    static constexpr G4double kSourcePositionZ = -1.0 * CLHEP::m;
    static constexpr G4double kBeamMeanX = 0.0 * CLHEP::mm;
    static constexpr G4double kBeamMeanY = -2.0 * CLHEP::mm;

    explicit PrimaryGeneratorAction(RunAction* runAct);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

private:
    G4ParticleGun* _particleGun;
    G4ParticleDefinition* _electron;
    G4double _kineticEnergy;
    RunAction* _runAction;
};

#endif
