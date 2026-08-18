#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "RunAction.hh"
#include "G4ParticleTable.hh"

class G4Event;

/**
 * Generates one 10-GeV-class electron per event for the demonstration.
 *
 * The transverse position is sampled from independent Gaussian distributions
 * with user-specified sigma_x and sigma_y.  The mean position is chosen to
 * match the scoring mesh / example geometry: x = 0 and y = -2 mm.
 */
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    explicit PrimaryGeneratorAction(RunAction*);
    ~PrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

private:
    G4ParticleGun* _particleGun;
    RunAction* _runAction;
};

#endif
