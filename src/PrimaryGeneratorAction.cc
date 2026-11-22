#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"

#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAct)
    : _particleGun(new G4ParticleGun(1)),
      _electron(G4Electron::Definition()),
      _runAction(runAct)
{
    // The command line gives the total energy; Geant4 expects the kinetic
    // energy, so the electron rest mass is subtracted.
    _kineticEnergy = _runAction->GetBeamEnergy() - _electron->GetPDGMass();
    if (_kineticEnergy <= 0)
    {
        G4ExceptionDescription msg;
        msg << "The beam energy (" << _runAction->GetBeamEnergy() / MeV
            << " MeV) must exceed the electron rest mass.";
        G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction()",
                    "CollTempRise003", FatalException, msg);
    }
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete _particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // Sample the incoming beam from a 2-D Gaussian distribution.
    const G4double xpos = G4RandGauss::shoot(kBeamMeanX, _runAction->GetBeamSigmaX());
    const G4double ypos = G4RandGauss::shoot(kBeamMeanY, _runAction->GetBeamSigmaY());
    const G4double zpos = kSourcePositionZ;

    // The gun is fully configured for every event, so the command-line
    // settings always take precedence over /gun/ commands.
    _particleGun->SetParticleDefinition(_electron);
    _particleGun->SetParticleEnergy(_kineticEnergy);
    _particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    _particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));

    // Save the generated primary for reproducibility checks.
    _runAction->FillPrimary(xpos, ypos, zpos, _kineticEnergy);

    _particleGun->GeneratePrimaryVertex(anEvent);
}
