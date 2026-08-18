#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(RunAction* runAct)
    : _particleGun(new G4ParticleGun(1)), _runAction(runAct)
{
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete _particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particleDef = particleTable->FindParticle("e-");

    // Sample the incoming beam from a 2-D Gaussian distribution.  The
    // longitudinal position is placed 1 m upstream of the collimator.
    const G4double xpos = G4RandGauss::shoot(
        0.0 * mm, _runAction->GetBeamSigmaX());

    const G4double ypos = G4RandGauss::shoot(
        -2.0 * mm, _runAction->GetBeamSigmaY());

    const G4double zpos = -1.0 * m;

    // Geant4 expects kinetic energy here.  Subtracting the electron rest mass
    // converts the user-supplied total energy into kinetic energy.
    const G4double mass = particleDef->GetPDGMass();
    const G4double Ekin = _runAction->GetBeamEnergy() - mass;

    _particleGun->SetParticleDefinition(particleDef);
    _particleGun->SetParticleEnergy(Ekin);
    _particleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    _particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));

    // Save the generated primary parameters for reproducibility checks.
    _runAction->beamGenPosX = xpos / mm;
    _runAction->beamGenPosY = ypos / mm;
    _runAction->beamGenPosZ = zpos / mm;
    _runAction->beamGenEkin = _particleGun->GetParticleEnergy() / GeV;
    _runAction->tree->Fill();
    ++_runAction->eventID;

    _particleGun->GeneratePrimaryVertex(anEvent);
}
