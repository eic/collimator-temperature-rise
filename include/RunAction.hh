#ifndef RunAction_h
#define RunAction_h 1

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TString.h"
#pragma GCC diagnostic pop

#include "G4Timer.hh"
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class TFile;
class TTree;

/**
 * Holds the run configuration given on the command line and writes it,
 * together with the generated primary particles, to a ROOT file.
 *
 * Output file (one per run; overwritten by every /run/beamOn):
 *   tree "meta" : one entry per run with the run configuration
 *                 collMat, collLength [mm], collTanTh, beamSigmaX [mm],
 *                 beamSigmaY [mm], beamEnergy [GeV], randomSeed, nEvents
 *   tree "tree" : one entry per event with the generated primary
 *                 eventID, beamGenPosX/Y/Z [mm], beamGenEkin [GeV]
 *
 * The deposited-energy mesh itself is produced by Geant4's scoring manager
 * and written to edep.csv by run.mac.
 *
 * The configuration getters return values in Geant4 internal units.
 */
class RunAction : public G4UserRunAction
{
public:
    RunAction(const G4String& fileName, const G4String& collMat,
              G4double collLength, G4double collTanTh,
              G4double beamSigmaX, G4double beamSigmaY, G4double beamEnergy);
    ~RunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    /// Record one generated primary (position, kinetic energy) in "tree".
    /// Called by PrimaryGeneratorAction once per event.
    void FillPrimary(G4double x, G4double y, G4double z, G4double ekin);

    /// Random seed used for the run (stored in "meta").
    void SetRandomSeed(long seed) { _randomSeed = seed; }

    const G4String& GetOutputFilename() const { return _fileName; }
    const G4String& GetCollMaterial() const { return _collMat; }
    G4double GetCollLength() const { return _collLength; }
    G4double GetCollTanTh() const { return _collTanTh; }
    G4double GetBeamSigmaX() const { return _beamSigmaX; }
    G4double GetBeamSigmaY() const { return _beamSigmaY; }
    G4double GetBeamEnergy() const { return _beamEnergy; }

private:
    G4Timer _timer;

    // Run configuration (Geant4 internal units).
    G4String _fileName;
    G4String _collMat;
    G4double _collLength;
    G4double _collTanTh;
    G4double _beamSigmaX;
    G4double _beamSigmaY;
    G4double _beamEnergy;
    long _randomSeed = 0;

    // ROOT output.  The trees are owned by the file.
    TFile* _file = nullptr;
    TTree* _tree = nullptr;
    TTree* _meta = nullptr;

    // Branch buffers of "tree" (one entry per event).
    G4int _eventID = 0;
    G4double _beamGenPosX = 0;   // [mm]
    G4double _beamGenPosY = 0;   // [mm]
    G4double _beamGenPosZ = 0;   // [mm]
    G4double _beamGenEkin = 0;   // [GeV]

    // Branch buffers of "meta" (one entry per run).
    TString _collMatStr;
    G4double _metaCollLength = 0;   // [mm]
    G4double _metaCollTanTh = 0;
    G4double _metaBeamSigmaX = 0;   // [mm]
    G4double _metaBeamSigmaY = 0;   // [mm]
    G4double _metaBeamEnergy = 0;   // [GeV]
    Long64_t _metaRandomSeed = 0;
    Long64_t _metaNEvents = 0;
};

#endif
