#ifndef RunAction_h
#define RunAction_h 1

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TString.h"
#pragma GCC diagnostic pop

#include "G4UserRunAction.hh"
#include "G4Timer.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

class G4Run;

/**
 * Stores the run configuration and a small amount of event-level beam
 * information in a ROOT file.
 *
 * The deposited-energy mesh itself is produced by Geant4's scoring manager
 * and written to edep.csv by run.mac.  The ROOT file written here is metadata
 * and event bookkeeping.
 */
class RunAction : public G4UserRunAction
{
public:
    RunAction(G4String, G4String, G4double, G4double,
              G4double, G4double, G4double);
    ~RunAction();

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*) override;

    G4String GetOutputFilename() { return _fileName; }
    G4String GetCollMaterial()   { return _collMat; }
    G4double GetCollLength()     { return _collLength; }
    G4double GetCollTanTh()      { return _collTanTh; }
    G4double GetBeamSigmaX()     { return _beamSigmaX; }
    G4double GetBeamSigmaY()     { return _beamSigmaY; }
    G4double GetBeamEnergy()     { return _beamEnergy; }

    TFile* hfile = nullptr;
    TTree* tree = nullptr;
    TTree* meta = nullptr;

    G4int eventID = 0;
    G4double beamGenPosX = 0;
    G4double beamGenPosY = 0;
    G4double beamGenPosZ = 0;
    G4double beamGenEkin = 0;
    TString collMatStr;

private:
    G4Timer* _timer;
    G4String _fileName;
    G4String _collMat;
    G4double _collLength;
    G4double _collTanTh;
    G4double _beamSigmaX;
    G4double _beamSigmaY;
    G4double _beamEnergy;
};

#endif
