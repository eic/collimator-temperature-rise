#include "RunAction.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include "TFile.h"
#include "TTree.h"
#pragma GCC diagnostic pop

#include "G4Run.hh"
#include "G4SystemOfUnits.hh"

RunAction::RunAction(const G4String& fileName, const G4String& collMat,
                     G4double collLength, G4double collTanTh,
                     G4double beamSigmaX, G4double beamSigmaY,
                     G4double beamEnergy)
    : _fileName(fileName),
      _collMat(collMat),
      _collLength(collLength),
      _collTanTh(collTanTh),
      _beamSigmaX(beamSigmaX),
      _beamSigmaY(beamSigmaY),
      _beamEnergy(beamEnergy)
{
    _collMatStr = _collMat.c_str();
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4cout << "[RunAction] BeginOfRunAction" << G4endl;
    _timer.Start();

    // Store the run configuration and the generated primaries in a compact
    // ROOT file.  Geant4 scoring data are written separately by the scoring
    // commands in run.mac.
    _file = TFile::Open(_fileName.c_str(), "RECREATE");
    if (!_file || _file->IsZombie())
    {
        G4ExceptionDescription msg;
        msg << "Could not create the ROOT output file " << _fileName;
        G4Exception("RunAction::BeginOfRunAction()", "CollTempRise010",
                    FatalException, msg);
        return;
    }

    _tree = new TTree("tree", "Generated primaries: positions [mm], kinetic energy [GeV]");
    _tree->Branch("eventID", &_eventID);
    _tree->Branch("beamGenPosX", &_beamGenPosX);
    _tree->Branch("beamGenPosY", &_beamGenPosY);
    _tree->Branch("beamGenPosZ", &_beamGenPosZ);
    _tree->Branch("beamGenEkin", &_beamGenEkin);

    // Explicit units for the stored configuration.
    _metaCollLength = _collLength / mm;
    _metaCollTanTh = _collTanTh;
    _metaBeamSigmaX = _beamSigmaX / mm;
    _metaBeamSigmaY = _beamSigmaY / mm;
    _metaBeamEnergy = _beamEnergy / GeV;
    _metaRandomSeed = _randomSeed;
    _metaNEvents = 0;

    _meta = new TTree("meta", "Run configuration: lengths [mm], beam energy [GeV]");
    _meta->Branch("collMat", &_collMatStr);
    _meta->Branch("collLength", &_metaCollLength);
    _meta->Branch("collTanTh", &_metaCollTanTh);
    _meta->Branch("beamSigmaX", &_metaBeamSigmaX);
    _meta->Branch("beamSigmaY", &_metaBeamSigmaY);
    _meta->Branch("beamEnergy", &_metaBeamEnergy);
    _meta->Branch("randomSeed", &_metaRandomSeed);
    _meta->Branch("nEvents", &_metaNEvents);

    _eventID = 0;
}

void RunAction::FillPrimary(G4double x, G4double y, G4double z, G4double ekin)
{
    if (!_tree) return;

    _beamGenPosX = x / mm;
    _beamGenPosY = y / mm;
    _beamGenPosZ = z / mm;
    _beamGenEkin = ekin / GeV;
    _tree->Fill();
    ++_eventID;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    if (!_file) return;

    _metaNEvents = run->GetNumberOfEvent();

    _file->cd();
    _meta->Fill();
    _meta->Write();
    _tree->Write();
    _file->Close();   // also deletes the trees attached to the file
    delete _file;
    _file = nullptr;
    _tree = nullptr;
    _meta = nullptr;

    _timer.Stop();
    G4cout << "[INFO] Output file = " << _fileName << " (" << _metaNEvents
           << " events)" << G4endl;
    G4cout << "[INFO] CPU time: " << _timer << G4endl;
}
