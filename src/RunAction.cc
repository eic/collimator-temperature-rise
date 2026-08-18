#include "RunAction.hh"
#include <cstdlib>

RunAction::RunAction(
    G4String fileName, G4String collMat, G4double collLength,
    G4double collTanTh, G4double beamSigmaX, G4double beamSigmaY,
    G4double beamEnergy)
    : _timer(new G4Timer),
      _fileName(fileName),
      _collMat(collMat),
      _collLength(collLength),
      _collTanTh(collTanTh),
      _beamSigmaX(beamSigmaX),
      _beamSigmaY(beamSigmaY),
      _beamEnergy(beamEnergy)
{
    collMatStr = _collMat.c_str();
}

RunAction::~RunAction()
{
    delete _timer;
}

void RunAction::BeginOfRunAction(const G4Run*)
{
    G4cout << "[RunAction] BeginOfRunAction" << G4endl;
    _timer->Start();

    // Store the input configuration and generated primary information in a
    // compact ROOT file.  Geant4 scoring data are written separately by the
    // scoring commands in run.mac.
    hfile = new TFile(_fileName, "RECREATE");
    if (!hfile || hfile->IsZombie())
    {
        G4cerr << "[ERROR] Could not create ROOT output file: "
               << _fileName << G4endl;
        std::exit(EXIT_FAILURE);
    }

    tree = new TTree("tree", "Generated primary particles");
    tree->Branch("eventID", &eventID);
    tree->Branch("beamGenPosX", &beamGenPosX);
    tree->Branch("beamGenPosY", &beamGenPosY);
    tree->Branch("beamGenPosZ", &beamGenPosZ);
    tree->Branch("beamGenEkin", &beamGenEkin);

    meta = new TTree("meta", "Simulation settings");
    meta->Branch("collMat", &collMatStr);
    meta->Branch("collLength", &_collLength);
    meta->Branch("collTanTh", &_collTanTh);
    meta->Branch("beamSigmaX", &_beamSigmaX);
    meta->Branch("beamSigmaY", &_beamSigmaY);
    meta->Branch("beamEnergy", &_beamEnergy);

    eventID = 0;
}

void RunAction::EndOfRunAction(const G4Run*)
{
    if (!hfile) return;

    hfile->cd();
    meta->Fill();
    meta->Write();
    tree->Write();
    hfile->Close();

    _timer->Stop();
    G4cout << "[INFO] Output file = " << _fileName << G4endl;
    G4cout << "[INFO] CPU time: " << *_timer << G4endl;
}
