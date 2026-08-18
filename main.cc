/*
 * Collimator Temperature-Rise Example
 *
 * Minimal Geant4 example accompanying the collimator robustness study.
 *
 * The program transports electrons through a simple tapered aluminum
 * collimator and uses a Geant4 box scoring mesh to record the deposited
 * energy density.  The thermal conversion and robustness analysis are
 * performed separately by ROOT macros in ana/.
 *
 * Command-line arguments:
 *   1. Geant4 macro file
 *   2. Collimator material
 *   3. Collimator tip length [mm]
 *   4. tan(theta) of the taper
 *   5. Beam sigma_x [mm]
 *   6. Beam sigma_y [mm]
 *   7. Beam energy [GeV]
 *   8. ROOT output file
 *
 * Example:
 *   ./collimator_temp_rise run.mac Al 50 0.1 1 1 10 output.root
 */

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4ScoringManager.hh"

#include <ctime>
#include <cstdlib>

int main(int argc, char** argv)
{
    G4String macFileName, collMat, fileName;
    G4double collLength = 0, collTanTh = 0;
    G4double beamSigmaX = 0, beamSigmaY = 0, beamEnergy = 0;

    // The example deliberately keeps the run configuration on the command
    // line so that the same executable can be used for different test cases.
    if (argc == 9)
    {
        macFileName = argv[1];
        collMat = argv[2];
        collLength = std::atof(argv[3]) * mm;
        collTanTh = std::atof(argv[4]);
        beamSigmaX = std::atof(argv[5]) * mm;
        beamSigmaY = std::atof(argv[6]) * mm;
        beamEnergy = std::atof(argv[7]) * GeV;
        fileName = argv[8];

        G4cout << "\n\n===========================================\n"
               << " Collimator temperature-rise example\n"
               << "===========================================\n"
               << " Macro          : " << macFileName << G4endl
               << " Material       : " << collMat << G4endl
               << " Tip length     : " << collLength/mm << " mm" << G4endl
               << " tan(theta)     : " << collTanTh << G4endl
               << " Beam sigma_x   : " << beamSigmaX/mm << " mm" << G4endl
               << " Beam sigma_y   : " << beamSigmaY/mm << " mm" << G4endl
               << " Beam energy    : " << beamEnergy/GeV << " GeV" << G4endl
               << " ROOT output    : " << fileName << G4endl
               << "===========================================\n"
               << G4endl;
    }
    else
    {
        G4cout << "\nUsage:\n  "
               << argv[0]
               << " <macro> <material> <length_mm> <tan_theta>"
               << " <sigma_x_mm> <sigma_y_mm> <beam_energy_GeV> <output.root>\n\n"
               << "Example:\n  "
               << argv[0]
               << " run.mac Al 50 0.1 1 1 10 output.root\n"
               << G4endl;
        return 1;
    }

    // CLHEP random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4Random::setTheSeed(static_cast<long>(std::time(nullptr)));
    G4Random::showEngineStatus();

    auto* runManager = new G4RunManager;

    auto* scoringManager = G4ScoringManager::GetScoringManager();
    scoringManager->SetVerboseLevel(1);

    G4PhysListFactory factory;
    const G4String physName = "FTFP_BERT_EMZ";
    auto* physics = factory.GetReferencePhysList(physName);
    runManager->SetUserInitialization(physics);

    auto* runAction = new RunAction(
        fileName, collMat, collLength, collTanTh,
        beamSigmaX, beamSigmaY, beamEnergy);

    runManager->SetUserInitialization(new DetectorConstruction(runAction));
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(new PrimaryGeneratorAction(runAction));

    runManager->Initialize();

    // Initialize visualization
    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    auto* uiManager = G4UImanager::GetUIpointer();

    if (G4StrUtil::contains(macFileName, "vis"))
    {
        auto* ui = new G4UIExecutive(argc, argv);
        uiManager->ApplyCommand("/control/execute " + macFileName);
        ui->SessionStart();
        delete ui;
    }
    else
    {
        uiManager->ApplyCommand("/control/execute " + macFileName);
    }

    delete visManager;
    delete runManager;
    return 0;
}
