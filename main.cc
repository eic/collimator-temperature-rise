/*
 * Collimator Temperature-Rise Example
 *
 * Minimal Geant4 example accompanying the collimator robustness study.
 *
 * The program transports electrons through a simple tapered aluminum
 * collimator.  A Geant4 command-based box scoring mesh (defined in run.mac)
 * records the deposited energy.  The thermal conversion and robustness
 * analysis are performed separately by the ROOT macros in ana/.
 *
 * Usage:
 *   collimator_temp_rise <macro> <material> <length_mm> <tan_theta>
 *                        <sigma_x_mm> <sigma_y_mm> <beam_energy_GeV>
 *                        <output.root> [seed]
 *
 *   macro            Geant4 macro file.  If its file name contains "vis"
 *                    (e.g. vis.mac), an interactive session is started.
 *   material         collimator material; only "Al" (G4_Al) is supported
 *   length_mm        collimator tip length L1 [mm]; 0 = no collimator
 *   tan_theta        tan(theta) of the collimator taper (> 0)
 *   sigma_x_mm       Gaussian beam size sigma_x [mm] (>= 0)
 *   sigma_y_mm       Gaussian beam size sigma_y [mm] (>= 0)
 *   beam_energy_GeV  total energy of the beam electrons [GeV]
 *   output.root      ROOT file for the run configuration and primaries
 *   seed             optional integer random seed (default: current time)
 *
 * Example:
 *   ./collimator_temp_rise run.mac Al 50 0.1 1 1 10 output.root
 */

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"

#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4ScoringManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>

namespace
{
void PrintUsage(const char* program)
{
    G4cout << "\nUsage:\n  " << program
           << " <macro> <material> <length_mm> <tan_theta>"
           << " <sigma_x_mm> <sigma_y_mm> <beam_energy_GeV> <output.root> [seed]\n\n"
           << "Example:\n  " << program << " run.mac Al 50 0.1 1 1 10 output.root\n"
           << G4endl;
}

/// Parse a floating-point command-line argument.  Returns false if the text
/// is not a complete, finite number.
bool ParseDouble(const char* text, G4double& value)
{
    char* end = nullptr;
    errno = 0;
    value = std::strtod(text, &end);
    return end != text && *end == '\0' && errno == 0 && std::isfinite(value);
}

/// Parse an integer command-line argument (the random seed).
bool ParseLong(const char* text, long& value)
{
    char* end = nullptr;
    errno = 0;
    value = std::strtol(text, &end, 10);
    return end != text && *end == '\0' && errno == 0;
}
}  // namespace

int main(int argc, char** argv)
{
    // The run configuration is given on the command line so that the same
    // executable can be used for different test cases.
    if (argc != 9 && argc != 10)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    const G4String macFileName = argv[1];
    const G4String collMat = argv[2];
    const G4String fileName = argv[8];

    G4double collLength = 0, collTanTh = 0;
    G4double beamSigmaX = 0, beamSigmaY = 0, beamEnergy = 0;
    if (!ParseDouble(argv[3], collLength) || !ParseDouble(argv[4], collTanTh) ||
        !ParseDouble(argv[5], beamSigmaX) || !ParseDouble(argv[6], beamSigmaY) ||
        !ParseDouble(argv[7], beamEnergy))
    {
        G4cerr << "ERROR: <length_mm>, <tan_theta>, <sigma_x_mm>, <sigma_y_mm> and "
                  "<beam_energy_GeV> must be numbers."
               << G4endl;
        PrintUsage(argv[0]);
        return 1;
    }
    if (collLength < 0 || collTanTh <= 0 || beamSigmaX < 0 || beamSigmaY < 0 ||
        beamEnergy <= 0)
    {
        G4cerr << "ERROR: require length_mm >= 0, tan_theta > 0, sigma_x_mm >= 0, "
                  "sigma_y_mm >= 0 and beam_energy_GeV > 0."
               << G4endl;
        return 1;
    }
    if (collMat != "Al" && collMat != "AL")
    {
        G4cerr << "ERROR: unsupported material \"" << collMat
               << "\". Only \"Al\" (G4_Al) is implemented." << G4endl;
        return 1;
    }
    if (!std::ifstream(macFileName))
    {
        G4cerr << "ERROR: cannot open the macro file \"" << macFileName << "\"." << G4endl;
        return 1;
    }

    // Convert to Geant4 internal units.
    collLength *= mm;
    beamSigmaX *= mm;
    beamSigmaY *= mm;
    beamEnergy *= GeV;

    // Random seed: taken from the command line for a reproducible run,
    // otherwise from the current time.  It is printed and stored in the
    // output ROOT file (tree "meta").
    long seed = static_cast<long>(std::time(nullptr));
    if (argc == 10 && !ParseLong(argv[9], seed))
    {
        G4cerr << "ERROR: the random seed must be an integer." << G4endl;
        return 1;
    }

    // Interactive mode: start a UI session (Qt if available) when the macro
    // file name, without its directory, contains "vis" (e.g. vis.mac).
    const std::string macBaseName =
        macFileName.substr(macFileName.find_last_of('/') + 1);
    G4UIExecutive* ui = nullptr;
    if (macBaseName.find("vis") != std::string::npos)
    {
        ui = new G4UIExecutive(argc, argv);
    }

    G4cout << "\n\n===========================================\n"
           << " Collimator temperature-rise example\n"
           << "===========================================\n"
           << " Macro          : " << macFileName << (ui ? " (interactive)" : " (batch)") << G4endl
           << " Material       : " << collMat << G4endl
           << " Tip length     : " << collLength / mm << " mm" << G4endl
           << " tan(theta)     : " << collTanTh << G4endl
           << " Beam sigma_x   : " << beamSigmaX / mm << " mm" << G4endl
           << " Beam sigma_y   : " << beamSigmaY / mm << " mm" << G4endl
           << " Beam energy    : " << beamEnergy / GeV << " GeV" << G4endl
           << " ROOT output    : " << fileName << G4endl
           << " Random seed    : " << seed << G4endl
           << "===========================================\n"
           << G4endl;

    // CLHEP random engine.
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    G4Random::setTheSeed(seed);
    G4Random::showEngineStatus();

    // Sequential run manager: RunAction writes a single ROOT file and is not
    // designed for multi-threaded event processing.
    auto* runManager = new G4RunManager;

    // Activate command-based scoring; the mesh itself is defined in run.mac.
    auto* scoringManager = G4ScoringManager::GetScoringManager();
    scoringManager->SetVerboseLevel(1);

    // Reference physics list with the most accurate standard EM option (EMZ).
    G4PhysListFactory factory;
    const G4String physName = "FTFP_BERT_EMZ";
    runManager->SetUserInitialization(factory.GetReferencePhysList(physName));

    // RunAction holds the run configuration and is shared with the
    // detector construction and the primary generator.
    auto* runAction = new RunAction(fileName, collMat, collLength, collTanTh,
                                    beamSigmaX, beamSigmaY, beamEnergy);
    runAction->SetRandomSeed(seed);

    runManager->SetUserInitialization(new DetectorConstruction(runAction));
    runManager->SetUserAction(runAction);
    runManager->SetUserAction(new PrimaryGeneratorAction(runAction));

    runManager->Initialize();

    // Visualization (used only by macros that issue /vis/ commands).
    auto* visManager = new G4VisExecutive;
    visManager->Initialize();

    auto* uiManager = G4UImanager::GetUIpointer();
    uiManager->ApplyCommand("/control/execute " + macFileName);
    if (ui)
    {
        ui->SessionStart();
        delete ui;
    }

    delete visManager;
    delete runManager;
    return 0;
}
