/*
 * readScoring.C
 *
 * Convert the Geant4 CSV scoring output into a ROOT THnSparse.
 *
 * Usage (from the results/ directory, see README):
 *   root -l -b -q '../ana/readScoring.C("../build/edep.csv","edep.root")'
 *
 * Input
 *   edep.csv written by  /score/dumpQuantityToFile boxMesh edep edep.csv.
 *   Geant4 writes one line per mesh cell (also for empty cells):
 *     iX, iY, iZ, total(value) [MeV], total(val^2), entry
 *   where, for each cell,
 *     total(value) = S1 = sum over events of the energy deposited per event
 *     total(val^2) = S2 = sum over events of (energy deposited per event)^2
 *     entry        = number of events that deposited energy in the cell.
 *
 * Output
 *   ROOT file containing the THnSparseF "h1" (axes x, y, z in cm):
 *     bin content = S1 [MeV]
 *     bin error^2 = S2 [MeV^2]  (used by drawScoring.C for the statistical
 *                                uncertainty)
 *   Empty cells are not stored.
 *
 * Mesh arguments
 *   The mesh arguments must describe the mesh defined in run.mac.  Note that
 *   /score/mesh/boxSize takes HALF-lengths, so the full box size is twice
 *   the run.mac values.  The defaults match run.mac:
 *     /score/mesh/boxSize 2 2 40 cm        -> 4 x 4 x 80 cm^3
 *     /score/mesh/nBin 40 40 800           -> 1 x 1 x 1 mm^3 cells
 *     /score/mesh/translate/xyz 0 -0.2 0 cm
 *   The macro stops with an error if the CSV file is inconsistent with the
 *   mesh arguments (wrong number of cells or out-of-range indices).
 */
#include <TAxis.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TString.h>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

void readScoring(
    TString fInName  = "../build/edep.csv",
    TString fOutName = "edep.root",
    Float_t boxSizeX = 4.0,     // [cm] FULL length = 2 x run.mac boxSize
    Float_t boxSizeY = 4.0,     // [cm] FULL length = 2 x run.mac boxSize
    Float_t boxSizeZ = 80.0,    // [cm] FULL length = 2 x run.mac boxSize
    Int_t   nBinX = 40,
    Int_t   nBinY = 40,
    Int_t   nBinZ = 800,
    Float_t translateX = 0.0,   // [cm] mesh center
    Float_t translateY = -0.2,  // [cm] mesh center
    Float_t translateZ = 0.0    // [cm] mesh center
)
{
    const Int_t ndim = 3;
    Int_t bins[ndim] = {nBinX, nBinY, nBinZ};
    Double_t xmin[ndim] = {translateX - boxSizeX / 2.0,
                           translateY - boxSizeY / 2.0,
                           translateZ - boxSizeZ / 2.0};
    Double_t xmax[ndim] = {translateX + boxSizeX / 2.0,
                           translateY + boxSizeY / 2.0,
                           translateZ + boxSizeZ / 2.0};

    std::cout << "Input file : " << fInName << std::endl;
    std::ifstream inputFile(fInName.Data());
    if (!inputFile.is_open())
    {
        std::cerr << "ERROR: Could not open " << fInName << std::endl;
        return;
    }

    auto* h1 = new THnSparseF("h1", "Deposited energy per mesh cell [MeV]",
                              ndim, bins, xmin, xmax);
    h1->GetAxis(0)->SetTitle("x [cm]");
    h1->GetAxis(1)->SetTitle("y [cm]");
    h1->GetAxis(2)->SetTitle("z [cm]");
    h1->Sumw2();   // store the per-cell sum of squares as the bin error^2

    std::string line;
    Long64_t nLines = 0;        // data lines read
    Long64_t nFilled = 0;       // non-empty cells stored
    Long64_t nOutOfRange = 0;   // indices outside the declared mesh

    while (std::getline(inputFile, line))
    {
        if (line.empty()) continue;

        // Header lines.  Check the unit of the scored quantity.
        if (line[0] == '#')
        {
            if (line.find("total(value)") != std::string::npos &&
                line.find("[MeV]") == std::string::npos)
            {
                std::cerr << "WARNING: the scored quantity is not in MeV; drawScoring.C "
                             "assumes MeV.\n  Header: " << line << std::endl;
            }
            continue;
        }

        int iX = 0, iY = 0, iZ = 0;
        double val = 0.0, val2 = 0.0;
        long long entry = 0;
        const int nRead = std::sscanf(line.c_str(), "%d,%d,%d,%lf,%lf,%lld",
                                      &iX, &iY, &iZ, &val, &val2, &entry);
        if (nRead != 6)
        {
            std::cerr << "WARNING: Skipping malformed line: " << line << std::endl;
            continue;
        }
        ++nLines;

        if (iX < 0 || iX >= nBinX || iY < 0 || iY >= nBinY || iZ < 0 || iZ >= nBinZ)
        {
            ++nOutOfRange;
            continue;
        }

        // Empty cells are written as "0,0,0" and are not stored.
        if (val == 0.0 && val2 == 0.0) continue;

        // Geant4 mesh indices start at 0, ROOT bins at 1.
        Int_t idx[ndim] = {iX + 1, iY + 1, iZ + 1};
        const Long64_t bin = h1->GetBin(idx, kTRUE);
        h1->SetBinContent(bin, val);
        h1->SetBinError2(bin, val2);
        ++nFilled;
    }
    inputFile.close();
    h1->SetEntries(nFilled);

    // Consistency checks between the CSV file and the mesh arguments.
    const Long64_t nExpected = static_cast<Long64_t>(nBinX) * nBinY * nBinZ;
    if (nOutOfRange > 0 || nLines != nExpected)
    {
        std::cerr << "ERROR: " << fInName << " does not match the mesh arguments ("
                  << nBinX << " x " << nBinY << " x " << nBinZ << " = " << nExpected
                  << " cells): read " << nLines << " cells, " << nOutOfRange
                  << " with out-of-range indices.\n"
                  << "       Pass the mesh used in run.mac to readScoring.C "
                     "(full box lengths = 2 x /score/mesh/boxSize)." << std::endl;
        delete h1;
        return;
    }

    auto* file = TFile::Open(fOutName.Data(), "RECREATE");
    if (!file || file->IsZombie())
    {
        std::cerr << "ERROR: Could not create " << fOutName << std::endl;
        delete file;
        delete h1;
        return;
    }

    file->cd();
    h1->Write();
    file->Close();

    std::cout << "Read " << nLines << " mesh cells (" << nBinX << " x " << nBinY
              << " x " << nBinZ << ", cell size " << boxSizeX / nBinX * 10 << " x "
              << boxSizeY / nBinY * 10 << " x " << boxSizeZ / nBinZ * 10
              << " mm^3); stored " << nFilled << " non-empty cells." << std::endl;
    std::cout << "Output file: " << fOutName << std::endl;

    delete file;
    delete h1;
}
