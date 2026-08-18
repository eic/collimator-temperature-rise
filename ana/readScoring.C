/*
 * Convert the Geant4 CSV scoring output into a ROOT THnSparse.
 *
 * Input:
 *   edep.csv written by /score/dumpQuantityToFile.
 *
 * Output:
 *   ROOT file containing THnSparseF "h1".
 *
 * The default mesh exactly matches run.mac:
 *   2 x 2 x 40 cm^3, 40 x 40 x 800 bins.
 */
#include <THnSparse.h>
#include <TString.h>
#include <TFile.h>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void readScoring(
    TString fInName = "../build/edep.csv",
    TString fOutName = "edep.root",
    Float_t boxSizeX = 4.0,   // [cm]
    Float_t boxSizeY = 4.0,   // [cm]
    Float_t boxSizeZ = 80.0,  // [cm]
    Int_t nBinX = 40,
    Int_t nBinY = 40,
    Int_t nBinZ = 800,
    Float_t translateX = 0.0,   // [cm]
    Float_t translateY = -0.2,  // [cm]
    Float_t translateZ = 0.0    // [cm]
)
{
    const Int_t ndim = 3;
    Int_t bins[ndim] = {nBinX, nBinY, nBinZ};
    Double_t xmin[ndim] = {
        translateX - boxSizeX/2.0,
        translateY - boxSizeY/2.0,
        translateZ - boxSizeZ/2.0
    };
    Double_t xmax[ndim] = {
        translateX + boxSizeX/2.0,
        translateY + boxSizeY/2.0,
        translateZ + boxSizeZ/2.0
    };

    auto* h1 = new THnSparseF("h1", "3-D deposited-energy density",
                              ndim, bins, xmin, xmax);

    std::cout << "Input file : " << fInName << std::endl;
    std::ifstream inputFile(fInName.Data());
    if (!inputFile.is_open())
    {
        std::cerr << "ERROR: Could not open " << fInName << std::endl;
        delete h1;
        return;
    }

    std::string line;
    Long64_t nLines = 0;
    Long64_t nFilled = 0;

    while (std::getline(inputFile, line))
    {
        // Ignore empty lines and Geant4 comment/header lines.
        if (line.empty() || line[0] == '#') continue;

        int iX = 0, iY = 0, iZ = 0, entry = 0;
        float val = 0.0f, val2 = 0.0f;

        // Geant4 writes: ix,iy,iz,value,error,entry.
        const int nRead = std::sscanf(
            line.c_str(), "%d,%d,%d,%f,%f,%d",
            &iX, &iY, &iZ, &val, &val2, &entry);

        if (nRead != 6)
        {
            std::cerr << "WARNING: Skipping malformed line: "
                      << line << std::endl;
            continue;
        }

        ++nLines;

        // Convert integer mesh indices to the physical bin-center position.
        const double x = translateX - boxSizeX/2.0
                       + (iX + 0.5) * boxSizeX/nBinX;
        const double y = translateY - boxSizeY/2.0
                       + (iY + 0.5) * boxSizeY/nBinY;
        const double z = translateZ - boxSizeZ/2.0
                       + (iZ + 0.5) * boxSizeZ/nBinZ;

        // val is the deposited energy reported by the Geant4 scorer.
        // val2 is used here as an occupancy/error indicator; empty entries
        // are not inserted into the sparse histogram.
        if (val2 > 0.0)
        {
            Double_t pos[3] = {x, y, z};
            h1->Fill(pos, val);
            ++nFilled;
        }
    }

    inputFile.close();

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

    std::cout << "Read " << nLines << " non-comment lines; filled "
              << nFilled << " mesh bins." << std::endl;
    std::cout << "Output file: " << fOutName << std::endl;

    delete file;
    delete h1;
}
