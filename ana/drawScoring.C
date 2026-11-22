/*
 * drawScoring.C
 *
 * Thermal analysis of the Geant4 deposited-energy mesh.
 *
 * Usage (from the results/ directory, see README):
 *   root -l -b -q '../ana/drawScoring.C("edep.root","Al",50)'
 *   root -l -b -q '../ana/drawScoring.C("edep.root","Al",50,1000,0.1,10,1,1)'
 *
 * Arguments (all optional; the defaults reproduce the README demonstration)
 *   fNameIn     ROOT file written by readScoring.C (THnSparse "h1")
 *   collMat     collimator material; only "Al" is implemented
 *   collLength  collimator tip length L1 [mm]          -> plot titles, outline
 *   nSim        number of simulated primary electrons  -> normalization;
 *               MUST equal /run/beamOn of the Geant4 run
 *   collTanTh   tan(theta) of the collimator taper     -> outline only
 *   beamEnergy  beam energy [GeV]                       -> plot titles only
 *   beamSigmaX  transverse beam size sigma_x [mm]       -> diffusion factor
 *   beamSigmaY  transverse beam size sigma_y [mm]       -> diffusion factor
 * Use the same values as on the collimator_temp_rise command line.
 *
 * Workflow
 *   1. Read the THnSparse "h1" produced by readScoring.C: energy deposited
 *      per mesh cell [MeV], summed over nSim primaries, with the per-cell
 *      sum of squares stored as the bin error^2.
 *   2. Apply the transverse thermal-diffusion factor.
 *   3. Divide by the cell volume -> energy density [MeV/cm^3].
 *   4. Normalize to the number of electrons delivered in one revolution
 *      period at beam current I, N_e = I / (e f_rev) -> [J/cm^3].
 *   5. Convert the energy density into an instantaneous temperature rise
 *      (temperature-dependent c_p and melting; see ThermalModel.h).
 *   6. Draw DeltaT maps/profiles and compute the volume above the
 *      thermal-shock and melting thresholds as a function of beam current.
 *
 * Outputs (written to the current directory; <mat> = collMat)
 *   cYZ_<mat>.png    DeltaT map in the Y-Z plane (two central X bins), Ibeam
 *   c3_<mat>.png     maximum DeltaT per Z slice for Ibeam, Ibeam1 and Ibeam2
 *   c4_<mat>.png     DeltaT(E_dep) curve used for cYZ (range of Ibeam)
 *   c5_<mat>.png     thermal-shock volume and melt volume (melting started) vs current
 *   c6_<mat>.png     DeltaT(E_dep) curve over the range of the current scan
 *   grShockVolume_<mat>.root, grMeltVolume_<mat>.root   the graphs of c5
 *
 * IMPORTANT
 *   The beam currents and the revolution frequency (USER SETTINGS below) are
 *   demonstration values, not machine parameters.  Change them for a
 *   physical study.
 */
#include "ThermalModel.h"

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <THnSparse.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPolyLine.h>
#include <TSpline.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

/// Tabulate DeltaT(E_dep) on [0, Emax] with nPoints equidistant points.
/// The graph is used as an interpolation table (TSpline3) and for plotting.
TGraph* buildDeltaTGraph(double Emax /*[J/cm^3]*/, int nPoints,
                         const MaterialModel& mat, const TString& collMat,
                         const TString& name)
{
    auto* gr = new TGraph(nPoints);
    gr->SetName(name);
    gr->SetTitle(Form("%s: #DeltaT vs E_{dep};E_{dep} [J/cm^{3}];#DeltaT [^{o}C]",
                      collMat.Data()));
    for (int i = 0; i < nPoints; ++i)
    {
        const double Edep = Emax * i / (nPoints - 1);   // [J/cm^3]
        gr->SetPoint(i, Edep, computeDeltaT(Edep, mat.props, mat.cp).deltaT);
    }
    return gr;
}

/// Draw the outline of the tapered collimator in the Z-Y plane [cm]: tip
/// length L1 at y = +H/2 and base length L2 = L1 + 2 H / tan(theta) at
/// y = -H/2 (same geometry as DetectorConstruction.cc).
void drawCollEdge(double L1 = 5.0 /*[cm]*/, double H = 3.0 /*[cm]*/, double tanTh = 0.1)
{
    H = std::abs(H);
    const int n = 5;
    const double L2 = L1 + 2.0 * H / tanTh;
    double z[n] = {-L2 / 2., -L1 / 2., L1 / 2., L2 / 2., -L2 / 2.};
    double y[n] = {-H / 2.,   H / 2.,   H / 2., -H / 2., -H / 2.};

    auto* coll = new TPolyLine(n, z, y);
    coll->SetLineColor(kBlack);
    coll->SetLineWidth(3);
    coll->SetFillStyle(0);
    coll->Draw("SAME");
}

/// Largest bin content of a THnSparse.
double getMaxTHnSparse(const THnSparse* h)
{
    double maxVal = 0.0;
    const Long64_t nbins = h->GetNbins();
    for (Long64_t i = 0; i < nbins; ++i)
        maxVal = std::max(maxVal, h->GetBinContent(i));
    return maxVal;
}

/// readScoring.C stores, for every cell, the sum over events of the squared
/// per-event deposit, S2 = sum_i x_i^2, as the bin error^2.  Replace it by the
/// variance of the cell sum S1 = sum_i x_i for a fixed number N of primaries:
///   Var(S1) = N / (N - 1) * (S2 - S1^2 / N).
/// Returns false if the histogram carries no per-bin error information.
bool setStatisticalErrors(THnSparse* h, Long64_t nSim)
{
    if (!h->GetCalculateErrors()) return false;
    if (nSim < 2) return true;   // variance not estimable; keep S2

    const Long64_t nbins = h->GetNbins();
    for (Long64_t i = 0; i < nbins; ++i)
    {
        const double s1  = h->GetBinContent(i);
        const double s2  = h->GetBinError2(i);
        const double var = (s2 - s1 * s1 / nSim) * nSim / (nSim - 1.0);
        h->SetBinError2(i, std::max(var, 0.0));
    }
    return true;
}

void drawScoring(
    TString  fNameIn    = "edep.root",
    TString  collMat    = "Al",
    Double_t collLength = 50,    // [mm] collimator tip length L1
    Long64_t nSim       = 1000,  // number of primaries (= /run/beamOn)
    Double_t collTanTh  = 0.1,   // tan(theta) of the taper
    Double_t beamEnergy = 10,    // [GeV]
    Double_t beamSigmaX = 1.0,   // [mm]
    Double_t beamSigmaY = 1.0    // [mm]
)
{
    using std::cout;
    using std::cerr;
    using std::endl;

    // =======================================================================
    // USER SETTINGS (demonstration values)
    // =======================================================================
    const double Ibeam    = 2e-3;    // [A] beam current for cYZ and c3
    const double Ibeam1   = 20e-3;   // [A] additional current shown in c3
    const double Ibeam2   = 150e-3;  // [A] additional current shown in c3
    const double revFreq  = 78e3;    // [Hz] revolution frequency
    const double IbeamMin = 1e-3;    // [A] lower end of the current scan (c5)
    const double IbeamMax = 2.5;     // [A] upper end of the current scan (c5)
    const int    nIbeam   = 50;      // number of currents in the (log) scan
    const int    nPointsDT = 1000;   // points of the DeltaT(E_dep) tables
    const double collH    = 3.0;     // [cm] collimator height H (DetectorConstruction.cc)
    const double minZ = -40, maxZ = 40;  // [cm] longitudinal range of the plots
    // =======================================================================

    cout << "drawScoring: input = " << fNameIn << ", material = " << collMat
         << ", L1 = " << collLength << " mm, nSim = " << nSim
         << ", tan(theta) = " << collTanTh << ", E = " << beamEnergy
         << " GeV, sigma_x/y = " << beamSigmaX << "/" << beamSigmaY << " mm"
         << endl;

    // ---- Input checks -----------------------------------------------------
    MaterialModel mat;
    if (!getMaterialModel(collMat, mat)) return;
    if (nSim <= 0)
    {
        cerr << "ERROR: nSim must be positive (number of simulated primaries)." << endl;
        return;
    }
    if (collTanTh <= 0 || beamSigmaX <= 0 || beamSigmaY <= 0)
    {
        cerr << "ERROR: collTanTh, beamSigmaX and beamSigmaY must be positive." << endl;
        return;
    }

    TFile* fIn = TFile::Open(fNameIn);
    if (!fIn || fIn->IsZombie())
    {
        cerr << "ERROR: cannot open " << fNameIn << endl;
        return;
    }
    auto* h = dynamic_cast<THnSparse*>(fIn->Get("h1"));
    if (!h || h->GetNdimensions() != 3)
    {
        cerr << "ERROR: " << fNameIn
             << " does not contain the 3-D THnSparse \"h1\" written by readScoring.C." << endl;
        return;
    }

    // ---- Statistical uncertainties -----------------------------------------
    if (!setStatisticalErrors(h, nSim))
    {
        cout << "WARNING: " << fNameIn << " has no per-cell variance (it was written by an\n"
             << "         older readScoring.C).  Statistical uncertainties are set to zero;\n"
             << "         re-run readScoring.C to obtain them." << endl;
        h->Sumw2();
        for (Long64_t i = 0; i < h->GetNbins(); ++i) h->SetBinError2(i, 0.0);
    }

    // Mesh cell size [cm] (uniform binning).
    const double binWidthX = h->GetAxis(0)->GetBinWidth(1);
    const double binWidthY = h->GetAxis(1)->GetBinWidth(1);
    const double binWidthZ = h->GetAxis(2)->GetBinWidth(1);
    cout << "Mesh cell size: " << binWidthX << " x " << binWidthY << " x "
         << binWidthZ << " cm^3" << endl;

    // ---- Transverse thermal-diffusion factor -------------------------------
    // During the characteristic time t = 1/f_rev (one revolution period) heat
    // diffuses over sigma_diff = sqrt(a t), with the thermal diffusivity
    // a = k / (rho c_p(T0)).  The peak energy density of a Gaussian beam with
    // sizes sigma_x, sigma_y is then reduced by the factor
    //   sigma_x sigma_y / sqrt((sigma_x^2 + sigma_diff^2)(sigma_y^2 + sigma_diff^2)),
    // which is applied uniformly to all mesh cells.  This phenomenological
    // correction is an analysis assumption; remove or replace it if the
    // analysis setup changes.  Reference retained from the original study:
    // https://inspirehep.net/files/f5414415b1051da61d51f427adb336be
    const double mat_k   = mat.props.k;       // [W/m/K]
    const double mat_rho = mat.props.rho;     // [g/cm^3]
    const double mat_cp  = mat.cp(T0);        // [J/g/K]
    const double sigma_x = beamSigmaX * 1e-3; // [m]
    const double sigma_y = beamSigmaY * 1e-3; // [m]
    const double t_diff  = 1.0 / revFreq;     // [s]
    const double a_diff  = mat_k / (mat_rho * 1e6 * mat_cp);   // [m^2/s]  (g/cm^3 -> g/m^3)
    const double sigma_diff = std::sqrt(a_diff * t_diff);       // [m]
    const double diffusionFactor =
        sigma_x * sigma_y /
        std::sqrt((sigma_x * sigma_x + sigma_diff * sigma_diff) *
                  (sigma_y * sigma_y + sigma_diff * sigma_diff));

    cout << "mat_rho = " << mat_rho << " [g/cm^3]" << endl;
    cout << "mat_cp(T0) = " << mat_cp << " [J/g/K]" << endl;
    cout << "mat_k = " << mat_k << " [W/m/K]" << endl;
    cout << "a_diff = " << a_diff << " [m^2/s]" << endl;
    cout << "sigma_diff = " << sigma_diff << " [m]" << endl;
    cout << "diffusionFactor = " << diffusionFactor << endl;

    h->Scale(diffusionFactor);

    // The scorer stores the energy deposited per mesh cell [MeV]; divide by
    // the cell volume to obtain an energy density [MeV/cm^3].
    h->Scale(1.0 / (binWidthX * binWidthY * binWidthZ));

    // ---- Normalization to the beam current ---------------------------------
    // The energy density is normalized to the number of electrons delivered
    // in one revolution period at beam current I (equivalently, the whole
    // stored beam lost in a single turn):  N_e = I / (e f_rev).
    //   E[J/cm^3] = E_MC[MeV/cm^3] / nSim * N_e * (1 MeV in J)
    //             = E_MC[MeV/cm^3] * 1e6 * I / (f_rev * nSim),
    // because (1 MeV in J) / e = 1e6 V exactly.
    auto* h1 = static_cast<THnSparse*>(h->Clone("hIbeam1"));
    auto* h2 = static_cast<THnSparse*>(h->Clone("hIbeam2"));
    h->Scale(1e6 * Ibeam / (revFreq * nSim));    // [J/cm^3] at Ibeam
    h1->Scale(1e6 * Ibeam1 / (revFreq * nSim));  // [J/cm^3] at Ibeam1
    h2->Scale(1e6 * Ibeam2 / (revFreq * nSim));  // [J/cm^3] at Ibeam2

    const double maxEdepDensity = getMaxTHnSparse(h);   // [J/cm^3] at Ibeam
    if (maxEdepDensity <= 0)
    {
        cerr << "ERROR: no energy deposited in the scoring mesh." << endl;
        return;
    }
    cout << "Max E_dep density at " << Ibeam * 1e3 << " mA = " << maxEdepDensity
         << " J/cm^3" << endl;

    // ---- Y-Z map: average of the two central X bins ------------------------
    const int nBinsX = h->GetAxis(0)->GetNbins();
    const int binCenterX1 = nBinsX / 2;
    const int binCenterX2 = binCenterX1 + 1;
    cout << "Y-Z map: average over X bins " << binCenterX1 << "-" << binCenterX2
         << " (" << h->GetAxis(0)->GetBinLowEdge(binCenterX1) << " < x < "
         << h->GetAxis(0)->GetBinUpEdge(binCenterX2) << " cm)" << endl;
    h->GetAxis(0)->SetRange(binCenterX1, binCenterX2);   // X-axis bin range cut
    TH2D* hYZ = h->Projection(1, 2);   // Projection(yDim, xDim): y = Y, x = Z
    hYZ->SetName("hYZ");
    h->GetAxis(0)->SetRange(0, 0);     // reset to the full range
    hYZ->Scale(0.5);                   // average of the two X bins

    // ---- DeltaT(E_dep) interpolation tables ---------------------------------
    TGraph* grDT  = buildDeltaTGraph(maxEdepDensity, nPointsDT, mat, collMat, "grDeltaT");
    TGraph* grDT1 = buildDeltaTGraph(getMaxTHnSparse(h1), nPointsDT, mat, collMat, "grDeltaT1");
    TGraph* grDT2 = buildDeltaTGraph(getMaxTHnSparse(h2), nPointsDT, mat, collMat, "grDeltaT2");

    auto* splDT  = new TSpline3("splDT", grDT);
    auto* splDT1 = new TSpline3("splDT1", grDT1);
    auto* splDT2 = new TSpline3("splDT2", grDT2);

    // ---- Plot style ----------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.07);
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t red[NRGBs]   = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t blue[NRGBs]  = {0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // ---- cYZ: DeltaT map at Ibeam -------------------------------------------
    auto* cYZ = new TCanvas("cYZ", "cYZ", 1800, 900);
    cYZ->cd();
    auto* pYZ = new TPad("pYZ", "", 0, 0, 1, 1);
    pYZ->SetTopMargin(0.1);
    pYZ->SetBottomMargin(0.1);
    pYZ->SetLeftMargin(0.1);
    pYZ->SetRightMargin(0.2);
    pYZ->Draw();
    pYZ->cd();

    for (int ix = 1; ix <= hYZ->GetNbinsX(); ++ix)
    {
        for (int iy = 1; iy <= hYZ->GetNbinsY(); ++iy)
        {
            const double Edep = hYZ->GetBinContent(ix, iy);   // [J/cm^3]
            hYZ->SetBinContent(ix, iy, splDT->Eval(Edep));    // [K]
        }
    }

    hYZ->SetTitle(Form("%g GeV e^{-} on %g mm %s @ %g mA / %g kHz",
                       beamEnergy, collLength, collMat.Data(), Ibeam * 1e3, revFreq * 1e-3));
    hYZ->GetXaxis()->SetTitle("Z [cm]");
    hYZ->GetYaxis()->SetTitle("Y [cm]");
    hYZ->GetZaxis()->SetTitle("#DeltaT [^{o}C]");
    hYZ->GetXaxis()->SetTitleSize(0.05);
    hYZ->GetYaxis()->SetTitleSize(0.05);
    hYZ->GetZaxis()->SetTitleSize(0.05);
    hYZ->GetXaxis()->SetLabelSize(0.05);
    hYZ->GetYaxis()->SetLabelSize(0.05);
    hYZ->GetZaxis()->SetLabelSize(0.05);
    hYZ->GetYaxis()->SetTitleOffset(1.0);
    hYZ->GetXaxis()->SetTitleOffset(1.0);
    hYZ->GetZaxis()->SetTitleOffset(1.0);
    hYZ->GetZaxis()->SetNdivisions(10);
    hYZ->SetMinimum(1e-2);
    hYZ->Draw("colz");
    gPad->SetLogz();
    gPad->SetGrid();
    gPad->Update();
    drawCollEdge(collLength * 0.1, collH, collTanTh);
    cYZ->SaveAs(Form("cYZ_%s.png", collMat.Data()));

    // ---- Maximum energy density per Z slice ----------------------------------
    const int iz = 2;   // Z axis index
    TAxis* zAxis = h->GetAxis(iz);
    const int nBinsZ = zAxis->GetNbins();
    std::vector<double> maxEdep(nBinsZ, 0.0), maxEdep1(nBinsZ, 0.0), maxEdep2(nBinsZ, 0.0);
    std::vector<double> maxRelEdepError(nBinsZ, 0.0), maxRelEdepError1(nBinsZ, 0.0),
                        maxRelEdepError2(nBinsZ, 0.0);

    // h, h1 and h2 are clones, so the bin index i refers to the same cell.
    const Long64_t nbins = h->GetNbins();
    Int_t coord[3];
    for (Long64_t i = 0; i < nbins; ++i)
    {
        const double edep = h->GetBinContent(i, coord);
        const int bz = coord[iz] - 1;   // ROOT bins are 1-based
        if (bz < 0 || bz >= nBinsZ) continue;   // under-/overflow (not expected)

        const double edep1 = h1->GetBinContent(i);
        const double edep2 = h2->GetBinContent(i);
        if (edep > maxEdep[bz])   { maxEdep[bz]  = edep;  maxRelEdepError[bz]  = h->GetBinError(i) / edep; }
        if (edep1 > maxEdep1[bz]) { maxEdep1[bz] = edep1; maxRelEdepError1[bz] = h1->GetBinError(i) / edep1; }
        if (edep2 > maxEdep2[bz]) { maxEdep2[bz] = edep2; maxRelEdepError2[bz] = h2->GetBinError(i) / edep2; }
    }

    // DeltaT of the hottest cell in each slice.  The uncertainty is the
    // relative statistical uncertainty of that cell's energy density.
    auto* grMaxTempZ  = new TGraphErrors();
    auto* grMaxTempZ1 = new TGraphErrors();
    auto* grMaxTempZ2 = new TGraphErrors();
    double maxVal = 0, maxValErr = 0;
    double maxVal1 = 0, maxValErr1 = 0;
    double maxVal2 = 0, maxValErr2 = 0;

    for (int bz = 0; bz < nBinsZ; ++bz)
    {
        const double z = zAxis->GetBinCenter(bz + 1);

        const double val  = splDT->Eval(maxEdep[bz]);
        const double val1 = splDT1->Eval(maxEdep1[bz]);
        const double val2 = splDT2->Eval(maxEdep2[bz]);

        const double err  = val  * maxRelEdepError[bz];
        const double err1 = val1 * maxRelEdepError1[bz];
        const double err2 = val2 * maxRelEdepError2[bz];

        grMaxTempZ->SetPoint(bz, z, val);   grMaxTempZ->SetPointError(bz, 0, err);
        grMaxTempZ1->SetPoint(bz, z, val1); grMaxTempZ1->SetPointError(bz, 0, err1);
        grMaxTempZ2->SetPoint(bz, z, val2); grMaxTempZ2->SetPointError(bz, 0, err2);

        if (maxVal < val)   { maxVal  = val;  maxValErr  = err; }
        if (maxVal1 < val1) { maxVal1 = val1; maxValErr1 = err1; }
        if (maxVal2 < val2) { maxVal2 = val2; maxValErr2 = err2; }
    }
    cout << "dT_max(" << Ibeam  * 1e3 << " mA) = " << maxVal  << " +/- " << maxValErr  << " K" << endl;
    cout << "dT_max(" << Ibeam1 * 1e3 << " mA) = " << maxVal1 << " +/- " << maxValErr1 << " K" << endl;
    cout << "dT_max(" << Ibeam2 * 1e3 << " mA) = " << maxVal2 << " +/- " << maxValErr2 << " K" << endl;

    // ---- c3: maximum DeltaT along Z -------------------------------------------
    auto* c3 = new TCanvas("c3", "c3", 900, 900);
    c3->cd();
    auto* p3 = new TPad("p3", "", 0, 0, 1, 1);
    p3->SetTopMargin(0.1);
    p3->SetBottomMargin(0.12);
    p3->SetLeftMargin(0.20);
    p3->SetRightMargin(0.05);
    p3->Draw();
    p3->cd();

    grMaxTempZ->SetMinimum(0);
    grMaxTempZ->SetMaximum(maxVal2 * 1.1);
    grMaxTempZ->GetXaxis()->SetRangeUser(minZ, maxZ);
    grMaxTempZ->SetTitle(Form("%g GeV e^{-} on %g mm %s", beamEnergy, collLength, collMat.Data()));
    grMaxTempZ->GetXaxis()->SetTitle("Z [cm]");
    grMaxTempZ->GetYaxis()->SetTitle("#DeltaT [^{o}C]");
    grMaxTempZ->GetXaxis()->SetTitleSize(0.05);
    grMaxTempZ->GetYaxis()->SetTitleSize(0.05);
    grMaxTempZ->GetXaxis()->SetLabelSize(0.05);
    grMaxTempZ->GetYaxis()->SetLabelSize(0.05);
    grMaxTempZ->GetYaxis()->SetTitleOffset(2.0);
    grMaxTempZ->GetXaxis()->SetTitleOffset(1.0);

    grMaxTempZ->SetLineWidth(3);
    grMaxTempZ1->SetLineWidth(3);
    grMaxTempZ2->SetLineWidth(3);

    grMaxTempZ->SetLineColor(kBlue);
    grMaxTempZ1->SetLineColor(kRed);
    grMaxTempZ2->SetLineColor(kBlack);

    grMaxTempZ->Draw("AL");
    grMaxTempZ1->Draw("SAME L");
    grMaxTempZ2->Draw("SAME L");

    // Reference lines: practical thermal-shock threshold and the temperature
    // interval from T0 to the melting point.
    const double dT_shock = mat.props.dT_shock;   // [K]
    const double Tm_T0    = mat.props.Tm - T0;    // [K]

    auto* shock_line = new TLine(minZ, dT_shock, maxZ, dT_shock);
    shock_line->SetLineWidth(4);
    shock_line->SetLineStyle(7);
    shock_line->SetLineColor(kMagenta);
    shock_line->Draw("SAME L");

    auto* Tm_T0_line = new TLine(minZ, Tm_T0, maxZ, Tm_T0);
    Tm_T0_line->SetLineWidth(3);
    Tm_T0_line->SetLineStyle(9);
    Tm_T0_line->SetLineColor(kOrange);
    Tm_T0_line->Draw("SAME L");

    gPad->SetGrid();

    auto* leg = new TLegend(0.2, 0.45, 0.5, 0.85);
    leg->AddEntry(grMaxTempZ, Form("%g mA", Ibeam * 1e3), "l");
    leg->AddEntry(grMaxTempZ1, Form("%g mA", Ibeam1 * 1e3), "l");
    leg->AddEntry(grMaxTempZ2, Form("%g mA", Ibeam2 * 1e3), "l");
    leg->AddEntry(shock_line, Form("#DeltaT_{shock} = %.f ^{o}C", dT_shock), "l");
    leg->AddEntry(Tm_T0_line, Form("T_{m}-T_{0} = %.f ^{o}C", Tm_T0), "l");
    leg->Draw();

    c3->SaveAs(Form("c3_%s.png", collMat.Data()));

    // ---- c4: DeltaT(E_dep) over the range of Ibeam ----------------------------
    auto* c4 = new TCanvas("c4", "c4", 900, 900);
    c4->cd();
    auto* p4 = new TPad("p4", "", 0, 0, 1, 1);
    p4->SetTopMargin(0.1);
    p4->SetBottomMargin(0.12);
    p4->SetLeftMargin(0.16);
    p4->SetRightMargin(0.01);
    p4->Draw();
    p4->cd();
    grDT->GetXaxis()->SetTitleSize(0.05);
    grDT->GetYaxis()->SetTitleSize(0.05);
    grDT->GetXaxis()->SetLabelSize(0.05);
    grDT->GetYaxis()->SetLabelSize(0.05);
    grDT->GetYaxis()->SetTitleOffset(1.6);
    grDT->GetXaxis()->SetTitleOffset(1.0);
    grDT->GetXaxis()->SetNdivisions(5);
    grDT->SetLineWidth(2);
    grDT->Draw("APL");

    c4->SaveAs(Form("c4_%s.png", collMat.Data()));

    // =======================================================================
    // Beam-current scan: volume above the thermal-shock and melting thresholds
    // =======================================================================
    // The scan is logarithmic to show the onset and the rapid growth of the
    // affected volume over several decades of beam current.
    //
    // DeltaT(E_dep) is monotonic, so both conditions are evaluated exactly
    // as "E_dep > E_threshold", where E_threshold follows from the inverse of
    // the DeltaT(E_dep) relation (no interpolation involved):
    //   shock volume : DeltaT > DeltaT_shock
    //   melt volume  : melting has started, i.e. the cell has reached Tm and
    //                  is partially or completely molten (E_dep above the
    //                  energy density needed to heat the material to Tm).
    const double eShock = energyDensityForDeltaT(dT_shock, mat.props, mat.cp);      // [J/cm^3]
    const double eMelt  = energyDensityForMeltingOnset(mat.props, mat.cp);          // [J/cm^3]
    const double eMeltComplete = energyDensityForDeltaT(Tm_T0, mat.props, mat.cp);  // [J/cm^3]
    cout << "Energy-density thresholds: DeltaT > " << dT_shock << " K for E_dep > "
         << eShock << " J/cm^3; melting starts for E_dep > " << eMelt
         << " J/cm^3 (completely molten above " << eMeltComplete << " J/cm^3)" << endl;

    // DeltaT(E_dep) over the full range of the scan (plotted in c6).
    TGraph* grDT_shock = buildDeltaTGraph(maxEdepDensity * IbeamMax / Ibeam, nPointsDT,
                                          mat, collMat, "grDeltaT_scan");

    std::vector<double> beamCurrent(nIbeam), IbeamArray(nIbeam);
    const double logIbeamMin = std::log10(IbeamMin);
    const double logIbeamMax = std::log10(IbeamMax);
    for (int j = 0; j < nIbeam; ++j)
    {
        IbeamArray[j]  = std::pow(10.0, logIbeamMin + j * (logIbeamMax - logIbeamMin) / (nIbeam - 1));
        beamCurrent[j] = IbeamArray[j] * 1e3;   // A -> mA
    }
    cout << "Current scan: " << IbeamArray.front() << " - " << IbeamArray.back()
         << " A (" << nIbeam << " points)" << endl;

    std::vector<double> shockVolume(nIbeam, 0.0);   // [mm^3]
    std::vector<double> meltVolume(nIbeam, 0.0);    // [mm^3]

    for (Long64_t i = 0; i < nbins; ++i)
    {
        const double edep = h->GetBinContent(i, coord);   // [J/cm^3] at Ibeam
        if (edep <= 0) continue;
        const double dV = h->GetAxis(0)->GetBinWidth(coord[0]) *
                          h->GetAxis(1)->GetBinWidth(coord[1]) *
                          h->GetAxis(2)->GetBinWidth(coord[2]);   // [cm^3]

        for (int j = 0; j < nIbeam; ++j)
        {
            const double edep_j = edep * IbeamArray[j] / Ibeam;   // scale to current j
            if (edep_j > eShock) shockVolume[j] += dV * 1e3;     // cm^3 -> mm^3
            if (edep_j > eMelt)  meltVolume[j]  += dV * 1e3;     // cm^3 -> mm^3
        }
    }

    // ---- c5: threshold volumes vs beam current ---------------------------------
    auto* c5 = new TCanvas("c5", "c5", 900, 900);
    c5->cd();
    auto* p5 = new TPad("p5", "", 0, 0, 1, 1);
    p5->SetTopMargin(0.1);
    p5->SetBottomMargin(0.12);
    p5->SetLeftMargin(0.22);
    p5->SetRightMargin(0.01);
    p5->Draw();
    p5->cd();

    auto* grShockVolume = new TGraph(nIbeam, beamCurrent.data(), shockVolume.data());
    auto* grMeltVolume  = new TGraph(nIbeam, beamCurrent.data(), meltVolume.data());

    grShockVolume->SetTitle(Form("%g GeV e^{-}", beamEnergy));
    grShockVolume->GetXaxis()->SetTitleSize(0.05);
    grShockVolume->GetYaxis()->SetTitleSize(0.05);
    grShockVolume->GetXaxis()->SetLabelSize(0.05);
    grShockVolume->GetYaxis()->SetLabelSize(0.05);
    grShockVolume->GetYaxis()->SetTitleOffset(2.0);
    grShockVolume->GetXaxis()->SetTitleOffset(1.0);
    grShockVolume->GetXaxis()->SetNdivisions(5);
    grShockVolume->GetXaxis()->SetTitle("I_{beam} [mA]");
    grShockVolume->GetYaxis()->SetTitle("V [mm^{3}]");

    grShockVolume->SetLineColor(kRed);
    grShockVolume->SetLineWidth(3);
    grShockVolume->Draw("AL");

    grMeltVolume->SetLineColor(kBlue);
    grMeltVolume->SetLineWidth(3);
    grMeltVolume->SetLineStyle(7);
    grMeltVolume->Draw("SAME L");

    auto* leg5 = new TLegend(0.25, 0.55, 0.9, 0.8);
    leg5->AddEntry(grShockVolume, "#DeltaT > #DeltaT^{practical}_{shock}", "pl");
    leg5->AddEntry(grMeltVolume, "#DeltaT #geq T_{m}-T_{0} (melting)", "pl");
    leg5->Draw();

    gPad->SetGrid();
    gPad->SetLogx();

    c5->SaveAs(Form("c5_%s.png", collMat.Data()));

    grShockVolume->SetName("grShockVolume");
    grShockVolume->SaveAs(Form("grShockVolume_%s.root", collMat.Data()));
    grMeltVolume->SetName("grMeltVolume");
    grMeltVolume->SaveAs(Form("grMeltVolume_%s.root", collMat.Data()));

    // ---- c6: DeltaT(E_dep) over the range of the current scan -----------------
    auto* c6 = new TCanvas("c6", "c6", 900, 900);
    c6->cd();
    auto* p6 = new TPad("p6", "", 0, 0, 1, 1);
    p6->SetTopMargin(0.1);
    p6->SetBottomMargin(0.12);
    p6->SetLeftMargin(0.16);
    p6->SetRightMargin(0.01);
    p6->Draw();
    p6->cd();
    grDT_shock->GetXaxis()->SetTitleSize(0.05);
    grDT_shock->GetYaxis()->SetTitleSize(0.05);
    grDT_shock->GetXaxis()->SetLabelSize(0.05);
    grDT_shock->GetYaxis()->SetLabelSize(0.05);
    grDT_shock->GetYaxis()->SetTitleOffset(1.6);
    grDT_shock->GetXaxis()->SetTitleOffset(1.0);
    grDT_shock->GetXaxis()->SetNdivisions(5);
    grDT_shock->SetLineWidth(2);
    grDT_shock->Draw("APL");

    c6->SaveAs(Form("c6_%s.png", collMat.Data()));

    // ---- Summary: maximum and total deposited energy (both scale with I) -----
    double totEdep = 0.0;   // [J] at Ibeam
    for (Long64_t i = 0; i < nbins; ++i)
    {
        const double edep = h->GetBinContent(i, coord);
        totEdep += edep * h->GetAxis(0)->GetBinWidth(coord[0]) *
                          h->GetAxis(1)->GetBinWidth(coord[1]) *
                          h->GetAxis(2)->GetBinWidth(coord[2]);
    }
    for (double I : {0.1, 0.2})   // [A]
    {
        cout << "I_beam = " << I * 1e3 << " mA --> Max E_dep = "
             << maxEdepDensity * I / Ibeam << " J/cm^3" << endl;
    }
    for (double I : {0.1, 0.2})   // [A]
    {
        cout << "I_beam = " << I * 1e3 << " mA --> Total E_dep = "
             << totEdep * I / Ibeam << " J" << endl;
    }
}
