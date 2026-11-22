/*
 * getCp.C
 *
 * Plot the temperature-dependent specific heat capacity of aluminum used by
 * drawScoring.C (cp_al() in ThermalModel.h).
 *
 * Usage (from any directory):
 *   root -l -b -q '../ana/getCp.C'
 *
 * Data source: NIST Chemistry WebBook, aluminum, condensed-phase
 * thermochemistry (Shomate equation).  The curve is drawn over the validity
 * range of the NIST parameterization (298-2790.8 K); the step at the melting
 * point is the solid -> liquid transition.
 *
 * Output: c1_heat_capacity.png
 */
#include "ThermalModel.h"

#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPad.h>

void getCp()
{
    const double tMin = 298.15;   // [K] lower validity limit of the solid fit
    const double tMax = 2790.0;   // [K] upper validity limit of the liquid fit
    const int nPoints = 500;

    auto* gr = new TGraph();
    for (int j = 0; j < nPoints; ++j)
    {
        const double t = tMin + j * (tMax - tMin) / (nPoints - 1);
        gr->SetPoint(j, t, cp_al(t));
    }

    auto* c1 = new TCanvas("c1", "c1", 1200, 600);
    c1->cd();
    auto* p1 = new TPad("p1", "", 0, 0, 1, 1);
    p1->SetTopMargin(0.1);
    p1->SetBottomMargin(0.12);
    p1->SetLeftMargin(0.16);
    p1->SetRightMargin(0.04);
    p1->Draw();
    p1->cd();

    gr->SetLineStyle(1);
    gr->SetLineColor(kRed);
    gr->SetLineWidth(3);

    gr->SetMinimum(0);
    gr->SetMaximum(3);

    gr->SetTitle("Heat Capacity for Al");
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetTitleOffset(1.6);
    gr->GetXaxis()->SetTitleOffset(1.0);
    gr->GetXaxis()->SetNdivisions(5);
    gr->GetXaxis()->SetLimits(0, 3000);
    gr->GetXaxis()->SetTitle("T [K]");
    gr->GetYaxis()->SetTitle("c_{p} [J/g/K]");

    gr->Draw("AL");

    gPad->SetGrid();

    c1->SaveAs("c1_heat_capacity.png");
}
