/*
 * Plot the temperature-dependent specific heat capacity used by drawScoring.C.
 *
 * Data source: NIST Chemistry WebBook, aluminum condensed-phase data.
 * cp_al(T) returns J/g/K.
 */
// Calculate molecular heat capacity
double cp_mol(const std::vector<double>& par, double t)
{
	double tt = t / 1000.0;
	double tt2 = tt * tt;
	double tt3 = tt2 * tt;
	return par[0] + par[1] * tt + par[2] * tt2 + par[3] * tt3 + par[4] / tt2;
}

//--- Heat capacity of [J/g*K] as a function of temperature [K]
// Al: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7429905&Units=SI&Mask=2#Thermo-Condensed
double cp_al(double t)
{
	std::vector<double> par;

	// Solid Phase Heat Capacity
	if (t < 933.0) {par = {28.08920,-5.414849,8.560423,3.427370,-0.277375};}
	// Liquid Phase Heat Capacity
	else {par = {31.75104,3.935826e-8,-1.786515e-8,2.694171e-9,5.480037e-9};}

	double cp = cp_mol(par, t);       // [J/mol/K]
	const double mol_mass = 26.98;    // [g/mol]
	return cp / mol_mass;             // [J/g/K]
}

// Generate and draw heat capacity graphs
void getCp()
{
	const int nPoints = 70;
	const double t0 = 273.0;
	const double dt = 60.0;

	TGraph* gr = new TGraph();

	for (int j = 0; j < nPoints; j++) 
	{
		double t = t0 + j * dt;
		double cp = cp_al(t);

		gr->SetPoint(j, t, cp);
	}

	TCanvas* c1 = new TCanvas("c1","c1",1200,600);
	c1->cd();
	TPad* p1 = new TPad("p1","",0,0,1,1);
   	p1->SetTopMargin(0.1);
	p1->SetBottomMargin(0.12);
	p1->SetLeftMargin(0.16);
	p1->SetRightMargin(0.01);
	p1->Draw();
	p1->cd();

	gr->SetLineStyle(1);
	gr->SetLineColor(kRed);
	gr->SetLineWidth(3);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColorAlpha(kBlack,0.6);

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
	gr->GetXaxis()->SetTitle("T [K]");
	gr->GetYaxis()->SetTitle("c_{p} [J/g/K]");

	gr->Draw("APL");

	gPad->SetGrid();

	c1->SaveAs("c1_heat_capacity.png");

	return;
}

