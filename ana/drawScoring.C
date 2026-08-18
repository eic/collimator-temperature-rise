/*
 * Thermal analysis of the Geant4 deposited-energy mesh.
 *
 * Workflow:
 *   1. Read THnSparseF "h1" produced by readScoring.C.
 *   2. Convert Monte Carlo energy density into a physical deposited-energy
 *      density for a chosen beam current and revolution frequency.
 *   3. Apply the transverse thermal-diffusion scaling used in the study.
 *   4. Convert deposited energy density [J/cm^3] to instantaneous DeltaT.
 *      The temperature calculation uses temperature-dependent aluminum
 *      heat capacity and includes the melting transition.
 *   5. Produce maps/profiles of DeltaT.
 *   6. Estimate the volume above the practical thermal-shock threshold and
 *      the melting threshold as a function of beam current.
 *
 * IMPORTANT:
 *   The beam-current and revolution-frequency values below are demonstration
 *   values, not universal machine parameters.  Change them for a physical
 *   study.  Likewise, nSim must match the number of primary electrons used
 *   to generate the input scoring file.
 */
// Calculate molecular heat capacity
double cp_mol(const std::vector<double>& par, double t)
{
	double tt = t / 1000.0;
	double tt2 = tt * tt;
	double tt3 = tt2 * tt;
	return par[0] + par[1] * tt + par[2] * tt2 + par[3] * tt3 + par[4] / tt2;
}

const double T0 = 300.0;

// Generic material definition
struct MaterialProps
{
	double rho;          // g/cm^3
	double Tm;           // K (melting or sublimation)
	double Lf;           // J/cm^3 (0 if none)
	double dT_shock;     // K
	double k;	     // W/m/K
};

// Generic ΔT solver (single engine)
struct DeltaTResult
{
	double deltaT;
	bool melted;
	bool shockFailed;
};

DeltaTResult computeDeltaT(
    double Edep,                    // J/cm^3
    const MaterialProps& mat,
    std::function<double(double)> cp_mat
)
{
	double T = T0;
	double E = Edep;
	double dE = 0.1;

	bool melted = false;

	/* ---- solid heating ---- */
	while (E > 0 && T < mat.Tm)
	{
		T += dE / (mat.rho * cp_mat(T));
		E -= dE;
	}

	/* ---- phase change ---- */
	if (T >= mat.Tm)
	{
		melted = true;
		if (mat.Lf > 0.0) // melt
		{
			if (E < mat.Lf)
				return { mat.Tm - T0, true, (mat.Tm - T0) > mat.dT_shock };
			E -= mat.Lf;
		}
		else // sublimate
		{
			return { mat.Tm - T0, false, (mat.Tm - T0) > mat.dT_shock };
		}
	}

	/* ---- liquid heating ---- */
	while (E > 0)
	{
		T += dE / (mat.rho * cp_mat(T));
		E -= dE;
	}

	double dT = T - T0;

	return { dT, melted, dT > mat.dT_shock };
}

// Material instantiations (precise & explicit)
MaterialProps Al 	= { 2.70, 934.0, 1045.0,  74.0, 210};
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

TGraph* buildDeltaTGraph(double Emax /*[J/cm^3]*/, int nPoints = 4000, TString collMat = "Al")
{
	cout<<__FUNCTION__<<endl;
	auto gr = new TGraph(nPoints);
	gr->SetName(Form("grDeltaT_%s",collMat.Data()));
	gr->SetTitle(Form("%s: #DeltaT vs E_{dep};E_{dep} [J/cm^{3}];#DeltaT [^{o}C]",collMat.Data()));

	for(int i = 0; i < nPoints; ++i)
	{
		double Edep = Emax * i / (nPoints - 1);  // J/cm^3
		double dT   = 0; // °C
		if(collMat == "Al")	{auto res = computeDeltaT(Edep, Al, cp_al); dT = res.deltaT;}
		else{cout<<"Error: Wrong coll. material: "<<collMat<<endl;}
		gr->SetPoint(i, Edep, dT);
		if(i % 100 == 0){cout<<i<<"/"<<nPoints<<endl;}
	}

	return gr;
}

void drawCollEdge(double L1 = 5.0 /*[cm]*/, double H = 3.0 /*[cm]*/, double tanTh = 0.1)
{
	H = abs(H);
	const int n = 5;
	double L2 = L1 + 2 * H / tanTh;
	double z[] = { -L2/2., -L1/2., L1/2.,  L2/2., -L2/2.};
	double x[] = { -H/2.,   H/2.,  H/2.,  -H/2.,  -H/2.};

	TPolyLine* coll = new TPolyLine(n, z, x);
	coll->SetLineColor(kBlack);
	coll->SetLineWidth(3);
	coll->SetFillStyle(0);
	coll->Draw("SAME");

	return;
}

double getMaxTHnSparse(THnSparse* h)
{
	double maxVal = 0.0;

	Long64_t nbins = h->GetNbins();
	for (Long64_t i = 0; i < nbins; ++i)
	{
		double v = h->GetBinContent(i);
		if (v > maxVal)
		maxVal = v;
	}
	return maxVal;
}

void drawScoring(
	TString fNameIn = "edep.root",
	TString collMat = "Al", 
	Double_t collLength = 50 /*[mm]*/
)
{
	cout<<__FUNCTION__<<endl;
	// number of electrons per run
	Int_t nSim = 1000; // Must match /run/beamOn in the Geant4 input macro.
	cout<<"nSim = "<<nSim<<" [e-]"<<endl;

	const Float_t minZ 	= -40, maxZ = 40; // [cm] - longitudinal range for drawings
	const Float_t Ibeam  	= 2e-3; // [A]
	const Float_t Ibeam1  	= 20e-3; // [A] Demonstration beam current.
	const Float_t Ibeam2  	= 150e-3; // [A] Demonstration beam current.
	const Float_t revFreq 	= 78e3; // [Hz] Demonstration revolution frequency.

	TFile *fIn = TFile::Open(fNameIn.Data());
	THnSparseF* h = dynamic_cast<THnSparseF*>(fIn->Get("h1"));
	h->Sumw2();
	// identify the two central X bins
	double binWidthX = h->GetAxis(0)->GetBinWidth(1);
	double binWidthY = h->GetAxis(1)->GetBinWidth(1);
	double binWidthZ = h->GetAxis(2)->GetBinWidth(1);
	cout<<"binWidthX = "<<binWidthX<<"; binWidthY = "<<binWidthY<<"; binWidthZ = "<<binWidthZ<<endl;

	// Thermal-diffusion scaling.  The reference and physical interpretation
	// should be documented in the paper when this factor is used quantitatively.
	// Reference retained from the original study:
	// https://inspirehep.net/files/f5414415b1051da61d51f427adb336be
	Double_t mat_k = 0, mat_rho = 0, mat_cp = 0;
	if(collMat == "Al")	{mat_k = Al.k; mat_rho = Al.rho; mat_cp = cp_al(T0);}
	else{cout<<"ERROR:: Wrong material name: "<<collMat<<endl;return;}
	
	// Termal diffuse scaling
	Double_t sigma_x = 1.0*1e-3; // [m]
	Double_t sigma_y = 1.0*1e-3; // [m]
	Double_t t_diff = 1.0/revFreq; // [sec]
	Double_t a_diff	= mat_k / (mat_rho * (1e6) * mat_cp); // g/cm^3 -> g/m^3
	Double_t sigma_diff = TMath::Sqrt(a_diff * t_diff);
	Double_t diffusionFactor = sigma_x * sigma_y / (TMath::Sqrt(sigma_x*sigma_x + sigma_diff*sigma_diff)*TMath::Sqrt(sigma_y*sigma_y + sigma_diff*sigma_diff));

	cout<<"mat_rho = "<<mat_rho<<" [kg/cm^3]"<<endl;
	cout<<"mat_cp = "<<mat_cp<<" [J/g/K]"<<endl;
	cout<<"mat_k = "<<mat_k<<" [W/m/K]"<<endl;
	cout<<"a_diff = "<<a_diff<<" [m^2/sec]"<<endl;
	cout<<"sigma_diff = "<<sigma_diff<<" [m]"<<endl;
	cout<<"diffusionFactor = "<<diffusionFactor<<endl;

	// Apply the phenomenological transverse diffusion correction used in the
	// original study.  Remove/replace this step if the analysis setup changes.
	h->Scale(diffusionFactor);

	// The scorer stores deposited energy per mesh cell. Divide by the cell
	// volume to obtain an energy density before converting MeV to joules.
	h->Scale(1./(binWidthX*binWidthY*binWidthZ)); // MeV -> MeV/cm^3

	int ndim = h->GetNdimensions();
	THnSparseF* h1 = (THnSparseF*)h->Clone("h1");
	THnSparseF* h2 = (THnSparseF*)h->Clone("h2");

	// Scale the per-primary energy density to continuous-beam conditions:
	//   E_density * I / (e * f_rev * N_sim)
	// where 1e6 converts MeV to eV and the numerical expression below uses
	// the electron-charge/current relation embedded in the original analysis.
	h->Scale(1e6*Ibeam/(revFreq*nSim));
	h1->Scale(1e6*Ibeam1/(revFreq*nSim)); // MeV/cm^3 -> J/(cm^3)
	h2->Scale(1e6*Ibeam2/(revFreq*nSim)); // MeV/cm^3 -> J/(cm^3)

	int nBinsX = h->GetAxis(0)->GetNbins();
	int binCenterX1 = nBinsX/2; int binCenterX2 = binCenterX1 + 1;
	cout<<"binCenterX1 = "<<binCenterX1<<" - binCenterX2 = "<<binCenterX2<<endl;
	h->GetAxis(0)->SetRange(binCenterX1, binCenterX2); // apply X-axis bin range cut
	TH2F* hYZ = (TH2F*) h->Projection(1,2); hYZ->SetName("hYZ");
	h->GetAxis(0)->SetRange(0, 0); // reset to full range
	hYZ->Scale(1/2.); // since we projected over two bins in X

	int nBinsY = h->GetAxis(1)->GetNbins();
	int binCenterY1 = nBinsY/2; int binCenterY2 = binCenterY1 + 1;
	cout<<"binCenterY1 = "<<binCenterY1<<" - binCenterY2 = "<<binCenterY2<<endl;
	h->GetAxis(1)->SetRange(binCenterY1, binCenterY2); // apply Y-axis bin range cut
	TH2F* hXZ = (TH2F*) h->Projection(0,2); hXZ->SetName("hXZ");
	h->GetAxis(1)->SetRange(0, 0); // reset to full range
	hXZ->Scale(1/2.); // since we projected over two bins in Y

	TGraph* grDT = buildDeltaTGraph(getMaxTHnSparse(h),1000,collMat);
	TGraph* grDT1 = buildDeltaTGraph(getMaxTHnSparse(h1),1000,collMat);
	TGraph* grDT2 = buildDeltaTGraph(getMaxTHnSparse(h2),1000,collMat);

	TSpline3* splDT = new TSpline3("splDT", grDT);
	TSpline3* splDT1 = new TSpline3("splDT1", grDT1);
	TSpline3* splDT2 = new TSpline3("splDT2", grDT2);

	gStyle->SetOptStat(0); 
	gStyle->SetTitleFontSize(0.07);
	const Int_t NRGBs = 5;
	const Int_t NCont = 255;
	Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
	Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
	Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
	Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
	TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
	gStyle->SetNumberContours(NCont);

	TCanvas* cYZ = new TCanvas("cYZ","cYZ",1800,900);
	cYZ->cd();
	TPad* pYZ = new TPad("pYZ","",0,0,1,1);
   	pYZ->SetTopMargin(0.1);
	pYZ->SetBottomMargin(0.1);
	pYZ->SetLeftMargin(0.1);
	pYZ->SetRightMargin(0.2);
	pYZ->Draw();
	pYZ->cd();

	for(int ix = 1; ix <= hYZ->GetNbinsX(); ++ix) 
	{
		for(int iy = 1; iy <= hYZ->GetNbinsY(); ++iy) 
		{
			double Edep = hYZ->GetBinContent(ix, iy);   // J/cm^3
			double dT   = splDT->Eval(Edep);
			hYZ->SetBinContent(ix, iy, dT);
		}
	}

	hYZ->SetTitle(
		Form("10 GeV e^{-} on %.3f mm %s @ %.f mA / %.f kHz",
			collLength,collMat.Data(),Ibeam*1e3,revFreq*1e-3));
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
	drawCollEdge(collLength*1e-1,3.0,0.1);
	cYZ->SaveAs(Form("cYZ_%s.png",collMat.Data()));

	int iz = 2;  // Z axis index
	TAxis* zAxis = h->GetAxis(iz);
	int nBinsZ = zAxis->GetNbins();
	std::vector<double> maxEdep(nBinsZ, 0.0);
	std::vector<double> maxEdep1(nBinsZ, 0.0);
	std::vector<double> maxEdep2(nBinsZ, 0.0);
	std::vector<double> maxRelEdepError(nBinsZ, 0.0);
	std::vector<double> maxRelEdepError1(nBinsZ, 0.0);
	std::vector<double> maxRelEdepError2(nBinsZ, 0.0);

	Long64_t nbins = h->GetNbins();
	int coord[10];
	for (Long64_t i = 0; i < nbins; ++i)
	{
		double edep = h->GetBinContent(i,coord);
		double edep_err = h->GetBinError(i);
		int bz1 = coord[iz] - 1;

		double edep1 = h1->GetBinContent(i,coord);
		double edep2 = h2->GetBinContent(i,coord);
		double edep1_err = h1->GetBinError(i);
		double edep2_err = h2->GetBinError(i);
		if (edep > maxEdep[bz1]){maxEdep[bz1] = edep; maxRelEdepError[bz1] = edep_err/edep;}
		if (edep1 > maxEdep1[bz1]){maxEdep1[bz1] = edep1; maxRelEdepError1[bz1] = edep1_err/edep1;}
		if (edep2 > maxEdep2[bz1]){maxEdep2[bz1] = edep2; maxRelEdepError2[bz1] = edep2_err/edep2;}
	}

	TGraphErrors* grMaxTempZ = new TGraphErrors();
	TGraphErrors* grMaxTempZ1 = new TGraphErrors();
	TGraphErrors* grMaxTempZ2 = new TGraphErrors();
	double maxVal = 0, maxValErr = 0;
	double maxVal1 = 0, maxValErr1 = 0;
	double maxVal2 = 0, maxValErr2 = 0;

	for (int bz = 0; bz < nBinsZ; ++bz)
	{
		double z = zAxis->GetBinCenter(bz + 1);

		double val = splDT->Eval(maxEdep[bz]);
		double val1 = splDT1->Eval(maxEdep1[bz]);
		double val2 = splDT2->Eval(maxEdep2[bz]);

		double err = splDT->Eval(maxEdep[bz]) * maxRelEdepError[bz];
		double err1 = splDT1->Eval(maxEdep1[bz]) * maxRelEdepError1[bz];
		double err2 = splDT2->Eval(maxEdep2[bz]) * maxRelEdepError2[bz];

		grMaxTempZ->SetPoint(bz, z, val); grMaxTempZ->SetPointError(bz, 0, err);
		grMaxTempZ1->SetPoint(bz, z, val1); grMaxTempZ1->SetPointError(bz, 0, err1);
		grMaxTempZ2->SetPoint(bz, z, val2); grMaxTempZ2->SetPointError(bz, 0, err2);

		if(maxVal < val){maxVal = val; maxValErr = err;}
		if(maxVal1 < val1){maxVal1 = val1; maxValErr1 = err1;}
		if(maxVal2 < val2){maxVal2 = val2; maxValErr2 = err2;}
	}
	cout<<"dT_max = "<<maxVal<<" +/- "<<maxValErr<<endl;

	TCanvas* c3 = new TCanvas("c3","c3",900,900);
	c3->cd();
	TPad* p3 = new TPad("p3","",0,0,1,1);
   	p3->SetTopMargin(0.1);
	p3->SetBottomMargin(0.12);
	p3->SetLeftMargin(0.20);
	p3->SetRightMargin(0.05);
	p3->Draw();
	p3->cd();

	grMaxTempZ->SetMinimum(0);
	grMaxTempZ->SetMaximum(maxVal2*1.1);
	grMaxTempZ->GetXaxis()->SetRangeUser(minZ,maxZ);
	grMaxTempZ->SetTitle(
		Form("10 GeV e^{-} on %.3f mm %s",
			collLength,collMat.Data()));
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

	// dT_shock is the practical instantaneous temperature-rise threshold
	// adopted for the demonstration aluminum robustness criterion.
	double dT_shock = 0, Tm_T0 = 0;
	if(collMat == "Al")	{dT_shock = Al.dT_shock; Tm_T0 = Al.Tm-T0;}
	else{cout<<"ERROR:: Wrong material name: "<<collMat<<endl;return;}

	TLine* shock_line = new TLine(minZ,dT_shock,maxZ,dT_shock);
	shock_line->SetLineWidth(4);
	shock_line->SetLineStyle(7);
	shock_line->SetLineColor(kMagenta);
	shock_line->Draw("SAME L");

	TLine* Tm_T0_line = new TLine(minZ,Tm_T0,maxZ,Tm_T0);
	Tm_T0_line->SetLineWidth(3);
	Tm_T0_line->SetLineStyle(9);
	Tm_T0_line->SetLineColor(kOrange);
	Tm_T0_line->Draw("SAME L");

	gPad->SetGrid();

	TLegend* leg = new TLegend(0.2,0.45,0.5,0.85);
	leg->AddEntry(grMaxTempZ,Form("%.f mA",Ibeam*1e3),"l");
	leg->AddEntry(grMaxTempZ1,Form("%.f mA",Ibeam1*1e3),"l");
	leg->AddEntry(grMaxTempZ2,Form("%.f mA",Ibeam2*1e3),"l");
	leg->AddEntry(shock_line,Form("#DeltaT_{shock} = %.f ^{o}C",dT_shock),"l");
	leg->AddEntry(Tm_T0_line,Form("T_{m}-T_{0} = %.f ^{o}C",Tm_T0),"l");
	leg->Draw();

	c3->SaveAs(Form("c3_%s.png",collMat.Data()));

	TCanvas* c4 = new TCanvas("c4","c4",900,900);
	c4->cd();
	TPad* p4 = new TPad("p4","",0,0,1,1);
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

	c4->SaveAs(Form("c4_%s.png",collMat.Data()));

	//______________________________________________________________________________________
	// Scan beam current and calculate the volume exceeding the thermal-shock
	// and melting thresholds.  The current scan is logarithmic to show the
	// onset and rapid growth of the affected volume over several decades.
	const double IbeamMax = 2.5; // [A]
	TGraph* grDT_shock = buildDeltaTGraph(getMaxTHnSparse(h) * IbeamMax / Ibeam,1000,collMat);
	TSpline3* splDT_shock = new TSpline3("splDT_shock", grDT_shock);

	const Int_t nIbeam = 50;
    	double IbeamArray[nIbeam];

    	double logIbeamMin = TMath::Log10(1e-3);
    	double logIbeamMax = TMath::Log10(IbeamMax);

    	for (int i = 0; i < nIbeam; i++) 
	{
        	IbeamArray[i] = TMath::Power(10, logIbeamMin + i * (logIbeamMax - logIbeamMin) / (nIbeam-1));
    	}
	cout<<"IbeamArray[first] = "<<IbeamArray[0]<<" - IbeamArray[last] = "<<IbeamArray[nIbeam-1]<<endl;

	double shockVolume[nIbeam] = {};
	double meltVolume[nIbeam] = {};
	double beamCurrent[nIbeam] = {};
	double maxEdepDensityVal[nIbeam] = {};
	double totEdepVal[nIbeam] = {};

	for (Long64_t i = 0; i < nbins; ++i)
	{
		Int_t coord[10];
		double edep = h->GetBinContent(i,coord);
		double edep_err = h->GetBinError(i);
		if(edep <= 0){continue;}
		double dV =
			h->GetAxis(0)->GetBinWidth(coord[0]) *
			h->GetAxis(1)->GetBinWidth(coord[1]) *
			h->GetAxis(2)->GetBinWidth(coord[2]);

		for(Int_t j = 0; j < nIbeam; j++)
		{
			double Ibeam_j = IbeamArray[j]; // [A]
			double edep_j = edep * Ibeam_j / Ibeam; 	// scale to the given beam current
			double edep_err_j = edep_err * Ibeam_j / Ibeam; // scale to the given beam current
			if(maxEdepDensityVal[j] < edep_j){maxEdepDensityVal[j] = edep_j;}
			totEdepVal[j] += edep_j * dV;
			double dT_j = splDT_shock->Eval(edep_j);

			if(dT_j > dT_shock)	{shockVolume[j] += dV * 1e3;} // cm^3 -> mm^3
			if(dT_j > Tm_T0)	{meltVolume[j] += dV * 1e3;} // cm^3 -> mm^3

			beamCurrent[j] = Ibeam_j * 1e3; // A -> mA
		}
	}


	TCanvas* c5 = new TCanvas("c5","c5",900,900);
	c5->cd();
	TPad* p5 = new TPad("p5","",0,0,1,1);
   	p5->SetTopMargin(0.1);
	p5->SetBottomMargin(0.12);
	p5->SetLeftMargin(0.22);
	p5->SetRightMargin(0.01);
	p5->Draw();
	p5->cd();

	TGraph* grShockVolume = new TGraph(nIbeam,beamCurrent,shockVolume);
	TGraph* grMeltVolume = new TGraph(nIbeam,beamCurrent,meltVolume);

	grShockVolume->SetTitle("10 GeV e^{-}");
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

	TLegend* leg5 = new TLegend(0.25,0.55,0.9,0.8);
	leg5->AddEntry(grShockVolume,"#DeltaT > #DeltaT^{practical}_{shock}","pl");
	leg5->AddEntry(grMeltVolume,"#DeltaT > T_{m}-T_{0}","pl");
	leg5->Draw();

	gPad->SetGrid();
	gPad->SetLogx();

	c5->SaveAs(Form("c5_%s.png",collMat.Data()));

	grShockVolume->SetName("grShockVolume");
	grShockVolume->SaveAs(Form("grShockVolume_%s.root",collMat.Data()));
	grMeltVolume->SetName("grMeltVolume");
	grMeltVolume->SaveAs(Form("grMeltVolume_%s.root",collMat.Data()));

	TCanvas* c6 = new TCanvas("c6","c6",900,900);
	c6->cd();
	TPad* p6 = new TPad("p6","",0,0,1,1);
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

	c6->SaveAs(Form("c6_%s.png",collMat.Data()));


	TGraph* grMaxEdepDensity = new TGraph(nIbeam,beamCurrent,maxEdepDensityVal);

	cout<<"I_beam = 100 mA --> Max E_dep = "<<grMaxEdepDensity->Eval(100)<<" J/cm^3"<<endl;
	cout<<"I_beam = 200 mA --> Max E_dep = "<<grMaxEdepDensity->Eval(200)<<" J/cm^3"<<endl;

	TGraph* grTotEdep = new TGraph(nIbeam,beamCurrent,totEdepVal);

	cout<<"I_beam = 100 mA --> Total E_dep = "<<grTotEdep->Eval(100)<<" J"<<endl;
	cout<<"I_beam = 200 mA --> Total E_dep = "<<grTotEdep->Eval(200)<<" J"<<endl;

	return;
}


