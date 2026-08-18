/*
 * Calculate the practical thermal-shock temperature-rise threshold for the
 * aluminum demonstration.
 *
 * The criterion is based on the elastic thermal stress relation used in the
 * original study:
 *
 *   DeltaT_shock = sigma_ultimate * (1 - nu) / (E * alpha)
 *
 * with consistent SI units.
 */
void getDeltaTempShock()
{
	//-- Young’s modulus = Modulus of Elasticity 
	// Aluminum: https://www.matweb.com/search/DataSheet.aspx?MatGUID=0cd1edf33ac145ee93a0aa6fc666c0e0
	Double_t E = 68; // GPa
	//-- Coefficient of thermal expansion (CTE)
	Double_t a = 24.0; // um/m/K
	//-- Poisson’s ratio
	// Al: https://alloysintl.com/inventory/aluminum-alloys-supplier-old/aluminum-1100/
	Double_t nu = 0.330; 
	//-- Thermal stress = Tensile Strength, Ultimate  
	// Al: https://alloysintl.com/inventory/aluminum-alloys-supplier-old/aluminum-1100/
	Double_t sigma_ultimate = 90; // MPa

	cout<<"Al : E = "<<E<<" [GPa] | alpha = "<<a<<" [um/m/K] | nu = "<<nu<<" | sigma = "<<sigma_ultimate<<" [MPa]";
	cout<<" ==> dT_shock = "<<sigma_ultimate*1e-3*(1.0-nu)/(E*a*1e-6)<<" [K]"<<endl;

	return;
}
