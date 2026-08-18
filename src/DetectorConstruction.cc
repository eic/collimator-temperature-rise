#include "DetectorConstruction.hh"

// This example intentionally keeps the geometry simple.  The production 
// collimator model can be substituted for the G4Trd below without changing
// the downstream scoring/thermal-analysis concept.

DetectorConstruction::DetectorConstruction(RunAction* runAct) : G4VUserDetectorConstruction(),
    _runAction(runAct)
{}

DetectorConstruction::~DetectorConstruction()
{}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    // The world is vacuum and is deliberately much larger than the simple
    // collimator.  The beam starts 1 m upstream of the collimator.
	// World
	//
	G4Box* world_solid = new G4Box("world",2*m,2*m,2*m);
	G4LogicalVolume* world_log = new G4LogicalVolume(world_solid, GetMaterial("Vacuum"), "world");
	G4VPhysicalVolume* world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),world_log,"world",0,false,0);
	world_log->SetVisAttributes(G4VisAttributes::GetInvisible());

	// Collimator
	//
	if( _runAction->GetCollLength() > 0)
	{
		G4cout<<"[INFO] DetectorConstruction::Construct ==> Building a simple collimator"<<G4endl;
		// *********************************
		// tanTh = 2 * H / (L2 - L1)
		// L2 - L1 = 2 * H / tanTh
		// L2 = L1 + 2 * H / tanTh
		// .................................
		// ....................|<-L1->|.....
		// .....................------......
		// ..Beam..=====>....../   ^  \.....
		// .................../    |H  \....
		// ................../)Th  v    \...
		// .................--------------..
		// ................|<------------>|.
		// .......................L2........
		// .................................
		// ********************************
		G4double L1 = _runAction->GetCollLength();
		// H is the transverse half-height parameter of the trapezoidal jaw.
		// W is the jaw width in the orthogonal transverse direction.
		G4double H = 30 * mm; // head height
		G4double L2 = L1 + 2. * H / _runAction->GetCollTanTh();
		G4double W = 30 * mm; // jaw width
		G4cout<<"L1 = "<<L1/mm<<" [mm]; L2 = "<<L2/mm<<" [mm]; W = "<<W/mm<<" [mm]; H = "<<H/mm<<" [mm]"<<G4endl;

		// G4Trd has different y half-lengths at the two z faces.  After the
		// 90-degree rotation about X, the trapezoid taper is along the beam (Z).
		G4Trd* coll_trd = new G4Trd("coll_trd", W/2., W/2., L2/2., L1/2., H/2.);
		G4LogicalVolume* coll_log = new G4LogicalVolume(
			coll_trd,GetMaterial(_runAction->GetCollMaterial()),"coll_trd");
		G4RotationMatrix* coll_rot = new G4RotationMatrix();
		coll_rot->rotateX(90 * deg);

		G4double posX = 0 * mm, posY = 0 * mm, posZ = 0 * mm;

		new G4PVPlacement(
			coll_rot,
			G4ThreeVector(posX,posY,posZ),
			coll_log,"coll_box",world_log,false,0);
		coll_log->SetVisAttributes(G4VisAttributes(G4Color::Magenta()));
	}

	return world_phys;
}

G4Material* DetectorConstruction::GetMaterial(G4String matName)
{
    // The public demonstration uses only "Al", mapped to Geant4's G4_Al.

	G4NistManager* nistManager = G4NistManager::Instance();
	G4Material* material = nistManager->FindOrBuildMaterial("G4_Galactic");
	
	if 	(matName == "Vacuum")	{material = new G4Material("Vacuum", 1., 1.01*g/mole, 1e-30*g/cm3);}
	else if (matName == "Al" || matName == "AL")	{material = nistManager->FindOrBuildMaterial("G4_Al");}
	else
	{
		G4cout<<"[WARNING] Wrong material name: "<<matName<<" --> return default: Galactic"<<G4endl;
	}
	G4cout<<"[MATERIAL] "<<*material<<G4endl;

	return material;
}

