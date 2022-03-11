#include "B4aDetectorConstruction.hh"
#include "B4aPrimaryGeneratorAction.hh"
#include "B4aDetectorMessenger.hh"

//#include "B4aKrattaParameterisation.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4PVParameterised.hh"
#include "G4ReflectedSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Transform3D.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

using namespace std;

G4ThreadLocal
G4GlobalMagFieldMessenger* B4aDetectorConstruction::fMagFieldMessenger = 0;



B4aDetectorConstruction::B4aDetectorConstruction()
 : G4VUserDetectorConstruction(),
  fStepLimit(NULL),
   //fGapPV(0),
   fCheckOverlaps(true),bt(0),fTargetRegion(0),fKrWinRegion(0),fPDRegion(0),fKrWrRegion(0)
{

detectorMessenger = new B4aDetectorMessenger(this);
}





B4aDetectorConstruction::~B4aDetectorConstruction()
{
delete fStepLimit ;
}



G4VPhysicalVolume* B4aDetectorConstruction::Construct()
{

  DefineMaterials();


  return DefineVolumes();
}



void B4aDetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  //G4NistManager* nistManager = G4NistManager::Instance();
  //nistManager->FindOrBuildMaterial("G4_Pb");
  G4NistManager* man = G4NistManager::Instance();

  G4String name, symbol;
  G4double z, a, density, ncomponents;
  G4int natoms;
  static const double     pi  = 3.14159265358979323846;
  static const double  twopi  = 2*pi;

  // Vacuum
  //G4Material* vaccum =new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
G4Material* vaccum = man->FindOrBuildMaterial("G4_Galactic");

  // Air
  //G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  //G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  //G4Material* air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  //air->AddElement(N, 70.*perCent);
  //air->AddElement(O, 30.*perCent);
  man->FindOrBuildMaterial("G4_AIR");

  //Aluminium
  man->FindOrBuildMaterial("G4_Al");

  //Ti
  man->FindOrBuildMaterial("G4_Ti");


  //CsI
  G4Material *CsI = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");


  //Tal
  G4Element *Tl = man->FindOrBuildElement("Tl");
    //new G4Element("Talium", "Tl", z=81, a= 11.72*g/mole);

  //CsI (Tl)
  G4Material *CsI_Tl = new G4Material("CsI_Tl",density= 4.51*g/cm3,ncomponents=2);
  CsI_Tl->AddMaterial(CsI,99.6*perCent);
  CsI_Tl->AddElement(Tl,0.4*perCent);
  CsI_Tl->GetIonisation()->SetMeanExcitationEnergy(553.1*eV);


  //Vikuiti ESR foil - polymer, non-metalic --  Polyethylene (??)
  man->FindOrBuildMaterial("G4_POLYETHYLENE");


  man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  //Si - photodiodes
  //G4Element* Si = new G4Element("Silicon", "Si", z=14., a=28.09*g/mole);
  man->FindOrBuildMaterial("G4_Si");

  //kratta window
  man->FindOrBuildMaterial("G4_Cu");

  // *********************

  //Target windows
  man->FindOrBuildMaterial("G4_MYLAR");

  //Target
  man->FindOrBuildMaterial("G4_Pyrex_Glass");

  //Target He3
  //G4Material* he3 = = new G4Material("Hel3", z=2., a=12.0078*g/mole,density=0.1650 *kg/m3);
  G4Material* he3 =new G4Material("Helium3", z=2., a=12.0078*g/mole,density=0.1650 *kg/m3, kStateGas, 293*kelvin, 1.*bar);

  //Mini Drift Chamber Gas
  G4Element* C = new G4Element("Carbon","C",z=6,a=12.01*g/mole);
  G4Element* F = new G4Element("Fluorine","F",z=9,a=18.9984032*g/mole);
  G4Material *CF4 = new G4Material("TetraFluoroMethane",density=3.72*mg/cm3,ncomponents=2,kStateGas,288.15*kelvin,1*atmosphere);
  CF4->AddElement(C,natoms=1);
  CF4->AddElement(F,natoms=2);

  //Graphite target
  //G4Element* C1 = new G4Element("Carbon","C",z=6,a=12.01*g/mole);
  //G4Material *CGr = new G4Material("Carbon",z=6.,a=12.01*g/mole,density=1.82*g/cm3);
  G4Material* CGr =  new G4Material("Graphite", density= 1.82*g/cm3, ncomponents=1);
  CGr->AddElement(C, natoms=1.);

  G4Element* D = new G4Element(name="Deuteron",symbol="D",z=1,a=2.0141018*g/mole);
  G4Material *CD2 = new G4Material(name="CD2",density=1.06*g/cm3,ncomponents=2);
  CD2->AddElement(C,natoms=2);
  CD2->AddElement(D,natoms=4);//lekko zmienione

  G4Material *CTarget = new G4Material(name="CTarget", density=2.1*g/cm3, ncomponents=1);
  CTarget->AddElement(C,natoms=1);
  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}



G4VPhysicalVolume* B4aDetectorConstruction::DefineVolumes()
{
  //scaling factor
  G4double aa=1.;
  //G4double aa=0.35;


  //Full Thicknesses
  G4double ESRfoilThick= 0.065*mm;
  G4double airLayer = 2.*mm;
  G4double pdThick = 0.5*mm;
  G4double housThick = 2.*mm;
  G4double krwinfoilThick = 0.01*mm;
  G4double frameWidth = 5.*mm;
  G4double frameThick = 2.*mm;
  G4double frameThickCoating = 3.*mm;

  G4double glassThick=2.*mm;
  G4double tarwinThick=0.01*mm;
  G4double graThick=1.*mm;
  G4double tiThick=0.013*mm;

  //C12target-detector full distance
  G4double ftar_det = 412.*mm;

  //Titarget-detector full distance
  G4double ftar_detTi = 600.*mm;

  //Half Sizes
   G4double fworldSizeXY = 800*mm;
   G4double fworldSizeZ  = 800*mm;

  //CsI - short
  G4double fCsiShort_xy1 =  14.*aa*mm;
  G4double fCsiShort_xy2 =  14.835*aa*mm;
  G4double fCsiShort_z   = 12.5 *mm;

  //CsI - long
  G4double fCsiLong_xy1 =  16.385*aa*mm;
  G4double fCsiLong_xy2 =  19.25*aa*mm;
  G4double fCsiLong_z   = 62.5 *mm;

  //PD
  G4double fPd_xy =  14.0*aa*mm;//28x28 mm^2
  G4double fPd_z    = 0.25*mm;

  //PD active layer
  G4double fPDac_z    = 0.23925*mm;

  //Front dead layer of PD
  G4double fPDdeadF_z    = 0.00075*mm;

  //Rear dead layer of PD
  G4double fPDdeadR_z    = 0.01*mm;

  //Kratta window
  G4double fKRwin_xy = 14.*aa*mm;


  //Kratta housing
  G4double fKRhouse_xy1 =  14.*aa*mm + frameWidth;
  G4double fKRhouse_xy2 =  23.*aa*mm + frameThick;
  G4double fKRhouse_z   =  88.*mm + frameThick;
  //G4double fKRhouse_z   =  118*mm+ frameThick;
  G4double dTh  =  1. *mm;


  //air gap 1
  G4double fKRairgap1_xy =  14.*aa*mm; //+ (frameThick-dTh);
  G4double fKRairgap1_z =  15*mm;
  //G4double fKRairgap1_z =  14.90*mm;

  //air gap 2
  G4double fKRairgap2_xy =  14.*aa*mm ;
  G4double fKRairgap2_z =  4.*mm;

  //air gap 3
  G4double fKRairgap3_xy =  14.*aa*mm;
  G4double fKRairgap3_z =  5.37*mm;

  //TARGET Ti-cell - flange
  G4double fRmax = 15.*mm;
  //G4double fRmin = fRmax-glassThick;
  G4double fh_z =  0.05*mm;

  //TARGET graphite-cell
  G4double fGra_xy = 15.*mm;
  G4double fGra_z =  200*mm;

  //TARGET Ti-cell - target -spallation
  G4double fRmax2 = 15.*mm;
  //G4double fRmin = fRmax-glassThick;
  G4double fh_z2 =  0.065*mm;

  //target side windows
  G4double ftarws_x = 8.*mm;
  //G4double ftarws_z = 2.*mm;
  G4double ftarws_y = 30.*mm;

  //target face windows
  G4double ftarwfR = 8.*mm;

  //Kratta Main
  G4double fKRRmin=ftar_det-100.*mm;
  G4double fKRRmax=600.*mm;

  // Get materials
  //******************************************
  G4Material* defaultMaterial = G4Material::GetMaterial("G4_AIR");
  G4Material* csi_mat = G4Material::GetMaterial("CsI_Tl");
  G4Material* pd_mat = G4Material::GetMaterial("G4_Si");
  G4Material* tgwin_mat = G4Material::GetMaterial("G4_MYLAR");
  G4Material* krwin_mat = G4Material::GetMaterial("G4_Cu");
  G4Material* cell_mat = G4Material::GetMaterial("G4_Pyrex_Glass");
  G4Material* target_mat = G4Material::GetMaterial("G4_Ti");
  G4Material* housing_mat = G4Material::GetMaterial("G4_Al");
  G4Material* csiwr_mat = G4Material::GetMaterial("G4_POLYETHYLENE");
  G4Material* gas_chamber = G4Material::GetMaterial("TetraFluoroMethane");
  G4Material* tarGra = G4Material::GetMaterial("Graphite");
  G4Material* plastic = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* vaccum = G4Material::GetMaterial("G4_Galactic");
  G4Material* CD2 = G4Material::GetMaterial("CD2");
  G4Material* CTarget = G4Material::GetMaterial("CTarget");
  G4Material* CGr = G4Material::GetMaterial("CGr");


  //if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
  //G4ExceptionDescription msg;
  //msg << "Cannot retrieve materials already defined.";
  //G4Exception("B4DetectorConstruction::DefineVolumes()", "MyCode0001", FatalException, msg);
  //}
  if (fTargetRegion){
    delete fTargetRegion;
    delete fGraRegion;
    delete fPDRegion;
    delete fKrWrRegion;
    delete fKrWinRegion;

}
  fTargetRegion= new G4Region("TgRegion");
  //fTargetRegion
  fPDRegion= new G4Region("PDRegion");
  fKrWrRegion= new G4Region("WrRegion");
  fKrWinRegion= new G4Region("WinRegion");
  fGraRegion= new G4Region("GraRegion");

  //
  // World
  //******

  G4VSolid* worldS = new G4Box("World",fworldSizeXY, fworldSizeXY, fworldSizeZ);

  G4LogicalVolume* worldLV = new G4LogicalVolume( worldS, vaccum, "World");

  G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps



  //
  // KRATTA housing
  //

 G4double housThick1=0.5*aa*housThick;

 G4VSolid* kr_housingS2a = new G4Trd("kr_housing2a",fKRhouse_xy1,fKRhouse_xy2,fKRhouse_xy1,fKRhouse_xy2,fKRhouse_z);
 G4VSolid* kr_housingS2b = new G4Trd("kr_housing2b",fKRhouse_xy1-housThick1,fKRhouse_xy2-housThick1,fKRhouse_xy1-housThick1,fKRhouse_xy2-housThick1,fKRhouse_z);
 G4VSolid* kr_housingS= new G4SubtractionSolid("kr_housingS", kr_housingS2a, kr_housingS2b,0, G4ThreeVector(0,0,0));

 G4VSolid* kr_housingS2a1 = new G4Trd("kr_housing2a",fKRhouse_xy1,fKRhouse_xy2,fKRhouse_xy1,fKRhouse_xy2,fKRhouse_z+housThick);
//G4VSolid* kr_housingS2b1 = new G4Trd("kr_housing2b1",fKRhouse_xy1-housThick1,fKRhouse_xy2-housThick1,fKRhouse_xy1-housThick1,fKRhouse_xy2-housThick1,fKRhouse_z);



 /*G4Box*  box =
  new G4Box("Box123",20*mm,30*mm,40*mm);

  G4LogicalVolume* box123LV =
    new G4LogicalVolume(box,            //its solid
                        defaultMaterial,             //its material
                        "box123");         //its name */

   G4VSolid* kr_obstacle = new G4Box("kr_obstacle",fKRwin_xy,fKRwin_xy,0.5*frameThick);
   G4VSolid* kr_coating = new G4Box("kr_coating",fKRwin_xy,fKRwin_xy,0.5*frameThickCoating);
   G4VSolid* kr_pipe = new G4Tubs("kr_pipe",0.*mm,3.*mm,500.*mm,0.*deg,360.*deg);

   G4LogicalVolume* kr_obstacleLV =
    new G4LogicalVolume(kr_obstacle,            //its solid
                        plastic,             //its material
                        "kr_obstacle");         //its name

    G4LogicalVolume* kr_coatingLV =
    new G4LogicalVolume(kr_coating,            //its solid
                        csiwr_mat,             //its material
                        "kr_coating");         //its name

    G4LogicalVolume* kr_pipeLV =
    new G4LogicalVolume(kr_pipe,            //its solid
                        csiwr_mat,             //its material
                        "kr_pipe");         //its name






  G4LogicalVolume* kr_housingLV
    = new G4LogicalVolume(
		 kr_housingS2a1,
		 //kr_housingS2,---->problem ??
                 housing_mat,
                 "kr_housing");

  G4LogicalVolume* kr_housingAirLV
    = new G4LogicalVolume(
		 kr_housingS2b,
                 defaultMaterial,
                 //housing_mat,
                 "kr_housingAir");


  //******************************************
  G4double pSPhi1=-45.*pi/180.;
  G4double pDPhi1=90*pi/180.;
  G4double rep_h=35*mm;
  //G4double pSTheta=45*pi/180.;
  //G4double pDTheta=49*pi/180.;
  G4double fOffset=0.;
  //


  //G4VSolid* krMainLS = new G4Sphere("KrMainL",300.*mm,600.*mm,(3./2.)*pi,pi,0,pi);
  G4VSolid* krMainS = new G4Sphere("KrMain",fKRRmin,fKRRmax,0,twopi,0,pi);
  //G4VSolid* krMainS = new G4Sphere("KrMain",0,fKRRmin,0,twopi,0,twopi);
  //Left
  // G4VSolid* krMainLS = new G4Tubs("KrMainL",fKRRmin,fKRRmax,rep_h,pSPhi1,pDPhi1);
  //G4LogicalVolume* krMainLLV = new G4LogicalVolume(krMainLS,defaultMaterial,"krMainLLV");
  G4LogicalVolume* krMainLV = new G4LogicalVolume(krMainS,vaccum,"krMainLV");


  G4RotationMatrix *rm0 = new G4RotationMatrix();
  //rm0->rotateX(90.*deg);
  rm0->rotateX(0.*deg);
  rm0->rotateY(0.*deg);
  rm0->rotateZ(0.*deg);
  //rm0->rotateZ(-90.*deg);

  G4PVPlacement *krMainPV = new G4PVPlacement(rm0,
					      G4ThreeVector(),
					      krMainLV,
					      "krMainPV",
					      worldLV,
					      false,
					      0,
					      fCheckOverlaps);


  G4ReflectX3D  Xreflection1;
  G4Transform3D transform1 =Xreflection1;
  G4PVPlacement *kr_housingPV ;
  G4double dist = fKRhouse_z+ftar_det;

  G4int noDet=32;

  G4double theta[noDet],phi[noDet];
  // G4double lat[noDet],lon[noDet];

  // theta[0]=15.*deg;   phi[0]=4.*deg;
  // theta[1]=38.*deg;   phi[1]=4.*deg;
  // theta[2]=43.*deg;   phi[2]=4.*deg;
  // theta[3]=48.*deg;   phi[3]=4.*deg;
  // theta[4]=53.*deg;   phi[4]=4.*deg;
  // theta[5]=58.*deg;   phi[5]=4.*deg;
  // theta[6]=63.*deg;   phi[6]=4.*deg;
  //
  // theta[7]=15.*deg;   phi[7]=-4.*deg;
  // theta[8]=38.*deg;   phi[8]=-4.*deg;
  // theta[9]=43.*deg;   phi[9]=-4.*deg;
  // theta[10]=48.*deg;  phi[10]=-4.*deg;
  // theta[11]=53.*deg;  phi[11]=-4.*deg;
  // theta[12]=58.*deg;  phi[12]=-4.*deg;
  // theta[13]=63.*deg;  phi[13]=-4.*deg;
  //
  // theta[14]=-15.*deg;   phi[14]=4.*deg;
  // theta[15]=-38.*deg;   phi[15]=4.*deg;
  // theta[16]=-43.*deg;   phi[16]=4.*deg;
  // theta[17]=-48.*deg;   phi[17]=4.*deg;
  // theta[18]=-53.*deg;   phi[18]=4.*deg;
  // theta[19]=-58.*deg;   phi[19]=4.*deg;
  // theta[20]=-63.*deg;   phi[20]=4.*deg;
  //
  // theta[21]=-15.*deg;   phi[21]=-4.*deg;
  // theta[22]=-38.*deg;   phi[22]=-4.*deg;
  // theta[23]=-43.*deg;   phi[23]=-4.*deg;
  // theta[24]=-48.*deg;  phi[24]=-4.*deg;
  // theta[25]=-53.*deg;  phi[25]=-4.*deg;
  // theta[26]=-58.*deg;  phi[26]=-4.*deg;
  // theta[27]=-63.*deg;  phi[27]=-4.*deg;


  theta[0]=15.*deg;   phi[0]=0.*deg;
  theta[1]=37.*deg;   phi[1]=0.*deg;
  theta[2]=43.*deg;   phi[2]=0.*deg;
  theta[3]=49.*deg;   phi[3]=0.*deg;
  theta[4]=55.*deg;   phi[4]=0.*deg;
  theta[5]=61.*deg;   phi[5]=0.*deg;

  theta[6]=37.*deg;   phi[6]=6.*deg;
  theta[7]=43.*deg;   phi[7]=6.*deg;
  theta[8]=49.*deg;   phi[8]=6.*deg;
  theta[9]=55.*deg;   phi[9]=6.*deg;
  theta[10]=61.*deg;   phi[10]=6.*deg;

  theta[11]=37.*deg;   phi[11]=-6.*deg;
  theta[12]=43.*deg;   phi[12]=-6.*deg;
  theta[13]=49.*deg;  phi[13]=-6.*deg;
  theta[14]=55.*deg;  phi[14]=-6.*deg;
  theta[15]=61.*deg;  phi[15]=-6.*deg;

  theta[16]=-15.*deg;   phi[16]=0.*deg;
  theta[17]=-37.*deg;   phi[17]=0.*deg;
  theta[18]=-43.*deg;   phi[18]=0.*deg;
  theta[19]=-49.*deg;   phi[19]=0.*deg;
  theta[20]=-55.*deg;   phi[20]=0.*deg;
  theta[21]=-61.*deg;   phi[21]=0.*deg;

  theta[22]=-37.*deg;   phi[22]=6.*deg;
  theta[23]=-43.*deg;   phi[23]=6.*deg;
  theta[24]=-49.*deg;   phi[24]=6.*deg;
  theta[25]=-55.*deg;   phi[25]=6.*deg;
  theta[26]=-61.*deg;   phi[26]=6.*deg;

  theta[27]=-37.*deg;   phi[27]=-6.*deg;
  theta[28]=-43.*deg;   phi[28]=-6.*deg;
  theta[29]=-49.*deg;  phi[29]=-6.*deg;
  theta[30]=-55.*deg;  phi[30]=-6.*deg;
  theta[31]=-61.*deg;  phi[31]=-6.*deg;




  // for(int i=0;i<28;i++){
  // 	theta[i]+=63.*deg;
  // }




  // lat[0]=2.43*deg;	lon[0]=-1.28*deg;
  // lat[1]=0.30*deg;	lon[1]=-1.28*deg;
  // lat[2]=0.98*deg;	lon[2]=-1.28*deg;
  // lat[3]=2.27*deg;	lon[3]=-1.28*deg;
  // lat[4]=2.73*deg;	lon[4]=-1.28*deg;
  // lat[5]=1.45*deg;	lon[5]=-1.28*deg;
  // lat[6]=0.17*deg;	lon[6]=-1.28*deg;

  // lat[7]=2.43*deg;	lon[7]=1.28*deg;
  // lat[8]=0.30*deg;	lon[8]=1.28*deg;
  // lat[9]=0.98*deg;	lon[9]=1.28*deg;
  // lat[10]=2.27*deg;	lon[10]=1.28*deg;
  // lat[11]=2.73*deg;	lon[11]=1.28*deg;
  // lat[12]=1.45*deg;	lon[12]=1.28*deg;
  // lat[13]=0.17*deg;	lon[13]=1.28*deg;

  // lat[14]=2.43*deg;	lon[14]=-1.28*deg;
  // lat[15]=0.30*deg;	lon[15]=-1.28*deg;
  // lat[16]=0.98*deg;	lon[16]=-1.28*deg;
  // lat[17]=2.27*deg;	lon[17]=-1.28*deg;
  // lat[18]=2.73*deg;	lon[18]=-1.28*deg;
  // lat[19]=1.45*deg;	lon[19]=-1.28*deg;
  // lat[20]=0.17*deg;	lon[20]=-1.28*deg;

  // lat[21]=2.43*deg;	lon[21]=1.28*deg;
  // lat[22]=0.30*deg;	lon[22]=1.28*deg;
  // lat[23]=0.98*deg;	lon[23]=1.28*deg;
  // lat[24]=2.27*deg;	lon[24]=1.28*deg;
  // lat[25]=2.73*deg;	lon[25]=1.28*deg;
  // lat[26]=1.45*deg;	lon[26]=1.28*deg;
  // lat[27]=0.17*deg;	lon[27]=1.28*deg;




  for (int i=0;i<noDet;i++){

    G4ThreeVector mv1, mv2;
    G4RotationMatrix rm1, rm2;
    G4Transform3D trv, trv2;

    G4double fx=dist*cos(phi[i])*sin(theta[i]);
    G4double fy;
    G4double fz=dist*cos(theta[i]);

    auto dist2 = dist - 100*mm;


   if(i<16){
   	rm1.rotateX(-phi[i]);
   	rm1.rotateY(theta[i]);
   	rm1.rotateZ(0.*deg);
   	rm2.rotateX(-phi[i]);
   	rm2.rotateY(theta[i]);
   	rm2.rotateZ(0.*deg);

   }
   else{
   	rm1.rotateX(-phi[i]);
   	rm1.rotateY(theta[i]);
   	rm1.rotateZ(0.*deg);
   	rm2.rotateX(-phi[i]);
   	rm2.rotateY(theta[i]);
   	rm2.rotateZ(0.*deg);
   }


    if(i<6){
      fy=0.*mm;
      // fy=29.*mm;
      mv2.setY(fy);

    }
    if(i>5&&i<11){
      fy=29.*mm*1.7;
      mv2.setY(fy-frameWidth-housThick);
    }
    if(i>10&&i<17){
      fy=-29.*mm*1.7;
      mv2.setY(fy+frameWidth+housThick);

    }
    if(i>15&&i<22)
    {
    fy=0.*mm;
    // fy=29.*mm;
    mv2.setY(fy);
    }
    if(i>21 && i<27){
    fy=29.*mm*1.7;
      mv2.setY(fy-frameWidth-housThick);
    }

    if(i>26){
      fy=-29.*mm*1.7;
      mv2.setY(fy+frameWidth+housThick);
    }


    mv1.setX(fx);
    mv1.setY(fy);
    mv1.setZ(fz);
    mv2.setX(dist2*cos(phi[i])*sin(theta[i]));
    mv2.setZ(dist2*cos(theta[i]));

    trv=G4Transform3D(rm1,mv1);
    trv2=G4Transform3D(rm2,mv2);

    kr_housingPV =  new G4PVPlacement(
				     trv,
				     kr_housingLV,
				     "kr_housingPV",
				     krMainLV,
				     //worldLV,
				     false,
				     i,
				     fCheckOverlaps);


    new G4PVPlacement(
                    trv2,
                    kr_obstacleLV,
                    "kr_obstaclePV",
                    krMainLV,
                    false,
                    i,
                    fCheckOverlaps);


    new G4PVPlacement(
                    trv2,
                    kr_coatingLV,
                    "kr_coatingPV",
                    krMainLV,
                    false,
                    i,
                    fCheckOverlaps); }
    /*new G4PVPlacement(
    				trv,

					kr_pipeLV,
					"kr_pipePV",
					//worldLV,
					krMainLV,
					false,
					i,
					fCheckOverlaps); */

  //box testowy
  /*new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    box123LV,                //its logical volume
                    "box123PV",              //its name
                    worldLV,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    fCheckOverlaps);          //overlaps checking*/






  // G4ThreeVector mv1a;
  // G4RotationMatrix rm1a;
  // G4Transform3D trva;

  //A3
  // mv1a.setX(132.*mm);
  // mv1a.setY(20.*mm);
  // mv1a.setZ(fKRhouse_z+ftar_det-20*mm);

  // rm1a.rotateX(-5.*deg);
  // rm1a.rotateY(16.*deg);
  // rm1a.rotateZ(0*deg);


  // trva=G4Transform3D(rm1a,mv1a);

  // kr_housingPV =  new G4PVPlacement(
		// 		     trva,
		// 		     kr_housingLV,
		// 		     "kr_housingPV",
		// 		     krMainLV,
		// 		     //worldLV,
		// 		     false,
		// 		     0,
		// 		     fCheckOverlaps);


    //****************************

    /*
   //Right-old

   G4double phi_r = phiMinR+ cN*widthKr1;

   G4double fx_r=dist*cos(phi_r)*sin(theta);
    G4double fy_r=dist*sin(phi_r)*sin(theta);
    G4double fz_r=dist*cos(theta);


    mv2.setX(fx_r);
    mv2.setY(fy_r);
    mv2.setZ(fz_r);

    rm2.rotateX(0.*deg);
    rm2.rotateY(90*deg+alpha+betha);
    rm2.rotateZ(phi_r);

    trv2=G4Transform3D(rm2,mv2);
    */


    /*
    kr_housingPV =  new G4PVPlacement(trv2,
				      kr_housingLV,
				      "kr_housingPV",
				      krMainLV,
				      //worldLV,
				      false,
				      cpNb_r1,
				      fCheckOverlaps);
    */

    //*******************************


 // G4double dist00=ftar_det+fKRhouse_z;
  G4double dist00=ftar_det;

    G4double fx_r0=dist00*cos(0.*pi/180.)*sin(89.*pi/180.);
    G4double fy_r0=dist00*sin(0.*pi/180.)*sin(89.*pi/180.);
    G4double fz_r0=dist00*cos(89.*pi/180.);


    G4ThreeVector mv00;
    G4RotationMatrix rm00;

    mv00.setX(fx_r0);
    mv00.setY(fy_r0);
    mv00.setZ(fz_r0);

    rm00.rotateX(0.*deg);
    rm00.rotateY(90.*deg);
    rm00.rotateZ(0*deg);



 G4Transform3D trv0=G4Transform3D(rm00,mv00);

  G4PVPlacement *kr_housingAirPV =  new G4PVPlacement(0,
						      G4ThreeVector(0, 0, 0),
						      kr_housingAirLV,
						      "kr_housingAirPV",
						      kr_housingLV,
						      false,
						      0,
						      fCheckOverlaps);


  /*
  kr_housingPV =  new G4PVPlacement(trv0,
				    kr_housingLV,
				    "kr_housingPV",
				    krMainLV,
				    false,
				    0,
				    fCheckOverlaps);
 */


  //*******************************************************
  //
  // kratta window
  //

  G4VSolid* kr_win = new G4Box("kr_win",fKRwin_xy,fKRwin_xy,0.5*frameThick);
  //volume which is a window in Kratta, filled with the Air
  G4LogicalVolume* kr_hwin1LV= new G4LogicalVolume(kr_win,defaultMaterial,"kr_winAir");
  //G4LogicalVolume* kr_hwin1LV= new G4LogicalVolume(kr_win,krwin_mat,"kr_win");
  //G4double box_z1=-90.*mm+0.5*(frameThick+krwinfoilThick) ;
  G4double box_z1=-90.*mm-0.5*frameThick ;

  G4PVPlacement *fkrwin1PV = new G4PVPlacement(
					       0,
					       G4ThreeVector(0., 0., box_z1),
					       kr_hwin1LV,
					       "kr_winAirPV",
					       kr_housingLV,
					       //worldLV,
					       false,
					       0,
					       fCheckOverlaps);


  G4VSolid* kr_hwinS = new G4Box("kr_hwinS",fKRwin_xy,fKRwin_xy,0.5*krwinfoilThick);
  G4LogicalVolume* kr_hwinLV= new G4LogicalVolume(kr_hwinS,tgwin_mat,"kr_win");

  G4double hentrWin1=0.5*(frameThick-krwinfoilThick);

G4RotationMatrix *rm02 = new G4RotationMatrix();
  rm02->rotateY(90.*deg);
  rm02->rotateX(0*deg);
  rm02->rotateZ(0.);



 G4PVPlacement *fkrwinPV = new G4PVPlacement(
                 0,
                 //rm02,
                 G4ThreeVector(0,0,0),
                 //G4ThreeVector(0., 0., hentrWin1),
                 kr_hwinLV,
                 "kr_winPV",
		 //worldLV,
		 kr_hwin1LV,
		 false,
                 0,
                 fCheckOverlaps);



  //*****************************************
  //
  // CsI-short-wrapping Volium
  //

  G4double wr=ESRfoilThick;

  G4VSolid* csi_wrS = new G4Trd("csi_wr",fCsiShort_xy2,fCsiLong_xy1,fCsiShort_xy2,fCsiLong_xy1,wr);
  G4LogicalVolume* csi_wrLV= new G4LogicalVolume(csi_wrS,csiwr_mat,"csiwr");

  G4double csi1wr_z=-44.935*mm;

  G4PVPlacement *fcsi_wrPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., csi1wr_z), // its position
                 csi_wrLV,       // its logical volume
                 "csi_wr",           // its name
                 kr_housingAirLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps





  //
  // CSI-short
  //

 G4double csi1_z=-57.5*mm;
 G4VSolid* csi_shortS = new G4Trd("csi_short",fCsiShort_xy1,fCsiShort_xy2,fCsiShort_xy1,fCsiShort_xy2,fCsiShort_z);

 //G4LogicalVolume* csi_shortLV= new G4LogicalVolume(csi_shortS,csi_mat,"csi_short");
 G4LogicalVolume* csi_shortLV= new G4LogicalVolume(csi_shortS,csi_mat,"csi_short");


 G4PVPlacement *fcsi_shortPV = new G4PVPlacement(
                 0,                // no rotation
                 //G4ThreeVector(0., 0., csi1_z-wr), // its position--old
                 G4ThreeVector(0., 0., csi1_z), // its position
                 csi_shortLV,       // its logical volume
                 "csi_short",           // its name
                 kr_housingAirLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // CSI-long
  //
  G4double csi2_z=17.63*mm;

  G4VSolid* csi_longS = new G4Trd("csi_long",fCsiLong_xy1,fCsiLong_xy2,fCsiLong_xy1,fCsiLong_xy2,fCsiLong_z);
  //G4LogicalVolume* csi_longLV= new G4LogicalVolume(csi_longS,csi_mat,"csi_long");
  G4LogicalVolume* csi_longLV= new G4LogicalVolume(csi_longS,csi_mat,"csi_long");

  G4PVPlacement *fcsi_longPV = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0., 0., csi2_z), // its position
                 csi_longLV,       // its logical volume
                 "csi_long",       // its name
                 //csi_longwrLV,   // its mother  volume---old
                 kr_housingAirLV,
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  //******************************************
  // PD
  //
  G4VSolid* PD_S = new G4Box("Pd",fPd_xy,fPd_xy,fPd_z);
  G4LogicalVolume* PD_LV = new G4LogicalVolume(PD_S, pd_mat,"Pd");
  //
  // PD-Active Layer
  //
  G4VSolid* pdacS = new G4Box("Pdac",fPd_xy,fPd_xy,fPDac_z);
  G4LogicalVolume* pdacLV = new G4LogicalVolume(pdacS, pd_mat,"Pdac");



  //
  // PD-Front Dead Layer
  //
  G4VSolid* pddfS = new G4Box("Pddf",fPd_xy,fPd_xy,fPDdeadF_z);
  G4LogicalVolume* pddfLV = new G4LogicalVolume(pddfS, pd_mat,"Pddf");


  //
  // PD-Rear Dead Layer
  //
  G4VSolid* pddrS = new G4Box("Pddr",fPd_xy,fPd_xy,fPDdeadR_z);
  G4LogicalVolume* pddrLV = new G4LogicalVolume(pddrS, pd_mat,"Pddr");




  G4double front_pd= (-1.)*(0.5*pdThick-fPDdeadF_z);
  G4double rear_pd= 0.5*pdThick-fPDdeadR_z;
  //G4double ac_pd= 0.5*pdThick-(2.*fPDdeadR_z+2.*fPDdeadF_z)/2.;
  G4double ac_pd= -1.*(0.5*pdThick-2.*fPDdeadF_z-fPDac_z);


  G4PVPlacement * fpddf_PV = new G4PVPlacement(0,
					       G4ThreeVector(0., 0.,front_pd),
					       pddfLV,
					       "Pdf",
					       PD_LV,
					       false,
					       0,
					       fCheckOverlaps);

  G4PVPlacement *fpddr_PV = new G4PVPlacement(0,
					      G4ThreeVector(0., 0.,rear_pd),
					      pddrLV,
					      "Pdr",
					      PD_LV,
					      false,
					      0,
					      fCheckOverlaps);



  G4PVPlacement *fpdac_PV = new G4PVPlacement(0,
					      G4ThreeVector(0., 0.,ac_pd),
					      pdacLV,
					      "PdAc",
					      PD_LV,
					      false,
					      0,
					      fCheckOverlaps);


  //********************************************

  G4PVPlacement *fPD_PV = new G4PVPlacement(0,
					    G4ThreeVector(0., 0., -74.75*mm),
					    PD_LV,
					    "PD",
					    kr_housingAirLV,
					    false,
					    0,
					    fCheckOverlaps);



  fPD_PV = new G4PVPlacement(0,
			     G4ThreeVector(0., 0., -70.25*mm),
			     PD_LV,
			     "PD",
			      kr_housingAirLV,
			     false,
			     1,
			     fCheckOverlaps);

  fPD_PV = new G4PVPlacement(0,
			     G4ThreeVector(0., 0., 80.38*mm),
			     PD_LV,
			     "PD",
			      kr_housingAirLV,
			     false,
			     2,
			     fCheckOverlaps);


//******************************************
// Air Layer1
//
  G4VSolid* airlayer1S = new G4Box("AirL1",fKRairgap1_xy,fKRairgap1_xy,0.5*fKRairgap1_z);
  G4LogicalVolume* airlayer1LV = new G4LogicalVolume(airlayer1S, defaultMaterial,"AirL1");

  G4PVPlacement *fairlayer1PV = new G4PVPlacement( 0,
						   //G4ThreeVector(0., 0., -82765.*mm),
						   G4ThreeVector(0., 0., -82.5*mm),
						    airlayer1LV,
						   "AirL1",
						   kr_housingAirLV,
						   false,
						   0,
						   fCheckOverlaps);

//******************************************
// Air Layer2
//
  G4VSolid* airlayer2S = new G4Box("AirL2",fKRairgap2_xy,fKRairgap2_xy,0.5*fKRairgap2_z);
  G4LogicalVolume* airlayer2LV = new G4LogicalVolume(airlayer2S, defaultMaterial,"AirL2");


  G4PVPlacement *fairlayer2PV = new G4PVPlacement( 0,
						   G4ThreeVector(0., 0.,-72.5*mm),
						   airlayer2LV,
						   "AirL2",
						   kr_housingAirLV,
						   false,
						   0,
						   fCheckOverlaps);

//******************************************
// Air Layer3
//
  G4VSolid* airlayer3S = new G4Box("AirL3",fKRairgap3_xy,fKRairgap3_xy,0.5*fKRairgap3_z);
  G4LogicalVolume* airlayer3LV = new G4LogicalVolume(airlayer3S, defaultMaterial,"AirL3");


  G4PVPlacement *fairlayer3PV = new G4PVPlacement( 0,
						   G4ThreeVector(0., 0.,83.315*mm),
						   airlayer3LV,
						   "AirL3",
						   kr_housingAirLV,
						   false,
						   0,
						   fCheckOverlaps);

  //*******************************************************
  /*

//for testing

  G4ReflectX3D  Xreflection;
  G4Transform3D transform =Xreflection;

  G4ReflectedSolid * ReflKrR = new G4ReflectedSolid("kr_housingR", kr_housingS, transform);
G4LogicalVolume* kr_housingRLV = new G4LogicalVolume(ReflKrR,defaultMaterial,"kr_housingRLV");

  G4PVPlacement *krMainRPV = new G4PVPlacement(0,
					      G4ThreeVector(),
					      kr_housingLV,
					      "krMainRPV",
					      krMainRLV,
					      false,
					      0,
					      fCheckOverlaps);
  */


  //******************************************
  // Ti-Target - flange
  //******************************************
  G4RotationMatrix *rm01 = new G4RotationMatrix();
  rm01->rotateX(0*deg);
  rm01->rotateY(0.*deg);
  rm01->rotateZ(0);


  //G4VSolid* targ_volS0 = new G4Box("graphiteS",fGra_xy,fGra_xy,0.5*graThick);
  G4VSolid* targ_volS0 = new G4Tubs("targ_vol0",0.,fRmax,0.5*fh_z,0.,twopi);
  G4LogicalVolume* targ_volLV0 = new G4LogicalVolume(targ_volS0,target_mat,"targ_volLV0");
  /*
  G4PVPlacement *ftarg_volPV0 = new G4PVPlacement(rm01,
						  G4ThreeVector((-1)*fGra_z,0.,0.),
						  //G4ThreeVector(0.,0.,0.),
						  targ_volLV0,
						  "targ_volPV0",
						  worldLV,
						  false,
						  0,
						  fCheckOverlaps);
  */

  //******************************************
  // Graphite-Target
  //******************************************
  /*
  G4VSolid* targetGraS = new G4Box("graphiteS",fGra_xy,fGra_xy,0.5*graThick);
  G4LogicalVolume* targetGraLV = new G4LogicalVolume(targetGraS,tarGra,"targetGraLV");

  G4PVPlacement *targetGraPV = new G4PVPlacement(rm01,
						 G4ThreeVector(0.,0.,0.),
						 targetGraLV,
						 "targetGraPV",
						 worldLV,
						 false,
						 0,
						 fCheckOverlaps);
  */


//******************************************
  // Ti-Target- 0.13 um
  //******************************************

  /*G4VSolid* tarTi2S = new G4Box("Ti2S",fGra_xy,fGra_xy,0.5*tiThick);
  G4LogicalVolume* tarTi2LV = new G4LogicalVolume(tarTi2S,target_mat,"tarTi2LV");

  G4PVPlacement *tarTi2PV = new G4PVPlacement(rm01,
						 G4ThreeVector(0.,0.,0.),
						 tarTi2LV,
						 "tarTi2PV",
						 worldLV,
						 false,
						 0,
						 fCheckOverlaps);*/

  //PRZEROBIĆ NA CYLINDER
  G4VSolid* tarTi2S = new G4Tubs("Ti2S",0.*mm,30.*mm,0.28*mm,0.*deg,360.*deg);    //CD2
  G4LogicalVolume* tarTi2LV = new G4LogicalVolume(tarTi2S,CD2,"tarTi2LV");

  // G4VSolid* tarTi2S = new G4Tubs("Ti2S",0.*mm,30.*mm,0.1*mm,0.*deg,360.*deg);    //C
  // G4LogicalVolume* tarTi2LV = new G4LogicalVolume(tarTi2S,CD2,"tarTi2LV");

  G4ThreeVector mvTarTi;
  G4RotationMatrix rmTarTi;
  G4Transform3D trvTarTi;

  mvTarTi.setX(0);
  mvTarTi.setY(0);
  mvTarTi.setZ(0);

  rmTarTi.rotateX(0.*degree);
  // rmTarTi.rotateY(63.*degree);
  rmTarTi.rotateY(0.*degree);
  rmTarTi.rotateZ(0.*degree);

  trvTarTi=G4Transform3D(rmTarTi,mvTarTi);


  G4PVPlacement *tarTi2PV = new G4PVPlacement(
  						 trvTarTi,
						 tarTi2LV,
						 "tarTi2PV",
						 worldLV,
						 false,
						 0,
						 fCheckOverlaps);

//  ***********************************************

    fTargetRegion->AddRootLogicalVolume(targ_volLV0);
    //fGraRegion->AddRootLogicalVolume(targetGraLV);
    //fTargetRegion->AddRootLogicalVolume(tarTi2LV);

    fKrWrRegion->AddRootLogicalVolume(csi_wrLV);
    fKrWinRegion->AddRootLogicalVolume(kr_hwin1LV);
    fKrWinRegion->AddRootLogicalVolume(kr_hwinLV);

    fPDRegion->AddRootLogicalVolume(pddfLV);
    fPDRegion->AddRootLogicalVolume(pdacLV);
    fPDRegion->AddRootLogicalVolume(pddrLV);


//  ***********************************************
  //
  // print parameters
  //

  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    //  << "---> The calorimeter is " << nofLayers << " layers of: [ "
    //<< absoThickness/mm << "mm of " << absorberMaterial->GetName()
    //<< " + "
    //<< gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;


  //
  // Visualization attributes
  //
  // worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  //G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //calorLV->SetVisAttributes(simpleBoxVisAtt);
  /*
  G4Colour  white   ()              ;  // white
     G4Colour  white   (1.0, 1.0, 1.0) ;  // white
     G4Colour  gray    (0.5, 0.5, 0.5) ;  // gray
     G4Colour  black   (0.0, 0.0, 0.0) ;  // black
     G4Colour  red     (1.0, 0.0, 0.0) ;  // red
     G4Colour  green   (0.0, 1.0, 0.0) ;  // green
     G4Colour  blue    (0.0, 0.0, 1.0) ;  // blue
     G4Colour  cyan    (0.0, 1.0, 1.0) ;  // cyan
     G4Colour  magenta (1.0, 0.0, 1.0) ;  // magenta
     G4Colour  yellow  (1.0, 1.0, 0.0) ;  // yellow

  */



  G4VisAttributes* wo= new G4VisAttributes(G4Colour(0,0,0));
  wo->SetVisibility(false);
  //wo->SetForceAuxEdgeVisible(true);
  //wo->SetForceSolid(true);
  //tg->SetForceWireframe(true);
  worldLV->SetVisAttributes(wo);

  G4VisAttributes* tg= new G4VisAttributes(G4Colour(0,1,0));
  tg->SetVisibility(true);
  tg->SetForceAuxEdgeVisible(true);
  tg->SetForceSolid(true);
  //tg->SetForceWireframe(true);
  targ_volLV0->SetVisAttributes(tg);

   G4VisAttributes* tgM= new G4VisAttributes(G4Colour(0.2,0.2,0.2));
    tgM->SetVisibility(true);
    tgM->SetForceAuxEdgeVisible(true);
  //tg->SetForceSolid(true);
  //tgM->SetForceWireframe(true);
  krMainLV->SetVisAttributes(tgM);


  //blue
  G4VisAttributes* win= new G4VisAttributes(G4Colour(0,0,1.0));
  win->SetVisibility(true);
  //win->SetForceAuxEdgeVisible (true);
  win->SetForceSolid (true);
  //kr_hwin1LV->SetVisAttributes(win);


  //red
  G4VisAttributes* win1= new G4VisAttributes(G4Colour(1.,0.,0));
  win1->SetVisibility(true);
  //win->SetForceAuxEdgeVisible (true);
  win1->SetForceSolid (true);
  kr_hwin1LV->SetVisAttributes(win1);
  //targetGraLV->SetVisAttributes(win1);
  tarTi2LV->SetVisAttributes(win1);

  //green
  G4VisAttributes* cryst1= new G4VisAttributes(G4Colour(0,1.0,0));
  cryst1->SetVisibility(true);
  cryst1->SetForceAuxEdgeVisible (true);
  // cryst1->SetForceSolid (true);
  //csi_shortLV->SetVisAttributes(cryst);
  csi_longLV->SetVisAttributes(cryst1);
  //kr_housingLV->SetVisAttributes(cryst1);

 //magenta
  G4VisAttributes* cryst2= new G4VisAttributes(G4Colour(1.,0,1.));
  cryst2->SetVisibility(true);
  cryst2->SetForceAuxEdgeVisible (true);
  //cryst2->SetForceSolid (true);
  csi_shortLV->SetVisAttributes(cryst2);
  ////csi_longLV->SetVisAttributes(cryst);

 //magenta
  G4VisAttributes* krh= new G4VisAttributes(G4Colour(1.,0,1.));
  //krh->SetVisibility(true);
  krh->SetForceAuxEdgeVisible (true);
  //siatka
  krh->SetForceSolid (true);
  kr_housingLV->SetVisAttributes(krh);
  //csi_longLV->SetVisAttributes(cryst);



  //yellow
   G4VisAttributes* wrap= new G4VisAttributes(G4Colour(0.,0.,1));
  wrap->SetVisibility(true);
  wrap->SetForceAuxEdgeVisible (true);
  //wrap->SetForceSolid(true);
  csi_wrLV->SetVisAttributes(wrap);
  //csi_shortwrLV->SetVisAttributes(wrap);



  G4VisAttributes* entr= new G4VisAttributes(G4Colour(1,1,0));
  entr->SetVisibility(true);
  //entr->SetForceAuxEdgeVisible (true);
  entr->SetForceSolid(true);
 // kr_entrWinLV->SetVisAttributes(entr);


  G4VisAttributes* entr1= new G4VisAttributes(G4Colour(1,1,0));
  entr1->SetVisibility(true);
  //entr1->SetForceAuxEdgeVisible (true);
  entr1->SetForceSolid(true);
  kr_hwin1LV->SetVisAttributes(entr1);




  G4VisAttributes* krair= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  krair->SetVisibility(true);
  krair->SetForceAuxEdgeVisible(false);
  //krair->SetForceSolid(true);
  //tg->SetForceWireframe(true);
  airlayer1LV->SetVisAttributes(krair);
  airlayer2LV->SetVisAttributes(krair);
  airlayer3LV->SetVisAttributes(krair);

  G4VisAttributes* pdd= new G4VisAttributes(G4Colour(1.0, 1.0, 0));
  pdd->SetVisibility(true);

  pdd->SetForceSolid(true);
  pddfLV->SetVisAttributes(pdd);
  pddrLV->SetVisAttributes(pdd);

  G4VisAttributes* pda= new G4VisAttributes(G4Colour(1.0, 0, 1.0));
  pda->SetVisibility(true);
  pda->SetForceSolid(true);
  pdacLV->SetVisAttributes(pda);

  // User Limits
  //
  // set tracking constraints in a given logical volume
  //
  // Sets a max step length in the tracker region, with G4StepLimiter

  //G4double maxStep = 0.01*mm;
  //fStepLimit = new G4UserLimits(maxStep);
  //csi_shortLV->SetUserLimits(fStepLimit);
  //kr_hwin1LV->SetUserLimits(fStepLimit);

 //targ_volLV0->SetUserLimits(fStepLimit);
  /// Set additional contraints on the track, with G4UserSpecialCuts
  ///
  /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
  /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
  ///                                           maxLength,
  ///                                           maxTime,
  ///                                           minEkin));



  //
  // Always return the physical World
  //
  return worldPV;
}
void B4aDetectorConstruction::SetBtEnergy(G4double BtEnergy)  {

bt=BtEnergy;

 G4RunManager::GetRunManager()->PhysicsHasBeenModified();


}
void B4aDetectorConstruction::SetParamUpdate(){

// G4RunManager::GetRunManager()->SetUserAction(new B4aPrimaryGeneratorAction());

}


void B4aDetectorConstruction::ConstructSDandField()
{
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.

  G4ThreeVector fieldValue = G4ThreeVector();
  // fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  //fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  //G4AutoDelete::Register(fMagFieldMessenger);
}
