#include "B4aKrattaParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
using namespace std;


B4aKrattaParameterisation::B4aKrattaParameterisation(  
        G4int    noDet, 
        G4double thMin,
	G4double thMax,
	G4double phiMin,
	G4double phiMax,
	G4double widthKr1, 
        G4double distX,
	G4double distY,
	G4double distZ,
	G4double spacing
	
)
 : G4VPVParameterisation()
{
  fNoDet  =  noDet; 
  fThMin =  thMin;
  fThMax =  thMax;
  fPhiMin = phiMin;
  fPhiMax = phiMax;
  fWidth   =  widthKr1;
  fDistX    = distX;
  fDistY    = distY;
  fDistZ    = distZ;
  fSpacing =  spacing;

   if( fNoDet > 0 ){
     /*     
for (G4int copyNo=0;copyNo<fNoDet;copyNo++)
       {
	 
	 fPhi[copyNo]=fPhiMin+copyNo*fWidth;
	 fTheta[copyNo]=fThMin + 0.5*fWidth;
	 
	 fPsi[copyNo]= copyNo*40*mm;
       }
     */


     //fRmaxIncr =  0.5 * (lengthFinal-lengthInitial)/(noChambers-1);
     //if (spacingZ < widthChamber) {
     //  G4Exception("B2bChamberParameterisation::B2bChamberParameterisation()",
     //              "InvalidSetup", FatalException,
     //              "Width>Spacing");
      }
   
}



B4aKrattaParameterisation::~B4aKrattaParameterisation()
{ }




void B4aKrattaParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{
  // Note: copyNo will start with zero!

  //G4double tg_th=fWidth/fDist;
  //G4double th=atan(tg_th);
  //G4double theta = fStartTh + copyNo*fWidth;
  //G4ThreeVector origin(fDistX,fDistY,fDistZ);
  //physVol->SetTranslation(origin);

  G4RotationMatrix *rm0 = new G4RotationMatrix();
  //rm0->rotateX(fPhi[copyNo]);
  rm0->rotateY(-90.*deg);
  //rm0->rotateZ(0.*deg);
  
  physVol->SetRotation(rm0);
  

  //G4ThreeVector mv1(fDistX,0.,0.);
  //physVol->SetTranslation(mv1);


  G4double fx,fy,fz;

  G4double phi= fPhiMin+copyNo*fWidth;  
  fx=(fDistX+91*mm)*cos(phi);
  fy=(fDistX+91*mm)*sin(phi);


 
 //fx=(fDistX+91*mm)*cos(fPhi[copyNo]);
  //fy=(fDistX+91*mm)*sin(fPhi[copyNo]);

  //fx=(fDistX)*cos((7./4.)*3.14);
  //fy=(fDistX)*sin((7./4.)*3.14);
 
  G4ThreeVector mv1(fx,fy,0.);
  physVol->SetTranslation(mv1);


  //fx=fPsi[copyNo]*cos(fPhi[copyNo])*sin(fTheta[copyNo]);
  //fy=fPsi[copyNo]*sin(fPhi[copyNo])*sin(fTheta[copyNo]);
  //fz=fDistX*cos(fTheta[copyNo]);//*fPsi[copyNo];

  //fx=fDistX*cos(fPhi[copyNo])*sin(fTheta[copyNo]);
  //fy=fDistX*sin(fPhi[copyNo])*sin(fTheta[copyNo]);
  //fz=fDistX*cos(fTheta[copyNo]);//*fPsi[copyNo];
  
  
  //cout<<fx<<" "<<fy<<endl;
  cout<<copyNo<<" "<<phi*180./3.14<<endl;
  //G4ThreeVector origin2(fx,fy,0.);

  //G4ThreeVector origin2(fx,fy,0.);
  //physVol->SetTranslation(origin2);

  

  
  G4RotationMatrix *rm=new G4RotationMatrix();
  //rm->set(fPhi[copyNo],fTheta[copyNo],fPsi[copyNo]) ;
  //rm->rotateX(fPhi[copyNo]);
  rm->rotateZ(phi);
  //rm->rotateZ(fPhi[copyNo]);
  //rm->rotateY(fTheta[copyNo]);
  //rm->rotateZ(fPhi[copyNo]);
 
  //physVol->SetRotation(rm);
  //  G4Transform3D *tr=new G4Transform3D(rm,mv1);
  //G4Transform3D m=
  //new G4RotateY3D(30.*deg);
//physVol->SetRotation(m);
}



//void B4aKrattaParameterisation::ComputeDimensions(){}


