#include "B4aPrimaryGeneratorAction.hh"
#include "B4aDetectorConstruction.hh"
#include "B4aActionInitialization.hh"

#include "B4aNucleonElasticXS.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "B4aHistoManager.hh"
#include "TF1.h"
#include "TGraph.h"
// #include "G4AutoLock.hh"
// namespace { G4Mutex PrimMutex = G4MUTEX_INITIALIZER;}

#ifndef min
#define min(x,y) ((x)<(y) ? (x):(y))
#endif
#ifndef max
#define max(x,y) ((x)>(y) ? (x):(y))
#endif


#include <iostream>
#include <fstream>
using namespace std;



G4float csmax;

B4aPrimaryGeneratorAction::B4aPrimaryGeneratorAction(B4aDetectorConstruction* myDC, HistoManager *histo):myDetector(myDC),bt(0),fHistoManager(histo)
{
  // rnd= new TRandom3();
	// rnd = new G4RandGauss();
  //#include "B4aMaterial.cfg"
  //generator_min=myDC->GetKinematicsMin();
  //generator_max=myDC->GetKinematicsMax();

  icros=myDC->GetNeumann();
  bfwhmx = myDC ->GetBfwhmX();
  bfwhmy = myDC->GetBfwhmY();
  bt = myDC->GetBtEnergy();
  pz = myDC->GetPz() ;
  pzz = myDC->GetPzz() ;
  themin = myDC->GetThetaMin()*180/M_PI;
  themax = myDC->GetThetaMax()*180/M_PI;
  themin2 = myDC->GetTheta2Min()*180/M_PI ;
  themax2 = myDC->GetTheta2Max()*180/M_PI ;
  fimin = myDC->GetPhiMin()*180/M_PI ;
  fimax = myDC->GetPhiMax()*180/M_PI ;
  //fimax = myDC->GetPhiMax()*180/M_PI ;

  tXplace = myDC->GetTargetXplace();
  tYplace = myDC->GetTargetYplace();
  tZplace = myDC->GetTargetZplace();
  thigh = myDC->GetTargetHigh();
  G4cout<<"\n";

  G4cout<<"--------------->PARAMS: "<<"thmin: "<<themin<<" thmax: "<<themax<<endl;
  G4cout<<"--------------->PARAMS: "<<"fimin: "<<fimin<<" fimax: "<<fimax<<endl;
  G4cout<<"--------------->PARAMS: "<<bfwhmx<<" "<<thigh<<" "<<bt<<endl;

  int n_particle = 1;
  bt1 = 0.;		bt2 = 0.;		bt3 = 0.;
  t1  = 0.;		t2  = 0.;		t3  = 0.;
  fi1 = 0.;             fi2 = 0.;		fi3 = 0.;
  bt /=1000.;



  G4int nofParticles = 1;
  particleGun1 = new G4ParticleGun();
  particleGun2 = new G4ParticleGun();

  // default particle kinematic
  //
  //G4ParticleDefinition* particleDefinition;
   particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  //   particleDefinition  = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
     particleGun1->SetParticleDefinition(particleDefinition);
     p_mass = particleDefinition->GetPDGMass()/1000.; //GeV

     particleGun2->SetParticleDefinition(particleDefinition);



   //G4GenericIon::GenericIon() ;
     Ti48_mass = 44.6523;//GeV


gr = new TGraph("../data/pd_200.dat","%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");

  G4cout<<"--------------->NUMBER OF POINTS: "<< gr->GetN()<<endl;
}



B4aPrimaryGeneratorAction::~B4aPrimaryGeneratorAction()
{
  delete particleGun1;
  delete particleGun2;

}

void B4aPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*
 G4double worldZHalfLength = 0;
  G4LogicalVolume* worlLV
    = G4LogicalVolumeStore::GetInstance()->GetVolume("World");
  G4Box* worldBox = 0;
  if ( worlLV) worldBox = dynamic_cast< G4Box*>(worlLV->GetSolid());
  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  }
  */
  // Set gun position
  //fParticleGun ->SetParticlePosition(G4ThreeVector(0., 0., -worldZHalfLength));

  //fParticleGun->GeneratePrimaryVertex(anEvent);
  //**************************************************************


// G4AutoLock lock(&PrimMutex);
  //ProcNb();
  // lock.lock();
  Pos();
G4ThreeVector dir = flatGen();

  // double theta, phi;
  // theta = 58*G4UniformRand()+11;
  // phi = 25*G4UniformRand()-12.5;
  //
  //
  //
  // G4ThreeVector dir(1.0, 1.0, 1.0);
  // dir.setPhi(phi*deg);
  // dir.setTheta(theta*deg);
  // dir.setMag(1.0);

  bt=gr->Eval(dir.theta())/1000;
  G4cout<<" theta="<<dir.theta()*180/M_PI<<" bt="<<bt<<G4endl;
//bt=0.1;

  //G4ThreeVector dir = flatGen();

// lock.unlock();
  // GetStartEnergy(bt1*1000.,bt2*1000.);
  // GetStartAngleTheta(t1,t2);
  // GetStartAnglePhi(fi1,fi2);
  // GetStartPosition(vertex[0],vertex[1],vertex[2]);


//G4cout<<"****************************************&^%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$theta="<<dir.getTheta()<<" phi="<<dir.getPhi()<<"\n";
  particleGun1->SetParticleMomentumDirection(dir);



    particleGun1->SetParticleEnergy(bt*GeV);
    //particleGun2->SetParticleEnergy(bt2*GeV);

    G4ThreeVector PPTV(vertex[0],vertex[1],vertex[2]);
    // PPTV.rotateY(63*deg);
    particleGun1->SetParticlePosition(PPTV);

  //  particleGun2->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));

  // G4ThreeVector vec(1.0, 1.0, 1.0);
  // vec.setPhi(phi*deg);
  // vec.setTheta(theta*deg);
  // vec.setMag(1.0);

    particleGun1->GeneratePrimaryVertex(anEvent);
  //  particleGun2->GeneratePrimaryVertex(anEvent);
    fHistoManager->FillNtuple(41, vertex[0]);
    fHistoManager->FillNtuple(42, vertex[1]);
    fHistoManager->FillNtuple(43, vertex[2]);
    fHistoManager->FillNtuple(44, dir.getTheta()*180./CLHEP::pi);
    fHistoManager->FillNtuple(45, dir.getPhi()*180./CLHEP::pi);
    fHistoManager->FillNtuple(46, bt*GeV);



}










void B4aPrimaryGeneratorAction::Pos(void)
{
  G4double bsgx,bsgy;
  thigh = 0.5*mm;


  bsgx = bfwhmx/(2.*sqrt(2.*log(2.)));
  bsgy = bfwhmy/(2.*sqrt(2.*log(2.)));


  vertex[0] = tXplace + G4RandGauss::shoot(0,1);
  vertex[1] = tYplace + G4RandGauss::shoot(0,1);
  vertex[2] = tZplace -thigh + G4UniformRand()*2*thigh;


 }




G4ThreeVector B4aPrimaryGeneratorAction::flatGen(){
  double theta, phi;
  // if(G4UniformRand()-0.5>0){
  //   theta = 52*G4UniformRand()+76;
  // }
  // else{
  //   theta = 52*G4UniformRand()-2;
  // }
  // phi = 14*G4UniformRand()-7;


  // if(G4UniformRand()-0.5>0){
  //   theta = 33*G4UniformRand()+36;
  // }
  // else{
  //   theta = -33*G4UniformRand()-36;
  // }

  theta = 58*G4UniformRand()+11;
  phi = 25*G4UniformRand()-12.5;


	//
	// if(G4UniformRand()-0.5>0){
	// 	theta = 56*G4UniformRand()+12;
	// }
	// else{
	// 	theta = -56*G4UniformRand()-12;
	// }
	// phi = 14*G4UniformRand()-7;

//  G4cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! thcos="<<thcos<<" phi="<<phi<<"\n";
  G4ThreeVector vec(1.0, 1.0, 1.0);
  vec.setPhi(phi*deg);
  vec.setTheta(theta*deg);
  vec.setMag(1.0);
 // G4cout<<"****************************************&^%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$theta="<<vec.getTheta()<<" phi="<<vec.getPhi()<<"\n";
  return vec;

}
