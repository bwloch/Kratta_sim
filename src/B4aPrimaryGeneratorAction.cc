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

  RandomInit();

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


   //cout<<"xxxxxxxx"<<p_mass<<endl;
   //myXS_test = new B4aNucleonElasticXS(particleDefinition);

   /*
  particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("Ti48");
     particleGun2->SetParticleDefinition(particleDefinition);
   Ti48_mass = particleDefinition->GetPDGMass()/1000.;
   */
   //G4GenericIon::GenericIon() ;
   //G4IonTable *ionTable=G4IonTable::GetIonTableda();
   //G4ParticleDefinition *Ti48=ionTable->GetIon(22,48,0.);
   //G4GenericIon::GenericIon() ;
     Ti48_mass = 44.6523;//GeV




     //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
     //fParticleGun->SetParticleEnergy(850.*MeV);

gr = new TGraph("../data/pp_200.dat","%lg %*lg %lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg");

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

  upunif(momentum);
// G4AutoLock lock(&PrimMutex);
  //ProcNb();
  // lock.lock();
  Pos();

  double theta, phi;
  theta = 58*G4UniformRand()+11;
  phi = 25*G4UniformRand()-12.5;



  G4ThreeVector dir(1.0, 1.0, 1.0);
  dir.setPhi(phi*deg);
  dir.setTheta(theta*deg);
  dir.setMag(1.0);

  bt=gr->Eval(theta)/1000;


  //G4ThreeVector dir = flatGen();

// lock.unlock();
  // GetStartEnergy(bt1*1000.,bt2*1000.);
  // GetStartAngleTheta(t1,t2);
  // GetStartAnglePhi(fi1,fi2);
  // GetStartPosition(vertex[0],vertex[1],vertex[2]);


//G4cout<<"****************************************&^%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$theta="<<dir.getTheta()<<" phi="<<dir.getPhi()<<"\n";
  particleGun1->SetParticleMomentumDirection(dir);
//  particleGun1->SetParticleMomentumDirection(G4ThreeVector(0.01,0.01,0.01));
  //particleGun2->SetParticleMomentumDirection(G4ThreeVector(momentum[3],momentum[4],momentum[5]).unit());
  //particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
  //particleGun2->SetParticleMomentumDirection(G4ThreeVector(momentum[3],momentum[4],momentum[5]));


// if(bt<1.) bt=0.03*G4UniformRand();

 // bt=0.05*G4UniformRand();
 // bt=0.1;
 // bt=G4UniformRand()*0.07+0.07;
    particleGun1->SetParticleEnergy(bt*GeV);
    //particleGun2->SetParticleEnergy(bt2*GeV);

//cout<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<endl;
//cout<<"momentum: "<<momentum[0]<<" "<<momentum[3]<<endl;

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

   //**************************************************************
  /*
  Pos();
  elastic(momentum);
  //upunif(momentum);


  GetStartEnergy(bt1*1000.);
  GetStartAngleTheta(t1);

  //myEventWeight=myXS_test->GetWeight(t1);
  //cout<<"-------------"<<t1<<" "<<myEventWeight<<endl;
  //myEA->fEventWeight=myEventWeight;

  GetStartAnglePhi(fi1);
  GetStartPosition(vertex[0],vertex[1],vertex[2]);
  //cout<<"-------------"<<t1<<" "<<fi1<<endl;

  particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
  //particleGun1->SetParticleMomentumDirection(G4ThreeVector(momentum[0],momentum[1],momentum[2]));
  particleGun1->SetParticlePosition(G4ThreeVector(vertex[0],vertex[1],vertex[2]));


 particleGun1->SetParticleEnergy(bt1*GeV);
//particleGun1->SetParticleEnergy((sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]+0.938*0.938)-0.938)*GeV);
//particleGun1->SetParticleEnergy(bt1*GeV);

//cout<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<endl;
//cout<<"momentum: "<<momentum[0]<<" "<<momentum[2]<<endl;

// cout<<bt<<" "<<sqrt(momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]+0.938*0.938)-0.938<<endl;

    particleGun1->GeneratePrimaryVertex(anEvent);
    //particleGun2->GeneratePrimaryVertex(anEvent);

    */


}

void B4aPrimaryGeneratorAction::RandomInit(G4int level )
{
  CLHEP::Ranlux64Engine *Engine1 = new CLHEP::Ranlux64Engine(aj1bx,level);
  CLHEP::Ranlux64Engine *Engine2 = new CLHEP::Ranlux64Engine(aj1by,level);

  //CLHEP::RandFlat *GD1 = new CLHEP::RandFlat(*Engine1);
  //CLHEP::RandFlat *GD2 = new CLHEP::RandFlat(*Engine2);
  //GaussDist[0] = GD1;
  //GaussDist[1] = GD2;

  CLHEP::RandGauss *GD1 = new CLHEP::RandGauss(*Engine1);
  CLHEP::RandGauss *GD2 = new CLHEP::RandGauss(*Engine2);
  GaussDist[0] = GD1;
  GaussDist[1] = GD2;


  CLHEP::Ranlux64Engine *Engine3 = new CLHEP::Ranlux64Engine(aj1phi,level);  //init generator with seed
  CLHEP::Ranlux64Engine *Engine4 = new CLHEP::Ranlux64Engine(aj1theta,level);
  CLHEP::Ranlux64Engine *Engine5 = new CLHEP::Ranlux64Engine(ajcs1,level);
  CLHEP::Ranlux64Engine *Engine6 = new CLHEP::Ranlux64Engine(aju,level);
  CLHEP::Ranlux64Engine *Engine7 = new CLHEP::Ranlux64Engine(aj2theta,level);
  CLHEP::Ranlux64Engine *Engine8 = new CLHEP::Ranlux64Engine(aj1ekin,level);
  CLHEP::Ranlux64Engine *Engine9 = new CLHEP::Ranlux64Engine(aj2phi,level);
  CLHEP::Ranlux64Engine *Engine10 = new CLHEP::Ranlux64Engine(ajcs1,level);
  CLHEP::Ranlux64Engine *Engine11 = new CLHEP::Ranlux64Engine(ajbz,level);

  CLHEP::RandFlat *GD3 = new CLHEP::RandFlat(*Engine3);
  CLHEP::RandFlat *GD4 = new CLHEP::RandFlat(*Engine4);
  CLHEP::RandFlat *GD5 = new CLHEP::RandFlat(*Engine5);
  CLHEP::RandFlat *GD6 = new CLHEP::RandFlat(*Engine6);
  CLHEP::RandFlat *GD7 = new CLHEP::RandFlat(*Engine7);
  CLHEP::RandFlat *GD8 = new CLHEP::RandFlat(*Engine8);
  CLHEP::RandFlat *GD9 = new CLHEP::RandFlat(*Engine9);
  CLHEP::RandFlat *GD10 = new CLHEP::RandFlat(*Engine10);
  CLHEP::RandFlat *GD11 = new CLHEP::RandFlat(*Engine11);

  FlatDist[0] = GD3;
  FlatDist[1] = GD4;
  FlatDist[2] = GD5;
  FlatDist[3] = GD6;
  FlatDist[4] = GD7;
  FlatDist[5] = GD8;
  FlatDist[6] = GD9;
  FlatDist[7] = GD10;
  FlatDist[8] = GD11;
}


 G4double B4aPrimaryGeneratorAction::RandomGauss(G4double seed, G4double mean , G4double deviation )
{
  G4double num;
  if (seed == aj1bx)
  {
    num = (*GaussDist[0]).fire(mean, deviation);
    (*GaussDist[0]).fire(mean, deviation);
    //num = (*FlatDist[0]).fire(mean, deviation);
    //(*FlatDist[0]).fire(mean, deviation);
   return num;
  }
  if (seed == aj1by)
  {
    num = (*GaussDist[1]).fire(mean, deviation);
    (*GaussDist[1]).fire(mean, deviation);
    //num = (*GD2).fire(mean, deviation);
    //(*GD2).fire(mean, deviation);
    return num;
  }
  G4cout <<"Error in choice random gauss distribution!! Seed = "<<seed<<G4endl;
  exit (1);
}



 G4double B4aPrimaryGeneratorAction::RandomFlat(G4double seed, G4double m , G4double n  )
{
  if (seed == aj1phi)   return (*FlatDist[0]).fire(m, n);
  if (seed == aj1theta) return (*FlatDist[1]).fire(m, n);// metod fire(m,n) - > return double ]m,n[
  if (seed == ajcs1) 	return (*FlatDist[2]).fire(m, n);
  if (seed == aju) 	return (*FlatDist[3]).fire(m, n);
  if (seed == aj2theta) return (*FlatDist[4]).fire(m, n);
  if (seed == aj1ekin) 	return (*FlatDist[5]).fire(m, n);
  if (seed == aj2phi) 	return (*FlatDist[6]).fire(m, n);
  if (seed == ajcs1) 	return (*FlatDist[7]).fire(m, n);
  if (seed == ajbz) 	return (*FlatDist[8]).fire(m, n);

  G4cout <<"Error in choice random : "<<m<<"-"<<n<<" !! Seed = "<<seed<<G4endl;
  exit (1);
}


void B4aPrimaryGeneratorAction::Pos(void)
{
  G4double bsgx,bsgy;
  thigh = 0.5*mm;


  bsgx = bfwhmx/(2.*sqrt(2.*log(2.)));
  bsgy = bfwhmy/(2.*sqrt(2.*log(2.)));

  // recalculated by Geant to mm; bfwhmx, bfwhmy read in default mm.
  //vertex[0] = tXplace + RandomGauss(aj1bx,0, bsgx);


  //vertex[0] = tXplace + RandomGauss(aj1bx,-bsgx, bsgx);
  //vertex[1] = tYplace + RandomGauss(aj1by,-bsgy, bsgy);
  //tZplace=-10;
  //vertex[2] = tZplace - thigh + RandomFlat(ajbz)*2*thigh;


  // vertex[0] = tXplace + RandomGauss(aj1bx, 0, bsgx);
  // vertex[1] = tYplace + RandomGauss(aj1by,0, bsgy);

  vertex[0] = tXplace + G4RandGauss::shoot(0,1);
  vertex[1] = tYplace + G4RandGauss::shoot(0,1);
  vertex[2] = tZplace -thigh + G4UniformRand()*2*thigh;

  // cout<<"@#!@!#!@"<<G4RandGauss::shoot(0,1)<<endl;


  // vertex[0] = tXplace + RandomGauss(aj1bx, 0, bsgx);
  // vertex[1] = tZplace - thigh + RandomFlat(ajbz)*2*thigh;
  // vertex[2] = tYplace + RandomGauss(aj1by, 0, bsgy);

  //cout<<thigh<<endl;

  //cout<<"vertexZ: "<<vertex[2]<<endl;
 }


G4double* B4aPrimaryGeneratorAction::elastic(double *ptot)
{
  G4double pm, thcos, pphi, bp1, bpproj;
  G4double thcosm,thcos0;
  //cout<<"xxxxxxxxxxx"<<themin<<" "<<themax<<endl;

  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
  thcos = thcosm + (thcos0 - thcosm)*RandomFlat(aj1theta);
  t1 = 180./CLHEP::pi*acos(thcos);
  //pphi = CLHEP::pi*(2.* RandomFlat(aj1phi));
  pphi =(fimin + (fimax-fimin)*CLHEP::RandFlat::shoot())*CLHEP::pi/180.;

  //fi1=pphi*180./CLHEP::pi;
  fi1=pphi*180./CLHEP::pi;

  //if   (npd_choice == -1)
  pm = p_mass;  		  //proton
  //else if(npd_choice==-2) pm = d_mass;
  //else pm=n_mass;

  bp1 = sqrt(bt*(bt + 2*pm));
  ptot[2] = bp1*thcos;  			//ptot1(3)
  bpproj = sqrt(bp1*bp1 - ptot[2]*ptot[2]);
  ptot[0] = bpproj*cos(pphi);			//ptot1(1)
  ptot[1] = bpproj*sin(pphi);			//ptot1(2)
  bt1 = bt;
  //rs  t1 = 180./pi*themin;
  return ptot;
}


G4double* B4aPrimaryGeneratorAction::upunif(double *ptot)
{
  G4double ptot1[3],ptot2[3],bp1=0.,bp2=0.,bpproj;	// auxillary variable
  G4double thcosm,thcos0,thcosm2,thcos02,thcos,thcos2,pphi; // auxillary variable
  G4double *w_bp1,*w_bp2,p1,p2, m1, m2;

  w_bp1 = &bp1;
  w_bp2 = &bp2;

 m1=p_mass;
 m2=Ti48_mass;

  for (int i=0;i<3;i++)
  {
    ptot1[i] = ptot[i];
    ptot2[i] = ptot[i+3];
  }
  //cout<<"min1: "<<themin<<" max1: "<<themax<<endl;
  //cout<<"min2: "<<themin2<<" max2: "<<themax2<<endl;

  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
  thcosm2 = cos(themax2*CLHEP::pi/180.);
  thcos02 = cos(themin2*CLHEP::pi/180.);

  thcos = thcosm + (thcos0 - thcosm)*RandomFlat(aj1theta);
  //rs pphi = pi*(2.*RandomFlat(aj1phi)-1.);

//prev
  // pphi = CLHEP::pi*(2.*RandomFlat(aj1phi));

  pphi = fimin + (fimax - fimin)*RandomFlat(aj2phi);

  //cout<<"yyyyyyyyyyyyyyyy"<<pphi<<endl;

  //G4cout <<"thcos : "<<thcos<<G4endl;
  thcos2 = gelkin(thcos,w_bp1,w_bp2);	//elastic kinematics

  bt1 = sqrt(m1*m1 + bp1*bp1) - m1;
  bt2 = sqrt(m2*m2 + bp2*bp2) - m2;

  ptot1[2] = bp1*thcos;
  bpproj = sqrt(bp1*bp1 - ptot1[2]*ptot1[2]);
  ptot1[0] = bpproj*cos(pphi*CLHEP::pi/180.);
  ptot1[1] = bpproj*sin(pphi*CLHEP::pi/180.);

  ptot2[2] = bp2*thcos2;
  bpproj = sqrt(bp2*bp2 - ptot2[2]*ptot2[2]);
  ptot2[0] = bpproj*cos(CLHEP::pi + pphi*CLHEP::pi/180.);
  ptot2[1] = bpproj*sin(CLHEP::pi + pphi*CLHEP::pi/180.);

  t1 = 180./CLHEP::pi*acos(thcos);
  t2 = 180./CLHEP::pi*acos(thcos2);

  //cout<<"t1: "<<t1<<endl;
  //cout<<"t2: "<<t2<<endl;

  //rs
  fi1 =   pphi*180./CLHEP::pi;
  fi2 = fi1+180.;
  if(fi2>360)fi2=fi2-360;

  //cout<<"xxxxxx "<<fi1<<endl;


  for (int j=0;j<3;j++)
  {
    ptot[j] = ptot1[j];
    ptot[j+3] = ptot2[j];
  }
  return ptot;
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

double* B4aPrimaryGeneratorAction::ugelast(double* ptot)
{
// elastic scattering - random generation of proton emmision angle
// according to cross section and analysing power
//
// gelkin for calculation of deuteron emmision angle and of particle energies
ugelast_read();
  double thcosm,thcos0,thcosm2,thcos02,bp1,bp2,thelab,thcos,thsin;
  double difft=0.,difft1=0.,csi,cstest,it11,t20,t22,pphi=0.,czvec,cztens1,cztens2,fphi=0.;
  double bt1_temp=0.,bt2_temp=0.,thcos2=0.,ptot1[3],ptot2[3],bpproj,u,pm,dm;
  double csmax = 15.;
  int ithet,it1,i,ok;

  double *w_bp1,*w_bp2;

  t1 = 0.;
  t2 = 0.;

  w_bp1 = &bp1;
  w_bp2 = &bp2;
 /*
  //Target ISCh
  if ((npd_choice == 4)||(npd_choice == 7))
    {
      pm = d_mass;
    }
  else
    {
      pm = p_mass;
    }
  //Beam/projectile
  if (npd_choice == 7)
    {
      dm = p_mass;
    }
  else
    {
      dm = d_mass;
    }
    */
  //    cout<<"part1 :"<<pm;
  //    cout<<"part2:"<<dm;
  thcosm = cos(themax*CLHEP::pi/180.);
  thcos0 = cos(themin*CLHEP::pi/180.);
  thcosm2 = cos(themax2*CLHEP::pi/180.);
  thcos02 = cos(themin2*CLHEP::pi/180.);
//cout<<"------> START"<<endl;

  do
  {
    ok = 1;
    thelab = CLHEP::pi*0.5*RandomFlat(aj1theta);
    thcos = cos(thelab);

    if((thcos < thcosm) || (thcos > thcos0)) {ok = 0;continue;}
    thsin = sin(thelab);
    thelab *= 180./CLHEP::pi;
    ithet = 0;
    for(i=0;i<70;i++)
    {
      if(thelab <= thl[i] && thelab > thl[i+1])
      {
        ithet = i;
        difft = thl[i+1] - thl[i];
        difft1 = thelab - thl[i];
        break;
      }
    }
    if(fabs(difft) < 0.0001) difft = 0.0001;

// linear interpolation of unpolarized cross section
    csi  = ds[ithet] + (ds[ithet+1] - ds[ithet])/difft*difft1;
    csi  = csi * thsin;
    cstest=csmax*RandomFlat(ajcs1);
    //cout<<RandomFlat(ajcs1)<<endl;

    //cout<<csi<<" "<<cstest<<endl;
    if (csi < cstest) {ok =0; continue;};
/*
    if(fabs(pz) >= 0.0001 || fabs(pzz) >= 0.0001)
    {

// linear interpolation of analysing powers
      it11 = ap1[ithet] + (ap1[ithet+1] - ap1[ithet])/difft*difft1;
      t20  = ap2[ithet] + (ap2[ithet+1] - ap2[ithet])/difft*difft1;
      t22  = ap3[ithet] + (ap3[ithet+1] - ap3[ithet])/difft*difft1;

// von Neumann's method applied to ph
      do
      {
        pphi = CLHEP::pi*(2.*RandomFlat(aj1phi) - 1.);
        u = 2.*RandomFlat(aju);
        czvec = sqrt(3.)*pz*it11*cos(pphi);
        cztens1 = 1./sqrt(8.)*pzz*t20;
        cztens2 = - sqrt(3.)/2.*pzz*t22*cos(2.*pphi);
        fphi = 1. + sqrt(3.)*pz*it11*cos(pphi) + 1./sqrt(8.)*pzz*t20 - sqrt(3.)/2.*pzz*t22*cos(2.*pphi);
      }
      while(u > fphi);
    }
    else
    {
      pphi = CLHEP::pi*(2.*RandomFlat(aj1phi) - 1.);
      fphi = CLHEP::pi+pphi;
    }
*/
    bp1 = 0.;
    bp2 = 0.;
    thcos2 = gelkin(thcos,w_bp1,w_bp2);		//elastic kinematics

    bt1_temp = sqrt(pm*pm + bp1*bp1) - pm;
    bt2_temp = sqrt(dm*dm + bp2*bp2) - dm;

  }
  //while(((thcos2 > thcos02) || (thcos2 < thcosm2))||(ok == 0));
while((ok == 0));


  bt1 = bt1_temp;	//G4cout <<"Energy bt1 : "<<bt1*GeV<<G4endl;
  bt2 = bt2_temp;	//G4cout <<"Energy bt2 : "<<bt2*GeV<<G4endl;
  //    cout<<"Energy1tuuu :"<<bt1;
  //  cout<<"energy2tuuu:"<<bt2;
  ptot1[2] = bp1*thcos;
  bpproj = sqrt(bp1*bp1 - ptot1[2]*ptot1[2]);
  ptot1[0] = bpproj*cos(pphi);
  ptot1[1] = bpproj*sin(pphi);
  ptot2[2] = bp2*thcos2;
  bpproj = sqrt(bp2*bp2 - ptot2[2]*ptot2[2]);
  ptot2[0] = bpproj*cos(CLHEP::pi+pphi);
  ptot2[1] = bpproj*sin(CLHEP::pi+pphi);

  t1 = 180./CLHEP::pi*acos(thcos);
  t2 = 180./CLHEP::pi*acos(thcos2);
  fi1 = pphi*180./CLHEP::pi;
  fi2 = fphi*180./CLHEP::pi;

  it1 = (int)thelab + 1;

  for (i=0;i<3;i++)
  {
    ptot[i] = ptot1[i];
    ptot[i+3] = ptot2[i];
  }

  return ptot;
}

G4double B4aPrimaryGeneratorAction::gelkin(G4double thcos, G4double *wbp1, G4double *wbp2)
{  // output : bp1,bp2,thcos2
  G4double et,etp,beta,gamma,betsq,gamsq,rang3,sinth,sinsq; 	  	//auxillary variable
  G4double bp0,x1,x2,x3,x4,ep3c,ep4c,p3c,p4c,delta,a,b,c,test,coscm,thecm;  //auxillary variable
  G4double ep3l,ep4l,p3l,e3,e4,p4l,cosang4,coslb;  				  //auxillary variable
  G4double bej, ang4;
  G4double thcos2,m1,m2;

  m1=p_mass;
  m2=Ti48_mass;

  bp0 = sqrt(bt*(bt + 2.*m1));
  x1 = m1;
  x2 = m2;
  x3 = m1;
  x4 = m2;

  et = x1 + x2 + bt;
  etp = sqrt(x1*x1 + x2*x2 + 2.*x2*(x1 + bt));
  beta = bp0/et;
  gamma = et/etp;
  betsq = beta*beta;
  gamsq = gamma*gamma;
  rang3 = acos(thcos);
  sinth = sin(rang3);
  coslb=cos(rang3);
  sinsq = sinth*sinth;

// Compute relativistic energy of outgoing particles
  ep3c = (etp*etp + x3*x3 - x4*x4)/(2.0*etp);
  ep4c = etp - ep3c;
  p3c = ep3c*ep3c - x3*x3;
  p4c = ep4c*ep4c - x4*x4;
  if (p4c <= 0.) p4c = 0.;
  else p4c = sqrt(p4c);
  if (p3c <= 0.) p3c = 0.;
  else p3c = sqrt(p3c);
  delta = beta*ep3c/p3c;
  a = betsq*gamsq*sinsq + 1.;
  b = delta*gamsq*sinsq;
  c = (gamsq*delta*delta+1.)*sinsq - 1.;
  test = b*b - a*c;
  if (test < -1e-6*b*b) return 0.;
  test = sqrt(max(test,0.));
  coscm = (-b+test)/a;
if (fabs (fabs(coscm)-1.)< 1e-6){
    if (coscm > 1.) coscm = 1.;
    else coscm = - 1.;
  }

  if (coscm > 1.|| coscm > coslb ||coslb<=0.0 )
    {

      coscm = (-b-test)/a;

 if (fabs(fabs(coscm)-1.)< 1e-6){
    if (coscm > 1.) coscm = 1.;
    else coscm = - 1.;
  }
 if (coscm < -1.){

   return 0;

 }

    }

  if (coscm < -1.){

   return 0;

 }







  thecm = acos(coscm);
  ep3l = gamma*(ep3c + beta*p3c*coscm);
  ep4l = et - ep3l;
  p3l = sqrt(ep3l*ep3l - x3*x3);
  e3 = ep3l - x3;
  e4 = bt - e3;
  bej = e3;
  p3l = sqrt(ep3l*ep3l - x3*x3);
  p4l = ep4l*ep4l - x4*x4;
  if (p4l <= 0.) p4l = 0.;
  else p4l = sqrt(p4l);
  *wbp1 = p3l;
  *wbp2 = p4l;
  if (p4l == 0.) cosang4 = 0.;
  else
  {
    cosang4 = gamma*(beta*ep4c - p3c*coscm)/p4l;
    if (fabs(cosang4) > 1.)
    {
      if (cosang4 > 1.) cosang4 = 1.;
      else cosang4 = - 1.;
    }
  }
  ang4 = acos(cosang4) * 180./CLHEP::pi;
  thcos2 = cosang4;

  //G4cout<<"thcos : "<<ang4<<endl;
  return thcos2;
}


void B4aPrimaryGeneratorAction::ugelast_read(void)
{
  ifstream file;
  file.open("/home/ciepal/Main/myGEANT/Kratta/src/sig_anpow.dat");
  if (!file)
  {
    G4cout <<"Cannot open the file : sig_anpow.dat !!!"<<G4endl;
    exit(1);
   }
  for (int i=0;i<70;i++)
  {
    file >> bt1_mev>>th>>thl[i]>>ds[i]>>ay[i]>>ap1[i]>>ap2[i]>>t21>>ap3[i];

    if(file.eof())
    {
      G4cout << " Unexpected end of file : /home/sworst/Breakup/data/sig_anpow.dat !!! Index : "<<i<<G4endl;
      exit(1);
    }
    if (file.fail())
    {
      G4cout << "Unexpected error in file : /home/sworst/Breakup/data/sig_anpow.dat !!! Index : "<<i<<G4endl;
      exit(1);
    }
  }
  file.close();
}
