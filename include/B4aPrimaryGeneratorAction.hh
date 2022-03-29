#ifndef B4aPrimaryGeneratorAction_h
#define B4aPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "B4aPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "TRandom3.h"
#include "TGraph.h"
#include <iostream>
#include <fstream>

#include "globals.hh"


class G4ParticleGun;
class G4Event;
class B4aDetectorConstruction;
class G4ParticleDefinition;
class Ranlux64Engine;
class RandFlat;
class B4aNucleonElasticXS;
class B4aEventAction;
class HistoManager;

class B4aPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  B4aPrimaryGeneratorAction(B4aDetectorConstruction*, HistoManager*);
  virtual ~B4aPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);

  G4double bt1,bt2,bt3,t1,t2,t3,fi1,fi2,fi3;
  G4float bt1_mev,th,thl[70],ds[70],ay[70],ap1[70],ap2[70],t21,ap3[70];
  G4float sig[20][20][15][200],ay0[20][20][15][200];
    //axx[20][20][15][200],ayy[20][20][15][200],axn[20][20][15][200],
    //	  axy[20][20][15][200],axz[20][20][15][200];

    G4int aj1bx, aj1by, ajbz, aj1theta, aj2theta, aj1phi, aj2phi, aj1ekin, aja10,
    	aja20, aju, ajcs1, ajcs02, aj1r1, aj1r2, npd_choice, icros ;
    G4double bfwhmx, bfwhmy, bt, pz, pzz, themin, themax, themin2, themax2, fimin,fimax;
    G4double tXplace, tYplace,tZplace,thigh;
    G4double generator_min, generator_max;
    TGraph *gr;
  //ifstream *file;
  G4double vertex[3];
  G4double GetVertexX(void){return vertex[0];} ;
  G4double GetVertexY(void){return vertex[1];} ;
  G4double GetVertexZ(void){return vertex[2];} ;

  B4aNucleonElasticXS *myXS_test;
  //B4aEventAction *myEA;

  //G4double myEventWeight;
  //virtual G4double GetWeight(void){return myEventWeight;};

  inline static G4double* GetStartEnergy (G4double en1 = -1., G4double en2 = -1., G4double en3 = -1.)
    {
      static G4double tes1[3];
      if (en1 != -1.)
      {
        tes1[0] = en1;
	tes1[1] = en2;
	tes1[2] = en3;
      }
      return tes1;
    };

    inline static G4double* GetStartAngleTheta (G4double ang1 = -1., G4double ang2 = -1., G4double ang3 = -1.)
    {
      static G4double tes2[3];
      if (ang1 != -1.)
      {
        tes2[0] = ang1;
	tes2[1] = ang2;
	tes2[2] = ang3;
      }
      return tes2;
    };

    inline static G4double* GetStartAnglePhi (G4double ang1 = -1., G4double ang2 = -1., G4double ang3 = -1.)
    {
      static G4double tes3[3];
      if (ang1 != -1.)
      {
        tes3[0] = ang1;
	tes3[1] = ang2;
	tes3[2] = ang3;
      }
      return tes3;
    };

    inline static G4double* GetStartPosition (double v1 = 0., double v2 = 0., double v3 = 0.)
    {
      static G4double tes3[3];
          if ((v1 != 0.)||(v2 != 0.)||(v3 != 0.))
      {
        tes3[0] = v1;	tes3[1] = v2;	tes3[2] = v3;
      }
      return tes3;
    }

    inline static G4int ProcNb(G4int num = 10)
    {
      static G4int temp;
      if (num != 10) temp = num;
      return temp;
    };

    inline static G4int GetChoice (G4int num = 10)
    {
      static G4int temp;
      if (num != 10) temp = num;
      return temp;
    };


  // set methods
  //void SetRandomFlag(G4bool value);

private:
  //G4ParticleGun*  fParticleGun; // G4 particle gun
  // G4RandGauss *rnd;
  G4ParticleGun* particleGun1;		//part 1
  G4ParticleGun* particleGun2;		//part 2

  G4ParticleDefinition *particleDefinition;
  B4aDetectorConstruction* myDetector;
  HistoManager *fHistoManager;


  void Pos(void);				//generate vertex position


    G4ThreeVector flatGen();


    G4double momentum[9];

    G4double p_mass, Ti48_mass, he3_mass, d_mass, he4_mass, t_mass;




};



#endif
