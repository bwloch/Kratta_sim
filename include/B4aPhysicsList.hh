#ifndef B4aPhysicsList_h
#define B4aPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "G4VModularPhysicsList.hh"

#include "globals.hh"



class B4aPhysicsList:public G4VUserPhysicsList
  //public G4VModularPhysicsList
//public G4VUserPhysicsList 
{
public:
  B4aPhysicsList();
  ~B4aPhysicsList();

public:
  virtual void SetCuts();
  //virtual void RegisterHadrons(G4String);

protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
    
  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  //virtual void ConstructHad();
  virtual void ConstructOp();


  virtual void AddTransportation();

private:
  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForProton;
  G4double cutForDeuteron;
  G4double cutForNeutron;


  // these methods Construct particles 
  void ConstructMyBosons();
  void ConstructMyLeptons();
  void ConstructMyHadrons();
  void ConstructMyShortLiveds();

};

#endif
