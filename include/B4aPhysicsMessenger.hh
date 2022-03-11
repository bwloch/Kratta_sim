#ifndef B4aPhysicsMessenger_h
#define B4aPhysicsMessenger_h 1

class B4aPhysicsList;
class B4aSteppingAction;

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIcmdWithoutParameter;
//class DetectorConstruction;
class G4UIdirectory;
//class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
//class G4UIcmdWithAString;
//class G4UIcmdWithoutParameter;

class B4aPhysicsMessenger: public G4UImessenger
{


  public:
   

B4aPhysicsMessenger(B4aPhysicsList* myPhysicsList);
   ~B4aPhysicsMessenger();

  void SetNewValue(G4UIcommand * command, G4String newValues);
  private:
  B4aPhysicsList* phyList;
  G4UIdirectory*  physicsDir;
  G4UIdirectory*  neutronsDir;
  G4UIdirectory*  fileOutDir;

//////////////////////////////
// Neutron Parameters
//////////////////////////////

  G4UIcmdWithAnInteger* NeutronTypeCmd;
  G4UIcmdWithAnInteger* NeutronModelCmd;
  G4UIcmdWithAnInteger* energyBroadeningCmd;
  G4UIcmdWithAnInteger* filesCmd;
};

#endif

