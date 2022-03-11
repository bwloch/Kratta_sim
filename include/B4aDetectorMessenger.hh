#ifndef B4aDetectorMessenger_h
#define B4aDetectorMessenger_h 1


#include "globals.hh"
#include "G4UImessenger.hh"

class B4aDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class B4aDetectorMessenger: public G4UImessenger
{


  public:
    B4aDetectorMessenger(B4aDetectorConstruction* Det);
   ~B4aDetectorMessenger();
   

    //void SetNewValue(G4UIcommand*, G4String);
  void SetNewValue(G4UIcommand * command, G4String newValues);
  
  private:
    B4aDetectorConstruction* myDetector;


  
  G4UIdirectory*  paramDir;
  G4UIdirectory	*generationDir, *TargetDir;


  //G4UIcmdWithAString*        MatCmd;

 // G4UIcmdWithAString*        FoilCmd;
//////////////////////////////
// Generation Parameters
//////////////////////////////

  //G4UIcmdWithAnInteger* NpdChoiceCmd;
  G4UIcmdWithAnInteger* IfNeumannCmd;
  G4UIcmdWithADouble* BfwhmX_Cmd;
  G4UIcmdWithADouble* BfwhmY_Cmd;
  G4UIcmdWithADoubleAndUnit* BtEnergyCmd;
  G4UIcmdWithADouble* Pz_Cmd;
  G4UIcmdWithADouble* Pzz_Cmd;
  G4UIcmdWithADouble* GenMinCmd;
  G4UIcmdWithADouble* GenMaxCmd;
  
  G4UIcmdWithADoubleAndUnit* ThetaMinCmd;
  G4UIcmdWithADoubleAndUnit* ThetaMaxCmd;
  G4UIcmdWithADoubleAndUnit* Theta2MinCmd;
  G4UIcmdWithADoubleAndUnit* Theta2MaxCmd;
  G4UIcmdWithADoubleAndUnit* PhiMinCmd;
  G4UIcmdWithADoubleAndUnit* PhiMaxCmd;
  G4UIcmdWithoutParameter* ParamUpdate;
  
  
//////////////////////////////
// Target
//////////////////////////////

  //G4UIcmdWithAnInteger* TargetIsCmd;
  G4UIcmdWithADoubleAndUnit*TframeOutRadiusCmd;
  G4UIcmdWithADoubleAndUnit*TargetOutRadiusCmd;
  G4UIcmdWithADoubleAndUnit* TargetHighCmd ;
// Target placement
  G4UIcmdWithADoubleAndUnit* TargetXplaceCmd ;
  G4UIcmdWithADoubleAndUnit* TargetYplaceCmd ;
  G4UIcmdWithADoubleAndUnit* TargetZplaceCmd ;



};

#endif

