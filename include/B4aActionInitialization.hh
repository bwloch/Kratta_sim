
#ifndef B4aActionInitialization_h
#define B4aActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class B4aDetectorConstruction;
class B4aPrimaryGeneratorAction;
class HistoManager ;
/// Action initialization class.
///

class B4aActionInitialization : public G4VUserActionInitialization
{
  public:
    //B4aActionInitialization(B4aDetectorConstruction*, B4aPrimaryGeneratorAction*, HistoManager*);
    B4aActionInitialization(B4aDetectorConstruction*);
    
    virtual ~B4aActionInitialization();

    virtual void BuildForMaster() const;
    virtual void Build() const;

  private:
    //B4aDetectorConstruction* fDetConstruction;
    B4aDetectorConstruction* myDetector;
    B4aPrimaryGeneratorAction* myPGA;
    HistoManager* myHisto;


};

#endif

    
