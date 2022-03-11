#include "B4aActionInitialization.hh"
#include "B4aPrimaryGeneratorAction.hh"
#include "B4aRunAction.hh"
#include "B4aEventAction.hh"
#include "B4aSteppingAction.hh"
#include "B4aDetectorConstruction.hh"
#include "B4aNucleonElasticXS.hh"
#include "B4aHistoManager.hh"



B4aActionInitialization::B4aActionInitialization
                            (B4aDetectorConstruction* detConstruction
			     //B4aPrimaryGeneratorAction* primGenAct
			     //HistoManager* histo
                            )
 : //G4VUserActionInitialization(),
   myDetector(detConstruction)
   //myPGA(primGenAct),
   //myHisto(histo)
{}



B4aActionInitialization::~B4aActionInitialization()
{}



void B4aActionInitialization::BuildForMaster() const
{

  HistoManager *myHisto= new HistoManager();
  SetUserAction(new B4aRunAction(myHisto));
}



void B4aActionInitialization::Build() const
{

  HistoManager *myHisto= new HistoManager();
  B4aPrimaryGeneratorAction *myPGA =  new B4aPrimaryGeneratorAction(myDetector,myHisto);
  //B4aNucleonElasticXS *myNEXS =  new B4aPrimaryGeneratorAction(myDetector);

  SetUserAction(myPGA);
  SetUserAction(new B4aRunAction(myHisto));

  B4aEventAction* eventAction = new B4aEventAction(myHisto);//, myPGA);
  SetUserAction(eventAction);
  //SetUserAction(new B4aSteppingAction(myDetector,eventAction));
  SetUserAction(new B4aSteppingAction(myDetector,eventAction,myHisto));

  //G4RunManager::GetRunManager()->SetUserAction(new B4aPrimaryGeneratorAction(myDetector));


}
