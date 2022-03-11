
#include "B4aDetectorConstruction.hh"
#include "B4aPrimaryGeneratorAction.hh"
#include "B4aActionInitialization.hh"
#include "B4aPhysicsList.hh"
#include "B4aHistoManager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#ifdef G4VIS_USE
//#include "B4aVisManager.hh"
#endif
#include <sys/stat.h>





namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " exampleB4a [-m macro ] [-u UIsession] [-t nThreads]" << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}



int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;
  G4String fileInp;
  G4String command = "/control/execute ";

#ifdef G4MULTITHREADED
  G4int nThreads = 6;



#endif
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) {
    runManager->SetNumberOfThreads(nThreads);
  }
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  // Set mandatory initialization classes
  //
  B4aDetectorConstruction* detConstruction = new B4aDetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  auto physicsList = new FTFP_BERT;
  physicsList->SetCutValue(0.001*CLHEP::mm,"proton");
  physicsList->SetCutValue(0.001*CLHEP::mm,"e-");
  physicsList->SetCutValue(0.001*CLHEP::mm,"neutron");
  runManager->SetUserInitialization(physicsList);

   // B4aPhysicsList* physicsList = new B4aPhysicsList();
   // runManager->SetUserInitialization(physicsList);


// HistoManager* myHisto = new HistoManager();
//
//  B4aPrimaryGeneratorAction* myPGA = new B4aPrimaryGeneratorAction(detConstruction,myHisto);



 // runManager->SetUserAction(myPGA);


  //runManager->SetUserInitialization(myPGA);


  B4aActionInitialization* actionInitialization
     //= new B4aActionInitialization(detConstruction,myPGA,myHisto);
     = new B4aActionInitialization(detConstruction);

  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  //
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  (argc>2) ? fileInp = argv[2]:fileInp = "myinit.mac";
  {
    struct stat fileData;
    if(stat(fileInp, &fileData)==-1)
      perror (fileInp), exit(1);
    G4cout<<"\n\t Dodano INIT plik MYINIT.MAC\n";
    G4cout<<" DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n DODANO INIT \n";
    UImanager->ApplyCommand(command+fileInp);
  }





  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    //G4String command = "/control/execute ";
  //  G4cout<<"MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n MAKRO COMMAND \n ";
  //  UImanager->ApplyCommand(command+macro);
  }
  else  {
    // interactive mode : define UI session
     UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      //UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}
