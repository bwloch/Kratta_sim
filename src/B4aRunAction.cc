#include "B4aRunAction.hh"
#include "B4aAnalysis.hh"
#include "B4aHistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



B4aRunAction::B4aRunAction(HistoManager * histo)
 : G4UserRunAction(), fHistoManager(histo)
{

  //fHistoManager=new HistoManager();
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  //analysisManager->SetVerboseLevel(1);
  //analysisManager->SetFirstHistoId(1);

  // Book histograms, ntuple
  //

  // Creating histograms
    //analysisManager->CreateH1("1","vx", 100, -6*mm, 6*mm);
  //analysisManager->CreateH1("2","vy", 100, -6*mm, 6*mm);
  //analysisManager->CreateH1("3","vz", 100, -52*mm, 52*mm);
  //analysisManager->CreateH1("4","trackL in gap", 100, 0., 50*cm);

  // Creating ntuple
  //
  //analysisManager->CreateNtuple("myNtupla", "myNtupla");
  //analysisManager->CreateNtupleDColumn("vertexX");
  //analysisManager->CreateNtupleDColumn("vertexY");
  //analysisManager->CreateNtupleDColumn("vertexZ");
  //analysisManager->CreateNtupleDColumn("Lgap");
  //analysisManager->FinishNtuple();
//fHistoManager->Book();
}



B4aRunAction::~B4aRunAction()
{
  delete G4AnalysisManager::Instance();
  //delete fHistoManager;
}



void B4aRunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  //if (analysisManager->IsActive() )analysisManager->OpenFile();
  // Open an output file
  //
   //G4String fileName = "polHe3";
   //analysisManager->OpenFile(fileName);

  fHistoManager->Book();

}



void B4aRunAction::EndOfRunAction(const G4Run* /*run*/)
{
fHistoManager->save();
  //save
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // //if (analysisManager->IsActive()){
  //   analysisManager->Write();
  //   analysisManager->CloseFile();

    //}



 // if ( analysisManager->GetH1(1) ) {
    //G4cout << G4endl << " ----> print histograms statistic ";
    //if(isMaster) {

      //fHistoManager->Book();
      //G4cout << "for the entire run " << G4endl << G4endl;
    //}
/*    else {
      G4cout << "for the local thread " << G4endl << G4endl;
    }
    */
    /*
    G4cout << " EAbs : mean = "
       << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

    G4cout << " EGap : mean = "
       << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy")
       << " rms = "
       << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;

    G4cout << " LAbs : mean = "
      << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

    G4cout << " LGap : mean = "
      << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length")
      << " rms = "
      << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
 */
  //fHistoManager->save();

  // save histograms & ntuple
  //
  //analysisManager->Write();
  //analysisManager->CloseFile();

}
