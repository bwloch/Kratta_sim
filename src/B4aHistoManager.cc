#include "B4aHistoManager.hh"
#include "B4aRunAction.hh"

#include "G4UnitsTable.hh"

using namespace std;



HistoManager::HistoManager()
: fFileName("steppingAction")
{
  MaxNtCol=50;
 // Book();

for (G4int i=0;i<MaxNtCol;i++){
fNtColId[i]=0;

}

}


HistoManager::~HistoManager()
{
 // delete G4AnalysisManager::Instance();
}

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    // analysisManager->SetFileName(fFileName);
	//G4bool fileOpen=
	analysisManager->OpenFile(fFileName);
	//analysisManager->SetVerboseLevel(7);

  //analysisManager->SetActivation(true);     //enable inactivation of histograms

  //Ntuple
  analysisManager->SetFirstNtupleId(1);

  analysisManager->CreateNtuple("T","myNtupla");

  fNtColId[0]=analysisManager->CreateNtupleDColumn(1,"edep_target"); //0

  //Kratta
  fNtColId[1]=analysisManager->CreateNtupleDColumn(1,"edep_krHouse");//1
  fNtColId[2]=analysisManager->CreateNtupleDColumn(1,"kr_nb");//2
  fNtColId[3]=analysisManager->CreateNtupleDColumn(1,"edep_krWin");//3
  fNtColId[4]=analysisManager->CreateNtupleDColumn(1,"edep_air1");//4
  //fNtColId[10]=analysisManager->CreateNtupleDColumn(1,"edep_PD");//

  fNtColId[5]=analysisManager->CreateNtupleDColumn(1,"edep_PDdeadF0");//5
  fNtColId[6]=analysisManager->CreateNtupleDColumn(1,"edep_PDdeadR0");//6
  fNtColId[7]=analysisManager->CreateNtupleDColumn(1,"edep_PDactive0");//7

  fNtColId[8]=analysisManager->CreateNtupleDColumn(1,"edep_PDdeadF1");//8
  fNtColId[9]=analysisManager->CreateNtupleDColumn(1,"edep_PDdeadR1");//9
  fNtColId[10]=analysisManager->CreateNtupleDColumn(1,"edep_PDactive1");//10


  fNtColId[11]=analysisManager->CreateNtupleDColumn(1,"edep_tg2");//11
  fNtColId[12]=analysisManager->CreateNtupleDColumn(1,"weights");//12
  fNtColId[13]=analysisManager->CreateNtupleDColumn(1,"edep_air2");//13


  fNtColId[14]=analysisManager->CreateNtupleDColumn(1,"edep_CsIshort");//14
  //fNtColId[17]=analysisManager->CreateNtupleDColumn(1,"edep_CsIshortwr");//14
  fNtColId[15]=analysisManager->CreateNtupleDColumn(1,"edep_csi_wr");//15
  fNtColId[16]=analysisManager->CreateNtupleDColumn(1,"edep_CsIlong");//16
  fNtColId[17]=analysisManager->CreateNtupleDColumn(1,"eventNum");//17

  fNtColId[18]=analysisManager->CreateNtupleDColumn(1,"edep_air3");//18
  fNtColId[19]=analysisManager->CreateNtupleDColumn(1,"particleID");//19

    fNtColId[20]=analysisManager->CreateNtupleDColumn(1,"primaryVx");//20
    fNtColId[21]=analysisManager->CreateNtupleDColumn(1,"primaryVy");//21
    fNtColId[22]=analysisManager->CreateNtupleDColumn(1,"primaryVz");//22

    fNtColId[23]=analysisManager->CreateNtupleDColumn(1,"momentumX");//23
    fNtColId[24]=analysisManager->CreateNtupleDColumn(1,"momentumY");//24
    fNtColId[25]=analysisManager->CreateNtupleDColumn(1,"momentumZ");//25
    fNtColId[26]=analysisManager->CreateNtupleDColumn(1,"momentumTot");//26

    fNtColId[27]=analysisManager->CreateNtupleDColumn(1,"energyTot");//27
    fNtColId[28]=analysisManager->CreateNtupleDColumn(1,"kinEnergyTot");//28
    fNtColId[29]=analysisManager->CreateNtupleDColumn(1,"theta");//29
    fNtColId[30]=analysisManager->CreateNtupleDColumn(1,"phi");//30

    fNtColId[31]=analysisManager->CreateNtupleDColumn(1,"world_air");//31
    analysisManager->CreateNtupleDColumn(1,"x_posKrHousTr");//32
    analysisManager->CreateNtupleDColumn(1,"y_posKrHousTr");//33
    analysisManager->CreateNtupleDColumn(1,"z_posKrHousTr");//34
    analysisManager->CreateNtupleDColumn(1,"x_momKrHousSt");//35
    analysisManager->CreateNtupleDColumn(1,"y_momKrHousSt");//36
    analysisManager->CreateNtupleDColumn(1,"z_momKrHousSt");//37
    fNtColId[38]=analysisManager->CreateNtupleDColumn(1,"edep_Plastic");//38


/*
    analysisManager->CreateNtupleDColumn(1,"Ekin_TgPost");//38
    analysisManager->CreateNtupleDColumn(1,"Ekin_AirPrev");//39
    //analysisManager->CreateNtupleDColumn(1,"Ekin_AirPost");
    analysisManager->CreateNtupleDColumn(1,"Ekin_KrWinPrev");//40
    //analysisManager->CreateNtupleDColumn(1,"Ekin_KrWinPost");
    analysisManager->CreateNtupleDColumn(1,"Ekin_Air1Prev");//41
    //analysisManager->CreateNtupleDColumn(1,"Ekin_Air1Post");
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd0fdPrev");//42
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd0acPrev");//43
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd0rdPrev");//44

    analysisManager->CreateNtupleDColumn(1,"Ekin_Air2Prev");//45
    //analysisManager->CreateNtupleDColumn(1,"Ekin_Air2Post");
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd1fdPrev");//46
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd1acPrev");//47
    analysisManager->CreateNtupleDColumn(1,"Ekin_Pd1rdPrev");//48


    analysisManager->CreateNtupleDColumn(1,"Ekin_CsIsPrev");//49
    //analysisManager->CreateNtupleDColumn(1,"Ekin_CsIsPost");
    analysisManager->CreateNtupleDColumn(1,"Ekin_WrPrev");//50
    //analysisManager->CreateNtupleDColumn(1,"Ekin_WrPost");

    analysisManager->CreateNtupleDColumn(1,"Ekin_CsIlPrev");//51
    //analysisManager->CreateNtupleDColumn(1,"Ekin_CsIlPost");
*/
	analysisManager->CreateNtupleDColumn(1,"Ekin_CsIsPrev");//39
	analysisManager->CreateNtupleDColumn(1,"Ekin_CsIsPost");//40
  analysisManager->CreateNtupleDColumn(1,"OurprimaryVx"); //41
  analysisManager->CreateNtupleDColumn(1,"OurprimaryVy");//42
  analysisManager->CreateNtupleDColumn(1,"OurprimaryVz");//43
  analysisManager->CreateNtupleDColumn(1,"Ourtheta");//44
  analysisManager->CreateNtupleDColumn(1,"Ourphi");//45
  analysisManager->CreateNtupleDColumn(1,"Ourenergy");//46
    analysisManager->CreateNtupleDColumn(1,"Ekin_KrHouse");//52

    analysisManager->FinishNtuple(1);




  // Define histograms start values
 /*
  const G4int kMaxHisto = 13;
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9",
                         "10","11","12"};
  const G4String title[] =
                { "dummy",             z                             //0
                  "kinetic energy of scattered primary particle",   //1
                  "kinetic energy of gamma",                        //2
                  "kinetic energy of neutrons",                     //3
                  "kinetic energy of protons",                      //4
                  "kinetic energy of deuterons",                    //5
                  "kinetic energy of alphas",                       //6
                  "kinetic energy of nuclei",                       //7
                  "kinetic energy of mesons",                       //8
                  "kinetic energy of baryons",                      //9
                  "Q = Ekin out - Ekin in",                         //10
                  "Pbalance = mag(P_out - P_in)",                   //11
                  "atomic mass of nuclei"                           //12
                 };
*/
  // Default values (to be reset via /analysis/h1/set command)
  //G4int nbins = 100;
  //G4double vmin = 0.;
  //G4double vmax = 100.;

  // Create all histograms as inactivated
  // as we have not yet set nbins, vmin, vmax
  //for (G4int k=0; k<kMaxHisto; k++) {
    //G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    //analysisManager->SetH1Activation(ih, false);
  //}
}
void HistoManager::save(){
  //if (fFactoryOn) {
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
    G4cout << "\n----> Histogram Tree is saved in " << fFileName << G4endl;

    delete G4AnalysisManager::Instance();
    //fFactoryOn = false;
    //}
}
//void HistoManager::FillNtuple(G4double energyAbs, G4double energyGap,
//                            G4double trackLAbs, G4double trackLGap)

void HistoManager::FillNtuple(G4int id, G4double value)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 1)


  //analysisManager->FillNtupleDColumn(1,fNtColId[id], value);
//  cout<<"------>> id: "<<id<<" value: "<<value<<endl;
  analysisManager->FillNtupleDColumn(1,id, value);


  //if(id==11)cout<<"------>>--------------------------->> "<<value<<" "<<id<<endl;
  //analysisManager->AddNtupleRow(1);

  // Fill 2nd ntuple ( id = 2)
  //  analysisManager->FillNtupleDColumn(2,fNtColId[2], trackLAbs);
  //analysisManager->FillNtupleDColumn(2,fNtColId[3], trackLGap);
  //analysisManager->AddNtupleRow(2);
}

void HistoManager::AddRow(G4int id)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Fill 1st ntuple ( id = 1)


  //analysisManager->FillNtupleDColumn(1,fNtColId[id], value);
  cout<<"------>> Add row"<<endl;
  analysisManager->AddNtupleRow(id);


  //if(id==11)cout<<"------>>--------------------------->> "<<value<<" "<<id<<endl;
  //analysisManager->AddNtupleRow(1);

  // Fill 2nd ntuple ( id = 2)
  //  analysisManager->FillNtupleDColumn(2,fNtColId[2], trackLAbs);
  //analysisManager->FillNtupleDColumn(2,fNtColId[3], trackLGap);
  //analysisManager->AddNtupleRow(2);
}
