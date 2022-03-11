#include <iomanip>  

#include "B4aPhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4RegionStore.hh"
#include "G4ios.hh"
#include "G4UserLimits.hh"
#include "B4aNucleonElasticXS.hh"

B4aPhysicsList::B4aPhysicsList() : G4VUserPhysicsList() 
{

  //defaultCutValue     = 10*km; 
  //defaultCutValue     = 0.07*mm; 
  defaultCutValue     = 0.0001*mm; 
  
  //defaultCutValue     = 4.*mm; 
  
  cutForGamma         = defaultCutValue;
  cutForElectron      = defaultCutValue;
  cutForDeuteron      = defaultCutValue;
  cutForNeutron      = defaultCutValue;
  cutForProton=defaultCutValue;  
  VerboseLevel = 1;
  OpVerbLevel = 0;

  //SetCutsWithDefault();
  SetVerboseLevel(VerboseLevel);
}



B4aPhysicsList::~B4aPhysicsList() 
{;}



void B4aPhysicsList::ConstructParticle() 
{

   ConstructMyBosons();
  ConstructMyLeptons();
  ConstructMyHadrons();
  ConstructMyShortLiveds();

}



void B4aPhysicsList::ConstructMyBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();

  //OpticalPhotons
  G4OpticalPhoton::OpticalPhotonDefinition();

}



void B4aPhysicsList::ConstructMyLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}


#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// construct Hadrons://///////////////////////////////////////////////////
void B4aPhysicsList::ConstructMyHadrons()
{
 //  mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

 //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

 //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

}

#include "G4ShortLivedConstructor.hh"


void B4aPhysicsList::ConstructMyShortLiveds()
{
  G4ShortLivedConstructor slConstructor;
  slConstructor.ConstructParticle();
}




// Construct Processes //////////////////////////////////////////////////////
void B4aPhysicsList::ConstructProcess() 
{

  AddTransportation();
  ConstructEM();


  //ConstructOp();
  //ConstructHad();

  ConstructGeneral();
 
}


// Transportation ///////////////////////////////////////////////////////////
//#include "DMXMaxTimeCuts.hh"
//#include "DMXMinEkineCuts.hh"
#include "G4StepLimiter.hh"

void B4aPhysicsList::AddTransportation() {
  
  G4VUserPhysicsList::AddTransportation();
  /*
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // time cuts for ONLY neutrons:
    if(particleName == "neutron") 
      pmanager->AddDiscreteProcess(new DMXMaxTimeCuts());
    // Energy cuts to kill charged (embedded in method) particles:
    pmanager->AddDiscreteProcess(new DMXMinEkineCuts());

    // Step limit applied to all particles:
    pmanager->AddProcess(new G4StepLimiter,       -1,-1,1);

  }		    
  */  
}


// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"


// e-
#include "G4eMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"


// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"


// alpha and GenericIon and deuterons, triton, He3:

//muon:
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCapture.hh"

//OTHERS:
#include "G4hIonisation.hh" 
#include "G4hMultipleScattering.hh"
#include "G4hBremsstrahlung.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
#include "G4EnergyLossTables.hh"

//em process options to allow msc step-limitation to be switched off
#include "G4EmProcessOptions.hh"

void B4aPhysicsList::ConstructEM() {
  
  //set a finer grid of the physic tables in order to improve precision
  //former LowEnergy models have 200 bins up to 100 GeV
  //G4EmProcessOptions opt;
  //opt.SetMaxEnergy(100*GeV);
  //opt.SetDEDXBinning(200);
  //opt.SetLambdaBinning(200);

  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();
    
    if (particleName == "gamma") 
      {
	//gamma
	G4RayleighScattering* theRayleigh = new G4RayleighScattering();
	theRayleigh->SetEmModel(new G4LivermoreRayleighModel());  //not strictly necessary
	pmanager->AddDiscreteProcess(theRayleigh);

	G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
	thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
	pmanager->AddDiscreteProcess(thePhotoElectricEffect);
	
	G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
	theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
	pmanager->AddDiscreteProcess(theComptonScattering);
	
	G4GammaConversion* theGammaConversion = new G4GammaConversion();
	theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
	pmanager->AddDiscreteProcess(theGammaConversion);

      } 
    else if (particleName == "e-") 
      {
	//electron
	// process ordering: AddProcess(name, at rest, along step, post step)
	// Multiple scattering
	G4eMultipleScattering* msc = new G4eMultipleScattering();
	msc->SetStepLimitType(fUseDistanceToBoundary);
	pmanager->AddProcess(msc,-1, 1, 1);

	// Ionisation
	G4eIonisation* eIonisation = new G4eIonisation();
	eIonisation->SetEmModel(new G4LivermoreIonisationModel());
	eIonisation->SetStepFunction(0.2, 100*um); //improved precision in tracking  
	pmanager->AddProcess(eIonisation,-1, 2, 2);
	
	// Bremsstrahlung
	G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
	eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel());
	pmanager->AddProcess(eBremsstrahlung, -1,-3, 3);
      } 
    else if (particleName == "e+") 
      {
	//positron	
	G4eMultipleScattering* msc = new G4eMultipleScattering();
	msc->SetStepLimitType(fUseDistanceToBoundary);
	pmanager->AddProcess(msc,-1, 1, 1);
	
	// Ionisation
	G4eIonisation* eIonisation = new G4eIonisation();
	eIonisation->SetStepFunction(0.2, 100*um);     
	pmanager->AddProcess(eIonisation,                 -1, 2, 2);

	//Bremsstrahlung (use default, no low-energy available)
	pmanager->AddProcess(new G4eBremsstrahlung(), -1,-1, 3);

	//Annihilation
	pmanager->AddProcess(new G4eplusAnnihilation(),0,-1, 4);      
      } 
    else if (particleName == "proton" || 
	     particleName == "pi+" || 
	     particleName == "pi-")
      {
	//multiple scattering
	pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      
	//ionisation
	G4hIonisation* hIonisation = new G4hIonisation();
	hIonisation->SetStepFunction(0.2, 50*um);
	pmanager->AddProcess(hIonisation,                     -1, 2, 2);      
	
	//bremmstrahlung
	pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      }
    else if( particleName == "alpha" || particleName == "He3")
      {
	//multiple scattering
	pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
	
	//ionisation
	G4ionIonisation* ionIoni = new G4ionIonisation();
	ionIoni->SetStepFunction(0.1, 20*um);
	pmanager->AddProcess(ionIoni,                   -1, 2, 2);
	pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);

      }
    else if(particleName == "proton" || particleName == "deuteron" || particleName == "triton")

      {
	//multiple scattering
	pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
	
	//ionisation
	//G4ionIonisation* ionIoni = new G4ionIonisation();
	//ionIoni->SetStepFunction(0.1, 20*um);
	//pmanager->AddProcess(ionIoni,                   -1, 2, 2);

	G4hIonisation* hIoni = new G4hIonisation();
	hIoni->SetStepFunction(0.2, 50*um);
	pmanager->AddProcess(hIoni,                     -1, 2, 2);      
	pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3); 



      }
    else if (particleName == "GenericIon")
      {
	// OBJECT may be dynamically created as either a GenericIon or nucleus
	// G4Nucleus exists and therefore has particle type nucleus
	// genericIon:
	
	//multiple scattering
	pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

	//ionisation
	G4ionIonisation* ionIoni = new G4ionIonisation();
	ionIoni->SetEmModel(new G4IonParametrisedLossModel());
	ionIoni->SetStepFunction(0.1, 20*um);
	pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      } 

    else if ((!particle->IsShortLived()) &&
	     (charge != 0.0) && 
	     (particle->GetParticleName() != "chargedgeantino")) 
      {
	//all others charged particles except geantino
        G4hMultipleScattering* aMultipleScattering = new G4hMultipleScattering();
        G4hIonisation* ahadronIon = new G4hIonisation();
	
	//multiple scattering
	pmanager->AddProcess(aMultipleScattering,-1,1,1);

	//ionisation
	pmanager->AddProcess(ahadronIon,       -1,2,2);      
      }
    
  }

  // turn off msc step-limitation - especially as electron cut 1nm
  //opt.SetMscStepLimitation(fMinimal);

  // switch on fluorescence, PIXE and Auger:
  //opt.SetFluo(true);
  //opt.SetPIXE(true);
  //opt.SetAuger(true);

}


// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
//#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

void B4aPhysicsList::ConstructOp() 
{

/*
  // default scintillation process
  G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
  // theScintProcessDef->DumpPhysicsTable();
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetScintillationYieldFactor(1.0); //
  theScintProcessDef->SetScintillationExcitationRatio(0.0); //
  theScintProcessDef->SetVerboseLevel(OpVerbLevel);

  // scintillation process for alpha:
  G4Scintillation* theScintProcessAlpha = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessAlpha->SetTrackSecondariesFirst(true);
  theScintProcessAlpha->SetScintillationYieldFactor(1.1);
  theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
  theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

  // scintillation process for heavy nuclei
  G4Scintillation* theScintProcessNuc = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessNuc->SetTrackSecondariesFirst(true);
  theScintProcessNuc->SetScintillationYieldFactor(0.2);
  theScintProcessNuc->SetScintillationExcitationRatio(1.0);
  theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

  // optical processes
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  //  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  //  theAbsorptionProcess->DumpPhysicsTable();
  //  theRayleighScatteringProcess->DumpPhysicsTable();
  theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
  // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
  theBoundaryProcess->SetVerboseLevel(OpVerbLevel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (theScintProcessDef->IsApplicable(*particle)) {
	//      if(particle->GetPDGMass() > 5.0*GeV) 
	if(particle->GetParticleName() == "GenericIon") {
	  pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxPostStep);
	}	  
	else if(particle->GetParticleName() == "alpha") {
	  pmanager->AddProcess(theScintProcessAlpha);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxPostStep);
	}
	else {
	  pmanager->AddProcess(theScintProcessDef);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	}	  
      }
      
      if (particleName == "opticalphoton") {
	pmanager->AddDiscreteProcess(theAbsorptionProcess);
	//	pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	pmanager->AddDiscreteProcess(theBoundaryProcess);
      }
    }
    */
}


// Hadronic processes ////////////////////////////////////////////////////////

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"

#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
//#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4EmExtraPhysics.hh"

#include "G4EmStandardPhysics_option4.hh"


// Inelastic processes:
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4He3InelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"


#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"

// High energy FTFP model and Bertini cascade
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
//#include "G4GGNuclNuclCrossSection.hh"

#include "G4HadronElastic.hh"
#include "G4HadronCaptureProcess.hh"


//#include "G4IonsShenCrossSection.hh"
//#include "G4TripathiCrossSection.hh"
//#include "G4TripathiLightCrossSection.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4NeutronInelasticXS.hh"



// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"


// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void B4aPhysicsList::ConstructGeneral() {

  //void B4aPhysicsList::ConstructHad() {
  //Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  elastic_he->SetMinEnergy( elastic_elimitPi );

  
  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin0 =   0.0*GeV;
  const G4double theBERTMin1 =  19.0*MeV;
  const G4double theBERTMax =    5.0*GeV;
  const G4double theHPMin =      0.0*GeV;
  const G4double theHPMax =     20.0*MeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  //G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());

  G4ComponentGGHadronNucleusXsc * ggHNXsec = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
  G4VCrossSectionDataSet * theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);

  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->AddDataSet(theGGHNEl);
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);
  
  auto theParticleIterator=GetParticleIterator();
  theParticleIterator->reset();

  while ((*theParticleIterator)()) 
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionPlusInelasticProcess* theInelasticProcess = 
	    new G4PionPlusInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( thePiData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
	  //Inelastic scattering
	  G4PionMinusInelasticProcess* theInelasticProcess = 
	    new G4PionMinusInelasticProcess("inelastic");
	  //theInelasticProcess->AddDataSet( thePiData );
	  //theInelasticProcess->RegisterMe( theFTFModel1 );
          //theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );	  
	  //Absorption
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
	}
      
      
      else if (particleName == "proton") 
	{
	//B4aNucleonElasticXS *myEP=new B4aNucleonElasticXS(particle);
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
				GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
          
          
          //pmanager->AddDiscreteProcess( theElasticProcess );
	  
	  // myElastic scattering
	 /* G4HadronElasticProcess* myElasticProcess = 
	    new G4HadronElasticProcess("elastic");
	  myElasticProcess->AddDataSet( new B4aNucleonElasticXS( G4Proton::Proton() ) );
	          
	  pmanager->AddDiscreteProcess(myElasticProcess );
	  */
	   
	  
 // Inelastic scattering
	  G4ProtonInelasticProcess* theInelasticProcess = 
	    new G4ProtonInelasticProcess("inelastic");
	  theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	 //theInelasticProcess->AddDataSet( new B4aNucleonElasticXS( G4Proton::Proton() ) );
	 
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
          //theInelasticProcess->RegisterMe( myEP );
          
	  pmanager->AddDiscreteProcess( theInelasticProcess );


	
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4NeutronHPElastic * theElasticNeutronHP = new G4NeutronHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4NeutronHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering		
	G4NeutronInelasticProcess* theInelasticProcess =
	  new G4NeutronInelasticProcess("inelastic");
	theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4NeutronHPInelastic * theNeutronInelasticHPModel = new G4NeutronHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4NeutronHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess =
	  new G4HadronCaptureProcess;
	G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
	theLENeutronCaptureModel->SetMinEnergy(theHPMin);
	theLENeutronCaptureModel->SetMaxEnergy(theHPMax);
	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
	theCaptureProcess->AddDataSet( new G4NeutronHPCaptureData);
	pmanager->AddDiscreteProcess(theCaptureProcess);

      }

      else if (particleName == "deuteron") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4DeuteronInelasticProcess* theInelasticProcess = 
	    new G4DeuteronInelasticProcess("inelastic");
	  //theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->AddDataSet(theGGHNInel);
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );


	  //*****************************
	  //test
	  //pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
	
	//ionisation
	//G4ionIonisation* ionIoni = new G4ionIonisation();
	//ionIoni->SetStepFunction(0.1, 20*um);
	//pmanager->AddProcess(ionIoni,                   -1, 2, 2);

	}
      
      else if (particleName == "triton") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4TritonInelasticProcess* theInelasticProcess = 
	    new G4TritonInelasticProcess("inelastic");
	  //theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->AddDataSet(theGGHNInel);
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      else if (particleName == "alpha") 
	{
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4AlphaInelasticProcess* theInelasticProcess = 
	    new G4AlphaInelasticProcess("inelastic");	 
          //theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->AddDataSet(theGGHNInel);
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
    else if (particleName == "He3") 
	    {
	  // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
	  //pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4He3InelasticProcess* theInelasticProcess = 
	    new G4He3InelasticProcess("inelastic");	 
          //theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->AddDataSet(theGGHNInel);
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

    }
  //}



  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      
      if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) 
	{ 
	  pmanager ->AddProcess(theDecayProcess);
	  // set ordering for PostStepDoIt and AtRestDoIt
	  pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	  pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
	}
    }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable = 
    G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();

  for (G4int i=0; i<theIonTable->Entries(); i++) 
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      G4String particleType = theIonTable->GetParticle(i)->GetParticleType();
      
      if (particleName == "GenericIon") 
	{
	  G4ProcessManager* pmanager = 
	    theIonTable->GetParticle(i)->GetProcessManager();
	  pmanager->SetVerboseLevel(VerboseLevel);
	  pmanager ->AddProcess(theRadioactiveDecay);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	} 
    }
}

// Cuts /////////////////////////////////////////////////////////////////////
void B4aPhysicsList::SetCuts() 
{
  /*
  if (verboseLevel >1)
    G4cout << "DMXPhysicsList::SetCuts:";
  
  if (verboseLevel>0){
    G4cout << "DMXPhysicsList::SetCuts:";
    G4cout << "CutLength : " 
	   << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
*/
  //special for low energy physics
  //G4double lowlimit=250*eV;  
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100.*GeV);




  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForDeuteron, "deuteron");
  SetCutValue(cutForNeutron, "neutron");
  SetCutValue(cutForProton, "proton");
  //if (verboseLevel>0) DumpCutValuesTable();


  //SetCutsWithDefault();

  // Production thresholds for detector regions
  G4Region* region;
  G4String regName;
  G4ProductionCuts* cuts;

  regName = "TgRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  //cuts->SetProductionCut(0.001*mm); 
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("e-"));
   cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("proton"));
  
  region->SetProductionCuts(cuts);

  regName = "WinRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  //cuts->SetProductionCut(0.001*mm); 
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("proton"));
  
  region->SetProductionCuts(cuts);

  regName = "PDRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("proton"));
  //cuts->SetProductionCut(0.001*mm,G4ProductionCuts::GetIndex("e-"));
  
  region->SetProductionCuts(cuts);

  regName = "WrRegion";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cuts = new G4ProductionCuts;
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("gamma"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("e-"));
  cuts->SetProductionCut(0.0001*mm,G4ProductionCuts::GetIndex("proton"));
  //cuts->SetProductionCut(0.001*mm,G4ProductionCuts::GetIndex("e-"));
  
  region->SetProductionCuts(cuts);



  //regName = "calorimeter";
  //region = G4RegionStore::GetInstance()->GetRegion(regName);
  //cuts = new G4ProductionCuts;
  //cuts->SetProductionCut(0.01*mm,G4ProductionCuts::GetIndex("gamma"));
  //cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e-"));
  //cuts->SetProductionCut(0.1*mm,G4ProductionCuts::GetIndex("e+"));
  //region->SetProductionCuts(cuts);

if (verboseLevel>0) DumpCutValuesTable();






}
