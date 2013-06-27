#include <string>

#include <TTree.h>

#include "CommonTools/include/Counters.hh"
#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Utils.hh"
#include "TopHiggsAnalysis/include/tHSelection.hh"
#include "CommonTools/include/EfficiencyEvaluator.hh"
#include "CommonTools/include/LeptonIdBits.h"
#include "CommonTools/include/PUWeight.h"

#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <TLorentzVector.h>

#include <TTree.h>

using namespace bits;

tHSelection::tHSelection(TTree *tree) 
  : TopHiggs(tree) {
  
  // choose the Higgs Mass
  std::string higgsConfigDir;
  std::string higgsConfigDirMass;
  std::ifstream setfile("config/higgs/higgsMass.txt");
  higgsConfigDir="config/";
  higgsConfigDirMass="config/tH125/";
  std::cout << "Reading configuration for Higgs mass = 125 GeV/c^2" << std::endl;
  
  // Selection efficiencies
  // The selection Configure method prepare the code for all the cuts defined in the txt file
  std::string fileCuts     = higgsConfigDirMass + "2e2nuCuts.txt";
  std::string fileSwitches = higgsConfigDir + "3l3n1b1qTHSwitches.txt";
  
  CutBasedHiggsSelection[eee].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER EEE"); 
  _selectionEEE = CutBasedHiggsSelection[eee].GetSelection();  

  CutBasedHiggsSelection[mmm].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER MMM"); 
  _selectionMMM = CutBasedHiggsSelection[mmm].GetSelection();

  CutBasedHiggsSelection[eem].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER EEM"); 
  _selectionEEM = CutBasedHiggsSelection[eem].GetSelection();

  CutBasedHiggsSelection[mme].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER MME"); 
  _selectionMME = CutBasedHiggsSelection[mme].GetSelection();

  // Z-peak cross-check
  CutBasedHiggsSelection[ee].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER EE"); 
  _selectionEE = CutBasedHiggsSelection[ee].GetSelection();

  CutBasedHiggsSelection[mm].Configure(fileCuts.c_str(),fileSwitches.c_str(),"FULL SELECTION EVENT COUNTER MM"); 
  _selectionMM = CutBasedHiggsSelection[mm].GetSelection();

  // Extra selection efficiencies - to be put here not to pass the full list of leptons to the preselection class
  _selectionEEE->addCut("etaElectronAcc");    
  _selectionEEE->addCut("ptElectronAcc");
  _selectionEEE->addCut("etaMuonAcc");
  _selectionEEE->addCut("ptMuonAcc");
  _selectionEEE->addCut("etUncorrJetAcc");  

  _selectionEEE->addSwitch("apply_kFactor");   
  _selectionEEE->addSwitch("isData");
  _selectionEEE->addSwitch("goodRunLS");
  _selectionEEE->addSwitch("applyBDTEleID");

  _selectionEEE->addStringParameter("electronIDType");
  _selectionEEE->addStringParameter("electronIDTypeLow");  
  _selectionEEE->addStringParameter("JESUncertainty");

  if(_selectionEEE->getSwitch("applyBDTEleID")) m_useBDTEleID = true;
  else m_useBDTEleID = false;

  // configuring the electron BDT
  fMVA = new ElectronIDMVA();
  fMVA->Initialize("BDTG method",
                   "elebdtweights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml",
                   "elebdtweights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml",
                   "elebdtweights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml",
                   "elebdtweights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml",
                   "elebdtweights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml",
                   "elebdtweights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",                
                   ElectronIDMVA::kWithIPInfo);

  std::cout << "estoy aqui" << std::endl;

  // configurating the new lepton BDT
  // Electrons
  // see src/TopHiggs.hh

  fMVAElectron = new IDForBsMVA(true);
  fMVAElectron->Initialize("leptonMVA/el_pteta_low_cb_BDTG.weights.xml",
			   "leptonMVA/el_pteta_low_fb_BDTG.weights.xml",
			   "leptonMVA/el_pteta_low_ec_BDTG.weights.xml",
			   "leptonMVA/el_pteta_high_cb_BDTG.weights.xml",  
			   "leptonMVA/el_pteta_high_fb_BDTG.weights.xml",
			   "leptonMVA/el_pteta_high_ec_BDTG.weights.xml",
			   true);

  // Muons
  // see src/TopHiggs.hh

  fMVAMuon = new IDForBsMVA(false);
  fMVAMuon->Initialize("leptonMVA/mu_pteta_low_b_BDTG.weights.xml",
		       "leptonMVA/mu_pteta_low_e_BDTG.weights.xml",
		       "leptonMVA/mu_pteta_high_b_BDTG.weights.xml",
		       "leptonMVA/mu_pteta_high_e_BDTG.weights.xml",  
		       // not used
		       "leptonMVA/mu_pteta_high_e_BDTG.weights.xml",
		       "leptonMVA/mu_pteta_high_e_BDTG.weights.xml",

		       false);

  // Reading GoodRUN LS
  std::cout << "[GoodRunLS]::goodRunLS is " << _selectionEEE->getSwitch("goodRunLS") << " isData is " 
	    <<  _selectionEEE->getSwitch("isData") << std::endl;

  // To read good run list!
  if (_selectionEEE->getSwitch("goodRunLS") && _selectionEEE->getSwitch("isData")) {
    std::string goodRunJsonFile = "config/json/goodCollisions2012.json";
    setJsonGoodRunList(goodRunJsonFile);
    fillRunLSMap();
  }

  // kinematics inicialize variables
  m_p3PFMET      = new TVector3(0.,0.,0.);
  m_metFromJets  = new TVector3(0.,0.,0.);
  m_pfMetJESUp   = new TVector3(0.,0.,0.);
  m_pfMetJESDown = new TVector3(0.,0.,0.);

  // Rochester correction for muons
  rmcor = new RochCor2012();

  for(int theChannel=0; theChannel<6; theChannel++) {
    m_p4Lepton1[theChannel] = new TLorentzVector(0.,0.,0.,0.);
    m_p4Lepton2[theChannel] = new TLorentzVector(0.,0.,0.,0.);
    m_p4Lepton3[theChannel] = new TLorentzVector(0.,0.,0.,0.);

    m_jetsSum[theChannel]       = new TLorentzVector(0.,0.,0.,0.);
    m_jetscbIDSum[theChannel]   = new TLorentzVector(0.,0.,0.,0.);
    m_uncorrJetsSum[theChannel] = new TLorentzVector(0.,0.,0.,0.);
  }

  // b-veto event variables
  m_maxDxyEvt = 0.0;
  m_maxDszEvt = 0.0;

  // histo to study jet/electron match
  H_deltaRuncorr = new TH1F("H_deltaRuncorr","uncorrected jets",100, 0.,2*TMath::Pi());
  H_deltaRcorr   = new TH1F("H_deltaRcorr",  "corrected jets",  100, 0.,2*TMath::Pi());

  // as defaults switch off jet ID studies 
  wantJetIdStuff = true;

  // counter
  nWWtoLLdecays   = 0;
  nWfromTtoLdecay = 0;
  nSSevent        = 0;
  nTotal          = 0;



}

tHSelection::~tHSelection(){

  for(int theChannel=0; theChannel<6; theChannel++) {  

    delete m_p4Lepton1[theChannel];
    delete m_p4Lepton2[theChannel];
    delete m_p4Lepton3[theChannel];

    delete m_jetsSum      [theChannel]; // MVA ID for jets
    delete m_jetscbIDSum  [theChannel]; // Cut Based ID for jets
    delete m_uncorrJetsSum[theChannel];

  }

  delete m_p3PFMET;  
  delete m_metFromJets;
  delete m_pfMetJESUp;
  delete m_pfMetJESDown;

  delete _selectionEEE;
  delete _selectionMMM;
  delete _selectionEEM;
  delete _selectionMME;

  delete _selectionEE;
  delete _selectionMM;

  myOutTree[eee]->save();
  myOutTree[mmm]->save();
  myOutTree[mme]->save();
  myOutTree[eem]->save();

  myOutTree[ee]->save();
  myOutTree[mm]->save();

  myGenLevelTree->save();
  
}

bool tHSelection::findMcTree(const char* processType) {

  _process = "UNDEFINED";
  _theGenEle = -1;
  _theGenPos = -1;
  
  // now we look for ee || mumu || emu
  // in the acceptance and with a loose pT threshold
  float etaEleAcc_  = 2.5;
  float ptEleAcc_   = 5.0; // GeV
  float etaMuonAcc_ = 2.4;
  float ptMuonAcc_  = 0.0; // GeV

  if(strcmp(processType,"WZ")==0 ){

    int nW = 0;
    int nZ = 0;

    //std::cout << " imc | st | id | mother id | motherid (motherid) " << std::endl;
    //std::cout << "-------------------------------------------------" << std::endl;
        
    for(int imc=0;imc<30;imc++) {

      if ( !statusMc[imc] == 3 ) continue; // I am only interested in hard-scattering products
      
      //std::cout << " " << imc << " | " << statusMc[imc] << " | " << idMc[imc] << " | " << idMc[mothMc[imc]] 
      //		<< " | " << idMc[mothMc[mothMc[imc]]]   << std::endl;
      
      if( idMc[imc] == 23 ) nZ++;
      if( fabs(idMc[imc]) == 24) nW++;

    
    }

    return ( nZ > 0 && nW > 0 );

  }

  // tH gen level
  if(strcmp(processType,"tH")==0) {

    bool tHevent = false;

    int index_Top              = 999; // top  quark
    int index_Higgs            = 999; // Higgs boson

    int index_WfromTop         = 999;
    int index_WplusfromH       = 999;
    int index_WminusfromH      = 999;

    int index_lminusfromWfromH = 999; // lepton(-) coming from W(-) from H
    int index_lplusfromWfromH  = 999; // lepton(+) coming from W(+) from H
    int index_lfromWfromT      = 999; // lepton(+) coming from W from T

    int index_nminusfromWfromH = 999; // neutrino(-) coming from W(-) from H
    int index_nplusfromWfromH  = 999; // neutrino(+) coming from W(+) from H
    int index_nfromWfromT      = 999; // neutrino(+) coming from W from T

    //std::cout << " imc | st | id | mother id | motherid (motherid) " << std::endl;
    //std::cout << "-------------------------------------------------" << std::endl;

    int nLeptonsfromWfromH = 0;
    int nLeptonsfromWfromT = 0;

    int nPlusLepton  = 0;
    int nMinusLepton = 0;

    int nWfromH            = 0;

    for(int imc=0;imc<30;imc++) {

      if ( !(statusMc[imc] == 3) ) continue; // I am only interested in hard-scattering products
      
      //std::cout << " " << imc << " | " << statusMc[imc] << " | " << idMc[imc] << " | " << idMc[mothMc[imc]] 
      //	<< " | " << idMc[mothMc[mothMc[imc]]]   << std::endl;

      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));

      if (imc == 7){
	_genForwardQuark_Pt  = pMc[imc]*fabs(sin(thetaMc[imc]));
	_genForwardQuark_Eta = etaMc[imc];
      }
      
      if( fabs(idMc[imc]) == 6  ) index_Top   = imc;
      if( fabs(idMc[imc]) == 25 ) index_Higgs = imc;

      if( fabs(idMc[imc]) == 24 && fabs(idMc[mothMc[imc]]) == 6  ) index_WfromTop = imc;
      
      if( idMc[imc] == -24 && fabs(idMc[mothMc[imc]]) == 25 ){
	index_WminusfromH = imc;
	nWfromH++;
      }
      if( idMc[imc] == +24 && fabs(idMc[mothMc[imc]]) == 25 ){
	index_WplusfromH  = imc;
	nWfromH++;
      }

      // Lepton(-) from W(-) decay
      if     ( (idMc[imc] == +11 || idMc[imc] == +13 || idMc[imc] == +15)  && idMc[mothMc[imc]] == -24 ){ 
	
	if      ( fabs(idMc[mothMc[mothMc[imc]]]) == 6  ){  // W from Top decay
	  index_lfromWfromT = imc;
	  nLeptonsfromWfromT++;
	  nMinusLepton++;
	}
	else if ( fabs(idMc[mothMc[mothMc[imc]]]) == 25 ){ // W from Higgs decay
	  index_lminusfromWfromH = imc; 
	  nLeptonsfromWfromH++;
	  nMinusLepton++;
	}

      }      

      // Neutrino(+) from W(-) decay
      else if( (idMc[imc] == -12 || idMc[imc] == -14 || idMc[imc] == -16)  && idMc[mothMc[imc]] == -24 ){ 

	if      ( fabs(idMc[mothMc[mothMc[imc]]]) == 6  ){ // W from Top decay
	  index_nfromWfromT = imc;
	  nLeptonsfromWfromT++;
	}
	else if ( fabs(idMc[mothMc[mothMc[imc]]]) == 25 ){ // W from Higgs decay
	  index_nplusfromWfromH = imc; 
	  nLeptonsfromWfromH++;
	}
	
      }

      // Lepton(+) from W(+) decay
      else if( (idMc[imc] == -11 || idMc[imc] == -13 || idMc[imc] == -15)  && idMc[mothMc[imc]] == +24 ){ 
	
	if      ( fabs(idMc[mothMc[mothMc[imc]]]) == 6  ){ // W from Top decay
	  index_lfromWfromT = imc;
	  nLeptonsfromWfromT++;
	  nPlusLepton++;
	}
	else if ( fabs(idMc[mothMc[mothMc[imc]]]) == 25 ){ // W from Higgs decay
	  index_lplusfromWfromH = imc; 
	  nLeptonsfromWfromH++;
	  nPlusLepton++;
	}

      }      

      // Neutrino(-) from W(+) decay
      else if( (idMc[imc] == 12 || idMc[imc] == +14 || idMc[imc] == +16)  && idMc[mothMc[imc]] == +24 ){ 

	if      ( fabs(idMc[mothMc[mothMc[imc]]]) == 6  ){ // W from Top decay
	  index_nfromWfromT = imc;
	  nLeptonsfromWfromT++;
	}
	else if ( fabs(idMc[mothMc[mothMc[imc]]]) == 25 ){ // W from Higgs decay
	  index_nminusfromWfromH = imc; 
	  nLeptonsfromWfromH++;
	}

      }

    }

    if ( nLeptonsfromWfromH == 4){
      nWWtoLLdecays++;
    }
    if ( nLeptonsfromWfromT == 2){ // (lepton + neutrino)
      nWfromTtoLdecay++;
    }

    if ( (nPlusLepton == 2 && nMinusLepton == 0) || ( nMinusLepton == 2 && nPlusLepton == 0 ) ){
      nSSevent++;
    }

    _genHiggsPt  = pMc[index_Higgs]*fabs(sin(thetaMc[index_Higgs]));
    _genHiggsEta = etaMc[index_Higgs]; 

    _genTopPt  = pMc[index_Top]*fabs(sin(thetaMc[index_Top]));
    _genTopEta = etaMc[index_Top]; 

    _genWpfromH_Pt  = pMc[index_WplusfromH]*fabs(sin(thetaMc[index_WplusfromH]));
    _genWpfromH_Eta = etaMc[index_WplusfromH];
    
    _genWmfromH_Pt  = pMc[index_WminusfromH]*fabs(sin(thetaMc[index_WminusfromH])) ;
    _genWmfromH_Eta = etaMc[index_WminusfromH];
    
    _genWfromT_Pt  = pMc[index_WfromTop]*fabs(sin(thetaMc[index_WfromTop]));
    _genWfromT_Eta = etaMc[index_WfromTop];

    if ( nLeptonsfromWfromH > 1 ){

      _genNeutrinoPlusfromWfromH_Pt  = pMc[index_nplusfromWfromH]*fabs(sin(thetaMc[index_nplusfromWfromH]));
      _genNeutrinoPlusfromWfromH_Eta = etaMc[index_nplusfromWfromH];
      
      _genNeutrinoMinusfromWfromH_Pt  = pMc[index_nminusfromWfromH]*fabs(sin(thetaMc[index_nminusfromWfromH]));
      _genNeutrinoMinusfromWfromH_Eta = etaMc[index_nminusfromWfromH];
      
      _genLeptonPlusfromWfromH_Pt  = pMc[index_lplusfromWfromH]*fabs(sin(thetaMc[index_lplusfromWfromH]));
      _genLeptonPlusfromWfromH_Eta = etaMc[index_lplusfromWfromH];
      
      _genLeptonMinusfromWfromH_Pt  = pMc[index_lminusfromWfromH]*fabs(sin(thetaMc[index_lminusfromWfromH]));
      _genLeptonMinusfromWfromH_Eta = etaMc[index_lminusfromWfromH];

    } 
     
    if(nLeptonsfromWfromT > 0){

      _genLeptonfromWfromT_Pt  = pMc[index_lfromWfromT]*fabs(sin(thetaMc[index_lfromWfromT]));
      _genLeptonfromWfromT_Eta = etaMc[index_lfromWfromT];

      _genNeutrinofromWfromT_Pt  = pMc[index_nfromWfromT]*fabs(sin(thetaMc[index_nfromWfromT]));
      _genNeutrinofromWfromT_Eta = etaMc[index_nfromWfromT]; 

    }
    
    if (index_Top < 30 && index_Higgs < 30 ) tHevent = true;

    return ( tHevent );

  }
  
  // signal: 2e2nu
  if(strcmp(processType,"HtoWWto2e2nu")==0) {
    int indeminus=999, indeplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc] == -11 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if( idMc[imc] ==  11 && fabs(etaMc[imc]) < etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    if( indeminus<25 && indeplus<25 ) {
      _theGenPos = indeplus;
      _theGenEle = indeminus;
    }
    return ( indeplus < 25 && indeminus < 25 );
  }

  // signal: 2m2nu
  if(strcmp(processType,"HtoWWto2m2nu")==0) {
    int indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc] == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc] ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
    }
    if( indmuminus<25 && indmuplus<25 ) {
      _theGenMuPlus  = indmuplus;
      _theGenMuMinus = indmuminus;
    }
    return ( indmuplus < 25 && indmuminus < 25 );
  }

  // signal: em2nu
  if(strcmp(processType,"HtoWWtoem2nu")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;

    bool isEM = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeplus = imc;
      if( idMc[imc]  ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  ==  11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos>ptMcMuMinus) isEM = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle>ptMcMuPlus) isEM = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isEM) || (indeminus<25 && indmuplus<25 && isEM) );
  }

  // signal: me2nu
  if(strcmp(processType,"HtoWWtome2nu")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;
    
    bool isME = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if( idMc[imc]  == -11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeplus = imc;
      if( idMc[imc]  ==  13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = imc;
      if( idMc[imc]  == -13 && fabs(etaMc[imc]) < etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if( idMc[imc]  ==  11 && fabs(etaMc[imc]) < etaEleAcc_  && ptMc > ptEleAcc_ )  indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos<ptMcMuMinus) isME = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle<ptMcMuPlus) isME = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isME) || (indeminus<25 && indmuplus<25 && isME) );
  }

  // signal ee excluding taus
  if(strcmp(processType,"HtoWWto2e2nu_prompt")==0) {
    int indeminus=999, indeplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if(idMc[imc] == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if(idMc[imc] == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    if( indeminus<25 && indeplus<25 ) {
      _theGenPos = indeplus;
      _theGenEle = indeminus;
    }
    return ( indeplus < 25 && indeminus < 25 );
  }

  // signal mm excluding taus
  if(strcmp(processType,"HtoWWto2m2nu_prompt")==0) {
    int indmuminus=999, indmuplus=999;
    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if(idMc[imc] == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if(idMc[imc] == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuminus = 11;
    }
    if( indmuminus<25 && indmuplus<25 ) {
      _theGenMuPlus = indmuplus;
      _theGenMuMinus = indmuminus;
    }
    return ( indmuplus < 25 && indmuminus < 25 );
  }

  // signal em excluding taus
  if(strcmp(processType,"HtoWWtoem2nu_prompt")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;

    bool isEM = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if(idMc[imc] == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if(idMc[imc] == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ )indmuminus = imc;
      if(idMc[imc] == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if(idMc[imc] == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }

    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos>ptMcMuMinus) isEM = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle>ptMcMuPlus) isEM = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isEM) || (indeminus<25 && indmuplus<25 && isEM) );
  }

  // signal me excluding taus
  if(strcmp(processType,"HtoWWtome2nu_prompt")==0) {
    int indeminus=999, indeplus=999, indmuminus=999, indmuplus=999;

    bool isME = false;

    for(int imc=6;imc<25;imc++) {
      float ptMc = pMc[imc]*fabs(sin(thetaMc[imc]));
      if(idMc[imc] == -11 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeplus = imc;
      if(idMc[imc] == 13 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ )indmuminus = imc;
      if(idMc[imc] == -13 && idMc[mothMc[imc]]==24 && fabs(etaMc[imc])<etaMuonAcc_ && ptMc > ptMuonAcc_ ) indmuplus = imc;
      if(idMc[imc] == 11 && idMc[mothMc[imc]]==-24 && fabs(etaMc[imc])<etaEleAcc_ && ptMc > ptEleAcc_ ) indeminus = imc;
    }
    
    if( indeplus<25 && indmuminus<25 ) {
      _theGenPos     = indeplus;
      _theGenMuMinus = indmuminus;
      float ptMcPos     = pMc[indeplus]*fabs(sin(thetaMc[indeplus]));      
      float ptMcMuMinus = pMc[indmuminus]*fabs(sin(thetaMc[indmuminus]));      
      if ( ptMcPos<ptMcMuMinus) isME = true;
    } else if( indeminus<25 && indmuplus<25 ) {
      _theGenEle = indeminus;
      _theGenMuPlus = indmuplus;
      float ptMcEle    = pMc[indeminus]*fabs(sin(thetaMc[indeminus]));      
      float ptMcMuPlus = pMc[indmuplus]*fabs(sin(thetaMc[indmuplus]));      
      if ( ptMcEle<ptMcMuPlus) isME = true;
    }
    
    return ( (indeplus<25 && indmuminus<25 && isME) || (indeminus<25 && indmuplus<25 && isME) );
  }
  
  // signal: 2l2nu
  if(strcmp(processType,"HtoWWto2l2nu")==0) {
    int indlminus=999, indlplus=999;
    for(int imc=6;imc<25;imc++) {
      if(idMc[imc]>10 && idMc[imc]<19 && idMc[mothMc[imc]]==-24)  indlminus=imc;
      if(idMc[imc]<-10 && idMc[imc]>-19 && idMc[mothMc[imc]]==24) indlplus=imc;
    }
    if(indlminus<25 && indlplus<25) {
      if( (idMc[indlminus]==-11) || (idMc[indlminus]==-13) || (idMc[indlminus]==-15))  
	_theGenEle = indlminus;
      if( (idMc[indlplus]==11) || (idMc[indlplus]==13) || (idMc[indlplus]==15) )
	_theGenPos = indlplus;
    }
    return (indlminus<25 && indlplus<25);
  }
  

  // WW: e / mu / tau
  else if(strcmp(processType,"WW")==0) {
    _process = "WW";
    TVector3 WminusP, WplusP;
    WminusP.SetMagThetaPhi(pMc[6],thetaMc[6],phiMc[6]);
    WplusP.SetMagThetaPhi(pMc[7],thetaMc[7],phiMc[7]);
    float pT = (WminusP+WplusP).Pt();
    _theGenEle = 6;
    _theGenPos = 7;
    _genHiggsPt = pT;
    return (
	    (abs(idMc[6])==24) && (abs(idMc[7])==24) &&
	    (abs(idMc[8])>10 && abs(idMc[8])<19 && abs(idMc[mothMc[8]])==24) &&
	    (abs(idMc[10])>10 && abs(idMc[10])<19 && abs(idMc[mothMc[10]])==24)
	    );
  }

  // w+jets: e / mu / tau
  else if(strcmp(processType,"Wjets")==0) {
    _process = "Wjets";
    return ( ((abs(idMc[8])==11) && abs(idMc[9])==12) || ((abs(idMc[8])==13) && abs(idMc[9])==14) || ((abs(idMc[8])==15) && abs(idMc[9])==16));
  }
  // ttbar: e / mu / tau
  else if(strcmp(processType,"ttbar")==0) {
    _process = "ttbar";
    _theGenEle = 13;
    _theGenPos = 15;
    return ( 
	    abs(idMc[9])==24 && abs(idMc[15])>10 && abs(idMc[15])<19 &&
	    abs(idMc[11])==24 && abs(idMc[13])>10 && abs(idMc[13])<19 &&
	    (idMc[13]*idMc[15]<0)
	    );
  }
  else if(strcmp(processType,"ZZleptonic")==0) {
    _process = "ZZleptonic";
    // 8,9; 10,11 are the daughters of the Z;Z
    return (fabs(idMc[8])>10 && fabs(idMc[8])<19 &&
	    fabs(idMc[9])>10 && fabs(idMc[9])<19 &&
	    fabs(idMc[10])>10 && fabs(idMc[10])<19 &&
	    fabs(idMc[11])>10 && fabs(idMc[11])<19);
  }
  else {
    std::cout << "This processType: " << processType << " is not expected, you should put MTtruth switch off" <<std::endl;
    return false;
  }
}

void tHSelection::Loop() {

  _verbose=false;
  if(fChain == 0) return;

  // kinematics reduced tree
  std::string reducedTreeNameEEE = _datasetName+"-datasetEEE.root";
  std::string reducedTreeNameMMM = _datasetName+"-datasetMMM.root";
  std::string reducedTreeNameEEM = _datasetName+"-datasetEEM.root";
  std::string reducedTreeNameMME = _datasetName+"-datasetMME.root";

  std::string reducedTreeNameEE = _datasetName+"-datasetEE.root";
  std::string reducedTreeNameMM = _datasetName+"-datasetMM.root";
 
  myOutTree[eee] = new RedTopHiggsTree(reducedTreeNameEEE.c_str());
  myOutTree[mmm] = new RedTopHiggsTree(reducedTreeNameMMM.c_str());
  myOutTree[eem] = new RedTopHiggsTree(reducedTreeNameEEM.c_str());
  myOutTree[mme] = new RedTopHiggsTree(reducedTreeNameMME.c_str());

  myOutTree [ee] = new RedTopHiggsTree(reducedTreeNameEE.c_str());
  myOutTree [mm] = new RedTopHiggsTree(reducedTreeNameMM.c_str());

  myGenLevelTree = new RedTopHiggsTree("genLevelInfo.root");

  if(!_selectionEEE->getSwitch("isData")) 
    myGenLevelTree->addtHMcTruthInfos();
  
  for(int theChannel=0; theChannel<6; theChannel++) {

    if(!_selectionEEE->getSwitch("isData")) {

      myOutTree[theChannel]->addMcTruthInfos();
      myOutTree[theChannel]->addPDFs();
      myOutTree[theChannel]->addtHMcTruthInfos();

    }
  
    myOutTree[theChannel]->addRunInfos();

    myOutTree[theChannel]->addMLVars();
    
    if (wantJetIdStuff) myOutTree[theChannel]->addJetsVars();

    myOutTree[theChannel]->addLatinos();
  
    myOutTree[theChannel]->addKinematics();

    myOutTree[theChannel]->addLeptonInfos();
    
    myOutTree[theChannel]->addSystematics();

    myOutTree[theChannel]->addMetStudies();

  }
  
  unsigned int lastLumi = 0;
  unsigned int lastRun  = 0;

  //PUWeight* fPUWeight = new PUWeight();

  Long64_t nbytes = 0, nb = 0;
  //Long64_t nentries = fChain->GetEntries();
  Long64_t nentries = 20000;
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%1000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    nTotal++;
    
    // for each event, reset the estimation of the variables
    resetKinematicsStart();

    // get the kFactor of the event (for signal)
    float weight = 1;
    
    float nputosave[3];
    if( !_selectionEEE->getSwitch("isData") ) {
      for(int i=0; i<3; i++) nputosave[i] = float(nPU[i]);
    } else {
      for(int i=0; i<3; i++) nputosave[i] = -1.;
    }

    // look to the MC truth decay tree 
    // bool decayEE = findMcTree("HtoWWto2e2nu");
    // bool decayMM = findMcTree("HtoWWto2m2nu");
    // bool decayEM = findMcTree("HtoWWtoem2nu");

    // MC Truth information
    bool promptEE, promptMM, promptEM, promptME;
    promptEE = promptMM = promptEM = promptME = false;

    _genmll = _genptll = _genyll = -1.;

    bool WZgenerated = false;
    bool tHgenerated = false;

    if( !_selectionEEE->getSwitch("isData") ) {
      tHgenerated = findMcTree("tH");  
      //WZgenerated = findMcTree("WZ");  
    }

    if (tHgenerated){
      myGenLevelTree -> filltHMcTruthInfos(_genHiggsPt,
					   _genHiggsEta,
					   
					   _genTopPt,
					   _genTopEta,
					   
					   _genWpfromH_Pt,
					   _genWpfromH_Eta,
					   
					   _genWmfromH_Pt,
					   _genWmfromH_Eta,
					   
					   _genWfromT_Pt,
					   _genWfromT_Eta,
					   
					   _genLeptonPlusfromWfromH_Pt,
					   _genLeptonPlusfromWfromH_Eta,
					   
					   _genLeptonMinusfromWfromH_Pt,
					   _genLeptonMinusfromWfromH_Eta,
					   
					   _genLeptonfromWfromT_Pt,
					   _genLeptonfromWfromT_Eta,
					   
					   _genNeutrinoPlusfromWfromH_Pt,
					   _genNeutrinoPlusfromWfromH_Eta,
					   
					   _genNeutrinoMinusfromWfromH_Pt,
					   _genNeutrinoMinusfromWfromH_Eta,
					   
					   _genNeutrinofromWfromT_Pt,
					   _genNeutrinofromWfromT_Eta,
					   
					   _genForwardQuark_Pt,
					   _genForwardQuark_Eta );
      myGenLevelTree->store();

    }    
    

    // Good Run selection
    if ( _selectionEEE->getSwitch("isData") && _selectionEEE->getSwitch("goodRunLS") && !isGoodRunLS() ) {
      if ( lastRun!= runNumber || lastLumi != lumiBlock) {
        lastRun = runNumber;
        lastLumi = lumiBlock;
        std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is rejected" << std::endl;
      }
      continue;
    }

    if ( _selectionEEE->getSwitch("isData") && _selectionEEE->getSwitch("goodRunLS") && 
	( lastRun!= runNumber || lastLumi != lumiBlock) ) {
      lastRun = runNumber;
      lastLumi = lumiBlock;
      std::cout << "[GoodRunLS]::Run " << lastRun << " LS " << lastLumi << " is OK" << std::endl;
    }
    
    // IMPORTANT: FOR DATA RELOAD THE TRIGGER MASK PER FILE WHICH IS SUPPOSED TO CONTAIN UNIFORM CONDITIONS X FILE
    reloadTriggerMask(runNumber);
    //bool passedHLT = (_selectionEEE->getSwitch("isData")) ? hasPassedHLT() : true;
    bool passedHLT[4];

    if (_selectionEEE->getSwitch("isData")){

      passedHLT[eee] = hasPassedHLT(eee);  
      passedHLT[mmm] = hasPassedHLT(mmm); 
      passedHLT[eem] = hasPassedHLT(eem);  
      passedHLT[mme] = hasPassedHLT(mme); 

      passedHLT[ee] = hasPassedHLT(ee);  
      passedHLT[mm] = hasPassedHLT(mm); 
    
    }else{

      passedHLT[eee] = true; 
      passedHLT[mmm] = true; 
      passedHLT[eem] = true;  
      passedHLT[mme] = true; 

      passedHLT[ee] = true;  
      passedHLT[mm] = true; 
      
    }

    // -------------------------------------------------------------
    // Vertex selection - we only consider the first vertex of the list ( = highest sumPT^2 )
    bool isGoodVertex = goodPV(0);
    
    // Count the number of good PVs (necessary for the MET cut)
    int nGoodPV=0;
    for(int v=0;v<nPV;++v) {
      if(goodPV(v)) nGoodPV++;
    }
    m_goodvertices = nGoodPV;

    // -------------------------------------------------------------
    
    // Get the best electrons and best muons for each channel

    /*
    std::cout<< " ************* " <<std::endl;
    std::cout<< " *** EVENT *** " <<std::endl;
    std::cout<< " ************* " <<std::endl;
    */


    // EEE Channel
    std::vector<int> theBestEEE = getBestEEE();

    // MMM Channel
    std::vector<int> theBestMMM = getBestMMM();

    // EEM Channel
    std::vector<int> theBestEEM = getBestEEM();

    // MME Channel
    std::vector<int> theBestMME = getBestMME();

    // EE Channel
    std::vector<int> theBestEE = getBestEE();
    std::cout<< " [Z-peak]: Checking EE " << std::endl;
    std::cout<< " (EE) Leptons index [0][1] = [" 
    	     << theBestEE.at(0) << "][" << theBestEE.at(1) << "]" << std::endl;

    // MM Channel
    std::vector<int> theBestMM = getBestMM();
    std::cout<< " [Z-peak]: Checking MM " << std::endl;
    std::cout<< " (MM) Leptons index [0][1] = [" 
    	     << theBestMM.at(0) << "][" << theBestMM.at(1) << "]" << std::endl;

    // Reconstructed channel
    m_channel[eee] = false;     
    m_channel[mmm] = false;
    m_channel[eem] = false;
    m_channel[mme] = false;

    m_channel[ee] = false;
    m_channel[mm] = false;

    if(isGoodVertex){
      
      // Give priority to the channels with muons
      if ( theBestMMM.at(0)> -1 && theBestMMM.at(1) > -1 && theBestMMM.at(2) > -1 ) {
	m_channel[mmm] = true;

	/*
	std::cout<< " event selected " << std::endl;

	int myMu1 = theBestMMM.at(0);
	int myMu2 = theBestMMM.at(1);
	int myMu3 = theBestMMM.at(2);

	float Muon1Pt = GetPt(pxMuon[myMu1],pyMuon[myMu1]);
	float Muon2Pt = GetPt(pxMuon[myMu2],pyMuon[myMu2]);
	float Muon3Pt = GetPt(pxMuon[myMu3],pyMuon[myMu3]);
	
	TLorentzVector Muon1, Muon2, Muon3;
	Muon1.SetPxPyPzE(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
	Muon2.SetPxPyPzE(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
	Muon3.SetPxPyPzE(pxMuon[myMu3], pyMuon[myMu3], pzMuon[myMu3], energyMuon[myMu3]);
	
	float ch1 = chargeMuon[theBestMMM.at(0)];
	float ch2 = chargeMuon[theBestMMM.at(1)];
	float ch3 = chargeMuon[theBestMMM.at(2)];

	std::cout << " ** Rochester Correction ** " << std::endl;
	std::cout << " [Before] pT[muon1][muon2][muon3] = [" 
		  << Muon1.Pt() << "][" << Muon2.Pt() << "][" << Muon3.Pt() << "]" << std::endl;
	
	rmcor->momcor_mc(Muon1, ch1, 0, 1.0);
	rmcor->momcor_mc(Muon2, ch2, 0, 1.0);
	rmcor->momcor_mc(Muon3, ch3, 0, 1.0);

	std::cout << " " << std::endl;
	std::cout << " [After] pT[muon1][muon2][muon3] = [" 
		  << Muon1.Pt() << "][" << Muon2.Pt() << "][" << Muon3.Pt() << "]" << std::endl;
	std::cout << " " << std::endl;
	*/

      }

      else if ( theBestMME.at(0)> -1 && theBestMME.at(1) > -1 && theBestMME.at(2) > -1 ) {
	m_channel[mme] = true;
      }

      else if ( theBestEEM.at(0)> -1 && theBestEEM.at(1) > -1 && theBestEEM.at(2) > -1 ) {
	m_channel[eem] = true;
      }

      else if ( theBestEEE.at(0)> -1 && theBestEEE.at(1) > -1 && theBestEEE.at(2) > -1 ) {
	m_channel[eee] = true;
      }

      // For Z-peak cross-check, select first MM
      if      ( theBestMM.at(0)> -1 && theBestMM.at(1) > -1 ) {
	m_channel[mm] = true;
      }else if( theBestEE.at(0)> -1 && theBestEE.at(1) > -1 ) {
	m_channel[ee] = true;
      }

      
    }

    /*
    std::cout << " This event has been selected as: " << std::endl;
    std::cout << " EEE = [" << m_channel[eee] << "]" << std::endl;
    std::cout << " MMM = [" << m_channel[mmm] << "]" << std::endl;
    std::cout << " EEM = [" << m_channel[eem] << "]" << std::endl;
    std::cout << " MME = [" << m_channel[mme] << "]" << std::endl;
    */

    //if( !m_channel[eee] && !m_channel[eem] && !m_channel[mme] && !m_channel[mmm] ) continue;
  
    //std::cout<< " we continue ..." << std::endl;

    // -------------------------------------------------------------
    // Set of kinematics: : now I've all the final leptons 
    resetKinematics();
    
    // MET is an event variable. Independent on the channel
    TVector3 pfmet(pxPFMet[0],pyPFMet[0],pzPFMet[0]);

    int sample = kMC2012;
    if(_selectionEEE->getSwitch("isData")) sample = kDATA2012;

    TVector3 corrpfmet = XYCorrectedMet(pfmet,sample);
    m_p3PFMET->SetXYZ(corrpfmet.x(),corrpfmet.y(),corrpfmet.z());
    m_theMET = m_p3PFMET->Pt();

    setKinematicsEEE(theBestEEE.at(0), theBestEEE.at(1), theBestEEE.at(2));
    setKinematicsMMM(theBestMMM.at(0), theBestMMM.at(1), theBestMMM.at(2));
    setKinematicsEEM(theBestEEM.at(0), theBestEEM.at(1), theBestEEM.at(2));
    setKinematicsMME(theBestMME.at(0), theBestMME.at(1), theBestMME.at(2));

    setKinematicsEE(theBestEE.at(0), theBestEE.at(1));
    setKinematicsMM(theBestMM.at(0), theBestMM.at(1));

    float dphiLLJ[6];
    float    btag[6];

    int       nsoftmu[6];
    int nsoftmunojets[6];
    int nextraleptons[6];

    float jesMtUp  [6];
    float jesMtDown[6];
    

    for(int ichan=0; ichan<6; ichan++) {

      // Initialize the btags for the leading and subleading jets to unphysical value
      leadJetBtag       [ichan] = -2000.;
      subleadJetBtag    [ichan] = -2000.;
      subLeadJetsMaxBtag[ichan] = -2000.;
      
      leadJetBtagBProb       [ichan] = -2000.;
      subleadJetBtagBProb    [ichan] = -2000.;
      subLeadJetsMaxBtagBProb[ichan] = -2000.;
      
      // Initialize the number of soft jets
      nsoftjets[ichan]     = 0; // MVA Id
      nsoftjetscbID[ichan] = 0; // Cut Based Id
      
      // Initialize the number of b-tagged of MVA Id jets
      m_numbtagCSVLmvaIDcentraljets[ichan] = 0; // loose
      m_numbtagCSVMmvaIDcentraljets[ichan] = 0; // medium
      m_numbtagCSVTmvaIDcentraljets[ichan] = 0; // tight

      // Initialize the number of b-tagged of Cut Based Id jets
      m_numbtagCSVLcbIDcentraljets[ichan] = 0; // loose
      m_numbtagCSVMcbIDcentraljets[ichan] = 0; // medium
      m_numbtagCSVTcbIDcentraljets[ichan] = 0; // tight

      m_nummvaIDforwardjets[ichan] = 0;
      m_numcbIDforwardjets [ichan] = 0; // loose

      // Initialize variables for jetId studies
      if (wantJetIdStuff) {
	
	leadJetPt[ichan]           = -2000.;
	leadJetEta[ichan]          = -2000.;
	leadJetLoosePFId[ichan]    = -2000.;
	leadJetMatchGen[ichan]     = -2000.;
	leadJetMvaJetId[ichan]     = -2000.;
	leadJetLooseId[ichan]      = -2000.;
	//
	subleadJetPt[ichan]        = -2000.;
	subleadJetEta[ichan]       = -2000.;
	subleadJetLoosePFId[ichan] = -2000.;
	subleadJetMatchGen[ichan]  = -2000.;
	subleadJetMvaJetId[ichan]  = -2000.;
	subleadJetLooseId[ichan]   = -2000.;
	
      }    
      
      // Initialize jet counters
      njets       [ichan] = numJets(eleCands[ichan],muCands[ichan],ichan);      // MVA Id jets
      ncbIDjets   [ichan] = numcbIDJets(eleCands[ichan],muCands[ichan],ichan);  // Cut Based Id jets
      nuncorrjets [ichan] = numUncorrJets(eleCands[ichan],muCands[ichan],ichan);// MVA Id jets 
      
      // if 1-jet bin, use deltaphi(ll-jet)
      // old HWW 
      dphiLLJ[ichan] = deltaPhiLLJet(ichan);   

      // b veto
      btag[ichan] = bVetoJets(eleCands[ichan],muCands[ichan],ichan);
      
      // Soft muon counter
      
      // w/o jet cleaning (for the selection)
      std::vector<int> emptyJets;
      emptyJets.clear();
      nsoftmu[ichan] = numSoftMuons(muCands[ichan],emptyJets);

      // w jet cleaning (for the top estimation)
      nsoftmunojets[ichan] = numSoftMuons(muCands[ichan],m_goodJets);

      // Extra lepton counter
      nextraleptons[ichan] = numExtraLeptons(eleCands[ichan],muCands[ichan]);

      // Calculate the smeared MET/MT // activate only when doing systematics, otherwise very slow
      //       JESPfMet(eleCands[ichan],muCands[ichan]);
      //       jesMtUp[ichan] = (transvMassJES(ichan)).first;
      //       jesMtDown[ichan] = (transvMassJES(ichan)).second;
    }

    float genPtHiggs = -1.;
    if ( !_selectionEEE->getSwitch("isData") ) {
      for(int imc=2;imc<10;imc++) {
        if(idMc[imc]==25 && statusMc[imc]==3) genPtHiggs = pMc[imc]*fabs(sin(thetaMc[imc]));
      }}

    // ---------------------------------------
    // filling counters for the different final states

    for(int theChannel=0; theChannel<6; theChannel++) {

      //std::cout<<"estoy en el channel theChannel == "<<theChannel<<std::endl;
      
      CutBasedHiggsSelection[theChannel].SetWeight       (weight);               
      //CutBasedHiggsSelection[theChannel].SetMcTruth    (promptEE); 
      CutBasedHiggsSelection[theChannel].SetHLT          (passedHLT[theChannel]);               
      CutBasedHiggsSelection[theChannel].SetIsChannel    (m_channel[theChannel]);     
      CutBasedHiggsSelection[theChannel].SetHighElePt    (hardestLeptonPt[theChannel]); 
      CutBasedHiggsSelection[theChannel].SetLowElePt     (slowestLeptonPt[theChannel]);  
      CutBasedHiggsSelection[theChannel].SetNJets        (njets[theChannel]);
      CutBasedHiggsSelection[theChannel].SetNcbIDJets    (ncbIDjets[theChannel]);
      CutBasedHiggsSelection[theChannel].SetNUncorrJets  (nuncorrjets[theChannel]);
      CutBasedHiggsSelection[theChannel].SetBTagJets     (btag[theChannel]);
      CutBasedHiggsSelection[theChannel].SetNSoftMuons   (nsoftmu[theChannel]);
      CutBasedHiggsSelection[theChannel].SetNExtraLeptons(nextraleptons[theChannel]);
      CutBasedHiggsSelection[theChannel].SetMet          (m_theMET);
      CutBasedHiggsSelection[theChannel].SetProjectedMet (m_projectedMet[theChannel]);
      CutBasedHiggsSelection[theChannel].SetNvtx         (nGoodPV);
      CutBasedHiggsSelection[theChannel].SetMetOverPtLL  (m_metOptll[theChannel]);
      CutBasedHiggsSelection[theChannel].SetDeltaPhiLLJet(dphiLLJ[theChannel]);   
      CutBasedHiggsSelection[theChannel].SetDeltaPhi     (m_deltaPhi[theChannel]);
      CutBasedHiggsSelection[theChannel].SetInvMass      (m_mll[theChannel]);
      CutBasedHiggsSelection[theChannel].SetDetaLeptons  (m_deltaEtaLeptons[theChannel]);
      CutBasedHiggsSelection[theChannel].SetPtll         (m_dilepPt[theChannel].Pt());
      CutBasedHiggsSelection[theChannel].SetWWInvMass    (m_transvMass[theChannel]);
      
      //if (m_channel[theChannel]) std::cout<<"m_chanel[theChannel] == true" <<std::endl;
      
      bool isSelected = true;
      
      // Call internally to some functions to obtain the bools for the steps
      bool dummybool = CutBasedHiggsSelection[theChannel].output();

      // latinos
      bool outputStep0  = CutBasedHiggsSelection[theChannel].outputStep0();
      bool outputStep1  = CutBasedHiggsSelection[theChannel].outputStep1();
      
      //if (outputStep0)   std::cout<< " outputStep0   = true"<<std::endl;
      //if (outputStep1)   std::cout<< " outputStep1   = true"<<std::endl;
      
      
      // Filling the tree

      /*
      // Not for the moment
      if(!_selectionEEE->getSwitch("isData")) {


	myOutTree[theChannel] -> fillMcTruth(promptEE,
	                                     _genmll ,
					     _genptll,
					     _genyll);
	myOutTree[theChannel] -> fillPDFs(wCTEQ66,wMRST2006NNLO,wNNPDF10100);
	
      }
      */ 

      //std::cout << "	run/lumi/event number " << runNumber << "/"<<lumiBlock<<"/"<<eventNumber<<std::endl;
      myOutTree[theChannel]->fillRunInfos(runNumber            , 
					  lumiBlock            , 
					  eventNumber          , 
					  weight               , 
					  passedHLT[theChannel], 
					  nputosave);
      
      // Index for the leading and subleading jets (obtained previously for each channel)
      int   theLJ = theLeadingJet[theChannel];
      int   theSJ = theSecondJet [theChannel];
      float ptLJ  = sqrt(  pxAK5PFPUcorrJet[theLJ]*pxAK5PFPUcorrJet[theLJ] 
			+ pyAK5PFPUcorrJet[theLJ]*pyAK5PFPUcorrJet[theLJ]);
      
      myOutTree[theChannel] -> fillAll(m_chMet          [theChannel], 
				       GetPt(m_p3PFMET->x(), m_p3PFMET->y()), GetPt(pxMet[0],pyMet[0]), 
				       m_projectedMet   [theChannel], 
				       m_deltaPhi       [theChannel], 
				       m_deltaErre      [theChannel], 
				       m_transvMass     [theChannel], 
				       m_mll            [theChannel], 
				       hardestLeptonPt  [theChannel], 
				       slowestLeptonPt  [theChannel], 
				       m_deltaEtaLeptons[theChannel], 
				       nGoodPV                      ,
				       true                         , 
				       true                         , 
				       true                         , 
				       true                         , 
				       isSelected);
      
      /// mmm channel
      if (theChannel == mmm){
	if(muCands[mmm].size()==3){
	  setLepIdVariables( (muCands[mmm])[0] , (muCands[mmm])[1] , (muCands[mmm])[2] , 13, 13, 13);
	}else{
	  setLepIdVariables(-1,-1,-1,13,13,13);
	}

      /// mme channel
      }else if (theChannel == mme){

	if(muCands[mme].size()==2 && eleCands[mme].size()==1){
	  setLepIdVariables( (muCands[mme])[0] , (muCands[mme])[1] , (eleCands[mme])[0] , 13, 13, 11);
	}else{
	  setLepIdVariables(-1,-1,-1,13,13,11);
	}

      /// eem channel
      }else if (theChannel == eem){

	if(muCands[eem].size()==1 && eleCands[eem].size()==2){
	  setLepIdVariables( (muCands[eem])[0] , (eleCands[eem])[0] , (eleCands[eem])[1] , 13, 11, 11);
	}else{
	  setLepIdVariables(-1,-1,-1,13,11,11);
	}
	
      /// eee channel
      }else if (theChannel == eee){

	if(eleCands[eee].size()==3){
	  setLepIdVariables( (eleCands[eee])[0] , (eleCands[eee])[1] , (eleCands[eee])[2] , 11, 11, 11);
	}else{
	  setLepIdVariables(-1,-1,-1,11,11,11);
	}

      }


      // mm channel, Z-peak
      if (theChannel == mm){
	if(muCands[mm].size()==2){
	  setLepIdVariables( (muCands[mm])[0] , (muCands[mm])[1] , -1 , 13, 13, 13 );
	}else{
	  setLepIdVariables(-1,-1,-1,13,13,13);
	}
      }

      // ee channel, Z-peak
      if (theChannel == ee){
	
	if(eleCands[ee].size()==2){
	  setLepIdVariables( (eleCands[ee])[0] , (eleCands[ee])[1] , -1 , 11, 11, 11);
	}else{
	  setLepIdVariables(-1,-1,-1,11,11,11);
	}
	
      } 


      myOutTree[theChannel] -> fillLeptons(myPt     , 
					   myEta    , 
					   myPhi    , 
					   myLepFlav, 
					   myLepId  ,
					   myLepIso , 
					   myConv);

      myOutTree[theChannel] -> fillMLVars(njets                      [theChannel], 
					  ncbIDjets                  [theChannel], 
					  nuncorrjets                [theChannel], 
					  m_maxDxyEvt                            , 
					  m_maxDszEvt                            , 
					  m_softbdisc                [theChannel],
					  m_hardbdisc                [theChannel], 
					  m_maxCombinedSecondaryVertexMVABJetTags, 
					  nsoftmu                    [theChannel], 
					  leadJetBtag                [theChannel], 
					  subleadJetBtag             [theChannel], 
					  subLeadJetsMaxBtag         [theChannel], 
					  leadJetBtagBProb           [theChannel], 
					  subleadJetBtagBProb        [theChannel], 
					  subLeadJetsMaxBtagBProb    [theChannel],
					  nextraleptons              [theChannel], 
					  nsoftmunojets              [theChannel],
					  m_numbtagjets              [theChannel], 
					  nsoftjets                  [theChannel],
					  m_numbtagCSVLcbIDcentraljets [theChannel],
					  m_numbtagCSVMcbIDcentraljets [theChannel],
					  m_numbtagCSVTcbIDcentraljets [theChannel],

					  m_numcbIDforwardjets [theChannel],

					  m_numbtagCSVLmvaIDcentraljets [theChannel],
					  m_numbtagCSVMmvaIDcentraljets [theChannel],
					  m_numbtagCSVTmvaIDcentraljets [theChannel],

					  m_nummvaIDforwardjets [theChannel]);
      
      if (wantJetIdStuff)     
	myOutTree[theChannel] -> fillJetsVars(leadJetPt          [theChannel],
					      leadJetEta         [theChannel], 
					      leadJetLoosePFId   [theChannel], 
					      leadJetMatchGen    [theChannel], 
					      leadJetMvaJetId    [theChannel], 
					      leadJetLooseId     [theChannel], 
					      subleadJetPt       [theChannel], 
					      subleadJetEta      [theChannel], 
					      subleadJetLoosePFId[theChannel], 
					      subleadJetMatchGen [theChannel], 
					      subleadJetMvaJetId [theChannel], 
					      subleadJetLooseId  [theChannel]);

      if (tHgenerated)
	myOutTree[theChannel] -> filltHMcTruthInfos(_genHiggsPt,
						    _genHiggsEta,
						    
						    _genTopPt,
						    _genTopEta,

						    _genWpfromH_Pt,
						    _genWpfromH_Eta,
						    
						    _genWmfromH_Pt,
						    _genWmfromH_Eta,
						    
						    _genWfromT_Pt,
						    _genWfromT_Eta,
						    
						    _genLeptonPlusfromWfromH_Pt,
						    _genLeptonPlusfromWfromH_Eta,
						    
						    _genLeptonMinusfromWfromH_Pt,
						    _genLeptonMinusfromWfromH_Eta,
						    
						    _genLeptonfromWfromT_Pt,
						    _genLeptonfromWfromT_Eta,
						    
						    _genNeutrinoPlusfromWfromH_Pt,
						    _genNeutrinoPlusfromWfromH_Eta,
						    
						    _genNeutrinoMinusfromWfromH_Pt,
						    _genNeutrinoMinusfromWfromH_Eta,
						    
						    _genNeutrinofromWfromT_Pt,
						    _genNeutrinofromWfromT_Eta,

						    _genForwardQuark_Pt,
						    _genForwardQuark_Eta );
      
      myOutTree[theChannel] -> fillLatinos ( outputStep0, 
					     outputStep1 ); 
      
      myOutTree[theChannel] -> fillRazor(m_MTR    [theChannel], 
					 m_MR     [theChannel], 
					 m_GammaMR[theChannel]);
      
      myOutTree[theChannel] -> fillMetStudies( m_projectedPFMet[theChannel], 
					       m_projectedTkMet[theChannel], 
					       significancePFMet        [0], 
					       significancePFChMet      [0], 
					       m_MTRcharged    [theChannel], 
					       -999.9                      ,  //dummy for dymva
					       rhoFastjet                  , 
					       rhoJetsFastJet); 
      
      std::vector<TLorentzVector> jesLJ = GetJetJesPcomponent(theLJ);
      std::vector<TLorentzVector> jesSJ = GetJetJesPcomponent(theSJ);

      float pxLJ[3] = { jesLJ[0].Px(), jesLJ[1].Px(), jesLJ[2].Px() };   
      float pyLJ[3] = { jesLJ[0].Py(), jesLJ[1].Py(), jesLJ[2].Py() };   
      float pzLJ[3] = { jesLJ[0].Pz(), jesLJ[1].Pz(), jesLJ[2].Pz() };
      float pxSJ[3] = { jesSJ[0].Px(), jesSJ[1].Px(), jesSJ[2].Px() };   
      float pySJ[3] = { jesSJ[0].Py(), jesSJ[1].Py(), jesSJ[2].Py() };   
      float pzSJ[3] = { jesSJ[0].Pz(), jesSJ[1].Pz(), jesSJ[2].Pz() };
      
      myOutTree[theChannel] -> fillKinematics( m_p3TKMET[theChannel].Px()   , 
					       m_p3TKMET[theChannel].Py()   , 
					       m_p3TKMET[theChannel].Pz()   , 
					       pxLJ                         , 
					       pyLJ                         , 
					       pzLJ                         , 
					       pxSJ                         ,  
					       pySJ                         , 
					       pzSJ                         ,
					       m_p4Lepton1[theChannel]->Px(), 
					       m_p4Lepton1[theChannel]->Py(), 
					       m_p4Lepton1[theChannel]->Pz(),
					       m_p4Lepton2[theChannel]->Px(), 
					       m_p4Lepton2[theChannel]->Py(), 
					       m_p4Lepton2[theChannel]->Pz(),
					       m_p4Lepton3[theChannel]->Px(), 
					       m_p4Lepton3[theChannel]->Py(), 
					       m_p4Lepton3[theChannel]->Pz(),
					       m_ch             [theChannel],
					       m_bdt            [theChannel],
					       m_jetsSum        [theChannel], 
					       m_uncorrJetsSum  [theChannel], 
					       m_p3PFMET); 
      
      // Dummy for the moment = -999., but to be changed to include the sc and r9 variables
      float dummyV[3];
      for (int ii=0; ii<3; ii++) dummyV[ii] = -999.;
      
      myOutTree[theChannel] -> fillSystematics(dummyV, dummyV, 
					       m_p4Lepton1Energy[theChannel], 
					       m_p4Lepton2Energy[theChannel],
					       m_p4Lepton3Energy[theChannel], 
					       m_p4Lepton1Type  [theChannel], 
					       m_p4Lepton2Type  [theChannel], 
					       m_p4Lepton3Type  [theChannel], 
					       m_metFromJets                , 
					       m_pfMetJESUp                 , 
					       m_pfMetJESDown               , 
					       jesMtUp          [theChannel], 
					       jesMtDown        [theChannel]);
      
      
      // dumping final tree, only if there are 3 final leptons
      if(m_channel[theChannel] && outputStep1) 
	myOutTree[theChannel] -> store();

    }

  }

  fMatch = new TFile("matching.root","RECREATE");
  fMatch->cd();
  H_deltaRuncorr->Write();
  H_deltaRcorr->Write();
  fMatch->Close();
}

void tHSelection::displayEfficiencies(std::string datasetName) {

  std::string::size_type loc = datasetName.find_first_of(".",0);
  if( loc != std::string::npos ) {
    datasetName.erase(loc);
  }
  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EEE selections: " << std::endl;
  CutBasedHiggsSelection[eee].displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full MMM selections: " << std::endl;
  CutBasedHiggsSelection[mmm].displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full EEM selections: " << std::endl;
  CutBasedHiggsSelection[eem].displayEfficiencies(datasetName);

  std::cout << "--------------------------------" << std::endl;
  std::cout << "Full MME selections: " << std::endl;
  CutBasedHiggsSelection[mme].displayEfficiencies(datasetName);

}

void tHSelection::counterEndJob(){
  
  std::cout << " " << std::endl;
  std::cout << " Total events:  " << std::endl;
  std::cout << nTotal << std::endl;
  std::cout << " Total events with H-->WW--> llnn:" << std::endl;
  std::cout << nWWtoLLdecays << std::endl;
  std::cout << " Total events with t-->Wb--> lnb:" << std::endl;
  std::cout << nWfromTtoLdecay << std::endl;
  std::cout << " Total events with SS leptons, coming from W from H and/or T:" << std::endl;
  std::cout << nSSevent << std::endl;

}

bool tHSelection::isSelectedMuon2012(int i){

  // HWW Selection : only not standalone MUON  
  Utils anaUtils;
  bool thisMuFlag1 = anaUtils.muonIdVal(muonIdMuon[i],AllGlobalMuons);
  bool thisMuFlag2 = anaUtils.muonIdVal(muonIdMuon[i],AllTrackerMuons);
  if (!thisMuFlag1 && !thisMuFlag2) return false;
  if (typeMuon[i]==8)               return false;  // to be used when new trees are available
  // Only not standalone muons latinos
  
  TLorentzVector Muon;
  Muon.SetPxPyPzE(double(pxMuon[i]), double(pyMuon[i]), double(pzMuon[i]), double(energyMuon[i]));

  RochCor2012 *tmprmcor;  
  tmprmcor = new RochCor2012();
  tmprmcor->momcor_mc(Muon, chargeMuon[i], 0, 1.0);
  
  //std::cout << " ** [GetBestMMM Function] :: Rochester Correction ** " << std::endl;
  //std::cout << " [Before] with TLorentz pT[muon] = [" << Muon.Pt() << "]" << std::endl;
  
  double MuonPt = double(GetPt(pxMuon[i],pyMuon[i]));
  //float MuonPt = Muon.Pt();
  /*
    std::cout << " Muon.Pt()  = " << Muon.Pt() << std::endl;
    std::cout << " GetPt(...) = " << MuonPt    << std::endl;
  */
  
  //std::cout << " [After] pT[muon] = [" << float(Muon.Pt()) << "]" << std::endl;
  
  //std::cout<< " muon [pT][eta] = ["<< MuonPt << "][" << etaMuon[i] << "]"<<std::endl;  
  
  if(_selectionEEE->getSwitch("etaMuonAcc") && !_selectionEEE->passCut("etaMuonAcc", etaMuon[i]) ) return false;
  if(_selectionEEE->getSwitch("ptMuonAcc")  && !_selectionEEE->passCut("ptMuonAcc" , MuonPt)     ) return false;
  
  bool theMuonID = true;
  isMuonID2012(i, &theMuonID);
  if (!theMuonID) return false;

  bool theMuonIso = true;
  isPFIsolatedMuon2012(i);
  if (!theMuonIso) return false;

  int ctfMuon   = trackIndexMuon[i]; 
  float dxyMuon = transvImpactParTrack[ctfMuon];
  float dzMuon  = muonDzPV(i,0);
  
  // Identification and Isolation MUON Selection
  if ( !theMuonID || !isPFIsolatedMuon2012(i) ) return false;      
  
  // Impact parameters MUON Selection
  if (MuonPt>20) {          // hardcoded
    if (_selectionEEE->getSwitch("muonIPhighPT") && (!_selectionEEE->passCut("muonIPhighPT",dxyMuon)) ) return false;   
  } else if (MuonPt<=20) {  // hardcoded
    if (_selectionEEE->getSwitch("muonIPlowPT")  && (!_selectionEEE->passCut("muonIPlowPT",dxyMuon)) ) return false;   
  }
  
  if (_selectionEEE->getSwitch("muonDz") && (!_selectionEEE->passCut("muonDz",dzMuon)) ) return false;         
  
  return true;

}

bool tHSelection::isSelectedMuon2013(int i){

  bool print = false;

  double MuonPt = double(GetPt(pxMuon[i],pyMuon[i]));

  // Kinematics
  if ( MuonPt < 5 || etaMuon[i] > 2.4 ) return false;

  if (print) std::cout<< " [Accep.]  = 1" << std::endl;

  if ( !pfmuonIdMuon[i] ) return false;

  int track      = trackIndexMuon[i];
  float theIp    = impactPar3DTrack[track];
  float theIpErr = impactPar3DErrorTrack[track];
  float sip3d    = theIp/theIpErr;
  if ( sip3d > 10 ) return false;

  int ctfMuon   = trackIndexMuon[i]; 
  float dxyMuon = transvImpactParTrack[ctfMuon];
  float dzMuon  = muonDzPV(i,0);
  if ( dxyMuon > 0.5 || dzMuon > 1.0 ) return false;

  // Isolation
  float relIso = (1.0/MuonPt)*(pfCandChargedIso04Muon[i] + TMath::Max(0.0,(pfCandNeutralIso04Muon[i]+pfCandPhotonIso04Muon[i]-0.5*pfCandChargedPUIso04Muon[i])));
  if ( relIso > 0.4 ) return false;

  // Lepton MVA ID
  float bdtValue = leptBDTForBs(fMVAMuon, i, false);
  if ( bdtValue < 0.7 ) return false;
  if (print) std::cout<< " [FinalBDT]= 1" << std::endl;  

  return true;

}

bool tHSelection::isSelectedElectron2012(int i){

  // HWW Selection

  bool print = false;
  
  float ElectronPt = GetPt(pxEle[i],pyEle[i]);
  
  if (print) std::cout<< " ELECTRON INFORMATION " << std::endl;
  if (print) std::cout<< " [pT][eta] = ["<< ElectronPt << "][" << etaEle[i] << "]"<<std::endl;  
  
  bool acceptanceEta = false;
  bool acceptancePt  = false;

  if(_selectionEEE->getSwitch("etaElectronAcc") && !_selectionEEE->passCut("etaElectronAcc",etaEle[i]) ) 
    return false;
  
  if(_selectionEEE->getSwitch("ptElectronAcc")  && !_selectionEEE->passCut("ptElectronAcc",ElectronPt) ) 
    return false;
  
  if (print) std::cout<< " [Accep.]  = 1" << std::endl;
  
  bool theElectronID, theElectronIsol, theElectronConvRej;
  theElectronID = theElectronIsol = theElectronConvRej = true;
  
  isEleID2012AndDenom(i,&theElectronID,&theElectronIsol,&theElectronConvRej);
  
  if (!theElectronID) return false;
  if (print) std::cout<< " [ID]      = " << theElectronID << std::endl;

  if (!theElectronIsol) return false;
  if (print) std::cout<< " [Iso]     = " << theElectronIsol << std::endl;

  if (!theElectronConvRej) return false;
  if (print) std::cout<< " [Con.Rej] = " << theElectronConvRej << std::endl;
  
  int   gsfTrack = gsfTrackIndexEle[i]; 
  float dxyEle = transvImpactParGsfTrack[gsfTrack];
  float dzEle  = eleDzPV(i,0);

  if (_selectionEEE->getSwitch("electronIP") && (!_selectionEEE->passCut("electronIP",dxyEle)) ) return false;
  if (_selectionEEE->getSwitch("electronDz") && (!_selectionEEE->passCut("electronDz",dzEle))  ) return false;
  
  return true;
  
}

bool tHSelection::isSelectedElectron2013(int i){

  // ttH Selection

  bool print = false;

  // Preselection requirements from AN-13-159

  float ElectronPt = GetPt(pxEle[i],pyEle[i]);
  
  if (print) std::cout<< " ELECTRON INFORMATION " << std::endl;
  if (print) std::cout<< " [pT][eta] = ["<< ElectronPt << "][" << etaEle[i] << "]"<<std::endl;    
  
  // Kinematics
  if ( ElectronPt < 7 || etaEle[i] > 2.5 ) return false;

  if (print) std::cout<< " [Accep.]  = 1" << std::endl;

  // BDT selection
  bool passBDT = false;
  if ( ElectronPt > 5 && ElectronPt < 10 ) {
    if       ( abs(etaEle[i])<0.8 ){
      if (eleBDT(fMVA,i) > 0.47) passBDT = true; 
    }else if ( abs(etaEle[i])>0.8 && abs(etaEle[i]) < 1.479 ) {
      if (eleBDT(fMVA,i) > 0.004) passBDT = true; 
    }else if ( abs(etaEle[i]) > 1.479 ) {
      if (eleBDT(fMVA,i) > 0.295) passBDT = true; 
    }
  }else if ( ElectronPt > 10 ){
    if       ( abs(etaEle[i])<0.8 ){
      if (eleBDT(fMVA,i) > 0.5) passBDT = true; 
    }else if ( abs(etaEle[i])>0.8 && abs(etaEle[i]) < 1.479 ) {
      if (eleBDT(fMVA,i) > 0.12) passBDT = true; 
    }else if ( abs(etaEle[i]) > 1.479 ) {
      if (eleBDT(fMVA,i) > 0.6) passBDT = true; 
    }
  }

  if (!passBDT) return false;
  if (print) std::cout<< " [BDTpre] = 1" << std::endl;  

  // Impact parameter cuts
  int gsftrack   = gsfTrackIndexEle[i];
  float theIp    = impactPar3DGsfTrack[gsftrack];
  float theIpErr = impactPar3DErrorGsfTrack[gsftrack]; 
  float sip3d    = theIp/theIpErr;
  if ( sip3d > 10 ) return false;
  
  float dxyEle = transvImpactParGsfTrack[gsftrack];
  float dzEle  = eleDzPV(i,0);      
  if ( dxyEle > 0.5 || dzEle > 1.0 ) return false;
  
  if (print) std::cout<< " [Imp.Pa.] = 1" << std::endl;  

  // Conversion rejection
  int missHits = expInnerLayersGsfTrack[gsftrack];
  if ( missHits > 1 ) return true;
  
  if (print) std::cout<< " [Conv.Re] = 1" << std::endl;  

  // Isolation
  float relIso = (1.0/ElectronPt)*(pfCandChargedIso04Ele[i] + TMath::Max(0.0,(pfCandNeutralIso04Ele[i]+pfCandPhotonIso04Ele[i]-0.5*pfCandChargedPUIso04Ele[i])));
  if ( relIso > 0.4 ) return false;

  if (print) std::cout<< " [Iso]     = 1" << std::endl;  

  // Lepton MVA ID
  float bdtValue = leptBDTForBs(fMVAElectron, i, true);
  if ( bdtValue < 0.7 ) return false;
  if (print) std::cout<< " [FinalBDT]= 1" << std::endl;  
  
  return true;

}

std::vector<int> tHSelection::getBestEEE(){

  selectedElectrons.clear();

  // Verbose
  bool check = false;

  // index for the three selected electrons
  int index1 = -1;
  int index2 = -1;
  int index3 = -1;
  
  // Maps for positive and negative electrons, separately
  std::map<float,std::pair<int,int> > PositiveElecIndexMap;
  std::map<float,std::pair<int,int> > NegativeElecIndexMap;
  
  // Map for all electrons (do not care about the charge)
  std::map<float,std::pair<int,int> > AllElecIndexMap;

  int NumberPositiveElec = 0;
  int NumberNegativeElec = 0;
  int NumberAllElec      = 0;

  if (check){
    std::cout<<""<<std::endl;
    std::cout<< " ... looping over the most general electron collection:" << std::endl;
  }

  for(int i=0;i<nEle;i++) {
    
    float ElectronPt     = GetPt(pxEle[i],pyEle[i]);
    float ElectronCharge = chargeEle[i];

    // Select the lepton ID
    bool tightElectronSelection     = true;
    bool LeptonMVAElectronSelection = false;

    if (tightElectronSelection){

      if (! isSelectedElectron2012(i))continue;

    }else if(LeptonMVAElectronSelection){

      if (! isSelectedElectron2013(i)) continue;

    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
    
    // Store this muon in the AllElecIndexMap
    AllElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i);
    NumberAllElec++;

    if (ElectronCharge < 0 ){ 
      NegativeElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i); 
      NumberNegativeElec++; 
    }else{ 
      PositiveElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i); 
      NumberPositiveElec++;
    }
    
  }
  
  // CHECK THE MAPS
  if(check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllElec;
    for(_iterAllElec = AllElecIndexMap.rbegin(); _iterAllElec != AllElecIndexMap.rend(); ++_iterAllElec) {
      std::cout << "AllElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllElec).first << " ][" 
		<< ((*_iterAllElec).second).first  << "," 
		<< ((*_iterAllElec).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterPosElec;
    for(_iterPosElec = PositiveElecIndexMap.rbegin(); _iterPosElec != PositiveElecIndexMap.rend(); ++_iterPosElec) {
      std::cout << "PositiveElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterPosElec).first << " ][" 
		<< ((*_iterPosElec).second).first  << "," 
		<< ((*_iterPosElec).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterNegElec;
    for(_iterNegElec = NegativeElecIndexMap.rbegin(); _iterNegElec != NegativeElecIndexMap.rend(); ++_iterNegElec) {
      std::cout << "NegativeElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterNegElec).first << " ][" 
		<< ((*_iterNegElec).second).first  << "," 
		<< ((*_iterNegElec).second).second << "]" << std::endl;
    }
  }

  double Pt1 = -999;
  double Pt2 = -999;
  double Pt3 = -999;

  double charge1 = 0;
  double charge2 = 0;
  double charge3 = 0;

  //SELECT THE HIGHEST PT ELECTRON
  
  if( NumberAllElec > 2 && NumberNegativeElec > 0 && NumberPositiveElec > 0 ){  

    int selectedElectrons = 0;

    std::map<float,std::pair<int,int> >::reverse_iterator it;
    for(it=AllElecIndexMap.rbegin();it!=AllElecIndexMap.rend();++it){

      if (selectedElectrons ==0){

	// The first electron we select is the one with highest pT
	index1  = ((*it).second).second;
	charge1 = ((*it).second).first;
	Pt1     =  (*it).first;

	selectedElectrons++;
	
      }else if(selectedElectrons == 1){
	
	// The second electron we select is the one with second highest pT
	index2  = ((*it).second).second;
	charge2 = ((*it).second).first;
	Pt2     =  (*it).first;

	selectedElectrons++;

      }else if(selectedElectrons == 2){

	int charge = ((*it).second).first;	

	if ( charge1*charge < 0 || charge2*charge < 0){
	  
	  index3  = ((*it).second).second;
	  charge3 = ((*it).second).first;
	  Pt3     = (*it).first;

	  selectedElectrons++;
	  break;
	  
	}else{
	  continue;
	}
	
      }
      
    }
    
    if (check){
    std::cout << "***********************" << std::endl;
    std::cout << "eee CHANNEL" << std::endl;
      std::cout << "the electron with 1st highest pT is electron[pT][ch,index] = ["
		<<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl; 
      std::cout << "the electron with 2nd highest pT is electron[pT][ch,index] = ["
		<<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl;     
      std::cout << "the electron with 3rd highest pT is electron[pT][ch,index] = ["
		<<Pt3<<"]["<<charge3<<","<<index3<<"]"<< std::endl; 
    }

  }

  PositiveElecIndexMap.clear();
  NegativeElecIndexMap.clear();
  AllElecIndexMap.clear();

  selectedElectrons.push_back(index1);
  selectedElectrons.push_back(index2);
  selectedElectrons.push_back(index3);

  return selectedElectrons;

}//getBestEEE()

std::vector<int> tHSelection::getBestMMM(){

  selectedMuons.clear();

  // Verbose
  bool check = false;
  
  // index for the three selected muons
  int index1 = -1;
  int index2 = -1;
  int index3 = -1;
  
  // Maps for positive and negative muons, separately
  std::map<float,std::pair<int,int> > PositiveMuonsIndexMap;
  std::map<float,std::pair<int,int> > NegativeMuonsIndexMap;
  
  // Map for all muons (do not care about the charge)
  std::map<float,std::pair<int,int> > AllMuonsIndexMap;

  int NumberPositiveMuons = 0;
  int NumberNegativeMuons = 0;
  int NumberAllMuons      = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general muon collection:" << std::endl;
  }

  RochCor2012 *tmprmcor;  

  for(int i=0;i<nMuon;i++) {
    
    float MuonPt     = GetPt(pxMuon[i],pyMuon[i]);
    float MuonCharge = chargeMuon[i];
    
    // Select the lepton ID
    bool tightMuonSelection     = true;
    bool LeptonMVAMuonSelection = false;
    
    if (tightMuonSelection){
      
      if (! isSelectedMuon2012(i))continue;
      
    }else if(LeptonMVAMuonSelection){
      
      if (! isSelectedMuon2013(i)) continue;
      
    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
   
    // Store this muon in the AllMuonsIndexMap
    AllMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i);
    NumberAllMuons++;

    if ( MuonCharge < 0 ){ 
      NegativeMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i); 
      NumberNegativeMuons++; 
    }else{ 
      PositiveMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i); 
      NumberPositiveMuons++;
    }
    
  }

  // CHECK THE MAPS
  if (check){

    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllMuons;
    for(_iterAllMuons = AllMuonsIndexMap.rbegin(); _iterAllMuons != AllMuonsIndexMap.rend(); ++_iterAllMuons) {
      std::cout << "AllMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllMuons).first << " ][" 
		<< ((*_iterAllMuons).second).first  << "," 
		<< ((*_iterAllMuons).second).second << "]" << std::endl;
    }
  
    std::map<float,std::pair<int,int> >::reverse_iterator _iterPosMuons;
    for(_iterPosMuons = PositiveMuonsIndexMap.rbegin(); _iterPosMuons != PositiveMuonsIndexMap.rend(); ++_iterPosMuons) {
      std::cout << "PositiveMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterPosMuons).first << " ][" 
		<< ((*_iterPosMuons).second).first  << "," 
		<< ((*_iterPosMuons).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterNegMuons;
    for(_iterNegMuons = NegativeMuonsIndexMap.rbegin(); _iterNegMuons != NegativeMuonsIndexMap.rend(); ++_iterNegMuons) {
      std::cout << "NegativeMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterNegMuons).first << " ][" 
		<< ((*_iterNegMuons).second).first  << "," 
		<< ((*_iterNegMuons).second).second << "]" << std::endl;
    }
  }

  float Pt1 = -999;
  float Pt2 = -999;
  float Pt3 = -999;

  double charge1 = 0;
  double charge2 = 0;
  double charge3 = 0;

  //SELECT THE HIGHEST PT MUON
  
  if( NumberAllMuons > 2 && NumberNegativeMuons > 0 && NumberPositiveMuons > 0 ){  

    int selectedMuons = 0;

    std::map<float,std::pair<int,int> >::reverse_iterator it;
    for(it=AllMuonsIndexMap.rbegin();it!=AllMuonsIndexMap.rend();++it){

      if (selectedMuons ==0){

	// The first muon we select is the one with highest pT
	index1  = ((*it).second).second;
	charge1 = ((*it).second).first;
	Pt1     =  (*it).first;

	selectedMuons++;
	
      }else if(selectedMuons == 1){

	// The second muon we select is the one with second highest pT
	index2  = ((*it).second).second;
	charge2 = ((*it).second).first;
	Pt2     =  (*it).first;

	selectedMuons++;

      }else if(selectedMuons == 2){

	// For the third muon selection, we take the muon with third highest pT that
	// also has opposite charge with respect to the second muon (we asumme the 
	// 2nd and 3rd muons come from the H --> W+W- --> l+l-nn decay

	int charge = ((*it).second).first;	

	if ( charge1*charge < 0 || charge2*charge < 0){
	  
	  index3  = ((*it).second).second;
	  charge3 = ((*it).second).first;
	  Pt3     = (*it).first;

	  selectedMuons++;
	  break;
	  
	}else{
	  continue;
	}
	
      }
      
    }

    if (check){
    std::cout << "***********************" << std::endl;
    std::cout << "mmm CHANNEL" << std::endl;
      std::cout << "the muon with 1st highest pT is muon[pT][ch,index] = ["
		<<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl; 
      std::cout << "the muon with 2nd highest pT is muon[pT][ch,index] = ["
		<<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl;     
      std::cout << "the muon with 3rd highest pT is muon[pT][ch,index] = ["
		<<Pt3<<"]["<<charge3<<","<<index3<<"]"<< std::endl; 
    }

  }

  PositiveMuonsIndexMap.clear();
  NegativeMuonsIndexMap.clear();
  AllMuonsIndexMap.clear();

  selectedMuons.push_back(index1);
  selectedMuons.push_back(index2);
  selectedMuons.push_back(index3);

  return selectedMuons;

}//getBestMMM()

std::vector<int> tHSelection::getBestMME(){

  selectedLeptons.clear();

  // Verbose
  bool check = false;

  // index for the three selected leptons
  int index1 = -1; // Muon
  int index2 = -1; // Muon
  int index3 = -1; // Electron
  
  // Map for muons
  // -First - element of the map: pT of the lepton
  // -Second- element of the map: a pair containing the charge and the index in the collection
  std::map<float,std::pair<int,int> > AllMuonsIndexMap;

  int NumberAllMuons      = 0;
  int NumberPositiveMuons = 0;
  int NumberNegativeMuons = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general muon collection:" << std::endl;
  }

  for(int i=0;i<nMuon;i++) {
    
    float MuonPt     = GetPt(pxMuon[i],pyMuon[i]);
    float MuonCharge = chargeMuon[i];
    
    // Select the lepton ID
    bool tightMuonSelection     = true;
    bool LeptonMVAMuonSelection = false;
    
    if (tightMuonSelection){
      
      if (! isSelectedMuon2012(i))continue;
      
    }else if(LeptonMVAMuonSelection){
      
      if (! isSelectedMuon2013(i)) continue;
      
    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }

    // Store this muon in the AllMuonsIndexMap
    AllMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i);
    NumberAllMuons++;

    // we also store the number of positive and negative muons
    if     (MuonCharge < 0 ) NumberNegativeMuons++; 
    else if(MuonCharge > 0 ) NumberPositiveMuons++;
    
  }
  
  // CHECK THE MAPS
  if (check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllMuons;
    for(_iterAllMuons = AllMuonsIndexMap.rbegin(); _iterAllMuons != AllMuonsIndexMap.rend(); ++_iterAllMuons) {
      std::cout << "AllMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllMuons).first << " ][" 
		<< ((*_iterAllMuons).second).first  << "," 
		<< ((*_iterAllMuons).second).second << "]" << std::endl;
      
    }
  }
  
  // Map for all electrons
  // -First - element of the map: pT of the lepton
  // -Second- element of the map: a pair containing the charge and the index in the collection
  std::map<float,std::pair<int,int> >AllElecIndexMap;

  int NumberAllElec      = 0;
  int NumberPositiveElec = 0;
  int NumberNegativeElec = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general electron collection:" << std::endl;
  }

  for(int i=0;i<nEle;i++) {

    float ElectronPt     = GetPt(pxEle[i],pyEle[i]);
    float ElectronCharge = chargeEle[i];
    
    // Select the lepton ID
    bool tightElectronSelection     = true;
    bool LeptonMVAElectronSelection = false;
    
    if (tightElectronSelection){
      
      if (! isSelectedElectron2012(i))continue;

    }else if(LeptonMVAElectronSelection){

      if (! isSelectedElectron2013(i)) continue;

    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
    
    AllElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i);
    NumberAllElec++;

    // we also store the number of positive and negative electrons
    if     (ElectronCharge < 0 ) NumberNegativeElec++; 
    else if(ElectronCharge > 0 ) NumberPositiveElec++;

  }
  
  // CHECK THE MAP
  if (check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllElec;
    for(_iterAllElec = AllElecIndexMap.rbegin(); _iterAllElec != AllElecIndexMap.rend(); ++_iterAllElec) {
      std::cout << "AllElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllElec).first << " ][" 
		<< ((*_iterAllElec).second).first  << "," 
		<< ((*_iterAllElec).second).second << "]" << std::endl;
      
    }
  }

  double Pt1 = -999;
  double Pt2 = -999;
  double Pt3 = -999;

  double charge1 = 0;
  double charge2 = 0;
  double charge3 = 0;

  // 2 muons + 1 electron
  
  // Possible combinations:
  //  ----------
  // | m+ m+ e- |
  // | m- m- e+ |
  // | m+ m- e- |
  // | m+ m- e+ |
  //  ----------

  if( NumberAllMuons > 1 && NumberAllElec > 0 ){

    // (m+ m- e+) OR (m+ m- e-)  
    if( NumberNegativeMuons > 0 && NumberPositiveMuons > 0 && (NumberNegativeElec > 0 || NumberPositiveElec > 0) ){

      if(check) std::cout << " --> We are in the case (m+ m- e+) OR (m+ m- e-)" << std::endl;
      
      // Select the muons (with opposite charges)
      
      // Loop over the muon collection for the first muon 
      std::map<float,std::pair<int,int> >::reverse_iterator it1;
      for(it1=AllMuonsIndexMap.rbegin();it1!=AllMuonsIndexMap.rend();++it1){
	
	// The first muon we select is the one with highest pT
	index1  = ((*it1).second).second;
	charge1 = ((*it1).second).first;
	Pt1     =  (*it1).first;
	
	break;

      }

      // Loop over the muon collection for the second muon selection
      std::map<float,std::pair<int,int> >::reverse_iterator it2;
      for(it2=AllMuonsIndexMap.rbegin();it2!=AllMuonsIndexMap.rend();++it2){
	
	int tmpcharge2 = ((*it2).second).first;	
	
	// The second muon we select is the one with second highest pT
	// also has opposite charge with respect to the first muon 
	if ( tmpcharge2*charge1 > 0 ) continue;
	
	index2  = ((*it2).second).second;
	charge2 = ((*it2).second).first;
	Pt2     = (*it2).first;
	
	break;
	
      }

      // Select the electron
      std::map<float,std::pair<int,int> >::reverse_iterator itEle;
      for(itEle=AllElecIndexMap.rbegin();itEle!=AllElecIndexMap.rend();++itEle){
	
	// We select the one with highest pT, in this channel, we do not care about the charge ...
	index3  = ((*itEle).second).second;
	charge3 = ((*itEle).second).first;
	Pt3     =  (*itEle).first;
	
	break;
	
      }
      
    }//end (m+ m- e+) OR (m+ m- e+)
    
    // (m+ m+ e-) OR (m- m- e+)
    else if ( (NumberPositiveMuons  > 1 && NumberNegativeMuons == 0 && NumberNegativeElec > 0)
	      ||
	      (NumberPositiveMuons == 0 && NumberNegativeMuons  > 1 && NumberPositiveElec > 0) ){
      
      if(check) std::cout << " --> We are in the case (m+ m+ e-) OR (m- m- e+) " << std::endl;
      
      int nselectedmuons = 0;
      
      // Select the muons 
      
      // Loop over the muon collection for the first muon 
      std::map<float,std::pair<int,int> >::reverse_iterator it1;
      for(it1=AllMuonsIndexMap.rbegin();it1!=AllMuonsIndexMap.rend();++it1){
	
	if (nselectedmuons == 0){
	  // The first muon we select is the one with highest pT
  
	  index1  = ((*it1).second).second;
	  charge1 = ((*it1).second).first;
	  Pt1     =  (*it1).first;
	  
	  nselectedmuons++;
	  
	}else if (nselectedmuons == 1){
	  // Once we have the first muon, we take the second one
	  // This is supposed to be of the same charge, as we have asked for that within the 'if'

	  index2  = ((*it1).second).second;
	  charge2 = ((*it1).second).first;
	  Pt2     =  (*it1).first;
	  
	  nselectedmuons++;
	  
	}else if (nselectedmuons > 1){
	  // Go out the loop, we have two muons already
	  
	  break;

	}
	
      }
      // Select the electron
      std::map<float,std::pair<int,int> >::reverse_iterator itEle;
      for(itEle=AllElecIndexMap.rbegin();itEle!=AllElecIndexMap.rend();++itEle){
	
	// We select the one with highest pT, in this channel, AND opposite charge wrt the muons
	
	int tmpcharge3 = ((*itEle).second).first;

	// check that the ch (muon) and the ch (electron) are different
	if (tmpcharge3 * charge1 > 0) continue;

	index3  = ((*itEle).second).second;
	charge3 = ((*itEle).second).first;
	Pt3     =  (*itEle).first;
	
	break;      
	
      }
      
    }// end (m+ m+ e-) OR (m- m- e+)
    
  }

  if(check){
    std::cout << "***********************" << std::endl;
    std::cout << "mme CHANNEL" << std::endl;
    std::cout << "the muon with highest pT is muon[pT][ch,index] = ["
	      <<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl;     
    std::cout << "the muon with second highest pT is muon[pT][ch,index] = ["
	      <<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl; 
    std::cout << "the selected electron is electron[pT][ch,index] = ["
	      <<Pt3<<"]["<<charge3<<","<<index3<<"]"<< std::endl; 
  }

  AllMuonsIndexMap.clear();
  AllElecIndexMap.clear();
  
  selectedLeptons.push_back(index1); // first muon
  selectedLeptons.push_back(index2); // second muon
  selectedLeptons.push_back(index3); // electron
  
  return selectedLeptons;
  
}//getBestMME

std::vector<int> tHSelection::getBestEEM(){

  selectedLeptons.clear();
  
  // Verbose
  bool check = false;

  // index for the three selected leptons
  int index1 = -1; // Electron
  int index2 = -1; // Electron
  int index3 = -1; // Muon
  
  // Map for muons
  // -First - element of the map: pT of the lepton
  // -Second- element of the map: a pair containing the charge and the index in the collection
  std::map<float,std::pair<int,int> > AllMuonsIndexMap;

  int NumberAllMuons      = 0;
  int NumberPositiveMuons = 0;
  int NumberNegativeMuons = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general muon collection:" << std::endl;
  }

  for(int i=0;i<nMuon;i++) {
    
   float MuonPt     = GetPt(pxMuon[i],pyMuon[i]);
    float MuonCharge = chargeMuon[i];
    
    // Select the lepton ID
    bool tightMuonSelection     = true;
    bool LeptonMVAMuonSelection = false;
    
    if (tightMuonSelection){
      
      if (! isSelectedMuon2012(i))continue;
      
    }else if(LeptonMVAMuonSelection){
      
      if (! isSelectedMuon2013(i)) continue;
      
    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }

    // Store this muon in the AllMuonsIndexMap
    AllMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i);
    NumberAllMuons++;

    // we also store the number of positive and negative muons
    if     (MuonCharge < 0 ) NumberNegativeMuons++; 
    else if(MuonCharge > 0 ) NumberPositiveMuons++;
    
  }
  
  // CHECK THE MAPS
  if(check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllMuons;
    for(_iterAllMuons = AllMuonsIndexMap.rbegin(); _iterAllMuons != AllMuonsIndexMap.rend(); ++_iterAllMuons) {
      
      std::cout << "AllMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllMuons).first << " ][" 
		<< ((*_iterAllMuons).second).first  << "," 
		<< ((*_iterAllMuons).second).second << "]" << std::endl;
      
    }
  }

  // Map for all electrons
  // -First - element of the map: pT of the lepton
  // -Second- element of the map: a pair containing the charge and the index in the collection
  std::map<float,std::pair<int,int> >AllElecIndexMap;

  int NumberAllElec      = 0;
  int NumberPositiveElec = 0;
  int NumberNegativeElec = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<<" ... Looping over the most general electron collection:" << std::endl;
  }

  for(int i=0;i<nEle;i++) {
    
    float ElectronPt     = GetPt(pxEle[i],pyEle[i]);
    float ElectronCharge = chargeEle[i];
    
    // Select the lepton ID
    bool tightElectronSelection     = true;
    bool LeptonMVAElectronSelection = false;
    
    if (tightElectronSelection){
      
      if (! isSelectedElectron2012(i))continue;

    }else if(LeptonMVAElectronSelection){

      if (! isSelectedElectron2013(i)) continue;

    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
    
    AllElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i);
    NumberAllElec++;

    // we also store the number of positive and negative electrons
    if     (ElectronCharge < 0 ) NumberNegativeElec++; 
    else if(ElectronCharge > 0 ) NumberPositiveElec++;

  }
  
  // CHECK THE MAP
  if(check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllElec;
    for(_iterAllElec = AllElecIndexMap.rbegin(); _iterAllElec != AllElecIndexMap.rend(); ++_iterAllElec) {
      std::cout << "AllElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllElec).first << " ][" 
		<< ((*_iterAllElec).second).first  << "," 
		<< ((*_iterAllElec).second).second << "]" << std::endl;
      
    }
  }

  double Pt1 = -999;
  double Pt2 = -999;
  double Pt3 = -999;

  double charge1 = 0;
  double charge2 = 0;
  double charge3 = 0;

  // 2 electrons + 1 muon
  
  // Possible combinations:
  //  ----------
  // | e+ e+ m- |
  // | e- e- m+ |
  // | e+ e- m- |
  // | e+ e- m+ |
  //  ----------

  if( NumberAllElec > 1 && NumberAllMuons > 0 ){

    // (e+ e- m+) OR (e+ e- m-)  
    if( NumberNegativeElec > 0 && NumberPositiveElec > 0 && (NumberNegativeMuons > 0 || NumberPositiveMuons > 0) ){

      if(check) std::cout << " --> We are in the case (e+ e- m+) OR (e+ e- m-)" << std::endl;
      
      // Select the electrons (with opposite charges)
      
      // Loop over the electron collection for the first electron
      std::map<float,std::pair<int,int> >::reverse_iterator it1;
      for(it1=AllElecIndexMap.rbegin();it1!=AllElecIndexMap.rend();++it1){
	
	// The first electron we select is the one with highest pT
	index1  = ((*it1).second).second;
	charge1 = ((*it1).second).first;
	Pt1     =  (*it1).first;
	
	break;

      }

      // Loop over the electron collection for the second electron selection
      std::map<float,std::pair<int,int> >::reverse_iterator it2;
      for(it2=AllElecIndexMap.rbegin();it2!=AllElecIndexMap.rend();++it2){
	
	int tmpcharge2 = ((*it2).second).first;	
	
	// The second electron we select is the one with second highest pT
	// also has opposite charge with respect to the first electron
	if ( tmpcharge2*charge1 > 0 ) continue;
	
	index2  = ((*it2).second).second;
	charge2 = ((*it2).second).first;
	Pt2     = (*it2).first;
	
	break;
	
      }

      // Select the muon
      std::map<float,std::pair<int,int> >::reverse_iterator itMuon;
      for(itMuon=AllMuonsIndexMap.rbegin();itMuon!=AllMuonsIndexMap.rend();++itMuon){
	
	// We select the one with highest pT, in this channel, we do not care about the charge ...
	index3  = ((*itMuon).second).second;
	charge3 = ((*itMuon).second).first;
	Pt3     =  (*itMuon).first;
	
	break;
	
      }
      
    }//end (e+ e- m+) OR (e+ e- m-)
    
    // (e+ e+ m-) OR (e- e- m+)
    else if ( (NumberPositiveElec  > 1 && NumberNegativeElec == 0 && NumberNegativeMuons > 0)
	      ||
	      (NumberPositiveElec == 0 && NumberNegativeElec  > 1 && NumberPositiveMuons > 0) ){
      
      if(check) std::cout << " --> We are in the case (e+ e+ m-) OR (e- e- m+) " << std::endl;
      
      int nselectedelectrons = 0;
      
      // Select the electrons
      
      // Loop over the electron collection for the first electron
      std::map<float,std::pair<int,int> >::reverse_iterator it1;
      for(it1=AllElecIndexMap.rbegin();it1!=AllElecIndexMap.rend();++it1){
	
	if (nselectedelectrons == 0){
	  // The first electron we select is the one with highest pT
  
	  index1  = ((*it1).second).second;
	  charge1 = ((*it1).second).first;
	  Pt1     =  (*it1).first;
	  
	  nselectedelectrons++;
	  
	}else if (nselectedelectrons == 1){
	  // Once we have the first electron, we take the second one
	  // This is supposed to be of the same charge, as we have asked for that within the 'if'

	  index2  = ((*it1).second).second;
	  charge2 = ((*it1).second).first;
	  Pt2     =  (*it1).first;
	  
	  nselectedelectrons++;
	  
	}else if (nselectedelectrons > 1){
	  // Go out the loop, we have two electrons already
	  
	  break;

	}
	
      }
      
      
      // Select the muon
      std::map<float,std::pair<int,int> >::reverse_iterator itMuon;
      for(itMuon=AllMuonsIndexMap.rbegin();itMuon!=AllMuonsIndexMap.rend();++itMuon){
	
	// We select the one with highest pT, in this channel, AND opposite charge wrt the electrons
	
	int tmpcharge3 = ((*itMuon).second).first;

	// check that the ch (electron) and the ch (muon) are different
	if (tmpcharge3 * charge1 > 0) continue;

	index3  = ((*itMuon).second).second;
	charge3 = ((*itMuon).second).first;
	Pt3     =  (*itMuon).first;
	
	break;      
	
      }
      
    }// end (e+ e+ m-) OR (e- e- m+)
    
  }
  
  if(check){
    std::cout << "***********************" << std::endl;
    std::cout << "emm CHANNEL" << std::endl;
    std::cout << "the electron with highest pT is electron[pT][ch,index] = ["
	      <<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl;     
    std::cout << "the electron with second highest pT is electron[pT][ch,index] = ["
	      <<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl; 
    std::cout << "the selected muon is muon[pT][ch,index] = ["
	      <<Pt3<<"]["<<charge3<<","<<index3<<"]"<< std::endl; 
  }

  AllMuonsIndexMap.clear();
  AllElecIndexMap.clear();
  
  selectedLeptons.push_back(index1); // first electron
  selectedLeptons.push_back(index2); // second electron
  selectedLeptons.push_back(index3); // muons
  
  return selectedLeptons;
  
}//getBestEEM

std::vector<int> tHSelection::getBestEE(){

  selectedElectrons.clear();

  // Verbose
  bool check = true;

  // index for the two selected electrons
  int index1 = -1;
  int index2 = -1;
  
  // Maps for positive and negative electrons, separately
  std::map<float,std::pair<int,int> > PositiveElecIndexMap;
  std::map<float,std::pair<int,int> > NegativeElecIndexMap;
  
  // Map for all electrons (do not care about the charge)
  std::map<float,std::pair<int,int> > AllElecIndexMap;

  int NumberPositiveElec = 0;
  int NumberNegativeElec = 0;
  int NumberAllElec      = 0;

  if (check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general electron collection:" << std::endl;
  }

  for(int i=0;i<nEle;i++) {
    
    float ElectronPt     = GetPt(pxEle[i],pyEle[i]);
    float ElectronCharge = chargeEle[i];

    // Select the lepton ID
    bool tightElectronSelection     = true;
    bool LeptonMVAElectronSelection = false;

    if (tightElectronSelection){

      if (! isSelectedElectron2012(i))continue;

    }else if(LeptonMVAElectronSelection){

      if (! isSelectedElectron2013(i)) continue;

    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
    
    // Store this muon in the AllElecIndexMap
    AllElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i);
    NumberAllElec++;

    if (ElectronCharge < 0 ){ 
      NegativeElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i); 
      NumberNegativeElec++; 
    }else{ 
      PositiveElecIndexMap[ElectronPt] = std::make_pair(ElectronCharge,i); 
      NumberPositiveElec++;
    }
    
  }
  
  // CHECK THE MAPS
  if(check){
    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllElec;
    for(_iterAllElec = AllElecIndexMap.rbegin(); _iterAllElec != AllElecIndexMap.rend(); ++_iterAllElec) {
      std::cout << "AllElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllElec).first << " ][" 
		<< ((*_iterAllElec).second).first  << "," 
		<< ((*_iterAllElec).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterPosElec;
    for(_iterPosElec = PositiveElecIndexMap.rbegin(); _iterPosElec != PositiveElecIndexMap.rend(); ++_iterPosElec) {
      std::cout << "PositiveElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterPosElec).first << " ][" 
		<< ((*_iterPosElec).second).first  << "," 
		<< ((*_iterPosElec).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterNegElec;
    for(_iterNegElec = NegativeElecIndexMap.rbegin(); _iterNegElec != NegativeElecIndexMap.rend(); ++_iterNegElec) {
      std::cout << "NegativeElecIndexMap[pT][ch,index] = [";
      std::cout << (*_iterNegElec).first << " ][" 
		<< ((*_iterNegElec).second).first  << "," 
		<< ((*_iterNegElec).second).second << "]" << std::endl;
    }
  }

  double Pt1 = -999;
  double Pt2 = -999;

  double charge1 = 0;
  double charge2 = 0;

  //SELECT THE HIGHEST PT ELECTRON
  
  if( NumberAllElec > 1 && NumberNegativeElec > 0 && NumberPositiveElec > 0 ){  

    int selectedElectrons = 0;

    std::map<float,std::pair<int,int> >::reverse_iterator it;
    for(it=AllElecIndexMap.rbegin();it!=AllElecIndexMap.rend();++it){

      if (selectedElectrons ==0){

	// The first electron we select is the one with highest pT
	index1  = ((*it).second).second;
	charge1 = ((*it).second).first;
	Pt1     =  (*it).first;

	selectedElectrons++;
	
      }else if(selectedElectrons == 1){

	// The second electron we select is the one with second highest pT and opposite sign wrt first electron	
	int charge = ((*it).second).first;
	if ( charge1*charge < 0 ){

	  index2  = ((*it).second).second;
	  charge2 = ((*it).second).first;
	  Pt2     =  (*it).first;
	  
	  selectedElectrons++;
	  break;
	}
	
      }else{
	continue;
      }
      
    }
    
    if (check){
    std::cout << "***********************" << std::endl;
    std::cout << "ee CHANNEL" << std::endl;
      std::cout << "the electron with 1st highest pT is electron[pT][ch,index] = ["
		<<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl; 
      std::cout << "the electron with 2nd highest pT is electron[pT][ch,index] = ["
		<<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl;     
    }

  }

  PositiveElecIndexMap.clear();
  NegativeElecIndexMap.clear();
  AllElecIndexMap.clear();

  selectedElectrons.push_back(index1);
  selectedElectrons.push_back(index2);

  return selectedElectrons;

}//getBestEE()

std::vector<int> tHSelection::getBestMM(){
  
  selectedMuons.clear();

  // Verbose
  bool check = true;
  
  // index for the two selected muons
  int index1 = -1;
  int index2 = -1;
  
  // Maps for positive and negative muons, separately
  std::map<float,std::pair<int,int> > PositiveMuonsIndexMap;
  std::map<float,std::pair<int,int> > NegativeMuonsIndexMap;
  
  // Map for all muons (do not care about the charge)
  std::map<float,std::pair<int,int> > AllMuonsIndexMap;

  int NumberPositiveMuons = 0;
  int NumberNegativeMuons = 0;
  int NumberAllMuons      = 0;

  if(check){
    std::cout<<""<<std::endl;
    std::cout<< " ... Looping over the most general muon collection:" << std::endl;
  }

  RochCor2012 *tmprmcor;  

  for(int i=0;i<nMuon;i++) {
    
    float MuonPt     = GetPt(pxMuon[i],pyMuon[i]);
    float MuonCharge = chargeMuon[i];
    
    // Select the lepton ID
    bool tightMuonSelection     = true;
    bool LeptonMVAMuonSelection = false;
    
    if (tightMuonSelection){
      
      if (! isSelectedMuon2012(i))continue;
      
    }else if(LeptonMVAMuonSelection){
      
      if (! isSelectedMuon2013(i)) continue;
      
    }else{
      std::cout << "ERROR: No ID Selected!!" << std::endl;
      break;
    }
   
    // Store this muon in the AllMuonsIndexMap
    AllMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i);
    NumberAllMuons++;

    if ( MuonCharge < 0 ){ 
      NegativeMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i); 
      NumberNegativeMuons++; 
    }else{ 
      PositiveMuonsIndexMap[MuonPt] = std::make_pair(MuonCharge,i); 
      NumberPositiveMuons++;
    }
    
  }

  // CHECK THE MAPS
  if (check){

    std::map<float,std::pair<int,int> >::reverse_iterator _iterAllMuons;
    for(_iterAllMuons = AllMuonsIndexMap.rbegin(); _iterAllMuons != AllMuonsIndexMap.rend(); ++_iterAllMuons) {
      std::cout << "AllMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterAllMuons).first << " ][" 
		<< ((*_iterAllMuons).second).first  << "," 
		<< ((*_iterAllMuons).second).second << "]" << std::endl;
    }
  
    std::map<float,std::pair<int,int> >::reverse_iterator _iterPosMuons;
    for(_iterPosMuons = PositiveMuonsIndexMap.rbegin(); _iterPosMuons != PositiveMuonsIndexMap.rend(); ++_iterPosMuons) {
      std::cout << "PositiveMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterPosMuons).first << " ][" 
		<< ((*_iterPosMuons).second).first  << "," 
		<< ((*_iterPosMuons).second).second << "]" << std::endl;
    }
    
    std::map<float,std::pair<int,int> >::reverse_iterator _iterNegMuons;
    for(_iterNegMuons = NegativeMuonsIndexMap.rbegin(); _iterNegMuons != NegativeMuonsIndexMap.rend(); ++_iterNegMuons) {
      std::cout << "NegativeMuonsIndexMap[pT][ch,index] = [";
      std::cout << (*_iterNegMuons).first << " ][" 
		<< ((*_iterNegMuons).second).first  << "," 
		<< ((*_iterNegMuons).second).second << "]" << std::endl;
    }
  }

  float Pt1 = -999;
  float Pt2 = -999;

  double charge1 = 0;
  double charge2 = 0;

  //SELECT THE HIGHEST PT MUON
  
  if( NumberAllMuons > 1 && NumberNegativeMuons > 0 && NumberPositiveMuons > 0 ){  

    int selectedMuons = 0;

    std::map<float,std::pair<int,int> >::reverse_iterator it;
    for(it=AllMuonsIndexMap.rbegin();it!=AllMuonsIndexMap.rend();++it){

      if (selectedMuons ==0){

	// The first muon we select is the one with highest pT
	index1  = ((*it).second).second;
	charge1 = ((*it).second).first;
	Pt1     =  (*it).first;

	selectedMuons++;
	
      }else if(selectedMuons == 1){

	int charge = ((*it).second).first;	

	if ( charge1*charge < 0 ){

	  // The second muon we select is the one with second highest pT and opposite charge wrt first muon
	  index2  = ((*it).second).second;
	  charge2 = ((*it).second).first;
	  Pt2     =  (*it).first;
	  
	  selectedMuons++;
	  break;
	}
	
      }else{
	continue;
      }
      
    }
    
    if (check){
    std::cout << "***********************" << std::endl;
    std::cout << "mm CHANNEL" << std::endl;
      std::cout << "the muon with 1st highest pT is muon[pT][ch,index] = ["
		<<Pt1<<"]["<<charge1<<","<<index1<<"]"<< std::endl; 
      std::cout << "the muon with 2nd highest pT is muon[pT][ch,index] = ["
		<<Pt2<<"]["<<charge2<<","<<index2<<"]"<< std::endl;     
    }

  }

  PositiveMuonsIndexMap.clear();
  NegativeMuonsIndexMap.clear();
  AllMuonsIndexMap.clear();

  selectedMuons.push_back(index1);
  selectedMuons.push_back(index2);

  return selectedMuons;

}//getBestMM()

void tHSelection::setKinematicsEEE(int myEle1, int myEle2, int myEle3){

  // Check that we have three selected electrons
  if ( myEle1 > -1 && myEle2 > -1 && myEle3 > -1 ) {

    // Definition of the variables to be used in this function
    float Electron1Charge = chargeEle[myEle1];  
    float Electron2Charge = chargeEle[myEle2];  
    float Electron3Charge = chargeEle[myEle3];  
    
    float Electron1Pt = GetPt(pxEle[myEle1],pyEle[myEle1]);
    float Electron2Pt = GetPt(pxEle[myEle2],pyEle[myEle2]);
    float Electron3Pt = GetPt(pxEle[myEle3],pyEle[myEle3]);
    
    TLorentzVector Electron1, Electron2, Electron3;
    Electron1.SetPxPyPzE(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    Electron2.SetPxPyPzE(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    Electron3.SetPxPyPzE(pxEle[myEle3], pyEle[myEle3], pzEle[myEle3], energyEle[myEle3]);
    
    // Store the selected electrons in eleCands[]. To be used later
    eleCands[eee].push_back(myEle1);
    eleCands[eee].push_back(myEle2);
    eleCands[eee].push_back(myEle3);
    
    // Look for the pair with closest m_ll to m_Z
    float m12, m13, m23;
    m12 = -999.9;
    m13 = -999.9;
    m23 = -999.9;
    // Must have opposite sign. Estimate the mass
    if(Electron1Charge*Electron2Charge < 0){
      m12 = (Electron1 + Electron2).M();
    }
    if(Electron1Charge*Electron3Charge < 0){
      m13 = (Electron1 + Electron3).M();
    }
    if(Electron2Charge*Electron3Charge < 0){
      m23 = (Electron2 + Electron3).M();
    }
    
    float closestmll = -999.9;
    float diffmll    = fabs(m12 - 91.1876); // start with m12
    // if m12 is negative, then diffml is going to be > 1091,
    // so in the second 'if' we will take m13 as the closest one. 
    if ( m12 > 0 ) {
      closestmll = m12; 
    }
    if ( m13 > 0 && ( fabs(m13 - 91.1876) < diffmll ) ) {
      closestmll = m13;
      diffmll = fabs(m13-91.1876);
    }
    if ( m23 > 0 && ( fabs(m23 - 91.1876) < diffmll ) ){
      closestmll = m23;
      diffmll = fabs(m23-91.1876);
    }
    
    hardestLeptonPt[eee] = TMath::Max( TMath::Max(Electron1Pt,Electron2Pt), Electron3Pt );
    slowestLeptonPt[eee] = TMath::Min( TMath::Min(Electron1Pt,Electron2Pt), Electron3Pt );
    
    m_p4Lepton1[eee] -> SetXYZT(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    m_p4Lepton2[eee] -> SetXYZT(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    m_p4Lepton3[eee] -> SetXYZT(pxEle[myEle3], pyEle[myEle3], pzEle[myEle3], energyEle[myEle3]);
    
    m_p4Lepton1Energy[eee] = energyEle[myEle1];
    m_p4Lepton2Energy[eee] = energyEle[myEle2];
    m_p4Lepton3Energy[eee] = energyEle[myEle3];
    
    m_p4Lepton1Type[eee] = 0;
    m_p4Lepton2Type[eee] = 0;
    m_p4Lepton3Type[eee] = 0;
    
    m_mll            [eee] = closestmll;
    m_deltaPhi       [eee] = -999.9;
    m_deltaErre      [eee] = -999.9;
    m_deltaEtaLeptons[eee] = -999.9;
    
    m_dilepPt [eee].SetXYZ( 0.0, 0.0, 0.0 );
    m_trilepPt[eee].SetXYZ( m_p4Lepton1[eee]->Vect().X()+m_p4Lepton2[eee]->Vect().X()+m_p4Lepton3[eee]->Vect().X(),
			    m_p4Lepton2[eee]->Vect().Y()+m_p4Lepton2[eee]->Vect().Y()+m_p4Lepton3[eee]->Vect().Y(), 0.0 );
    // usual definition
    m_transvMass[eee] = -999.9;
    // chris' variable
    m_GammaMR[eee] = -999.9;
    
    m_MR            [eee] = -999.9;
    m_MTR           [eee] = -999.9;
    m_metOptll      [eee] = -999.9;
    m_mT2           [eee] = -999.9;
    m_p3TKMET       [eee] = -999.9;
    m_chMet         [eee] = -999.9;
    m_projectedMet  [eee] = -999.9;
    m_projectedPFMet[eee] = -999.9;
    m_projectedTkMet[eee] = -999.9;
    m_MTRcharged    [eee] = -999.9;

    m_ch[eee][0] = chargeEle[myEle1];
    m_ch[eee][1] = chargeEle[myEle2];
    m_ch[eee][2] = chargeEle[myEle3];
    m_iso[eee][0] = pfCombinedIsoEle[myEle1] / Electron1Pt;
    m_iso[eee][1] = pfCombinedIsoEle[myEle2] / Electron2Pt;
    m_iso[eee][2] = pfCombinedIsoEle[myEle3] / Electron3Pt;
    m_lh[eee][0] = eleIdLikelihoodEle[myEle1];
    m_lh[eee][1] = eleIdLikelihoodEle[myEle2];    
    m_lh[eee][2] = eleIdLikelihoodEle[myEle3];    
    /*
    m_bdt[eee][0] = eleBDT(fMVA,myEle1);
    m_bdt[eee][1] = eleBDT(fMVA,myEle2);
    m_bdt[eee][2] = eleBDT(fMVA,myEle3);
    */
    m_bdt[eee][0] = leptBDTForBs(fMVAElectron,myEle1,true);
    m_bdt[eee][1] = leptBDTForBs(fMVAElectron,myEle2,true);
    m_bdt[eee][2] = leptBDTForBs(fMVAElectron,myEle3,true);

    m_chmaj[eee][0] = eleChargeMajority(myEle1);
    m_chmaj[eee][1] = eleChargeMajority(myEle2);    
    m_chmaj[eee][2] = eleChargeMajority(myEle3);

    Electron1.Clear();
    Electron2.Clear();
    Electron3.Clear();

  }

}//setKinematicsEEE()

void tHSelection::setKinematicsMMM(int myMu1, int myMu2, int myMu3){

  // Check that we have three selected electrons
  if ( myMu1 > -1 && myMu2 > -1 && myMu3 > -1 ) {

    // Definition of the variables to be used in this function
    float Muon1Charge = chargeMuon[myMu1];  
    float Muon2Charge = chargeMuon[myMu2];  
    float Muon3Charge = chargeMuon[myMu3];  
    
    float Muon1Pt = GetPt(pxMuon[myMu1],pyMuon[myMu1]);
    float Muon2Pt = GetPt(pxMuon[myMu2],pyMuon[myMu2]);
    float Muon3Pt = GetPt(pxMuon[myMu3],pyMuon[myMu3]);
    
    TLorentzVector Muon1, Muon2, Muon3;
    Muon1.SetPxPyPzE(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    Muon2.SetPxPyPzE(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    Muon3.SetPxPyPzE(pxMuon[myMu3], pyMuon[myMu3], pzMuon[myMu3], energyMuon[myMu3]);
    
    // Store the selected muons in muCands[]. To be used later
    muCands[mmm].push_back(myMu1);
    muCands[mmm].push_back(myMu2);
    muCands[mmm].push_back(myMu3);
    
    // Look for the pair with closest m_ll to m_Z
    float m12, m13, m23;
    m12 = -999.9;
    m13 = -999.9;
    m23 = -999.9;
    // Must have opposite sign. Estimate the mass
    if(Muon1Charge*Muon2Charge < 0){
      m12 = (Muon1 + Muon2).M();
    }
    if(Muon1Charge*Muon3Charge < 0){
      m13 = (Muon1 + Muon3).M();
    }
    if(Muon2Charge*Muon3Charge < 0){
      m23 = (Muon2 + Muon3).M();
    }

    //std::cout << "m12 = " << m12 << std::endl;
    //std::cout << "m13 = " << m13 << std::endl;
    //std::cout << "m23 = " << m23 << std::endl;
    
    float closestmll = -999.9;
    float diffmll    = fabs(m12 - 91.1876); // start with m12
    // if m12 is negative, then diffml is going to be > 1091,
    // so in the second 'if' we will take m13 as the closest one. 
    if ( m12 > 0 ) {
      closestmll = m12; 
    }
    if ( m13 > 0 && ( fabs(m13 - 91.1876) < diffmll ) ) {
      closestmll = m13;
      diffmll = fabs(m13-91.1876);
    }
    if ( m23 > 0 && ( fabs(m23 - 91.1876) < diffmll ) ){
      closestmll = m23;
      diffmll = fabs(m23-91.1876);
    }
    
    //std::cout << " The selected minv is " << closestmll << std::endl;

    hardestLeptonPt[mmm] = TMath::Max( TMath::Max(Muon1Pt,Muon2Pt), Muon3Pt );
    slowestLeptonPt[mmm] = TMath::Min( TMath::Min(Muon1Pt,Muon2Pt), Muon3Pt );
    
    m_p4Lepton1[mmm] -> SetXYZT(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    m_p4Lepton2[mmm] -> SetXYZT(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    m_p4Lepton3[mmm] -> SetXYZT(pxMuon[myMu3], pyMuon[myMu3], pzMuon[myMu3], energyMuon[myMu3]);
    
    m_p4Lepton1Energy[mmm] = energyMuon[myMu1];
    m_p4Lepton2Energy[mmm] = energyMuon[myMu2];
    m_p4Lepton3Energy[mmm] = energyMuon[myMu3];
    
    m_p4Lepton1Type[mmm] = 1;
    m_p4Lepton2Type[mmm] = 1;
    m_p4Lepton3Type[mmm] = 1;
    
    m_mll            [mmm] = closestmll;
    m_deltaPhi       [mmm] = -999.9;
    m_deltaErre      [mmm] = -999.9;
    m_deltaEtaLeptons[mmm] = -999.9;
    
    m_dilepPt [mmm].SetXYZ( 0.0, 0.0, 0.0 );
    m_trilepPt[mmm].SetXYZ( m_p4Lepton1[mmm]->Vect().X()+m_p4Lepton2[mmm]->Vect().X()+m_p4Lepton3[mmm]->Vect().X(),
			    m_p4Lepton2[mmm]->Vect().Y()+m_p4Lepton2[mmm]->Vect().Y()+m_p4Lepton3[mmm]->Vect().Y(), 0.0 );
    // usual definition
    m_transvMass[mmm] = -999.9;
    // chris' variable
    m_GammaMR[mmm] = -999.9;
    
    m_MR            [mmm] = -999.9;
    m_MTR           [mmm] = -999.9;
    m_metOptll      [mmm] = -999.9;
    m_mT2           [mmm] = -999.9;
    m_p3TKMET       [mmm] = -999.9;
    m_chMet         [mmm] = -999.9;
    m_projectedMet  [mmm] = -999.9;
    m_projectedPFMet[mmm] = -999.9;
    m_projectedTkMet[mmm] = -999.9;
    m_MTRcharged    [mmm] = -999.9;

    m_ch[mmm][0] = chargeMuon[myMu1];
    m_ch[mmm][1] = chargeMuon[myMu2];
    m_ch[mmm][2] = chargeMuon[myMu3];
    m_iso[mmm][0] = pfCombinedIsoMuon[myMu1] / Muon1Pt;
    m_iso[mmm][1] = pfCombinedIsoMuon[myMu2] / Muon2Pt;
    m_iso[mmm][2] = pfCombinedIsoMuon[myMu3] / Muon3Pt;
    m_lh[mmm][0] = -999.9;
    m_lh[mmm][1] = -999.9;    
    m_lh[mmm][2] = -999.9;    
    /*
    m_bdt[mmm][0] = -999.9;
    m_bdt[mmm][1] = -999.9;
    m_bdt[mmm][2] = -999.9;
    */
    m_bdt[mmm][0] = leptBDTForBs(fMVAMuon,myMu1,false);
    m_bdt[mmm][1] = leptBDTForBs(fMVAMuon,myMu2,false);
    m_bdt[mmm][2] = leptBDTForBs(fMVAMuon,myMu3,false);

    m_chmaj[mmm][0] = -999.9;
    m_chmaj[mmm][1] = -999.9;    
    m_chmaj[mmm][2] = -999.9;

    Muon1.Clear();
    Muon2.Clear();
    Muon3.Clear();

  }

}

void tHSelection::setKinematicsMME(int myMu1 , int myMu2 , int myEle1){

  // Check that we have three selected leptons
  if ( myMu1 > -1 && myMu2 > -1 && myEle1 > -1 ){

    // Definition of the variables to be used in this function
    float Muon1Charge     = chargeMuon[myMu1];  
    float Muon2Charge     = chargeMuon[myMu2];  
    float Electron1Charge = chargeEle[myEle1];  
    
    float Muon1Pt     = GetPt(pxMuon[myMu1],pyMuon[myMu1]);
    float Muon2Pt     = GetPt(pxMuon[myMu2],pyMuon[myMu2]);
    float Electron1Pt = GetPt(pxEle[myEle1],pyEle[myEle1]);
    
    TLorentzVector Muon1, Muon2, Electron1;
    Muon1.SetPxPyPzE    (pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    Muon2.SetPxPyPzE    (pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    Electron1.SetPxPyPzE(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
        
    // Store the selected muons in muCands[]
    muCands[mme].push_back(myMu1);
    muCands[mme].push_back(myMu2);
    
    // Store the selected electrons in eleCands[]
    eleCands[mme].push_back(myEle1);
    
    // Look for the pair with closest m_ll to m_Z
    float m12, m13, m23;
    m12 = -999.9;
    m13 = -999.9;
    m23 = -999.9;
    // Must have opposite sign. Estimate the mass
    if(Muon1Charge*Muon2Charge < 0){
      m12 = (Muon1 + Muon2).M();
    }
    if(Muon1Charge*Electron1Charge < 0){
      m13 = (Muon1 + Electron1).M();
    }
    if(Muon2Charge*Electron1Charge < 0){
      m23 = (Muon2 + Electron1).M();
    }

    //std::cout << "m12 = " << m12 << std::endl;
    //std::cout << "m13 = " << m13 << std::endl;
    //std::cout << "m23 = " << m23 << std::endl;

    
    float closestmll = -999.9;
    float diffmll    = fabs(m12 - 91.1876); // start with m12
    // if m12 is negative, then diffml is going to be > 1091,
    // so in the second 'if' we will take m13 as the closest one. 
    if ( m12 > 0 ) {
      closestmll = m12; 
    }
    if ( m13 > 0 && ( fabs(m13 - 91.1876) < diffmll ) ) {
      closestmll = m13;
      diffmll = fabs(m13-91.1876);
    }
    if ( m23 > 0 && ( fabs(m23 - 91.1876) < diffmll ) ){
      closestmll = m23;
      diffmll = fabs(m23-91.1876);
    }
    
    //std::cout << " The selected minv is " << closestmll << std::endl;

    hardestLeptonPt[mme] = TMath::Max( TMath::Max(Muon1Pt,Muon2Pt), Electron1Pt );
    slowestLeptonPt[mme] = TMath::Min( TMath::Min(Muon1Pt,Muon2Pt), Electron1Pt );
    
    m_p4Lepton1[mme] -> SetXYZT(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    m_p4Lepton2[mme] -> SetXYZT(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    m_p4Lepton3[mme] -> SetXYZT(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    
    m_p4Lepton1Energy[mme] = energyMuon[myMu1];
    m_p4Lepton2Energy[mme] = energyMuon[myMu2];
    m_p4Lepton3Energy[mme] = energyEle[myEle1];
    
    m_p4Lepton1Type[mme] = 1;
    m_p4Lepton2Type[mme] = 1;
    m_p4Lepton3Type[mme] = 0;
    
    m_mll            [mme] = closestmll;
    m_deltaPhi       [mme] = -999.9;
    m_deltaErre      [mme] = -999.9;
    m_deltaEtaLeptons[mme] = -999.9;
    
    m_dilepPt [mme].SetXYZ( 0.0, 0.0, 0.0 );
    m_trilepPt[mme].SetXYZ( m_p4Lepton1[mme]->Vect().X()+m_p4Lepton2[mme]->Vect().X()+m_p4Lepton3[mme]->Vect().X(),
			    m_p4Lepton1[mme]->Vect().Y()+m_p4Lepton2[mme]->Vect().Y()+m_p4Lepton3[mme]->Vect().Y(), 0.0 );
    // usual definition
    m_transvMass[mme] = -999.9;
    // chris' variable
    m_GammaMR[mme] = -999.9;
    
    m_MR            [mme] = -999.9;
    m_MTR           [mme] = -999.9;
    m_metOptll      [mme] = -999.9;
    m_mT2           [mme] = -999.9;
    m_p3TKMET       [mme] = -999.9;
    m_chMet         [mme] = -999.9;
    m_projectedMet  [mme] = -999.9;
    m_projectedPFMet[mme] = -999.9;
    m_projectedTkMet[mme] = -999.9;
    m_MTRcharged    [mme] = -999.9;

    m_ch[mme][0] = chargeMuon[myMu1];
    m_ch[mme][1] = chargeMuon[myMu2];
    m_ch[mme][2] = chargeEle[myEle1];
    m_iso[mme][0] = pfCombinedIsoMuon[myMu1] / Muon1Pt;
    m_iso[mme][1] = pfCombinedIsoMuon[myMu2] / Muon2Pt;
    m_iso[mme][2] = pfCombinedIsoEle[myEle1] / Electron1Pt;
    m_lh[mme][0] = -999.9;
    m_lh[mme][1] = -999.9;    
    m_lh[mme][2] = eleIdLikelihoodEle[myEle1];    
    /*
    m_bdt[mme][0] = -999.9;
    m_bdt[mme][1] = -999.9;
    m_bdt[mme][2] = eleBDT(fMVA,myEle1);
    */
    m_bdt[mme][0] = leptBDTForBs(fMVAMuon,myMu1,false);
    m_bdt[mme][1] = leptBDTForBs(fMVAMuon,myMu2,false);
    m_bdt[mme][2] = leptBDTForBs(fMVAElectron,myEle1,true);
    m_chmaj[mme][0] = -999.9;
    m_chmaj[mme][1] = -999.9;    
    m_chmaj[mme][2] = eleChargeMajority(myEle1);

    Muon1.Clear();
    Muon2.Clear();
    Electron1.Clear();

  }

}

void tHSelection::setKinematicsEEM(int myEle1 , int myEle2 , int myMu1){

  // Check that we have three selected leptons
  if ( myEle1 > -1 && myEle2 > -1 && myMu1 > -1 ){

    // Definition of the variables to be used in this function
    float Electron1Charge = chargeEle[myEle1];  
    float Electron2Charge = chargeEle[myEle2];  
    float Muon1Charge     = chargeMuon[myMu1];  

    float Electron1Pt = GetPt(pxEle[myEle1],pyEle[myEle1]);    
    float Electron2Pt = GetPt(pxEle[myEle2],pyEle[myEle2]);    
    float Muon1Pt     = GetPt(pxMuon[myMu1],pyMuon[myMu1]);
    
    TLorentzVector Electron1, Electron2, Muon1;
    Electron1.SetPxPyPzE(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    Electron2.SetPxPyPzE(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    Muon1.SetPxPyPzE    (pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    
    // Store the selected muons in muCands[]
    muCands[eem].push_back(myMu1);
    
    // Store the selected electrons in eleCands[]
    eleCands[eem].push_back(myEle1);
    eleCands[eem].push_back(myEle2);
    
    // Look for the pair with closest m_ll to m_Z
    float m12, m13, m23;
    m12 = -999.9;
    m13 = -999.9;
    m23 = -999.9;
    // Must have opposite sign. Estimate the mass
    if(Electron1Charge*Electron2Charge < 0){
      m12 = (Electron1 + Electron2).M();
    }
    if(Electron1Charge*Muon1Charge < 0){
      m13 = (Electron1 + Muon1).M();
    }
    if(Electron2Charge*Muon1Charge < 0){
      m23 = (Electron2 + Muon1).M();
    }

    //std::cout << "m12 = " << m12 << std::endl;
    //std::cout << "m13 = " << m13 << std::endl;
    //std::cout << "m23 = " << m23 << std::endl;

    
    float closestmll = -999.9;
    float diffmll    = fabs(m12 - 91.1876); // start with m12
    // if m12 is negative, then diffml is going to be > 1091,
    // so in the second 'if' we will take m13 as the closest one. 
    if ( m12 > 0 ) {
      closestmll = m12; 
    }
    if ( m13 > 0 && ( fabs(m13 - 91.1876) < diffmll ) ) {
      closestmll = m13;
      diffmll = fabs(m13-91.1876);
    }
    if ( m23 > 0 && ( fabs(m23 - 91.1876) < diffmll ) ){
      closestmll = m23;
      diffmll = fabs(m23-91.1876);
    }
    
    //std::cout << " The selected minv is " << closestmll << std::endl;

    hardestLeptonPt[eem] = TMath::Max( TMath::Max(Electron1Pt,Electron2Pt), Muon1Pt );
    slowestLeptonPt[eem] = TMath::Min( TMath::Min(Electron1Pt,Electron2Pt), Muon1Pt );
    
    m_p4Lepton1[eem] -> SetXYZT(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    m_p4Lepton2[eem] -> SetXYZT(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    m_p4Lepton3[eem] -> SetXYZT(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    
    m_p4Lepton1Energy[eem] = energyEle[myEle1];
    m_p4Lepton2Energy[eem] = energyEle[myEle2];
    m_p4Lepton3Energy[eem] = energyMuon[myMu1];
    
    m_p4Lepton1Type[eem] = 0; // electron
    m_p4Lepton2Type[eem] = 0; // electron
    m_p4Lepton3Type[eem] = 1; // muon
    
    m_mll            [eem] = closestmll;
    m_deltaPhi       [eem] = -999.9;
    m_deltaErre      [eem] = -999.9;
    m_deltaEtaLeptons[eem] = -999.9;
    
    m_dilepPt [eem].SetXYZ( 0.0, 0.0, 0.0 );
    m_trilepPt[eem].SetXYZ( m_p4Lepton1[eem]->Vect().X()+m_p4Lepton2[eem]->Vect().X()+m_p4Lepton3[eem]->Vect().X(),
			    m_p4Lepton1[eem]->Vect().Y()+m_p4Lepton2[eem]->Vect().Y()+m_p4Lepton3[eem]->Vect().Y(), 0.0 );
    // usual definition
    m_transvMass[eem] = -999.9;
    // chris' variable
    m_GammaMR[eem] = -999.9;
    
    m_MR            [eem] = -999.9;
    m_MTR           [eem] = -999.9;
    m_metOptll      [eem] = -999.9;
    m_mT2           [eem] = -999.9;
    m_p3TKMET       [eem] = -999.9;
    m_chMet         [eem] = -999.9;
    m_projectedMet  [eem] = -999.9;
    m_projectedPFMet[eem] = -999.9;
    m_projectedTkMet[eem] = -999.9;
    m_MTRcharged    [eem] = -999.9;

    m_ch[eem][0] = chargeEle[myEle1];
    m_ch[eem][1] = chargeEle[myEle2];
    m_ch[eem][2] = chargeMuon[myMu1];
    m_iso[eem][0] = pfCombinedIsoEle[myEle1] / Electron1Pt;
    m_iso[eem][1] = pfCombinedIsoEle[myEle2] / Electron2Pt;
    m_iso[eem][2] = pfCombinedIsoMuon[myMu1] / Muon1Pt;
    m_lh[eem][0] = eleIdLikelihoodEle[myEle1];
    m_lh[eem][1] = eleIdLikelihoodEle[myEle2];    
    m_lh[eem][2] = -999.9;    
    /*
    m_bdt[eem][0] = eleBDT(fMVA,myEle1);
    m_bdt[eem][1] = eleBDT(fMVA,myEle2);
    m_bdt[eem][2] = -999.9;
    */
    m_bdt[eem][0] = leptBDTForBs(fMVAElectron,myEle1,true);
    m_bdt[eem][2] = leptBDTForBs(fMVAElectron,myEle2,true);
    m_bdt[eem][3] = leptBDTForBs(fMVAMuon,myMu1,false);
    m_chmaj[eem][0] = eleChargeMajority(myEle1);
    m_chmaj[eem][1] = eleChargeMajority(myEle2);    
    m_chmaj[eem][2] = -999.9;

    Electron1.Clear();
    Electron2.Clear();
    Muon1.Clear();
    
  }

}

void tHSelection::setKinematicsEE(int myEle1, int myEle2){

  // Check that we have two selected electrons
  if ( myEle1 > -1 && myEle2 > -1 ) {

    // Definition of the variables to be used in this function
    float Electron1Charge = chargeEle[myEle1];  
    float Electron2Charge = chargeEle[myEle2];  
    
    float Electron1Pt = GetPt(pxEle[myEle1],pyEle[myEle1]);
    float Electron2Pt = GetPt(pxEle[myEle2],pyEle[myEle2]);
    
    TLorentzVector Electron1, Electron2, Electron3;
    Electron1.SetPxPyPzE(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    Electron2.SetPxPyPzE(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    
    // Store the selected electrons in eleCands[]. To be used later
    eleCands[ee].push_back(myEle1);
    eleCands[ee].push_back(myEle2);
    
    // Look for the pair with closest m_ll to m_Z
    float mll;
    mll = (Electron1 + Electron2).M();
    
    hardestLeptonPt[ee] = TMath::Max( Electron1Pt,Electron2Pt );
    slowestLeptonPt[ee] = TMath::Min( Electron1Pt,Electron2Pt );
    
    m_p4Lepton1[ee] -> SetXYZT(pxEle[myEle1], pyEle[myEle1], pzEle[myEle1], energyEle[myEle1]);
    m_p4Lepton2[ee] -> SetXYZT(pxEle[myEle2], pyEle[myEle2], pzEle[myEle2], energyEle[myEle2]);
    m_p4Lepton3[ee] -> SetXYZT(0.0, 0.0, 0.0, 0.0); // dummy
    
    m_p4Lepton1Energy[ee] = energyEle[myEle1];
    m_p4Lepton2Energy[ee] = energyEle[myEle2];
    m_p4Lepton3Energy[ee] = -999.9;
    
    m_p4Lepton1Type[ee] = 0;
    m_p4Lepton2Type[ee] = 0;
    m_p4Lepton3Type[ee] = 0; // don't care
    
    m_mll            [ee] = mll;
    m_deltaPhi       [ee] = -999.9;
    m_deltaErre      [ee] = -999.9;
    m_deltaEtaLeptons[ee] = -999.9;
    
    m_dilepPt [ee].SetXYZ( m_p4Lepton1[ee]->Vect().X()+m_p4Lepton2[ee]->Vect().X(), 
			   m_p4Lepton2[ee]->Vect().Y()+m_p4Lepton2[ee]->Vect().Y(), 0.0 );

    m_trilepPt[ee].SetXYZ( 0.0, 0.0, 0.0 );

    // usual definition
    m_transvMass[ee] = -999.9;
    // chris' variable
    m_GammaMR   [ee] = -999.9;
    
    m_MR            [ee] = -999.9;
    m_MTR           [ee] = -999.9;
    m_metOptll      [ee] = -999.9;
    m_mT2           [ee] = -999.9;
    m_p3TKMET       [ee] = -999.9;
    m_chMet         [ee] = -999.9;
    m_projectedMet  [ee] = -999.9;
    m_projectedPFMet[ee] = -999.9;
    m_projectedTkMet[ee] = -999.9;
    m_MTRcharged    [ee] = -999.9;

    m_ch[ee][0] = chargeEle[myEle1];
    m_ch[ee][1] = chargeEle[myEle2];
    m_ch[ee][2] = -999.9;
    m_iso[ee][0] = pfCombinedIsoEle[myEle1] / Electron1Pt;
    m_iso[ee][1] = pfCombinedIsoEle[myEle2] / Electron2Pt;
    m_iso[ee][2] = -999.9;
    m_lh[ee][0] = eleIdLikelihoodEle[myEle1];
    m_lh[ee][1] = eleIdLikelihoodEle[myEle2];    
    m_lh[ee][2] = -999.9;    
    /*
    m_bdt[eee][0] = eleBDT(fMVA,myEle1);
    m_bdt[eee][1] = eleBDT(fMVA,myEle2);
    m_bdt[eee][2] = eleBDT(fMVA,myEle3);
    */
    m_bdt[ee][0] = leptBDTForBs(fMVAElectron,myEle1,true);
    m_bdt[ee][1] = leptBDTForBs(fMVAElectron,myEle2,true);
    m_bdt[ee][2] = -999.9;

    m_chmaj[ee][0] = eleChargeMajority(myEle1);
    m_chmaj[ee][1] = eleChargeMajority(myEle2);    
    m_chmaj[ee][2] = -999.9;

    Electron1.Clear();
    Electron2.Clear();

  }

}//setKinematicsEE()

void tHSelection::setKinematicsMM(int myMu1, int myMu2){

  // Check that we have two selected muons
  if ( myMu1 > -1 && myMu2 > -1 ) {

    // Definition of the variables to be used in this function
    float Muon1Charge = chargeMuon[myMu1];  
    float Muon2Charge = chargeMuon[myMu2];  
    
    float Muon1Pt = GetPt(pxMuon[myMu1],pyMuon[myMu1]);
    float Muon2Pt = GetPt(pxMuon[myMu2],pyMuon[myMu2]);
    
    TLorentzVector Muon1, Muon2, Muon3;
    Muon1.SetPxPyPzE(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    Muon2.SetPxPyPzE(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    
    // Store the selected muons in muCands[]. To be used later
    muCands[mm].push_back(myMu1);
    muCands[mm].push_back(myMu2);
    
    // Look for the pair with closest m_ll to m_Z
    float mll = (Muon1 + Muon2).M();

    hardestLeptonPt[mm] = TMath::Max( Muon1Pt,Muon2Pt );
    slowestLeptonPt[mm] = TMath::Min( Muon1Pt,Muon2Pt );
    
    m_p4Lepton1[mm] -> SetXYZT(pxMuon[myMu1], pyMuon[myMu1], pzMuon[myMu1], energyMuon[myMu1]);
    m_p4Lepton2[mm] -> SetXYZT(pxMuon[myMu2], pyMuon[myMu2], pzMuon[myMu2], energyMuon[myMu2]);
    m_p4Lepton3[mm] -> SetXYZT(0.0, 0.0, 0.0, 0.0);
    
    m_p4Lepton1Energy[mm] = energyMuon[myMu1];
    m_p4Lepton2Energy[mm] = energyMuon[myMu2];
    m_p4Lepton3Energy[mm] = -999.9;
    
    m_p4Lepton1Type[mm] = 1;
    m_p4Lepton2Type[mm] = 1;
    m_p4Lepton3Type[mm] = 1;
    
    m_mll            [mm] = mll;
    m_deltaPhi       [mm] = -999.9;
    m_deltaErre      [mm] = -999.9;
    m_deltaEtaLeptons[mm] = -999.9;
    
    m_dilepPt [mm].SetXYZ( m_p4Lepton1[mmm]->Vect().X()+m_p4Lepton2[mmm]->Vect().X(), 
			   m_p4Lepton2[mmm]->Vect().Y()+m_p4Lepton2[mmm]->Vect().Y(), 0.0 );
    m_trilepPt[mm].SetXYZ( 0.0, 0.0, 0.0 );

    // usual definition
    m_transvMass[mm] = -999.9;
    // chris' variable
    m_GammaMR   [mm] = -999.9;
    
    m_MR            [mm] = -999.9;
    m_MTR           [mm] = -999.9;
    m_metOptll      [mm] = -999.9;
    m_mT2           [mm] = -999.9;
    m_p3TKMET       [mm] = -999.9;
    m_chMet         [mm] = -999.9;
    m_projectedMet  [mm] = -999.9;
    m_projectedPFMet[mm] = -999.9;
    m_projectedTkMet[mm] = -999.9;
    m_MTRcharged    [mm] = -999.9;

    m_ch[mm][0] = chargeMuon[myMu1];
    m_ch[mm][1] = chargeMuon[myMu2];
    m_ch[mm][2] = -999.9;
    m_iso[mm][0] = pfCombinedIsoMuon[myMu1] / Muon1Pt;
    m_iso[mm][1] = pfCombinedIsoMuon[myMu2] / Muon2Pt;
    m_iso[mm][2] = -999.9;
    m_lh[mm][0] = -999.9;
    m_lh[mm][1] = -999.9;    
    m_lh[mm][2] = -999.9;    
    /*
    m_bdt[mm][0] = -999.9;
    m_bdt[mm][1] = -999.9;
    m_bdt[mm][2] = -999.9;
    */
    m_bdt[mm][0] = leptBDTForBs(fMVAMuon,myMu1,false);
    m_bdt[mm][1] = leptBDTForBs(fMVAMuon,myMu2,false);
    m_bdt[mm][2] = -999.9;

    m_chmaj[mm][0] = -999.9;
    m_chmaj[mm][1] = -999.9;    
    m_chmaj[mm][2] = -999.9;

    Muon1.Clear();
    Muon2.Clear();

  }

}// setKinematicsMM()

void tHSelection::resetKinematicsStart() {


  // for tH we don't care about this, to be checked and removed
  theElectron  = -1;
  thePositron  = -1;
  theMuonMinus = -1;
  theMuonPlus  = -1;

  thePreElectron  = -1;
  thePrePositron  = -1;
  thePreMuonMinus = -1;
  thePreMuonPlus  = -1;

}

void tHSelection::resetKinematics() {

  m_p3PFMET->SetXYZ(0,0,0);

  selectedLeptons.clear();
  selectedMuons.clear();
  selectedElectrons.clear();

  m_goodJets.clear();
  m_goodcbIDJets.clear();

  _genHiggsPt = -999.9;
  _genHiggsEta = -999.9;

  _genTopPt = -999.9;
  _genTopEta = -999.9;

  _genWpfromH_Pt = -999.9;
  _genWpfromH_Eta = -999.9;

  _genWmfromH_Pt = -999.9;
  _genWmfromH_Eta = -999.9;

  _genWfromT_Pt = -999.9;
  _genWfromT_Eta = -999.9;

  _genLeptonPlusfromWfromH_Pt = -999.9;
  _genLeptonPlusfromWfromH_Eta = -999.9;

  _genLeptonMinusfromWfromH_Pt = -999.9;
  _genLeptonMinusfromWfromH_Eta = -999.9;

  _genLeptonfromWfromT_Pt = -999.9;
  _genLeptonfromWfromT_Eta = -999.9;

  _genNeutrinoPlusfromWfromH_Pt = -999.9;
  _genNeutrinoPlusfromWfromH_Eta = -999.9;

  _genNeutrinoMinusfromWfromH_Pt = -999.9;
  _genNeutrinoMinusfromWfromH_Eta = -999.9;
  
  _genNeutrinofromWfromT_Pt = -999.9;
  _genNeutrinofromWfromT_Eta = -999.9;
  
  _genForwardQuark_Pt = -999.9;
  _genForwardQuark_Eta = -999.9;
  
  for(int theChannel=0; theChannel<6; theChannel++) {

    eleCands[theChannel].clear();
    muCands[theChannel].clear();

    m_p4Lepton1[theChannel]->SetXYZT(0,0,0,0);                                                        
    m_p4Lepton2[theChannel]->SetXYZT(0,0,0,0);
    m_p4Lepton3[theChannel]->SetXYZT(0,0,0,0);

    hardestLeptonPt[theChannel]   = 0.;
    slowestLeptonPt[theChannel]   = 0.;
    m_mll[theChannel]             = 0.;
    m_deltaPhi[theChannel]        = 0.;
    m_deltaErre[theChannel]       = 0.;
    m_deltaEtaLeptons[theChannel] = 0.; 
    m_dilepPt[theChannel]         = 0.;
    m_transvMass[theChannel]      = 0.;
    m_metOptll[theChannel]        = 0.;
    m_mT2[theChannel]             = 0.;
    m_projectedMet[theChannel]    = 0.;
    m_chMet[theChannel]           = 0.;

    hardestLeptonPt[theChannel] = -999.9; 
    slowestLeptonPt[theChannel] = -999.9;
    leadJetBtag    [theChannel] = -999.9; 
    subleadJetBtag [theChannel] = -999.9; 
    subLeadJetsMaxBtag[theChannel] = 999.9;
    leadJetBtagBProb[theChannel] = -999.9; 
    subleadJetBtagBProb[theChannel] = -999.9;
    subLeadJetsMaxBtagBProb[theChannel] = -999,9;
    m_softbdisc[theChannel] = -999.9;
    m_hardbdisc[theChannel] = -999.9;
    njets[theChannel] = 0; 
    ncbIDjets[theChannel] = 0;
    nuncorrjets[theChannel] = 0;

  }
}




int tHSelection::numJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  //std::cout << "calling to numJet en el channel = " << theChannel << std::endl;

  m_goodJets.clear();

  int num = 0;
  float ETMax  = 0.;
  float ETMax2 = 0.;

  nsoftjets                    [theChannel] = 0;
  m_numbtagCSVMmvaIDcentraljets[theChannel] = 0;//MVA ID -- btag for central jets
  m_numbtagCSVLmvaIDcentraljets[theChannel] = 0;
  m_numbtagCSVTmvaIDcentraljets[theChannel] = 0;
  m_nummvaIDforwardjets        [theChannel] = 0;//MVA ID for forward jets
  
  // Initialize jet index
  theLeadingJet[theChannel] = -1; // Jet with highest pT
  theSecondJet [theChannel] = -1; // Second jet with highest pT  

  m_jetsSum[theChannel]->SetXYZT(0.,0.,0.,0);

  TString JESUncertainty(_selectionEEE->getStringParameter("JESUncertainty"));

  //std::cout << "!!! empiezo loop de jets" << std::endl;

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    //std::cout<<"JET j= "<<j<<std::endl;

    TVector3       p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    float pt    = GetPt(      pxAK5PFPUcorrJet[j],      pyAK5PFPUcorrJet[j]);
    float rawpt = GetPt(uncorrpxAK5PFPUcorrJet[j],uncorrpyAK5PFPUcorrJet[j]);

    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) 
      pt = ( GetJESCorrected(p4Jet,JESUncertainty.Data()) ).Pt();

    // PF jet ID variables
    float neutralHadFrac     = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    float neutralEmFraction  = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int nConstituents        = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] 
      + photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] 
      + muonMultiplicityAK5PFPUcorrJet[j] + HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int chargedMultiplicity  = chargedHadronMultiplicityAK5PFPUcorrJet[j] 
      + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction  = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    
    //if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
    //chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    // Loose MvaId
    // if (!isLooseJetMva(pt,etaAK5PFPUcorrJet[j],jetIdMvaPhilV1AK5PFPUcorrJet[j])) continue;

    // Boolean for the matching with the selected leptons on the event
    bool foundMatch = false;

    // Check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {

      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        H_deltaRcorr -> Fill(deltaR);
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Ele.Clear();

      }
    }
    if(foundMatch) continue;

    // Check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {

      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Muon.Clear();

      }
    }
    if(foundMatch) continue;
    

    if( pt > 5.0 ) (*m_jetsSum[theChannel]) += p4Jet;
    
    // Jet eta acceptance cut
    if(_selectionEEE->getSwitch("etaJetAcc") && !_selectionEEE->passCut("etaJetAcc",fabs(etaAK5PFPUcorrJet[j]))) continue;

    if( fabs(etaAK5PFPUcorrJet[j]==5) ) continue; 
    // Remove the place with wrong JECs 
    // (https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1259/1.html)
    
    // Classify the jet as lead or sublead, if it is the case
    if ( pt>ETMax2 && pt>ETMax ) {
      
      theSecondJet[theChannel] = theLeadingJet[theChannel];
      ETMax2 = ETMax;
      
      theLeadingJet   [theChannel] = j;
      leadJetBtag     [theChannel] = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];
      leadJetBtagBProb[theChannel] = jetBProbabilityBJetTagsAK5PFPUcorrJet[j];
      ETMax = pt;

    } else if ( pt>ETMax2 && pt<ETMax ) {

      theSecondJet       [theChannel] = j;
      subleadJetBtag     [theChannel] = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];
      subleadJetBtagBProb[theChannel] = jetBProbabilityBJetTagsAK5PFPUcorrJet[j];
      ETMax2 = pt;

    }
  
    if( pt > 10. ) nsoftjets[theChannel]++;

    // Jet Et (pT) acceptance cut
    if(_selectionEEE->getSwitch("etJetAcc") && !_selectionEEE->passCut("etJetAcc", pt)) continue;

    // b-tagging algorithms
    float CHE = trackCountingHighEffBJetTagsAK5PFPUcorrJet   [j];     
    float JBP = jetBProbabilityBJetTagsAK5PFPUcorrJet        [j];
    float CVS = combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[j];
    
    // CVS Working points, as in BTV Twiki
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP    

    if      ( fabs(etaAK5PFPUcorrJet[j]) < 2.4 ){ // Central Jets

      if     (CVS > 0.898) m_numbtagCSVTmvaIDcentraljets[theChannel]++;
      else if(CVS > 0.679) m_numbtagCSVMmvaIDcentraljets[theChannel]++;
      else if(CVS > 0.244) m_numbtagCSVLmvaIDcentraljets[theChannel]++;

    } else {
      m_nummvaIDforwardjets[theChannel]++;
    }

    if ( fabs(etaAK5PFPUcorrJet[j]) >= 2.4 ) continue; // for the numJets, only central

    m_goodJets.push_back(j);
    num++;

    p3Jet.Clear();
    p4Jet.Clear();
    
  }

  // Extra info for jet Id studies
  if( wantJetIdStuff ) {
    int theFirst  = theLeadingJet[theChannel];
    int theSecond = theSecondJet[theChannel];
    
    // Leading jet
    if (theFirst>-1) {
      float firstNeutralHadFrac    = neutralHadronEnergyAK5PFPUcorrJet[theFirst]/uncorrenergyAK5PFPUcorrJet[theFirst];
      float firstNeutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[theFirst]/uncorrenergyAK5PFPUcorrJet[theFirst];

      int firstNConstituents       = chargedHadronMultiplicityAK5PFPUcorrJet[theFirst] 
	+ neutralHadronMultiplicityAK5PFPUcorrJet[theFirst] + photonMultiplicityAK5PFPUcorrJet[theFirst] 
	+ electronMultiplicityAK5PFPUcorrJet[theFirst] + muonMultiplicityAK5PFPUcorrJet[theFirst] 
	+ HFHadronMultiplicityAK5PFPUcorrJet[theFirst] + HFEMMultiplicityAK5PFPUcorrJet[theFirst];

      float firstChargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[theFirst]/uncorrenergyAK5PFPUcorrJet[theFirst];

      int firstChargedMultiplicity  = chargedHadronMultiplicityAK5PFPUcorrJet[theFirst] 
	+ electronMultiplicityAK5PFPUcorrJet[theFirst] + muonMultiplicityAK5PFPUcorrJet[theFirst];

      float firstChargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[theFirst]/uncorrenergyAK5PFPUcorrJet[theFirst];

      float firstPt = sqrt( pxAK5PFPUcorrJet[theFirst]*pxAK5PFPUcorrJet[theFirst] +
			    pyAK5PFPUcorrJet[theFirst]*pyAK5PFPUcorrJet[theFirst] );
      
      leadJetPt       [theChannel] = firstPt;
      leadJetEta      [theChannel] = etaAK5PFPUcorrJet[theFirst];
      leadJetLoosePFId[theChannel] = isPFJetID(fabs(etaAK5PFPUcorrJet[theFirst]),firstNeutralHadFrac,
					       firstNeutralEmFraction,firstNConstituents,firstChargedHadFraction,
					       firstChargedMultiplicity,firstChargedEmFraction, TopHiggs::loose);

      leadJetMvaJetId [theChannel] = -999.9; //jetIdMvaPhilV1AK5PFPUcorrJet[theFirst];
      leadJetLooseId  [theChannel] = isLooseJetMva(leadJetPt[theChannel],etaAK5PFPUcorrJet[theFirst],
						   -999.9); // <-- dummy value -- jetIdMvaPhilV1AK5PFPUcorrJet[theFirst]);
    } else {

      leadJetPt       [theChannel] = -999.;
      leadJetEta      [theChannel] = -999.;
      leadJetLoosePFId[theChannel] = -999.;
      leadJetMvaJetId [theChannel] = -999.;
      leadJetLooseId  [theChannel] = -999.;

    }
    
    // Subleading jet
    if (theSecond>-1) {
      float secondNeutralHadFrac    = neutralHadronEnergyAK5PFPUcorrJet[theSecond]/uncorrenergyAK5PFPUcorrJet[theSecond];
      float secondNeutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[theSecond]/uncorrenergyAK5PFPUcorrJet[theSecond];

      int secondNConstituents       = chargedHadronMultiplicityAK5PFPUcorrJet[theSecond] 
	+ neutralHadronMultiplicityAK5PFPUcorrJet[theSecond] + photonMultiplicityAK5PFPUcorrJet[theSecond] 
	+ electronMultiplicityAK5PFPUcorrJet[theSecond] + muonMultiplicityAK5PFPUcorrJet[theSecond] 
	+ HFHadronMultiplicityAK5PFPUcorrJet[theSecond] + HFEMMultiplicityAK5PFPUcorrJet[theSecond];

      float secondChargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[theSecond]/uncorrenergyAK5PFPUcorrJet[theSecond];
      int secondChargedMultiplicity  = chargedHadronMultiplicityAK5PFPUcorrJet[theSecond] 
	+ electronMultiplicityAK5PFPUcorrJet[theSecond] + muonMultiplicityAK5PFPUcorrJet[theSecond];
      float secondChargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[theSecond]/uncorrenergyAK5PFPUcorrJet[theSecond];
      float secondPt = sqrt( pxAK5PFPUcorrJet[theSecond]*pxAK5PFPUcorrJet[theSecond] 
			     + pyAK5PFPUcorrJet[theSecond]*pyAK5PFPUcorrJet[theSecond] );

      subleadJetPt       [theChannel] = secondPt;
      subleadJetEta      [theChannel] = etaAK5PFPUcorrJet[theSecond];
      subleadJetLoosePFId[theChannel] = isPFJetID(fabs(etaAK5PFPUcorrJet[theSecond]),secondNeutralHadFrac,
						  secondNeutralEmFraction,secondNConstituents,secondChargedHadFraction,
						  secondChargedMultiplicity,secondChargedEmFraction,TopHiggs::loose);

      subleadJetMvaJetId [theChannel] = -999.9; //jetIdMvaPhilV1AK5PFPUcorrJet[theSecond];
      subleadJetLooseId  [theChannel] = isLooseJetMva(subleadJetPt[theChannel],etaAK5PFPUcorrJet[theSecond],
						      -999.9); // <-dummy value--jetIdMvaPhilV1AK5PFPUcorrJet[theSecond]);
    } else {

      subleadJetPt       [theChannel] = -999.;
      subleadJetEta      [theChannel] = -999.;
      subleadJetLoosePFId[theChannel] = -999.;
      subleadJetMvaJetId [theChannel] = -999.;
      subleadJetLooseId  [theChannel] = -999.;
    }
    
    // Match with gen jets
    int firstAss  = -999;
    int secondAss = -999;

    float firstDRmin  = 999.;
    float secondDRmin = 999.;

    if (theFirst>-1) {

      for (int iGen=0; iGen<nAK5GenJet; iGen++) {
	TVector3 t3GenJet(pxAK5GenJet[iGen],pyAK5GenJet[iGen],pzAK5GenJet[iGen]);
	TVector3 t3FirstRecoJet(pxAK5PFPUcorrJet[theFirst],pyAK5PFPUcorrJet[theFirst],pzAK5PFPUcorrJet[theFirst]);
	float genPt   = t3GenJet.Perp();
	float firstPt = sqrt( pxAK5PFPUcorrJet[theFirst]*pxAK5PFPUcorrJet[theFirst] 
			      + pyAK5PFPUcorrJet[theFirst]*pyAK5PFPUcorrJet[theFirst] );
	float firstDR = t3FirstRecoJet.DeltaR(t3GenJet);
	double firstExpres  = ErrEt(firstPt,t3FirstRecoJet.Eta());
	if ( (firstDR<firstDRmin) && ((fabs(firstPt-genPt))/genPt)<0.5 ) {
	  firstAss = iGen;
	  firstDRmin = firstDR;
	}

	t3GenJet.Clear();
	t3FirstRecoJet.Clear();

      }
      if (firstAss>-999) {
	float firstGenAssPt = sqrt(pxAK5GenJet[firstAss]*pxAK5GenJet[firstAss] 
				   + pyAK5GenJet[firstAss]*pyAK5GenJet[firstAss]);
	if (firstDRmin > 0.1 + 0.3 * exp(-0.05*(firstGenAssPt-10))) firstAss = -999;
      }
      if (firstAss>-1) leadJetMatchGen[theChannel] = 1;  
      if (firstAss<0)  leadJetMatchGen[theChannel] = 0;  

    } else {

      leadJetMatchGen[theChannel] = -999;

    }
    
    if (theSecond>-1) {

      for (int iGen=0; iGen<nAK5GenJet; iGen++) {
	TVector3 t3GenJet(pxAK5GenJet[iGen],pyAK5GenJet[iGen],pzAK5GenJet[iGen]);
	TVector3 t3SecondRecoJet(pxAK5PFPUcorrJet[theSecond],pyAK5PFPUcorrJet[theSecond],pzAK5PFPUcorrJet[theSecond]);
	float genPt = t3GenJet.Perp();
	float secondPt = sqrt( pxAK5PFPUcorrJet[theSecond]*pxAK5PFPUcorrJet[theSecond] 
			       + pyAK5PFPUcorrJet[theSecond]*pyAK5PFPUcorrJet[theSecond] );
	float secondDR = t3SecondRecoJet.DeltaR(t3GenJet);
	double secondExpres = ErrEt(secondPt,t3SecondRecoJet.Eta());
	if ( (secondDR<secondDRmin) && ((fabs(secondPt-genPt))/genPt)<0.5 ) {
	  secondAss = iGen;
	  secondDRmin = secondDR;
	}
	
	t3GenJet.Clear();
	t3SecondRecoJet.Clear();

      }
      if (secondAss>-999) {
	float secondGenAssPt = sqrt(pxAK5GenJet[secondAss]*pxAK5GenJet[secondAss] 
				    + pyAK5GenJet[secondAss]*pyAK5GenJet[secondAss]);
	if (secondDRmin > 0.1 + 0.3 * exp(-0.05*(secondGenAssPt-10))) secondAss = -999;
      }
      if (secondAss>-1) subleadJetMatchGen[theChannel] = 1;  
      if (secondAss<0)  subleadJetMatchGen[theChannel] = 0;  

    } else {

      subleadJetMatchGen[theChannel] =  -999;

    }
  }
 
  return num;
}

int tHSelection::numcbIDJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel) {

  int num=0;
  m_goodcbIDJets.clear();

  nsoftjetscbID               [theChannel]=0;
  m_numbtagCSVMcbIDcentraljets[theChannel] = 0;//MVA ID for jets
  m_numbtagCSVLcbIDcentraljets[theChannel] = 0;
  m_numbtagCSVTcbIDcentraljets[theChannel] = 0;
  m_numcbIDforwardjets        [theChannel] = 0;//MVA ID for jets

  m_jetscbIDSum[theChannel]->SetXYZT(0.,0.,0.,0);

  TString JESUncertainty(_selectionEEE->getStringParameter("JESUncertainty"));

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3       p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    float pt    = GetPt(      pxAK5PFPUcorrJet[j],      pyAK5PFPUcorrJet[j]);
    float rawpt = GetPt(uncorrpxAK5PFPUcorrJet[j],uncorrpyAK5PFPUcorrJet[j]);

    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) 
      pt = ( GetJESCorrected(p4Jet,JESUncertainty.Data()) ).Pt();

    float beta   = betastarclassicIdMvaAK5PFPUcorrJet[j];  // To be changed
    float betaTh = 0.64; // To be changed
    float rms    = dR2MeanIdMvaAK5PFPUcorrJet[j];  // To be changed
    float eta    = etaAK5PFPUcorrJet[j];

    bool passID = true;
    if( abs(eta) < 2.5 ){
      //if ( beta < 0 ) passID = false;
      if ( beta > 0.2*TMath::Log( m_goodvertices - betaTh ) ) passID = false;
      if ( rms > 0.06 ) passID = false;
    } else if ( abs(eta) < 3.0 ){
      if ( rms > 0.05 ) passID = false;
    }else{
      if ( rms > 0.055 ) passID = false;
    }

    if ( !passID ) continue;

    bool foundMatch = false;

    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {

      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Ele.Clear();

      }
    }
    if(foundMatch) continue;

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {

      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Muon.Clear();

      }
    }
    if(foundMatch) continue;
    
    if(pt>5.0) (*m_jetscbIDSum[theChannel]) += p4Jet;
    
    if(_selectionEEE->getSwitch("etaJetAcc") && !_selectionEEE->passCut("etaJetAcc",fabs(etaAK5PFPUcorrJet[j]))) continue;
    if(fabs(etaAK5PFPUcorrJet[j]==5)) continue; 
    // this is to remove the place with wrong JECs (https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1259/1.html)
    
    if(pt>10.) nsoftjetscbID[theChannel]++;
    if(_selectionEEE->getSwitch("etJetAcc") && !_selectionEEE->passCut("etJetAcc", pt)) continue;

    // th
    float CVS = combinedSecondaryVertexBJetTagsAK5PFPUcorrJet[j];

    // CHANGE CUTS for b-tag CVS
    
    if      ( fabs(etaAK5PFPUcorrJet[j]) < 2.4 ){ // Central Jets

      if     (CVS > 0.898) m_numbtagCSVTcbIDcentraljets[theChannel]++;
      else if(CVS > 0.679) m_numbtagCSVMcbIDcentraljets[theChannel]++;
      else if(CVS > 0.244) m_numbtagCSVLcbIDcentraljets[theChannel]++;

    }else{
      m_numcbIDforwardjets[theChannel]++;
    }

    if (fabs(etaAK5PFPUcorrJet[j]) >= 2.4 ) continue; // for the numJets, only central

    m_goodcbIDJets.push_back(j);
    num++;

    p3Jet.Clear();
    p4Jet.Clear();
    
  }

  return num;

}


int tHSelection::numUncorrJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel ) {

  int num=0;

  TString JESUncertainty(_selectionEEE->getStringParameter("JESUncertainty"));

  m_uncorrJetsSum[theChannel]->SetXYZT(0.,0.,0.,0.);

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    float uncorrEt = uncorrenergyAK5PFPUcorrJet[j]*fabs(sin(thetaAK5PFPUcorrJet[j]));
    TLorentzVector p4Jet;
    p4Jet.SetPtEtaPhiE(uncorrEt,etaAK5PFPUcorrJet[j],phiAK5PFPUcorrJet[j],uncorrenergyAK5PFPUcorrJet[j]);
    TVector3 p3Jet = p4Jet.Vect();

    TLorentzVector p4JESJet(p3Jet, uncorrenergyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) uncorrEt = (GetJESCorrected(p4JESJet,JESUncertainty.Data())).Pt();

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] 
      + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    
    //if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
    //chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    bool foundMatch=false;
    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {

      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Ele.Clear();
	
      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {

      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;

	p3Muon.Clear();

      }
    }
    if(foundMatch) continue;

    if(uncorrEt>5.0) (*m_uncorrJetsSum[theChannel]) += p4Jet;

    if(_selectionEEE->getSwitch("etaJetAcc")     &&!_selectionEEE->passCut("etaJetAcc",fabs(etaAK5PFPUcorrJet[j])))continue;
    if(_selectionEEE->getSwitch("etUncorrJetAcc")&&!_selectionEEE->passCut("etUncorrJetAcc", uncorrEt))   continue;
    
    num++;

    p3Jet.Clear();
    p4Jet.Clear();
    p4JESJet.Clear();

  }

  return num;
}

float tHSelection::bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel ) {

  TString JESUncertainty(_selectionEEE->getStringParameter("JESUncertainty"));

  float maxTCHE=-999;
  float maxJetBProb=-999;
  float outputSubLeadJets = -999;
  float outputSubLeadJetsBProb = -999;
  m_numbtagjets[theChannel]=0;

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    // no threshold is applied here on pt. Not affected by JES uncertainties
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);

    float pt = GetPt(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j]);
    if(JESUncertainty == TString("Up") || JESUncertainty == TString("Down")) 
      pt = (GetJESCorrected(p4Jet,JESUncertainty.Data())).Pt();

    if(_selectionEEE->getSwitch("etaJetAcc") && !_selectionEEE->passCut("etaJetAcc",fabs(etaAK5PFPUcorrJet[j]))) continue;
    
    // if(theSecondJet[theChannel]>-1 && fabs(etaAK5PFPUcorrJet[j])>4.5) continue;

    // hardcoded
    float rawpt = GetPt(uncorrpxAK5PFPUcorrJet[j],uncorrpyAK5PFPUcorrJet[j]);
    if(rawpt < 10.0) continue;

    //    if(weightedDz1AK5PFPUcorrJet[j] >= 2) continue;

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] 
      + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    
    //if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
    //chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;

    // Loose MvaId 
    //if (!isLooseJetMva(pt,etaAK5PFPUcorrJet[j],jetIdMvaPhilV1AK5PFPUcorrJet[j])) continue;

    bool foundMatch=false;
    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {

      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
	
	p3Ele.Clear();

      }
    }

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {

      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
	
	p3Muon.Clear();

      }
    }
    if(foundMatch) continue;

    float tmpTCHE = trackCountingHighEffBJetTagsAK5PFPUcorrJet[j];     
    float tmpJBP  = jetBProbabilityBJetTagsAK5PFPUcorrJet[j];

    if(tmpTCHE > maxTCHE) maxTCHE = tmpTCHE; 
    //     if(pt<=30 && tmpTCHE > maxTCHE) maxTCHE = tmpTCHE;
    //     if(pt>30 && tmpJBP > maxJetBProb) maxJetBProb = tmpJBP;

    if(j != theLeadingJet[theChannel] && tmpTCHE > outputSubLeadJets) outputSubLeadJets = tmpTCHE;
    if(j != theLeadingJet[theChannel] && tmpJBP > outputSubLeadJetsBProb) outputSubLeadJetsBProb = tmpJBP;

    if(tmpTCHE>1.6) m_numbtagjets[theChannel]++;

    p3Jet.Clear();
    p4Jet.Clear();
    
  }

  subLeadJetsMaxBtag[theChannel] = outputSubLeadJets;
  subLeadJetsMaxBtagBProb[theChannel] = outputSubLeadJetsBProb;

  int lj = theLeadingJet[theChannel];
  float ptleadjet = GetPt(pxAK5PFPUcorrJet[lj],pyAK5PFPUcorrJet[lj]);

  // hardcode the cuts
  float bitval = 1;
  if(maxTCHE>=2.1) bitval=0.;
  // give data-MC discrepancies, revert to 2011 TCHE only
  // if(maxJetBProb>=1.05) bitval=0.;
  
  m_softbdisc[theChannel]=maxTCHE;
  m_hardbdisc[theChannel]=maxJetBProb;

  return bitval;

}

float tHSelection::deltaPhiLLJet(int ichan) {   
  
  int myLeadingJet = theLeadingJet[ichan];
  int mySecondJet  = theSecondJet[ichan];

  // if the event has >=2 jet, the cut is on dphi(ll-jj). But the acceptance on jets has to be restricted from 5.0 to 4.5
  if( mySecondJet > -1 && myLeadingJet > -1 && m_dilepPt[ichan].Pt()>0 
      && fabs(etaAK5PFPUcorrJet[myLeadingJet])<4.5 && fabs(etaAK5PFPUcorrJet[mySecondJet])<4.5 ) {
    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);
    TVector3 subleadingJetP3(pxAK5PFPUcorrJet[mySecondJet],pyAK5PFPUcorrJet[mySecondJet],pzAK5PFPUcorrJet[mySecondJet]);
    TVector3 jjP3 = leadingJetP3+subleadingJetP3;
    
    float valueDeltaPhiLLJet = fabs(180./TMath::Pi() * jjP3.DeltaPhi(m_dilepPt[ichan]));

    leadingJetP3.Clear();
    subleadingJetP3.Clear();
    jjP3.Clear();

    return valueDeltaPhiLLJet;

    //return fabs(180./TMath::Pi() * jjP3.DeltaPhi(m_dilepPt[ichan]));
  }


  // if the event has <=1 jet, the cut is on dphi(ll-j)
  if(myLeadingJet > -1 && m_dilepPt[ichan].Pt()>0) {

    TVector3 leadingJetP3(pxAK5PFPUcorrJet[myLeadingJet],pyAK5PFPUcorrJet[myLeadingJet],pzAK5PFPUcorrJet[myLeadingJet]);  
    if(leadingJetP3.Pt()>15.0){
      
      float valueDeltaPhiLLJet = fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));

      leadingJetP3.Clear();

      return valueDeltaPhiLLJet;

      //return fabs(180./TMath::Pi() * leadingJetP3.DeltaPhi(m_dilepPt[ichan]));                  
    
    }else return 0.1;
  } else return 0.1;
}

int tHSelection::numSoftMuons(std::vector<int> muonToRemove, std::vector<int> jetsToRemove) {

  int num = 0;
  for(int i=0; i<nMuon; ++i) {

    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) { 
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;

    bool isAroundGoodJet=false;
    for(int jet=0; jet<(int)jetsToRemove.size(); jet++) {
      int j = jetsToRemove[jet];
      TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
      TVector3 p3Muon(pxMuon[i],pyMuon[i],pzMuon[i]);
      if(p3Jet.DeltaR(p3Muon)<0.5) isAroundGoodJet=true;
    }
    if(isAroundGoodJet) continue;

    float pt = GetPt(pxMuon[i],pyMuon[i]);
    if(pt < 3.0) continue;

    Utils anaUtils;
    if(!anaUtils.muonIdVal(muonIdMuon[i],AllTrackerMuons)||!anaUtils.muonIdVal(muonIdMuon[i],TMLastStationTight))continue;
       
    int track = trackIndexMuon[i];
    if(trackValidHitsTrack[track]<=10) continue;

    float dxyMuon= transvImpactParTrack[track];
    float dzMuon = fabs(muonDzPV(i,0));
    if(dxyMuon > 0.200) continue;     // hardcoded  
    if(dzMuon  > 0.200) continue;     // hardcoded  

    float isoSumRel = (sumPt03Muon[i] + emEt03Muon[i] + hadEt03Muon[i]) / pt;
    //    float isoSumRel = pfCombinedIsoMuon[i] / pt;
    if(pt>20 && isoSumRel<0.1) continue;  
    
    num++;
  }
  return num;
}

int tHSelection::numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove  ) {

  int numEle = 0;
  for(int i=0; i<nEle; ++i) {
    
    bool isSelEle=false;
    for(int eleSel=0; eleSel<(int)eleToRemove.size(); eleSel++) {
      if(i==eleToRemove[eleSel]) isSelEle=true;
    }
    if(isSelEle) continue;

    if(_selectionEEE->getSwitch("etaElectronAcc")&&!_selectionEEE->passCut("etaElectronAcc",etaEle[i]) ) continue;
    if(_selectionEEE->getSwitch("ptElectronAcc")&&!_selectionEEE->passCut("ptElectronAcc",GetPt(pxEle[i],pyEle[i])))continue;
    bool theId, theIso, theConvRej;
    theId = theIso = theConvRej = true;
    isEleID2012AndDenom(i,&theId,&theIso,&theConvRej);
    float pt = GetPt(pxEle[i],pyEle[i]);	

    if(!theId || !theIso || !theConvRej) continue;
    
    int track = gsfTrackIndexEle[i];
    float dxyEle = transvImpactParGsfTrack[track];
    float dzEle  = eleDzPV(i,0);
    if (_selectionEEE->getSwitch("electronIP") && (!_selectionEEE->passCut("electronIP",dxyEle)) ) continue;
    if (_selectionEEE->getSwitch("electronDz") && (!_selectionEEE->passCut("electronDz",dzEle)) ) continue;

    numEle++;
  }

  int numMu = 0;
  for(int i=0; i<nMuon; ++i) {
    
    bool isSelMuon=false;
    for(int muSel=0; muSel<(int)muonToRemove.size(); muSel++) {
      if(i==muonToRemove[muSel]) isSelMuon=true;
    }
    if(isSelMuon) continue;
    
    float ptMu = GetPt(pxMuon[i],pyMuon[i]);
    if(_selectionEEE->getSwitch("etaMuonAcc") && !_selectionEEE->passCut("etaMuonAcc",etaMuon[i]) ) continue;
    if(_selectionEEE->getSwitch("ptMuonAcc") && !_selectionEEE->passCut("ptMuonAcc",ptMu) ) continue;

    bool theId = true;
    isMuonID2012(i,&theId);
    if(!theId) continue;
    if( ! isPFIsolatedMuon2012(i) ) continue; 

    int track = trackIndexMuon[i];
    float dxy = transvImpactParTrack[track];
    float dz  = muonDzPV(i,0);  

    if (ptMu>20) {   // hardcoded
      if (_selectionEEE->getSwitch("muonIPhighPT") && (!_selectionEEE->passCut("muonIPhighPT",dxy)) ) continue;   
    } 
    if (ptMu<20) {   // hardcoded
      if (_selectionEEE->getSwitch("muonIPlowPT")  && (!_selectionEEE->passCut("muonIPlowPT",dxy)) ) continue;   
    }
    if (_selectionEEE->getSwitch("muonDz") && (!_selectionEEE->passCut("muonDz",dz)) )  continue;   

    numMu++;
  }
  
  return numEle + numMu;
}

void tHSelection::setLepIdVariables(int index1, int index2, int index3, int pdgid1, int pdgid2, int pdgid3) {

  Utils anaUtils;

  int selectedLeptons[3];
  selectedLeptons[0] = index1;
  selectedLeptons[1] = index2; 
  selectedLeptons[2] = index3; 
  int pdgids[3];
  pdgids[0] = pdgid1;
  pdgids[1] = pdgid2;
  pdgids[2] = pdgid3;

  for(int i = 0; i < 3; i++) {
    if(selectedLeptons[i] > -1) {

      int lepIndex = selectedLeptons[i];

      if(pdgids[i]==11) {

        myPt     [i] = GetPt(pxEle[lepIndex],pyEle[lepIndex]);
        myEta    [i] = etaEle[lepIndex];
        myPhi    [i] = phiEle[lepIndex];
	myLepFlav[i] = 11;
        myLepId  [i] = mvaidtrigEle[lepIndex];
        myLepIso [i] = corrEleIso2012(lepIndex);

        int gsf = gsfTrackIndexEle[lepIndex];  
        myConv   [i] = (!hasMatchedConversionEle[lepIndex] && expInnerLayersGsfTrack[gsf]==0);

      } else if(pdgids[i]==13) {

        myPt     [i] = GetPt(pxMuon[lepIndex],pyMuon[lepIndex]);
        myEta    [i] = etaMuon[lepIndex];
        myPhi    [i] = phiMuon[lepIndex];
        bool muid;
        isMuonID2012(lepIndex,&muid);
	myLepFlav[i] = 13;
        myLepId  [i] = (muid) ? 1. : 0.;
        myLepIso [i] = mvaisoMuon[lepIndex];      
        myConv   [i] = 1.;

      } else {
        cout << "MISTAKE! WRONG PDG TYPE!" << endl;
      }

    } else {

      myPt    [i] = -999.;
      myEta   [i] = -999.;
      myPhi   [i] = -999.;
      myLepId [i] = -999.;
      myLepIso[i] = -999.;
      myConv  [i] = -999.;

    }
  }
}

int tHSelection::getPV() {

  int hardestPV = -1;
  float sumPtMax = 0.0;
  for(int v=0; v<nPV; v++) {
    if(SumPtPV[v] > sumPtMax) {
      sumPtMax = SumPtPV[v];
      hardestPV = v;
    }
  }

  return hardestPV;

}

bool tHSelection::isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits) {
  TVector3 p3Track(pxTrack[iTrack],pyTrack[iTrack],pzTrack[iTrack]);
  double pt = p3Track.Pt();
  if(pt < ptMin) return false;
  if(pt > ptMax) return false;
  if(trackNormalizedChi2Track[iTrack] > chi2) return false; 
  if(fabs(p3Track.Eta()) > etaMax) return false;
  if(trackValidHitsTrack[iTrack] < nHits) return false;
  return true;
}

double tHSelection::mT(TVector3 plep, TVector3 pneu) {

  TVector3 pTlep(plep.X(),plep.Y(),0.0);
  TVector3 pTneu(pneu.X(),pneu.Y(),0.0);

  return sqrt(2 * (pTlep.Mag()*pTneu.Mag() - pTlep*pTneu));

}

double tHSelection::mT2(TVector3 plep1, TVector3 plep2, TVector3 ptmiss) {

  // need an external dependency: removed right now
  //   Mt2::TwoVector pTlep1(plep1.X(),plep1.Y());
  //   Mt2::TwoVector pTlep2(plep2.X(),plep2.Y());
  //   Mt2::TwoVector pTmiss(ptmiss.X(),ptmiss.Y());
  //   double invis_mass = 0.0;
  
  //   Mt2::SUSYPhys_Mt2_222_Calculator mt2Calculator;
  
  //   // Could tell the MT2 calculating object to be verbose, and print out
  //   // debug messages while it is thinking ... but we won't:
  
  //   mt2Calculator.setDebug(false);
  
  //   // Now we can actually calculate MT2:
  //   double mt2 = mt2Calculator.mt2_222( pTlep1, pTlep2, pTmiss, invis_mass);
  
  double mt2 = 0.0;
  return mt2;

}

float tHSelection::GetProjectedMet(TVector3 p1, TVector3 p2) {

  // calculate with PF met
  float projMET_pf = 0.0;
  float deltaPhi1_pf = fabs(p1.DeltaPhi(*m_p3PFMET));
  float deltaPhi2_pf = fabs(p2.DeltaPhi(*m_p3PFMET));
  float deltaphi_pf = TMath::Min(deltaPhi1_pf,deltaPhi2_pf);
  if(deltaphi_pf<TMath::Pi()/2.) projMET_pf = m_p3PFMET->Mag() * fabs(sin(deltaphi_pf));
  else projMET_pf = m_p3PFMET->Mag();

  // calculate with TKMET
  TVector3 p3tkMet = pfChargedMet(p1,p2);
  float projMET_tk = 0.0;
  float deltaPhi1_tk = fabs(p1.DeltaPhi(p3tkMet));
  float deltaPhi2_tk = fabs(p2.DeltaPhi(p3tkMet));
  float deltaphi_tk  = TMath::Min(deltaPhi1_tk,deltaPhi2_tk);
  if(deltaphi_tk<TMath::Pi()/2.) projMET_tk = p3tkMet.Mag() * fabs(sin(deltaphi_tk));
  else projMET_tk = p3tkMet.Mag();

  return TMath::Min(projMET_pf,projMET_tk);
}

// calculate with PF met
float tHSelection::GetProjectedPFMet(TVector3 p1, TVector3 p2) {  

  float projMET_pf = 0.0;
  float deltaPhi1_pf = fabs(p1.DeltaPhi(*m_p3PFMET));
  float deltaPhi2_pf = fabs(p2.DeltaPhi(*m_p3PFMET));
  float deltaphi_pf = TMath::Min(deltaPhi1_pf,deltaPhi2_pf);
  if(deltaphi_pf<TMath::Pi()/2.) projMET_pf = m_p3PFMET->Mag() * fabs(sin(deltaphi_pf));
  else projMET_pf = m_p3PFMET->Mag();

  return projMET_pf;
}

// calculate with TK met  
float tHSelection::GetProjectedTkMet(TVector3 p1, TVector3 p2) {  

  TVector3 p3tkMet = pfChargedMet(p1,p2);
  float projMET_tk = 0.0;
  float deltaPhi1_tk = fabs(p1.DeltaPhi(p3tkMet));
  float deltaPhi2_tk = fabs(p2.DeltaPhi(p3tkMet));
  float deltaphi_tk  = TMath::Min(deltaPhi1_tk,deltaPhi2_tk);
  if(deltaphi_tk<TMath::Pi()/2.) projMET_tk = p3tkMet.Mag() * fabs(sin(deltaphi_tk));
  else projMET_tk = p3tkMet.Mag();

  return projMET_tk;
}

/// specific for HWW that has multiple channels with different HLT requirements
bool tHSelection::reloadTriggerMask(int runN)
{
  std::vector<int> triggerMask;

  // load the triggers required for EEE
  for (std::vector<std::string>::const_iterator fIter=requiredTriggersEEE.begin();fIter!=requiredTriggersEEE.end();++fIter)
    {   
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          //if( !strcmp ((*fIter).c_str(), nameHLT->at(i).c_str() ) )
          // nameHLT[i] has ..._vXXX
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersEEE = triggerMask;
  
  // load the triggers NOT required for EEE                                                                           
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersEEE.begin();fIter!=notRequiredTriggersEEE.end();++fIter) {
    std::string pathName = getHLTPathForRun(runN,*fIter);
    for(unsigned int i=0; i<nameHLT->size(); i++) {
      if(nameHLT->at(i).find(pathName) != string::npos) {
        triggerMask.push_back( indexHLT[i] ) ;
        break;
      }
    }
  }
  m_notRequiredTriggersEEE = triggerMask;

  // load the triggers required for MMM
  triggerMask.clear();
  for (std::vector<std::string>::const_iterator fIter=requiredTriggersMMM.begin();fIter!=requiredTriggersMMM.end();++fIter)
    {   
      //      std::cout << "For MM required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersMMM = triggerMask;

  // load the triggers NOT required for MMM
  triggerMask.clear();
  for (std::vector<std::string>::const_iterator fIter=notRequiredTriggersMMM.begin();fIter!=notRequiredTriggersMMM.end();++fIter)
    {   
      //      std::cout << "For MM not required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_notRequiredTriggersMMM = triggerMask;

  // load the triggers required for EEM
  triggerMask.clear();
  for (std::vector<std::string>::const_iterator fIter=requiredTriggersEEM.begin();fIter!=requiredTriggersEEM.end();++fIter)
    {   
      //      std::cout << "For EM required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersEEM = triggerMask;

  // load the triggers NOT required for EEM
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersEEM.begin();fIter!=notRequiredTriggersEEM.end();++fIter)
    {   
      //      std::cout << "For EM not required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_notRequiredTriggersEEM = triggerMask;

  // load the triggers required for MME
  triggerMask.clear();
  for (std::vector<std::string>::const_iterator fIter=requiredTriggersMME.begin();fIter!=requiredTriggersMME.end();++fIter)
    {   
      //      std::cout << "For EM required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_requiredTriggersMME = triggerMask;

  // load the triggers NOT required for MME
  triggerMask.clear();
  for (std::vector< std::string >::const_iterator fIter=notRequiredTriggersMME.begin();fIter!=notRequiredTriggersMME.end();++fIter)
    {   
      //      std::cout << "For EM not required: " << *fIter << std::endl;
      std::string pathName = getHLTPathForRun(runN,*fIter);
      for(unsigned int i=0; i<nameHLT->size(); i++)
        {
          if(nameHLT->at(i).find(pathName) != string::npos)
            {
              triggerMask.push_back( indexHLT[i] ) ;
              break;
            }
        }
    }
  m_notRequiredTriggersMME = triggerMask;

}

bool tHSelection::hasPassedHLT(int channel) {

  Utils anaUtils;

  if(channel==eee) { 

    bool required    = anaUtils.getTriggersOR(m_requiredTriggersEEE   , firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersEEE, firedTrg);
    return (required && !notRequired);

  } else if(channel==mmm) {

    bool required    = anaUtils.getTriggersOR(m_requiredTriggersMMM   , firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersMMM, firedTrg);
    return (required && !notRequired);

  } else if(channel==eem) {

    bool required    = anaUtils.getTriggersOR(m_requiredTriggersEEM   , firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersEEM, firedTrg);
    return (required && !notRequired);

  } else if(channel==mme) {

    bool required    = anaUtils.getTriggersOR(m_requiredTriggersMME   , firedTrg);
    bool notRequired = anaUtils.getTriggersOR(m_notRequiredTriggersMME, firedTrg);
    return (required && !notRequired);

  }

  return true;

}

void tHSelection::setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {

  if     (channel==eee) requiredTriggersEEE=reqTriggers;
  else if(channel==mmm) requiredTriggersMMM=reqTriggers;
  else if(channel==eem) requiredTriggersEEM=reqTriggers;
  else if(channel==mme) requiredTriggersMME=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;

}

void tHSelection::setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel) {

  if     (channel==eee) notRequiredTriggersEEE=reqTriggers;
  else if(channel==mmm) notRequiredTriggersMMM=reqTriggers;
  else if(channel==eem) notRequiredTriggersEEM=reqTriggers;
  else if(channel==mme) notRequiredTriggersMME=reqTriggers;
  else std::cout << "WARNING: triggers are set for an unknown channel!" << std::endl;

}

void tHSelection::JESPfMet( std::vector<int> eleToRemove, std::vector<int> muonToRemove) {

  TLorentzVector jetSumNom, jetSumUp, jetSumDown;
  jetSumUp.SetXYZT(0.,0.,0.,0);
  jetSumDown.SetXYZT(0.,0.,0.,0);

  for(int j=0;j<nAK5PFPUcorrJet;j++) {

    TVector3 p3Jet(pxAK5PFPUcorrJet[j],pyAK5PFPUcorrJet[j],pzAK5PFPUcorrJet[j]);
    TLorentzVector p4Jet(p3Jet, energyAK5PFPUcorrJet[j]);
    float pt = p4Jet.Pt();

    if(_selectionEEE->getSwitch("etaJetAcc") && !_selectionEEE->passCut("etaJetAcc",fabs(etaAK5PFPUcorrJet[j]))) continue;


    TLorentzVector p4JetUp = GetJESCorrected(p4Jet,"Up");
    float ptUp = p4JetUp.Pt();
    float energyUp = p4JetUp.E();

    TLorentzVector p4JetDown = GetJESCorrected(p4Jet,"Down");
    float ptDown = p4JetDown.Pt();
    float energyDown = p4JetDown.E();

    TLorentzVector p4JetJesUp, p4JetJesDown;
    p4JetJesUp.SetPtEtaPhiE(ptUp,p4Jet.Eta(),p4Jet.Phi(),energyUp);
    p4JetJesDown.SetPtEtaPhiE(ptDown,p4Jet.Eta(),p4Jet.Phi(),energyDown);

    // PF jet ID variables
    float neutralHadFrac = neutralHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    float neutralEmFraction = neutralEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int nConstituents = chargedHadronMultiplicityAK5PFPUcorrJet[j] + neutralHadronMultiplicityAK5PFPUcorrJet[j] +
      photonMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] + muonMultiplicityAK5PFPUcorrJet[j] +
      HFHadronMultiplicityAK5PFPUcorrJet[j] + HFEMMultiplicityAK5PFPUcorrJet[j];
    float chargedHadFraction = chargedHadronEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    int chargedMultiplicity = chargedHadronMultiplicityAK5PFPUcorrJet[j] + electronMultiplicityAK5PFPUcorrJet[j] 
      + muonMultiplicityAK5PFPUcorrJet[j];
    float chargedEmFraction = chargedEmEnergyAK5PFPUcorrJet[j]/uncorrenergyAK5PFPUcorrJet[j];
    
    // if(!isPFJetID(fabs(etaAK5PFPUcorrJet[j]),neutralHadFrac,neutralEmFraction,nConstituents,
    //               chargedHadFraction,chargedMultiplicity,chargedEmFraction, Higgs::loose)) continue;
    
    bool foundMatch = false;

    // check if the electrons falls into the jet
    for(int i=0; i<(int)eleToRemove.size(); i++) {
      int ele = eleToRemove[i];
      if ( ele > -1 ) {
        TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Ele ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    // check if the muons falls into the jet
    for(int i=0; i<(int)muonToRemove.size(); i++) {
      int mu = muonToRemove[i];
      if ( mu > -1 ) {
        TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
        float deltaR =  fabs( p3Jet.DeltaR( p3Muon ) );
        H_deltaRcorr -> Fill(deltaR);
        // taking from ee config file, but jets veto is the same for all the channels
        if(_selectionEEE->getSwitch("jetConeWidth") && _selectionEEE->passCut("jetConeWidth",deltaR)) foundMatch=true;
      }
    }
    if(foundMatch) continue;

    if(pt>5.0) {
      jetSumNom += p4Jet;
      jetSumUp += p4JetJesUp;
      jetSumDown += p4JetJesDown;
    }

  }

  // add back the electron and muon candidates, with their calibrations
  TVector3 electronSum, muonSum;
  electronSum.SetXYZ(0,0,0);
  muonSum.SetXYZ(0,0,0);
  for(int i=0; i<(int)eleToRemove.size(); i++) {
    int ele = eleToRemove[i];
    if ( ele > -1 ) {
      TVector3 p3Ele(pxEle[ele],pyEle[ele],pzEle[ele]);
      electronSum += p3Ele;
    }
  }

  // check if the muons falls into the jet
  for(int i=0; i<(int)muonToRemove.size(); i++) {
    int mu = muonToRemove[i];
    if ( mu > -1 ) {
      TVector3 p3Muon(pxMuon[mu],pyMuon[mu],pzMuon[mu]);
      muonSum += p3Muon;
    }
  }

  TVector3 metFromJetsNom  = jetSumNom.Vect()  + electronSum + muonSum; 
  TVector3 metFromJetsUp   = jetSumUp.Vect()   + electronSum + muonSum; 
  TVector3 metFromJetsDown = jetSumDown.Vect() + electronSum + muonSum; 

  TVector3 diffUp   = metFromJetsUp   - metFromJetsNom;
  TVector3 diffDown = metFromJetsDown - metFromJetsNom;

  *m_metFromJets  = metFromJetsNom;
  *m_pfMetJESUp   = *m_p3PFMET + diffUp;
  *m_pfMetJESDown = *m_p3PFMET + diffDown;

}

std::pair<float,float> tHSelection::transvMassJES(int theChannel) {

  float mTUp   = sqrt( 2.*(m_dilepPt[theChannel].Pt())*(m_pfMetJESUp->Pt())*
		       (1- cos(m_pfMetJESUp->DeltaPhi(m_dilepPt[theChannel]))) );
  float mTDown = sqrt( 2.*(m_dilepPt[theChannel].Pt())*(m_pfMetJESDown->Pt())*
		       (1- cos(m_pfMetJESDown->DeltaPhi(m_dilepPt[theChannel]))) );

  return std::make_pair(mTUp,mTDown);

}

std::vector<TLorentzVector> tHSelection::GetJetJesPcomponent(int jet) {

  // [nom/+1s/-1s]

  if(jet<0) {
    TLorentzVector up(0,0,0,0);
    std::vector<TLorentzVector> zero;
    zero.push_back(up);
    zero.push_back(up);
    zero.push_back(up);    
    return zero;
  }

  TLorentzVector JP4(pxAK5PFPUcorrJet[jet],pyAK5PFPUcorrJet[jet],pzAK5PFPUcorrJet[jet],energyAK5PFPUcorrJet[jet]);
  TLorentzVector LJP4JesUp = GetJESCorrected(JP4,"Up"); // <--- crash here
  TLorentzVector LJP4JesDown = GetJESCorrected(JP4,"Down");

  std::vector<TLorentzVector> jes;
  jes.push_back(JP4);
  jes.push_back(LJP4JesUp);
  jes.push_back(LJP4JesDown);

  return jes;

}

void tHSelection::getDYGeneratorKinematics(int lepType) {

  TLorentzVector p1, p2;

  for(int imc=2;imc<40;imc++) {
    if(fabs(idMc[mothMc[imc]])==23) {
      float pt = pMc[imc]*fabs(sin(thetaMc[imc]));
      if(idMc[imc]==lepType) p1.SetPtEtaPhiE(pt,etaMc[imc],phiMc[imc],energyMc[imc]);
      if(idMc[imc]==-lepType) p2.SetPtEtaPhiE(pt,etaMc[imc],phiMc[imc],energyMc[imc]);
    }
  }
  
  if(p1.Pt()>0 && p2.Pt()>0) {
    _genmll = (p1+p2).M();
    _genptll = (p1+p2).Pt();
    _genyll = (p1+p2).Rapidity();
  } 

}

TVector3 tHSelection::getLeadingJet(int index, float ptThr) {
  TVector3 p3(0,0,0);
  if(index>-1) {
    p3.SetXYZ(pxAK5PFPUcorrJet[index],pyAK5PFPUcorrJet[index],pzAK5PFPUcorrJet[index]);
    if(p3.Pt()>ptThr) return p3;
  }
  return p3;
}

double tHSelection::ErrEt( double Et, double Eta) {
  
  double InvPerr2;
  
  double N, S, C, m;
  if(fabs(Eta) < 0.5 ) {
    N = 3.96859;
    S = 0.18348;
    C = 0.;
    m = 0.62627;
  } else if( fabs(Eta) < 1. ) {
    N = 3.55226;
    S = 0.24026;
    C = 0.;
    m = 0.52571;
  } else if( fabs(Eta) < 1.5 ) {
    N = 4.54826;
    S = 0.22652;
    C = 0.;
    m = 0.58963;
  } else if( fabs(Eta) < 2. ) {
    N = 4.62622;
    S = 0.23664;
    C = 0.;
    m = 0.48738;
  } else if( fabs(Eta) < 3. ) {
    N = 2.53324;
    S = 0.34306;
    C = 0.;
    m = 0.28662;
  } else if( fabs(Eta) < 5. ) {
    N = 2.95397;
    S = 0.11619;
    C = 0.;
    m = 0.96086;
  }
  
  // this is the absolute resolution (squared), not sigma(pt)/pt	
  // so have to multiply by pt^2, thats why m+1 instead of m-1	
  InvPerr2 =  (N * fabs(N) ) + (S * S) * pow(Et, m+1) + (C * C) * Et * Et ;
  
  return sqrt(InvPerr2)/Et;
}

bool tHSelection::isLooseJetMva(float pt, float eta, float id) {

  bool isOk = true;

  if (pt<10) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.0) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.2)  isOk = false;
  }

  if (pt<20 && pt>=10) {
    if (fabs(eta)<=2.5 && id<-0.4)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<-0.4) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<-0.4) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.4)   isOk = false;
  }

  if (pt<30 && pt>=20) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.2) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.6)  isOk = false;
  }

  if (pt<50 && pt>=30) {
    if (fabs(eta)<=2.5 && id<0.0)                   isOk = false;
    if (fabs(eta)>2.5 && fabs(eta)<=2.75 && id<0.0) isOk = false;
    if (fabs(eta)>2.75 && fabs(eta)<=3.0 && id<0.6) isOk = false;
    if (fabs(eta)>3.0 && fabs(eta)<=5.0 && id<0.2)  isOk = false;
  }

  return isOk;
}

