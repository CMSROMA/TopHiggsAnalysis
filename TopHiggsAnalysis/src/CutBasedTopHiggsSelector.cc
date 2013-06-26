#include "TopHiggsAnalysis/include/CutBasedTopHiggsSelector.hh"
#include <iostream>
#include <stdio.h>
#include <math.h>

CutBasedTopHiggsSelector::CutBasedTopHiggsSelector() {

  // latinos
  m_step0 = false;
  m_step1 = false;

  m_processID = -1;
}

CutBasedTopHiggsSelector::CutBasedTopHiggsSelector( const CutBasedTopHiggsSelector& selector ) {

  m_weight = selector.m_weight;
  m_foundMcTree = selector.m_foundMcTree;
  m_passedHLT = selector.m_passedHLT;
  m_isThisChannel = selector.m_isThisChannel;
  m_highPt = selector.m_highPt;
  m_isElectronId = selector.m_isElectronId;
  m_isPositronId = selector.m_isPositronId;
  m_isElectronIsol = selector.m_isElectronIsol;
  m_isPositronIsol = selector.m_isPositronIsol;
  m_isElectronConvRej = selector.m_isElectronConvRej;
  m_isPositronConvRej = selector.m_isPositronConvRej;
  m_isElectronIp = selector.m_isElectronIp;
  m_isPositronIp = selector.m_isPositronIp;
  m_invMass = selector.m_invMass;
  m_nJets            = selector.m_nJets;
  m_nUncorrJets      = selector.m_nUncorrJets;
  m_btagJets         = selector.m_btagJets;
  m_nSoftMuons       = selector.m_nSoftMuons;
  m_nExtraLeptons    = selector.m_nExtraLeptons;
  m_met = selector.m_met;
  m_projectedMet = selector.m_projectedMet;
  m_metOverPtLL = selector.m_metOverPtLL;
  m_dymva = selector.m_dymva;
  m_deltaPhiLLJet = selector.m_deltaPhiLLJet;
  m_deltaPhi = selector.m_deltaPhi;
  m_detaLeptons = selector.m_detaLeptons;
  m_maxPtElectron = selector.m_maxPtElectron;
  m_minPtElectron = selector.m_minPtElectron;
  m_ptll = selector.m_ptll;
  m_WWInvMass = selector.m_WWInvMass;
  m_nVtx = selector.m_nVtx;
  m_extraSlowLeptonPTMin = selector.m_extraSlowLeptonPTMin;
  m_processID = selector.m_processID;
  *_selection = *selector._selection;
  *globalCounter = *selector.globalCounter;
  multiProcessCounter = selector.multiProcessCounter;

  // latinos
  m_step0  = selector.m_step0;
  m_step1  = selector.m_step1;

}

CutBasedTopHiggsSelector::~CutBasedTopHiggsSelector() {}

void CutBasedTopHiggsSelector::Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle) {

  _selection = new Selection(std::string(fileCuts),std::string(fileSwitches));

  // these cuts are applied in the HiggsSelection class, but are configured here

  // General switches
  _selection->addSwitch("MCtruth");
  _selection->addSwitch("trigger");

  // Lepton selection switches
  /*
  _selection->addSwitch("leptonId");
  _selection->addSwitch("leptonIso");
  _selection->addSwitch("leptonD0"); 
  _selection->addSwitch("convRej");
  */

  _selection->addCut("electronIP");
  _selection->addCut("electronDz");
  _selection->addCut("muonIPhighPT");
  _selection->addCut("muonIPlowPT");
  _selection->addCut("muonDz");

  // Jet selection switches 
  _selection->addCut("etaJetAcc");
  _selection->addCut("etJetAcc");
  _selection->addCut("jetConeWidth");

  _selection->summary();

  globalCounter = new Counters();
  globalCounter->SetTitle(theTitle);
  globalCounter->AddVar("event"); 
  globalCounter->AddVar("MCtruth"); 
  globalCounter->AddVar("trigger");     // 0
  globalCounter->AddVar("preselected"); // 1

}

bool CutBasedTopHiggsSelector::output() {

  Counters *theCounter=0;

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter = multiProcessCounter.find(m_processID);

    if ( iter == multiProcessCounter.end() ) {
      
      std::cout << "First time I get process " << m_processID 
		<< ": adding a counter" << std::endl;

      char buffer[200];
      sprintf(buffer,"Event counter for process %d", m_processID);
      
      Counters *processCounter = new Counters();
      processCounter->SetTitle(buffer);
      processCounter->AddVar("event");
      processCounter->AddVar("MCtruth");
      processCounter->AddVar("trigger");
      processCounter->AddVar("preselected");

      multiProcessCounter.insert( std::make_pair(m_processID,processCounter) );      
    }

    theCounter = multiProcessCounter[m_processID];

  }
  
  else theCounter = globalCounter;

  // latinos
  m_step0 = false;
  m_step1 = false;

  theCounter->IncrVar("event",m_weight);
  
  if (_selection->getSwitch("MCtruth") && !m_foundMcTree) return false;
  theCounter->IncrVar("MCtruth",m_weight);

  if(_selection->getSwitch("trigger") && !m_passedHLT ) return false;
  theCounter->IncrVar("trigger",m_weight); 
  m_step0 = true;

  if(!m_isThisChannel) return false;
  theCounter->IncrVar("preselected",m_weight);
  m_step1 = true;
 
  return false; //?

}


void CutBasedTopHiggsSelector::displayEfficiencies(std::string datasetName) {

  if( m_processID > -1 ) {

    std::map<int, Counters*>::const_iterator iter;
    for( iter=multiProcessCounter.begin(); iter!=multiProcessCounter.end(); ++iter ) {

      Counters *theCounter = iter->second;

      theCounter->Draw();
      theCounter->Draw("trigger","MCtruth");
      theCounter->Draw("preselected","trigger");
    }
  }

  else {

    char namefile[500];
    sprintf(namefile,"%s-Counters.root",datasetName.c_str());
    
    globalCounter->Draw();
    globalCounter->Draw("trigger","MCtruth");
    globalCounter->Draw("preselected","trigger");
   
    globalCounter->Save(namefile,"update");
  }

}
