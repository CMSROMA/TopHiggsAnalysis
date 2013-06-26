//-------------------------------------------------------
// Description:
//    Class for selection of reconstructed qtH
// Original HiggsMLSelection. Authors:
//    Chiara Rovelli & Emanuele Di Marco
//    Universita' di Roma "La Sapienza" & INFN Roma
// Modified for tH analysis by C. Jorda
//-------------------------------------------------------

#ifndef tHSelection_h
#define tHSelection_h

#include <vector>
#include "CommonTools/include/Monitor.hh"
#include "EgammaAnalysisTools/include/CutBasedEleIDSelector.hh"
#include "EgammaAnalysisTools/include/ElectronLikelihood.h"
#include "TopHiggsAnalysis/include/TopHiggs.hh"
#include "HiggsAnalysisTools/include/CommonHiggsPreselector.hh"
#include "HiggsAnalysisTools/include/CutBasedHiggsSelector.hh"
#include "TopHiggsAnalysis/include/CutBasedTopHiggsSelector.hh"
#include "HiggsAnalysisTools/include/RedHiggsTree.h"
#include "TopHiggsAnalysis/include/RedTopHiggsTree.h"
#include "HiggsAnalysisTools/include/RedTriggerTree.hh"
#include "HiggsAnalysisTools/include/RedEleIDTree.h"
#include "HiggsAnalysisTools/include/kFactorEvaluator.hh"
#include <TVector3.h>
#include <TLorentzVector.h>

// Rochester correction for muons
#include "TopHiggsAnalysis/include/RochCor2012.h"

class tHSelection : public TopHiggs{
public:
  
  //! constructor
  tHSelection(TTree *tree=0);
  //! destructor
  virtual ~tHSelection();
  //! loop over events
  void Loop();
  //! set the name for dataset in output
  void SetDatasetName(std::string filename) {_datasetName=filename;};
  //! display the efficiency table
  void displayEfficiencies(std::string filename);
  //! set the required triggers masks (one per channel)
  void setRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  //! set the not-required triggers masks (one per channel)
  void setNotRequiredTriggers(const std::vector<std::string>& reqTriggers, int channel);
  //! print some counters at the end of the job
  void counterEndJob();

private:

  bool findMcTree(const char* processType);

  //! get the best MMM
  std::vector<int> getBestMMM(); // MMM channel

  //! get the best EEE
  std::vector<int> getBestEEE(); // EEE channel

  //! get the best MME
  std::vector<int> getBestMME(); // MME channel

  //! get the best EEM
  std::vector<int> getBestEEM(); // EEM channel

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsMMM(int myMu1 , int myMu2 , int myMu3 );
  void setKinematicsEEE(int myEle1, int myEle2, int myEle3);
  void setKinematicsMME(int myMu1 , int myMu2 , int myEle1);
  void setKinematicsEEM(int myEle1, int myEle2, int myMu1 );

  //! reset the kinematic quantities at the beginning of event and after the selection if needed
  void resetKinematicsStart();
  void resetKinematics();

  //! count jet multiplicity
  int numJets       ( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel );
  int numcbIDJets   ( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel );
  int numUncorrJets ( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel );

  //! calculate the Up/Down smeared met applying JES uncertainties
  void JESPfMet( std::vector<int> eleToRemove, std::vector<int> muonToRemove);

  //! calculate the Up/Down smeared MT
  std::pair<float,float> transvMassJES(int theChannel);

  //! calculate each component of a jet 3-momentum for: 0=nominal,1=JES up,2=JES down. Components are: 0/1/2 = x/y/z
  std::vector<TLorentzVector> GetJetJesPcomponent(int jet);

  //! give the highest b-tag of calojets in the event
  float bVetoJets( std::vector<int> eleToRemove, std::vector<int> muonToRemove, int theChannel );

  //! in the 1-jet bin, deltaphi between ll system and leading jet
  float deltaPhiLLJet(int ichan);

  //! count the soft muons
  int numSoftMuons(std::vector<int> muonToRemove, std::vector<int> jetsToRemove);

  //! count the extra leptons (id, iso, d0,acceptance etc) with pt>10 GeV
  int numExtraLeptons( std::vector<int> eleToRemove, std::vector<int> muonToRemove );

  //! get the kFactor of the event
  float getkFactor(std::string process);

  //! get the generator level quantities for DY->ll (ee/mumu only)
  void getDYGeneratorKinematics(int lepType);

  //! set the electron/muon ID variables to dump
  void setLepIdVariables(int index1, int index2, int index3, int pdgid1, int pdgid2, int pdgid3);

  //! search for the hardest lepton vertex
  int getPV();

  //! dxy parameter with respect to PV for electron tracks
  double trackDxyPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);

  //! dsz parameter with respect to PV for electron tracks
  double trackDszPV(float PVx, float PVy, float PVz, float eleVx, float eleVy, float eleVz, float elePx, float elePy, float elePz);

  //! methods for the jet veto: track quality
  bool isGoodTrack(int iTrack, float ptMin, float ptMax, float chi2, float etaMax, float nHits);

  //! methods for the jet veto: per jet variable
  std::vector<float> jetBTagVariables(int jetIndex);

  //! methods for the jet veto: event based variable
  void calcEventBVetoVariables(std::vector<int> jets);

  //! method to evaluate Mt from the lepton and neutrino pT's
  double mT(TVector3 plep, TVector3 pneu);

  //! method to evaluate Mt2 from the 2 leptons pT's and pTmiss
  double mT2(TVector3 plep1, TVector3 plep2, TVector3 ptmiss);

  //! Get pt given x/y coordinates
  float GetPt(float px, float py) { return TMath::Sqrt(px*px + py*py); }

  //! for the two leptons, find the closest one to MET in phi. if(deltaPhi(lepton-MET)<pi/2) projectedMET = MET * sin(deltaPhi(lepton-MET)); else projectedMET = MET
  float GetProjectedMet(TVector3 p1, TVector3 p2);
  float GetProjectedPFMet(TVector3 p1, TVector3 p2);
  float GetProjectedTkMet(TVector3 p1, TVector3 p2);

  //! reload the trigger mask_s_ (one per channel)
  bool reloadTriggerMask(int runN);

  //! get the trigger answer depending on the channel
  bool hasPassedHLT(int channel);

  //! get the leading jet three-momentum; 0 if it is below threshold in pT
  TVector3 getLeadingJet(int index, float ptThr=15.0);

  //! for jetId studies
  double ErrEt( double Et, double Eta);
  bool isLooseJetMva(float pt, float eta, float id);

  //! to evaluate eleID
  bool m_useBDTEleID;

  //! to evaluate preselection efficiency
  Selection *_preselection;
  CommonHiggsPreselector CommonHiggsPreselection;

  //! to evaluate full selection efficiency
  Selection *_selectionEEE, *_selectionMMM, *_selectionMME, *_selectionEEM;

  CutBasedTopHiggsSelector CutBasedHiggsSelection[4];

  //! be verbose during runtime
  bool _verbose;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  //! an integer defining the sub-channel
  enum { eee = 0, mmm = 1, eem = 2, mme = 3 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[4]; // channels defined above
  bool isOk[4];
  
  //! trigger masks
  std::vector<int> m_requiredTriggersEEE, m_notRequiredTriggersEEE;
  std::vector<int> m_requiredTriggersMMM, m_notRequiredTriggersMMM;
  std::vector<int> m_requiredTriggersEEM, m_notRequiredTriggersEEM;
  std::vector<int> m_requiredTriggersMME, m_notRequiredTriggersMME;

  std::vector<std::string>    requiredTriggersEEE,    requiredTriggersMMM,    requiredTriggersEEM,    requiredTriggersMME;
  std::vector<std::string> notRequiredTriggersEEE, notRequiredTriggersMMM, notRequiredTriggersEEM, notRequiredTriggersMME;

  //! kinematics of the event
  int theElectron,  thePositron;
  int theMuonMinus, theMuonPlus;
  int thePreElectron,  thePrePositron, thePreElectronME, thePreElectronEM;
  int thePreMuonMinus, thePreMuonPlus, thePreMuonME, thePreMuonEM;
  
  // -- note: the index [4] refers to the channel
  int theLeadingJet[4];
  int theSecondJet [4];
  std::vector<int> eleCands[4], muCands[4];
  TLorentzVector *m_p4Lepton1      [4], *m_p4Lepton2     [4], *m_p4Lepton3     [4];
  float           m_p4Lepton1Energy[4], m_p4Lepton2Energy[4], m_p4Lepton3Energy[4];
  int             m_p4Lepton1Type  [4], m_p4Lepton2Type  [4], m_p4Lepton3Type  [4];

  TVector3 *m_p3PFMET;
  TVector3 *m_metFromJets, *m_pfMetJESUp, *m_pfMetJESDown;
  TVector3 m_p3TKMET[4];
  float m_theMET;
  TLorentzVector *m_jetsSum[4], *m_jetscbIDSum[4], *m_uncorrJetsSum[4];

  TVector3 m_dilepPt[4];
  TVector3 m_trilepPt[4];
  float m_deltaPhi[4];
  float m_deltaErre[4];
  float m_deltaEtaLeptons[4];
  float m_mll[4];
  float m_transvMass[4];
  float m_MTR[4], m_MTRcharged[4], m_MR[4], m_GammaMR[4];
  float m_mT2[4];
  float m_projectedMet[4], m_chMet[4];
  float m_projectedPFMet[4], m_projectedTkMet[4];
  float m_metOptll[4];
  float hardestLeptonPt[4], slowestLeptonPt[4];
  float leadJetBtag[4], subleadJetBtag[4], subLeadJetsMaxBtag[4];
  float leadJetBtagBProb[4], subleadJetBtagBProb[4], subLeadJetsMaxBtagBProb[4];
  float m_softbdisc[4], m_hardbdisc[4];
  int   njets[4], ncbIDjets[4], nuncorrjets[4];
  int   m_goodvertices;

  // for jetId studies
  bool  wantJetIdStuff;
  float leadJetPt          [4], leadJetEta        [4], leadJetMvaJetId   [4];
  int   leadJetLoosePFId   [4], leadJetMatchGen   [4], leadJetLooseId    [4];
  float subleadJetPt       [4], subleadJetEta     [4], subleadJetMvaJetId[4];
  int   subleadJetLoosePFId[4], subleadJetMatchGen[4], subleadJetLooseId [4];

  // 4 channels, 3 leptons
  int   m_ch   [4][3];
  float m_iso  [4][3];
  float m_lh   [4][3];
  float m_bdt  [4][3];
  int   m_chmaj[4][3];

  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  int m_closestPV;
  int nsoftjets[4], nsoftjetscbID[4];

  // jet studies
  int m_numbtagjets[4];

  int m_numbtagCSVMmvaIDcentraljets[4];//MVA ID for jets
  int m_numbtagCSVLmvaIDcentraljets[4];
  int m_numbtagCSVTmvaIDcentraljets[4];

  int m_numbtagCSVMcbIDcentraljets[4];//Cut Based ID for jets
  int m_numbtagCSVLcbIDcentraljets[4];
  int m_numbtagCSVTcbIDcentraljets[4];

  int m_nummvaIDforwardjets[4];//MVA ID for jets
  int m_numcbIDforwardjets [4];//Cut Based ID for jets
  
  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  std::vector<int> selectedLeptons;
  std::vector<int> selectedMuons;
  std::vector<int> selectedElectrons;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;
  std::vector<int> m_goodcbIDJets;// Cut Based ID

  //! reduced tree for event selection (on at least 3 leptons events)
  RedTopHiggsTree *myOutTree[4];

  //! reduced tree for trigger studies (on all events)
  RedTriggerTree *myTriggerTree;

  //! variables for EleID
  int   myLepFlav[3];
  float myPt[3], myEta[3], myPhi[3];
  float myLepId[3], myLepIso[3], myConv[3];
  float myR9[3], mySCEnergy[3];

  //! new variables
  float m_eOverP[100];

  //! tH gen level

  float _genHiggsPt;
  float _genHiggsEta;

  float _genTopPt;
  float _genTopEta;

  float _genWpfromH_Pt;
  float _genWpfromH_Eta;

  float _genWmfromH_Pt;
  float _genWmfromH_Eta;

  float _genWfromT_Pt;
  float _genWfromT_Eta;
  
  float _genLeptonPlusfromWfromH_Pt;
  float _genLeptonPlusfromWfromH_Eta;

  float _genLeptonMinusfromWfromH_Pt;
  float _genLeptonMinusfromWfromH_Eta;

  float _genLeptonfromWfromT_Pt;
  float _genLeptonfromWfromT_Eta;

  float _genNeutrinoPlusfromWfromH_Pt;
  float _genNeutrinoPlusfromWfromH_Eta;

  float _genNeutrinoMinusfromWfromH_Pt;
  float _genNeutrinoMinusfromWfromH_Eta;

  float _genNeutrinofromWfromT_Pt;
  float _genNeutrinofromWfromT_Eta; 

  float _highestPtGen[1], _lowestPtGen[1];
  float _nGenJet[1];
  float _emFracEle[50], _hadFracEle[50];


  //! counters
  int nWWtoLLdecays;
  int nWfromTtoLdecay;
  int nSSevent;
  int nTotal;

  //! name of rootfile with dataset
  std::string _datasetName;

  //! to check the electron/jet matching
  TFile *fMatch;
  TH1F *H_deltaRcorr;
  TH1F *H_deltaRuncorr;

  //! kfactor evaluator offline
  kFactorEvaluator *calculator_;

  //! DY generator level quantities
  float _genmll, _genptll, _genyll;

  //! Rochester correction for muons
  RochCor2012 *rmcor;

};
#endif
