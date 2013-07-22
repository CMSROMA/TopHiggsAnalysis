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

  //! ID for leptons
  bool isSelectedMuon2012(int muonIndex);
  bool isSelectedMuon2013(int muonIndex);

  bool isSelectedElectron2012(int electronIndex);
  bool isSelectedElectron2013(int electronIndex);

  //! get the best MMM
  std::vector<int> getBestMMM(); // MMM channel

  //! get the best EEE
  std::vector<int> getBestEEE(); // EEE channel

  //! get the best MME
  std::vector<int> getBestMME(); // MME channel

  //! get the best EEM
  std::vector<int> getBestEEM(); // EEM channel

  //! get the best MM -- for Z-peak cross-checks
  std::vector<int> getBestMM();   

  //! get the best EE -- for Z-peak cross-checks
  std::vector<int> getBestEE();    

  //! set the 4 vectors, invariant mass, etc. after preselections and full selection
  void setKinematicsMMM(int myMu1 , int myMu2 , int myMu3 );
  void setKinematicsEEE(int myEle1, int myEle2, int myEle3);
  void setKinematicsMME(int myMu1 , int myMu2 , int myEle1);
  void setKinematicsEEM(int myEle1, int myEle2, int myMu1 );

  //! set kinematics for Z-peak  cross-check
  void setKinematicsMM(int myMu1 , int myMu2 );
  void setKinematicsEE(int myEle1, int myEle2 );

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

  //! for calculating lepton mva variables
  void estimateLeptonMVAvariables(int channel, int lepton1, int lepton2, int lepton3);

  //! to evaluate eleID
  bool m_useBDTEleID;

  //! to evaluate preselection efficiency
  Selection *_preselection;
  CommonHiggsPreselector CommonHiggsPreselection;

  //! to evaluate full selection efficiency
  Selection *_selectionEEE, *_selectionMMM, *_selectionMME, *_selectionEEM, *_selectionEE, *_selectionMM;

  CutBasedTopHiggsSelector CutBasedHiggsSelection[6];

  //! be verbose during runtime
  bool _verbose;

  //! process variables to initialize kFactors
  int _massVal;
  std::string _process;

  //! an integer defining the sub-channel
  enum { eee = 0, mmm = 1, eem = 2, mme = 3, ee = 4, mm = 5 };

  //! array containing the possibility of having reconstructed a certain sub-channel
  bool m_channel[6]; // channels defined above
  bool isOk[6];
  
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
  
  // -- note: the index [6] refers to the channel
  int theLeadingJet[6];
  int theSecondJet [6];
  std::vector<int> eleCands[6], muCands[6];
  TLorentzVector *m_p4Lepton1      [6], *m_p4Lepton2     [6], *m_p4Lepton3     [6];
  float           m_p4Lepton1Energy[6], m_p4Lepton2Energy[6], m_p4Lepton3Energy[6];
  int             m_p4Lepton1Type  [6], m_p4Lepton2Type  [6], m_p4Lepton3Type  [6];

  TVector3 *m_p3PFMET;
  TVector3 *m_metFromJets, *m_pfMetJESUp, *m_pfMetJESDown;
  TVector3 m_p3TKMET[6];
  float m_theMET;
  TLorentzVector *m_jetsSum[6], *m_jetscbIDSum[6], *m_uncorrJetsSum[6];

  TVector3 m_dilepPt[6];
  TVector3 m_trilepPt[6];
  float m_deltaPhi[6];
  float m_deltaErre[6];
  float m_deltaEtaLeptons[6];
  float m_mll[6];
  float m_transvMass[6];
  float m_MTR[6], m_MTRcharged[6], m_MR[6], m_GammaMR[6];
  float m_mT2[6];
  float m_projectedMet[6], m_chMet[6];
  float m_projectedPFMet[6], m_projectedTkMet[6];
  float m_metOptll[6];
  float hardestLeptonPt[6], slowestLeptonPt[6];
  float leadJetBtag[6], subleadJetBtag[6], subLeadJetsMaxBtag[6];
  float leadJetBtagBProb[6], subleadJetBtagBProb[6], subLeadJetsMaxBtagBProb[6];
  float m_softbdisc[6], m_hardbdisc[6];
  int   njets[6], ncbIDjets[6], nuncorrjets[6];
  int   m_goodvertices;

  // for jetId studies
  bool  wantJetIdStuff;
  float leadJetPt          [6], leadJetEta        [6], leadJetMvaJetId   [6];
  int   leadJetLoosePFId   [6], leadJetMatchGen   [6], leadJetLooseId    [6];
  float subleadJetPt       [6], subleadJetEta     [6], subleadJetMvaJetId[6];
  int   subleadJetLoosePFId[6], subleadJetMatchGen[6], subleadJetLooseId [6];

  // 4 channels, 3 leptons
  int   m_ch   [6][3];
  float m_iso  [6][3];
  float m_lh   [6][3];
  float m_bdt  [6][3];
  int   m_chmaj[6][3];

  //! B-Veto event variables
  float m_maxDxyEvt, m_maxDszEvt;
  float m_maxTrackCountingHighEffBJetTags, m_maxImpactParameterMVABJetTags, m_maxCombinedSecondaryVertexMVABJetTags; 
  int m_closestPV;
  int nsoftjets[6], nsoftjetscbID[6];

  // jet studies
  int m_numbtagjets[6];

  int m_numbtagCSVMcbIDaccepjets[6];//Cut Based ID for jets
  int m_numbtagCSVLcbIDaccepjets[6];
  int m_numbtagCSVTcbIDaccepjets[6];

  int m_numcbIDcentralLjets [6];//Cut Based ID for jets
  int m_numcbIDforwardLjets [6];//Cut Based ID for jets

  int m_numcbIDcentralMjets [6];//Cut Based ID for jets
  int m_numcbIDforwardMjets [6];//Cut Based ID for jets

  int m_numbtagCSVMmvaIDaccepjets[6];//MVA ID for jets
  int m_numbtagCSVLmvaIDaccepjets[6];
  int m_numbtagCSVTmvaIDaccepjets[6];

  int m_numbtagCSVMmvaIDcentraljets[6];//MVA ID for jets
  int m_numbtagCSVLmvaIDcentraljets[6];
  int m_numbtagCSVTmvaIDcentraljets[6];

  int m_nummvaIDcentraljets[6];//MVA ID for jets
  int m_nummvaIDforwardjets[6];//MVA ID for jets

  int m_nummvaIDaccepOUTjets[6];//MVA ID for jets
  int m_nummvaIDaccepINjets [6];//MVA ID for jets

  float ptCVStaggedM [6][5];
  float etaCVStaggedM[6][5];
  float phiCVStaggedM[6][5];
  float eneCVStaggedM[6][5];
  float cvsCVStaggedM[6][5];

  float ptForwardM [6][5];
  float etaForwardM[6][5];
  float phiForwardM[6][5];
  float eneForwardM[6][5];
  float cvsForwardM[6][5];

  float ptCentralM [6][5];
  float etaCentralM[6][5];
  float phiCentralM[6][5];
  float eneCentralM[6][5];
  float cvsCentralM[6][5];

  float ptCVStaggedL [6][5];
  float etaCVStaggedL[6][5];
  float phiCVStaggedL[6][5];
  float eneCVStaggedL[6][5];
  float cvsCVStaggedL[6][5];

  float ptForwardL [6][5];
  float etaForwardL[6][5];
  float phiForwardL[6][5];
  float eneForwardL[6][5];
  float cvsForwardL[6][5];

  float ptCentralL [6][5];
  float etaCentralL[6][5];
  float phiCentralL[6][5];
  float eneCentralL[6][5];
  float cvsCentralL[6][5];

  int _theGenEle, _theGenPos;
  int _theGenMuMinus, _theGenMuPlus;

  std::vector<int> selectedLeptons;
  std::vector<int> selectedMuons;
  std::vector<int> selectedElectrons;

  //! vector to store indices of candidate to include / exclude
  std::vector<int> m_goodJets;
  std::vector<int> m_goodcbIDJets;// Cut Based ID

  //! reduced tree for event selection (on at least 3 leptons events)
  RedTopHiggsTree *myOutTree[6];

  //! reduced tree for trigger studies (on all events)
  RedTriggerTree *myTriggerTree;

  //! reduced tree for gen level variables
  RedTopHiggsTree *myGenLevelTree;

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

  float _genForwardQuark_Pt;
  float _genForwardQuark_Eta;

  float _genbQuark_Pt;
  float _genbQuark_Eta;

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

  //! Lepton MVA input variables
  float leptPt[6][3];
  float leptEta[6][3];
  float neuRelIs[6][3];
  float chRelIso[6][3];
  float jetDR_in[6][3];
  float jetPtRatio_in[6][3];
  float jetBTagCSV_in[6][3];
  float sip3d[6][3];
  float mvaId[6][3];
  int   innerHits[6][3];
  float logdxy[6][3];
  float logdz[6][3];

};
#endif
