 
#ifndef RedTopHiggsTree_h
#define RedTopHiggsTree_h

#include "TLorentzVector.h"
#include "TVector3.h"

class TFile;
class TTree;

class G3EventProxy;

class RedTopHiggsTree {
public:
   RedTopHiggsTree(const char * filename = "eleID.root");
  ~RedTopHiggsTree();

  //! add more informations for analysis not cut based
  void addMLVars();
  //! add infos for jetId studies
  void addJetsVars();
  //! add the electron ID+iso variables for the selected best electrons
  void addLeptonInfos();
  //! add the CSA07 processID and weight block
  void addCSA07Infos();
  //! add the k-Factor (used for signal only)
  void addKFactor();
  //! add the MC truth informations
  void addMcTruthInfos();
  //! add the HLT electron triggers informations
  void addHLTElectronsInfos();
  //! add the HLT muon triggers informations
  void addHLTMuonsInfos();
  //! add run,lumi, event number (for data)
  void addRunInfos();
  //! add latinos
  void addLatinos();
  //! add kinematics
  void addKinematics();
  //! add razor variables
  void addRazor();
  //! add variables for W+jets
  void addFake();
  //! add variables to study shape systematics
  void addSystematics();
  //! for met studies
  void addMetStudies();
  //! add PDFs
  void addPDFs();
  //! add MCTruth for tH MC sampe
  void addtHMcTruthInfos();
  

  //! event by event final dataset fill
  void fillAll(float met, float pfmet, float cmet, float projmet, 
	       float dphi, float derre, float tmass, float mee, float max, float min, float deta, int nvtx,
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);

  void fillAll(float met, float pfmet, float cmet, float projmet, 
	       float dphi, float derre, float tmass, float mee, 
	       float max, float min, float maxEta, float minEta,
	       float deta, int nvtx,
	       bool finalLeptons, bool jetVeto, bool uncorrJetVeto, bool preDeltaPhi, bool finalSelection);

  void fillKinematics(float pxTkMet, float pyTkMet, float pzTkMet,
                      float pxLeadJet[3], float pyLeadJet[3], float pzLeadJet[3],
                      float pxSecJet[3], float pySecJet[3], float pzSecJet[3],
                      float pxL1, float pyL1, float pzL1,
                      float pxL2, float pyL2, float pzL2,
                      float pxL3, float pyL3, float pzL3,
                      int ch[3],
		      float bdt[3],
                      TLorentzVector *jetSum, TLorentzVector *uncorrJetSum, TVector3 *pfmet);
  
  void fillSystematics(float scE[3], float r9[3], float ene1, float ene2, float ene3, int ty1, int ty2, int ty3,
                       TVector3 *metFromJets, TVector3 *pfMetUp, TVector3 *pfMetDown, float mtUp, float mtDown); 
                  
  void fillRazor(float MTR, float mR, float gammaMR);

  void fillFake(int ntigh, float wfp, float wsfp, 
		float wfp15, float wsfp15, float wff15, float wsff15, float wpp15, float wspp15,
		float wfp35, float wsfp35, float wff35, float wsff35, float wpp35, float wspp35,
		float wfp50, float wsfp50, float wff50, float wsff50, float wpp50, float wspp50,
		float wfpQCD, float wsfpQCD, float wffQCD, float wsffQCD, float wppQCD, float wsppQCD);
  //! fill more informations for analysis not cut based
  void fillMLVars(int njets, int nuncorrjets, int   ncbIDjets, float dxyEVT, float dszEVT,
                  float softbdisc, float hardbdisc, float bTagSecVertex, int nSoftMu, 
                  float leadJetBTagTrackCount, float subleadJetBTagTrackCount, float subleadJetsMaxBTagTrackCount, 
                  float leadJetBTagBProb, float subleadJetBTagBProb, float subleadJetsMaxBTagBProb, 
                  int numExtraLep, int nSoftMuNoJets, int nSoftBJets, int nSoftJets,
		  int numbtagCSVLcbIDcentraljets,  int numbtagCSVMcbIDcentraljets,  int numbtagCSVTcbIDcentraljets,
		  int numcbIDforwardjets,  
		  int numbtagCSVLmvaIDcentraljets, int numbtagCSVMmvaIDcentraljets, int numbtagCSVTmvaIDcentraljets,
		  int nummvaIDforwardjets);

  //! fill infos for jetId studies
  void fillJetsVars(float ljpt, float ljeta, int ljpfid, int ljmatch, float ljmva, int ljl, float sljpt, float sljeta, int sljpfid, int sljmatch, float sljmva, int sljl);

  //! fill lepton ID variables
  void fillLeptons(float pt[3], float eta[3], float phi[3], int flavour[3],
                   float lepid[3], float lepiso[3], float lepconv[3]);
  //! fill the CSA07 processID and weight and lumi (in pb-1)
  void fillCSA07(double weight, double processId, float lumi=1000.);
  //! fill with the k-Factor (used for signal only)
  void fillKFactor(float kfactor, float genh, float ptlj );
  //! fill the MC truth informations
  void fillMcTruth(bool prompt, float genmll, float genptll, float genyll);
  //! fill the HLT electron triggers informations
  void fillHLTElectrons(bool singleEle, bool singleEleRelaxed, bool singleEleOR);
  //! fill the HLT muons triggers informations
  void fillHLTMuons(bool singleMuon, bool singleMuonRelaxed, bool singleMuonOR);
  //! fill the run,lumi, event number
  void fillRunInfos(int run, int lumi, long event, float puweight, bool HLT, float npu[3]);   
  //! latinos 
  void fillLatinos(bool s0, bool s1);
  //! met studies
  void fillMetStudies(float projPF, float projTk, float signPFMet, float signChMet, float m_MTRcha, float DYMVA, float rho, float rhojet );
  //! PDFs
  void fillPDFs(double cteq66[45], double mstw[31], double nnpdf[101]);
  //! MCTruth for tH
  void filltHMcTruthInfos(float genHiggsPt,
			  float genHiggsEta,
			  
			  float genTopPt, 
			  float genTopEta,
			  
			  float genWpfromH_Pt, 
			  float genWpfromH_Eta,
			  
			  float genWmfromH_Pt,
			  float genWmfromH_Eta,
			  
			  float genWfromT_Pt,
			  float genWfromT_Eta,
			  
			  float genLeptonPlusfromWfromH_Pt,
			  float genLeptonPlusfromWfromH_Eta,
			  
			  float genLeptonMinusfromWfromH_Pt,
			  float genLeptonMinusfromWfromH_Eta,
			  
			  float genLeptonfromWfromT_Pt,
			  float genLeptonfromWfromT_Eta,
			  
			  float genNeutrinoPlusfromWfromH_Pt,
			  float genNeutrinoPlusfromWfromH_Eta,
			  
			  float genNeutrinoMinusfromWfromH_Pt,
			  float genNeutrinoMinusfromWfromH_Eta,
			  
			  float genNeutrinofromWfromT_Pt,
			  float genNeutrinofromWfromT_Eta,

			  float genForwardQuark_Pt,
			  float genForwardQuark_Eta);
  

  //! effectively store the events in the tree
  void store();
  //! save in the ROOT file
  void save();

private:
  double myCTEQ66[45], myMSTW[31], myNNPDF[101];
  bool myHLT;
  bool myPromptDecay;
  float myGenptll, myGenyll, myGenmll;
  bool myHLTSingleElectron;
  bool myHLTSingleElectronRelaxed;
  bool myHLTSingleElectronOR;
  bool myHLTSingleMuon;
  bool myHLTSingleMuonRelaxed;
  bool myHLTSingleMuonOR;
  float myMet;       
  float myPFMet;       
  float myCaloMet;
  float myProjectedMet;
  float myDeltaPhi;  
  float myDeltaR;  
  float myTransvMass;
  float myEleInvMass;
  float maxPtEle;  
  float minPtEle;  
  float maxEtaEle;  
  float minEtaEle;  
  float myDetaLeptons;
  int myNVtx;
  float myNpu[3];
  int myNjets;
  int mycbIDNjets;
  int myNuncorrjets;
  float myDxyEVT;
  float myDszEVT;
  float mySoftBdisc;
  float myHardBdisc;
  float myBTagSecVertex;
  int myNSoftMu, myNSoftMuNoJets;
  float myLeadingJetBTagTrackCount, mySubleadingJetBTagTrackCount, mySubleadingJetsMaxBTagTrackCount;
  float myLeadingJetBTagJetBProb, mySubleadingJetBTagJetBProb, mySubleadingJetsMaxBTagJetBProb;
  int  myNumExtraLep, myNumSoftBJets, myNumSoftJets;

  int myNumbtagCSVLcbIDcentraljets , myNumbtagCSVMcbIDcentraljets , myNumbtagCSVTcbIDcentraljets;
  int myNumcbIDforwardjets;
  int myNumbtagCSVLmvaIDcentraljets, myNumbtagCSVMmvaIDcentraljets, myNumbtagCSVTmvaIDcentraljets;
  int myNummvaIDforwardjets;

  bool myFinalLeptons;
  bool myJetVeto;
  bool myUncorrJetVeto;
  bool myPreDeltaPhi;
  bool myFinalSelection;
  double myWeight;
  double myProcesId;
  float myLumi;
  float myKFactor, myPUWeight;
  float myGenHPt;
  float myLeadingJetPt;
  int myRun, myLS;
  long myEvent;
  float myPxTkMet, myPyTkMet, myPzTkMet;
  float myPxLeadJet[3], myPyLeadJet[3], myPzLeadJet[3];
  float myPxSecondJet[3], myPySecondJet[3], myPzSecondJet[3];
  float myPxL1, myPyL1, myPzL1;
  float myPxL2, myPyL2, myPzL2;
  float myPxL3, myPyL3, myPzL3;
  float myEneL1, myEneL2, myEneL3;
  int myTypeL1, myTypeL2, myTypeL3;

  float myMTR, myMR, myGammaMR, myDYMVA;
  float myRho, myRhoJet;

  float myProjPFMet, myProjPFChargedMet;
  float mySignPFMet, mySignPFChargedMet;
  float myMTRchargedMet;
  
  TLorentzVector *myJetsSum, *myUncorrJetsSum;
  TVector3 *myPfMet;

  //! for W+jets
  int myTight;
  float myWeightFP, myWeightStatFP;
  float myWeightFP15, myWeightStatFP15;
  float myWeightFF15, myWeightStatFF15;
  float myWeightPP15, myWeightStatPP15;
  float myWeightFP35, myWeightStatFP35;
  float myWeightFF35, myWeightStatFF35;
  float myWeightPP35, myWeightStatPP35;
  float myWeightFP50, myWeightStatFP50;
  float myWeightFF50, myWeightStatFF50;
  float myWeightPP50, myWeightStatPP50;
  float myWeightFPQCD, myWeightStatFPQCD;
  float myWeightFFQCD, myWeightStatFFQCD;
  float myWeightPPQCD, myWeightStatPPQCD;

  //! gen variables for tH
  float myGenHiggsPt;
  float myGenHiggsEta;

  float myGenTopPt;
  float myGenTopEta;

  float myGenWpfromH_Pt;
  float myGenWpfromH_Eta;

  float myGenWmfromH_Pt;
  float myGenWmfromH_Eta;

  float myGenWfromT_Pt;
  float myGenWfromT_Eta;
  
  float myGenLeptonPlusfromWfromH_Pt;
  float myGenLeptonPlusfromWfromH_Eta;

  float myGenLeptonMinusfromWfromH_Pt;
  float myGenLeptonMinusfromWfromH_Eta;

  float myGenLeptonfromWfromT_Pt;
  float myGenLeptonfromWfromT_Eta;

  float myGenNeutrinoPlusfromWfromH_Pt;
  float myGenNeutrinoPlusfromWfromH_Eta;

  float myGenNeutrinoMinusfromWfromH_Pt;
  float myGenNeutrinoMinusfromWfromH_Eta;

  float myGenNeutrinofromWfromT_Pt;
  float myGenNeutrinofromWfromT_Eta;

  float myGenForwardQuark_Pt;
  float myGenForwardQuark_Eta;

  // latinos
  bool mySteps[2];

  // lepton variables
  int   myLepFlav[3];
  float myPt[3], myEta[3], myPhi[3];
  float myLepId[3], myLepIso[3], myLepConv[3];
  float myLepBDT[3];
  TVector3 *myMetFromJets, *myPfMetUp, *myPfMetDown;
  float myMtUp, myMtDown;

  //! jet Id variables
  float myLeadJetPt,      myLeadJetEta,      myLeadJetIdMva;
  int   myLeadJetLooseId, myLeadJetGenMatch, myLeadJetPassLooseId;
  float mySubleadJetPt,      mySubleadJetEta,      mySubleadJetIdMva;
  int   mySubleadJetLooseId, mySubleadJetGenMatch, mySubleadJetPassLooseId;

  // lepton variables
  int myLepCharge[3];
  float myScEnergy[3], myR9[3];

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedTopHiggsTree_h
