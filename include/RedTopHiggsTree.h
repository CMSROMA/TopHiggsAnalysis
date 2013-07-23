
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
  //! add infos for more general jetId studies (2013)                                                                                                                                                         
  void add2013IDJetsVars();
  //! add more infos for jetId studies
  void addMoreJetsVars();
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
  //! add lepton mva variables
  void addLeptonMVAidVariables();

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
		  // cb
		  int numbtagCSVLcbIDaccepjets,  
		  int numbtagCSVMcbIDaccepjets,  
		  int numbtagCSVTcbIDaccepjets,

		  int numcbIDcentralLjets,  
		  int numcbIDforwardLjets,  

		  int numcbIDcentralMjets,  
		  int numcbIDforwardMjets,  

		  // mva
		  int numbtagCSVLmvaIDaccepjets,  
		  int numbtagCSVMmvaIDaccepjets,  
		  int numbtagCSVTmvaIDaccepjets,

		  int numbtagCSVLmvaIDcentraljets, 
		  int numbtagCSVMmvaIDcentraljets, 
		  int numbtagCSVTmvaIDcentraljets,

		  int nummvaIDcentraljets,
		  int nummvaIDforwardjets,

		  int nummvaIDaccepINjets  ,  
		  int nummvaIDaccepOUTjets);

  //! fill infos for jetId studies
  void fillJetsVars(float ljpt, float ljeta, int ljpfid, int ljmatch, float ljmva, int ljl, float sljpt, float sljeta, int sljpfid, int sljmatch, float sljmva, int sljl);

  //! fill more infos for a more general jetId study (2013)                                                                                                                                                   
  void fill2013ExtraJetsVars(float jetPuId_eta[20], float jetPuId_phi[20], float jetPuId_energy[20], float jetPuId_pt[20], float jetPuId_csv[20], float jetPuId_betastar[20], float jetPuId_rms[20], 
			     int jetPuId_cutBased[20], float jetPuId_mvaBased[20], int jetPuId_associated[20]);


  //! fill more infos for jetId studies
  void fillMoreJetsVars(float ptCVStaggedM[5], float etaCVStaggedM[5], float phiCVStaggedM[5], float eneCVStaggedM[5], float cvsCVStaggedM[5],
			float   ptForwardM[5], float   etaForwardM[5], float   phiForwardM[5], float   eneForwardM[5], float cvsForwardM[5],
			float   ptCentralM[5], float   etaCentralM[5], float   phiCentralM[5], float   eneCentralM[5], float cvsCentralM[5],

			float ptCVStaggedL[5], float etaCVStaggedL[5], float phiCVStaggedL[5], float eneCVStaggedL[5], float cvsCVStaggedL[5],
			float   ptForwardL[5], float   etaForwardL[5], float   phiForwardL[5], float   eneForwardL[5], float   cvsForwardL[5],
			float   ptCentralL[5], float   etaCentralL[5], float   phiCentralL[5], float   eneCentralL[5], float   cvsCentralL[5]);
			
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
			  float genForwardQuark_Eta,

			  float genbQuark_Pt,
			  float genbQuark_Eta);
  
  //! add lepton mva variables
  void fillLeptonMVAidVariables(float leptPt[3],
				float leptEta[3],
				float neuRelIs[3],
				float chRelIso[3],
				float jetDR_in[3],
				float jetPtRatio_in[3],
				float jetBTagCSV_in[3],
				float sip3d[3],
				float mvaId[3],
				int   innerHits[3],
				float logdxy[3],
				float logdz[3]);
  
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

  int myNumbtagCSVLcbIDaccepjets; 
  int myNumbtagCSVMcbIDaccepjets; 
  int myNumbtagCSVTcbIDaccepjets;

  int myNumcbIDcentralLjets;
  int myNumcbIDforwardLjets;

  int myNumcbIDcentralMjets;
  int myNumcbIDforwardMjets;

  int myNumbtagCSVLmvaIDaccepjets; 
  int myNumbtagCSVMmvaIDaccepjets; 
  int myNumbtagCSVTmvaIDaccepjets;

  int myNumbtagCSVLmvaIDcentraljets; 
  int myNumbtagCSVMmvaIDcentraljets; 
  int myNumbtagCSVTmvaIDcentraljets;

  int myNummvaIDcentraljets;
  int myNummvaIDforwardjets;

  int myNummvaIDaccepINjets;
  int myNummvaIDaccepOUTjets;

  float myptCVStaggedM [5];
  float myetaCVStaggedM[5];
  float myphiCVStaggedM[5];
  float myeneCVStaggedM[5];
  float mycvsCVStaggedM[5];

  float myptForwardM [5];
  float myetaForwardM[5];
  float myphiForwardM[5];
  float myeneForwardM[5];
  float mycvsForwardM[5];

  float myptCentralM [5];
  float myetaCentralM[5];
  float myphiCentralM[5];
  float myeneCentralM[5];
  float mycvsCentralM[5];

  float myptCVStaggedL [5];
  float myetaCVStaggedL[5];
  float myphiCVStaggedL[5];
  float myeneCVStaggedL[5];
  float mycvsCVStaggedL[5];

  float myptForwardL [5];
  float myetaForwardL[5];
  float myphiForwardL[5];
  float myeneForwardL[5];
  float mycvsForwardL[5];

  float myptCentralL [5];
  float myetaCentralL[5];
  float myphiCentralL[5];
  float myeneCentralL[5];
  float mycvsCentralL[5];

  float myJetPuId_eta[20];
  float myJetPuId_phi[20];
  float myJetPuId_energy[20];
  float myJetPuId_pt[20];
  float myJetPuId_csv[20];
  float myJetPuId_betastar[20];
  float myJetPuId_rms[20];
  int myJetPuId_cutBased[20];
  float myJetPuId_mvaBased[20];
  int myJetPuId_associated[20];

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

  float myGenbQuark_Pt;
  float myGenbQuark_Eta;

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

  // lepton mva variables
    //! Lepton MVA input variables
  float myleptPt[3];
  float myleptEta[3];
  float myneuRelIs[3];
  float mychRelIso[3];
  float myjetDR_in[3];
  float myjetPtRatio_in[3];
  float myjetBTagCSV_in[3];
  float mysip3d[3];
  float mymvaId[3];
  int   myinnerHits[3];
  float mylogdxy[3];
  float mylogdz[3];

  TFile* myFile;
  TTree* myTree;
 
};

#endif // RedTopHiggsTree_h
