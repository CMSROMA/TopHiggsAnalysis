//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 22 21:45:43 2013 by ROOT version 5.30/02
// from TTree latino/latino tree
// found on file: tH125q_blvu_Yt1_H126toWW_datasetEE.root
//////////////////////////////////////////////////////////

#ifndef jetIdStudy_h
#define jetIdStudy_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class jetIdStudy {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         puweight;
   Bool_t          hlt;
   Float_t         met;
   Float_t         pfMet;
   Float_t         caloMet;
   Float_t         projMet;
   Float_t         deltaPhi;
   Float_t         deltaR;
   Float_t         transvMass;
   Float_t         eleInvMass;
   Float_t         maxPtEle;
   Float_t         minPtEle;
   Float_t         maxEtaEle;
   Float_t         minEtaEle;
   Float_t         detaLeptons;
   Int_t           nVtx;
   Bool_t          finalLeptons;
   Bool_t          jetVeto;
   Bool_t          uncorrJetVeto;
   Bool_t          preDeltaPhi;
   Bool_t          finalSelection;
   Bool_t          promptDecay;
   Float_t         genmll;
   Float_t         genptll;
   Float_t         genyll;
   Double_t        cteq66W[100];
   Double_t        mstwW[100];
   Double_t        nnpdfW[101];
   Float_t         genHiggsPt;
   Float_t         genHiggsEta;
   Float_t         genTopPt;
   Float_t         genTopEta;
   Float_t         genWpfromH_Pt;
   Float_t         genWpfromH_Eta;
   Float_t         genWmfromH_Pt;
   Float_t         genWmfromH_Eta;
   Float_t         genWfromT_Pt;
   Float_t         genWfromT_Eta;
   Float_t         genLeptonPlusfromWfromH_Pt;
   Float_t         genLeptonPlusfromWfromH_Eta;
   Float_t         genLeptonMinusfromWfromH_Pt;
   Float_t         genLeptonMinusfromWfromH_Eta;
   Float_t         genLeptonfromWfromT_Pt;
   Float_t         genLeptonfromWfromT_Eta;
   Float_t         genNeutrinoPlusfromWfromH_Pt;
   Float_t         genNeutrinoPlusfromWfromH_Eta;
   Float_t         genNeutrinoMinusfromWfromH_Pt;
   Float_t         genNeutrinoMinusfromWfromH_Eta;
   Float_t         genNeutrinofromWfromT_Pt;
   Float_t         genNeutrinofromWfromT_Eta;
   Float_t         genForwardQuark_Pt;
   Float_t         genForwardQuark_Eta;
   Float_t         genbQuark_Pt;
   Float_t         genbQuark_Eta;
   Int_t           run;
   Int_t           lumi;
   ULong64_t       event;
   Float_t         npu[3];
   Int_t           njets;
   Int_t           ncbIDjets;
   Int_t           nuncorrjets;
   Float_t         dxyEVT;
   Float_t         dszEVT;
   Float_t         softtche;
   Float_t         hardbjpb;
   Float_t         bTagSecVertex;
   Int_t           nSoftMu;
   Int_t           nSoftMuNoJets;
   Float_t         leadingJetBTagTrackCount;
   Float_t         subleadingJetBTagTrackCount;
   Float_t         subleadingJetsMaxBTagTrackCount;
   Float_t         leadingJetBTagJetBProb;
   Float_t         subleadingJetBTagJetBProb;
   Float_t         subleadingJetsMaxBTagJetBProb;
   Int_t           numExtraLep;
   Int_t           nsoftjet;
   Int_t           nsoftbjet;
   Int_t           numbtagCSVLcbIDaccepjets;
   Int_t           numbtagCSVMcbIDaccepjets;
   Int_t           numbtagCSVTcbIDaccepjets;
   Int_t           numcbIDcentralLjets;
   Int_t           numcbIDforwardLjets;
   Int_t           numcbIDcentralMjets;
   Int_t           numcbIDforwardMjets;
   Int_t           numbtagCSVLmvaIDaccepjets;
   Int_t           numbtagCSVMmvaIDaccepjets;
   Int_t           numbtagCSVTmvaIDaccepjets;
   Int_t           numbtagCSVLmvaIDcentraljets;
   Int_t           numbtagCSVMmvaIDcentraljets;
   Int_t           numbtagCSVTmvaIDcentraljets;
   Int_t           nummvaIDcentraljets;
   Int_t           nummvaIDforwardjets;
   Int_t           nummvaIDaccepINjets;
   Int_t           nummvaIDaccepOUTjets;
   Float_t         jetPuId_eta[20];
   Float_t         jetPuId_phi[20];
   Float_t         jetPuId_energy[20];
   Float_t         jetPuId_pt[20];
   Float_t         jetPuId_csv[20];
   Float_t         jetPuId_betastar[20];
   Float_t         jetPuId_rms[20];
   Int_t           jetPuId_cutBased[20];
   Float_t         jetPuId_mvaBased[20];
   Int_t           jetPuId_associated[20];
   Float_t         ptCVStaggedL[5];
   Float_t         etaCVStaggedL[5];
   Float_t         phiCVStaggedL[5];
   Float_t         eneCVStaggedL[5];
   Float_t         cvsCVStaggedL[5];
   Float_t         ptForwardL[5];
   Float_t         etaForwardL[5];
   Float_t         phiForwardL[5];
   Float_t         eneForwardL[5];
   Float_t         cvsForwardL[5];
   Float_t         ptCentralL[5];
   Float_t         etaCentralL[5];
   Float_t         phiCentralL[5];
   Float_t         eneCentralL[5];
   Float_t         cvsCentralL[5];
   Float_t         ptCVStaggedM[5];
   Float_t         etaCVStaggedM[5];
   Float_t         phiCVStaggedM[5];
   Float_t         eneCVStaggedM[5];
   Float_t         cvsCVStaggedM[5];
   Float_t         ptForwardM[5];
   Float_t         etaForwardM[5];
   Float_t         phiForwardM[5];
   Float_t         eneForwardM[5];
   Float_t         cvsForwardM[5];
   Float_t         ptCentralM[5];
   Float_t         etaCentralM[5];
   Float_t         phiCentralM[5];
   Float_t         eneCentralM[5];
   Float_t         cvsCentralM[5];
   Bool_t          step[29];
   Float_t         pxChMet;
   Float_t         pyChMet;
   Float_t         pzChMet;
   Float_t         pxLeadJet[3];
   Float_t         pyLeadJet[3];
   Float_t         pzLeadJet[3];
   Float_t         pxSecondJet[3];
   Float_t         pySecondJet[3];
   Float_t         pzSecondJet[3];
   Float_t         pxL1;
   Float_t         pyL1;
   Float_t         pzL1;
   Float_t         pxL2;
   Float_t         pyL2;
   Float_t         pzL2;
   Float_t         pxL3;
   Float_t         pyL3;
   Float_t         pzL3;
   Int_t           ch[3];
   Float_t         bdt[3];
   Float_t         pt[3];
   Float_t         eta[3];
   Float_t         phi[3];
   Int_t           flavour[3];
   Float_t         lepid[3];
   Float_t         lepiso[3];
   Float_t         lepconv[3];
   Float_t         scEnergy[3];
   Float_t         R9[3];
   Float_t         eneL1;
   Float_t         eneL2;
   Float_t         eneL3;
   Int_t           typeL1;
   Int_t           typeL2;
   Int_t           typeL3;
   Float_t         transvMassUp;
   Float_t         transvMassDown;
   Float_t         projPFMet;
   Float_t         projPFChargedMet;
   Float_t         signPFMet;
   Float_t         signPFChargedMet;
   Float_t         mtrchargedMet;
   Float_t         dymva1;
   Float_t         rho;
   Float_t         rhoJet;
   Float_t         leptPt[3];
   Float_t         leptEta[3];
   Float_t         neuRelIs[3];
   Float_t         chRelIso[3];
   Float_t         jetDR_in[3];
   Float_t         jetPtRatio_in[3];
   Float_t         jetBTagCSV_in[3];
   Float_t         sip3d[3];
   Float_t         mvaId[3];
   Int_t           innerHits[3];
   Float_t         logdxy[3];
   Float_t         logdz[3];

   // List of branches
   TBranch        *b_puweight;   //!
   TBranch        *b_hlt;   //!
   TBranch        *b_met;   //!
   TBranch        *b_pfMet;   //!
   TBranch        *b_caloMet;   //!
   TBranch        *b_projMet;   //!
   TBranch        *b_deltaPhi;   //!
   TBranch        *b_deltaR;   //!
   TBranch        *b_transvMass;   //!
   TBranch        *b_eleInvMass;   //!
   TBranch        *b_maxPtEle;   //!
   TBranch        *b_minPtEle;   //!
   TBranch        *b_maxEtaEle;   //!
   TBranch        *b_minEtaEle;   //!
   TBranch        *b_detaLeptons;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_finalLeptons;   //!
   TBranch        *b_jetVeto;   //!
   TBranch        *b_uncorrJetVeto;   //!
   TBranch        *b_preDeltaPhi;   //!
   TBranch        *b_finalSelection;   //!
   TBranch        *b_promptDecay;   //!
   TBranch        *b_genmll;   //!
   TBranch        *b_genptll;   //!
   TBranch        *b_genyll;   //!
   TBranch        *b_cteq66W;   //!
   TBranch        *b_mstwW;   //!
   TBranch        *b_nnpdfW;   //!
   TBranch        *b_genHiggsPt;   //!
   TBranch        *b_genHiggsEta;   //!
   TBranch        *b_genTopPt;   //!
   TBranch        *b_genTopEta;   //!
   TBranch        *b_genWpfromH_Pt;   //!
   TBranch        *b_genWpfromH_Eta;   //!
   TBranch        *b_genWmfromH_Pt;   //!
   TBranch        *b_genWmfromH_Eta;   //!
   TBranch        *b_genWfromT_Pt;   //!
   TBranch        *b_genWfromT_Eta;   //!
   TBranch        *b_genLeptonPlusfromWfromH_Pt;   //!
   TBranch        *b_genLeptonPlusfromWfromH_Eta;   //!
   TBranch        *b_genLeptonMinusfromWfromH_Pt;   //!
   TBranch        *b_genLeptonMinusfromWfromH_Eta;   //!
   TBranch        *b_genLeptonfromWfromT_Pt;   //!
   TBranch        *b_genLeptonfromWfromT_Eta;   //!
   TBranch        *b_genNeutrinoPlusfromWfromH_Pt;   //!
   TBranch        *b_genNeutrinoPlusfromWfromH_Eta;   //!
   TBranch        *b_genNeutrinoMinusfromWfromH_Pt;   //!
   TBranch        *b_genNeutrinoMinusfromWfromH_Eta;   //!
   TBranch        *b_genNeutrinofromWfromT_Pt;   //!
   TBranch        *b_genNeutrinofromWfromT_Eta;   //!
   TBranch        *b_genForwardQuark_Pt;   //!
   TBranch        *b_genForwardQuark_Eta;   //!
   TBranch        *b_genbQuark_Pt;   //!
   TBranch        *b_genbQuark_Eta;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_ncbIDjets;   //!
   TBranch        *b_nuncorrjets;   //!
   TBranch        *b_dxyEVT;   //!
   TBranch        *b_dszEVT;   //!
   TBranch        *b_softtcche;   //!
   TBranch        *b_hardbjbp;   //!
   TBranch        *b_bTagSecVertex;   //!
   TBranch        *b_nSoftMu;   //!
   TBranch        *b_nSoftMuNoJets;   //!
   TBranch        *b_leadingJetBTagTrackCount;   //!
   TBranch        *b_subleadingJetBTagTrackCount;   //!
   TBranch        *b_subleadingJetsMaxBTagTrackCount;   //!
   TBranch        *b_leadingJetBTagJetBProb;   //!
   TBranch        *b_subleadingJetBTagJetBProb;   //!
   TBranch        *b_subleadingJetsMaxBTagJetBProb;   //!
   TBranch        *b_numExtraLep;   //!
   TBranch        *b_nsoftjet;   //!
   TBranch        *b_nsoftbjet;   //!
   TBranch        *b_numbtagCSVLcbIDaccepjets;   //!
   TBranch        *b_numbtagCSVMcbIDaccepjets;   //!
   TBranch        *b_numbtagCSVTcbIDaccepjets;   //!
   TBranch        *b_numcbIDcentralLjets;   //!
   TBranch        *b_numcbIDforwardLjets;   //!
   TBranch        *b_numcbIDcentralMjets;   //!
   TBranch        *b_numcbIDforwardMjets;   //!
   TBranch        *b_numbtagCSVLmvaIDaccepjets;   //!
   TBranch        *b_numbtagCSVMmvaIDaccepjets;   //!
   TBranch        *b_numbtagCSVTmvaIDaccepjets;   //!
   TBranch        *b_numbtagCSVLmvaIDcentraljets;   //!
   TBranch        *b_numbtagCSVMmvaIDcentraljets;   //!
   TBranch        *b_numbtagCSVTmvaIDcentraljets;   //!
   TBranch        *b_nummvaIDcentraljets;   //!
   TBranch        *b_nummvaIDforwardjets;   //!
   TBranch        *b_nummvaIDaccepINjets;   //!
   TBranch        *b_nummvaIDaccepOUTjets;   //!
   TBranch        *b_jetPuId_eta;   //!
   TBranch        *b_jetPuId_phi;   //!
   TBranch        *b_jetPuId_energy;   //!
   TBranch        *b_jetPuId_pt;   //!
   TBranch        *b_jetPuId_csv;   //!
   TBranch        *b_jetPuId_betastar;   //!
   TBranch        *b_jetPuId_rms;   //!
   TBranch        *b_jetPuId_cutBased;   //!
   TBranch        *b_jetPuId_mvaBased;   //!
   TBranch        *b_jetPuId_associated;   //!
   TBranch        *b_ptCVStaggedL;   //!
   TBranch        *b_etaCVStaggedL;   //!
   TBranch        *b_phiCVStaggedL;   //!
   TBranch        *b_eneCVStaggedL;   //!
   TBranch        *b_cvsCVStaggedL;   //!
   TBranch        *b_ptForwardL;   //!
   TBranch        *b_etaForwardL;   //!
   TBranch        *b_phiForwardL;   //!
   TBranch        *b_eneForwardL;   //!
   TBranch        *b_cvsForwardL;   //!
   TBranch        *b_ptCentralL;   //!
   TBranch        *b_etaCentralL;   //!
   TBranch        *b_phiCentralL;   //!
   TBranch        *b_eneCentralL;   //!
   TBranch        *b_cvsCentralL;   //!
   TBranch        *b_ptCVStaggedM;   //!
   TBranch        *b_etaCVStaggedM;   //!
   TBranch        *b_phiCVStaggedM;   //!
   TBranch        *b_eneCVStaggedM;   //!
   TBranch        *b_cvsCVStaggedM;   //!
   TBranch        *b_ptForwardM;   //!
   TBranch        *b_etaForwardM;   //!
   TBranch        *b_phiForwardM;   //!
   TBranch        *b_eneForwardM;   //!
   TBranch        *b_cvsForwardM;   //!
   TBranch        *b_ptCentralM;   //!
   TBranch        *b_etaCentralM;   //!
   TBranch        *b_phiCentralM;   //!
   TBranch        *b_eneCentralM;   //!
   TBranch        *b_cvsCentralM;   //!
   TBranch        *b_step;   //!
   TBranch        *b_pxChMet;   //!
   TBranch        *b_pyChMet;   //!
   TBranch        *b_pzChMet;   //!
   TBranch        *b_pxLeadJet;   //!
   TBranch        *b_pyLeadJet;   //!
   TBranch        *b_pzLeadJet;   //!
   TBranch        *b_pxSecondJet;   //!
   TBranch        *b_pySecondJet;   //!
   TBranch        *b_pzSecondJet;   //!
   TBranch        *b_pxL1;   //!
   TBranch        *b_pyL1;   //!
   TBranch        *b_pzL1;   //!
   TBranch        *b_pxL2;   //!
   TBranch        *b_pyL2;   //!
   TBranch        *b_pzL2;   //!
   TBranch        *b_pxL3;   //!
   TBranch        *b_pyL3;   //!
   TBranch        *b_pzL3;   //!
   TBranch        *b_ch;   //!
   TBranch        *b_bdt;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_eta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_flavour;   //!
   TBranch        *b_lepid;   //!
   TBranch        *b_lepiso;   //!
   TBranch        *b_lepconv;   //!
   TBranch        *b_scEnergy;   //!
   TBranch        *b_R9;   //!
   TBranch        *b_eneL1;   //!
   TBranch        *b_eneL2;   //!
   TBranch        *b_eneL3;   //!
   TBranch        *b_typeL1;   //!
   TBranch        *b_typeL2;   //!
   TBranch        *b_typeL3;   //!
   TBranch        *b_transvMassUp;   //!
   TBranch        *b_transvMassDown;   //!
   TBranch        *b_projPFMet;   //!
   TBranch        *b_projPFChargedMet;   //!
   TBranch        *b_signPFMet;   //!
   TBranch        *b_signPFChargedMet;   //!
   TBranch        *b_mtrchargedMet;   //!
   TBranch        *b_dymva1;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoJet;   //!
   TBranch        *b_leptPt;   //!
   TBranch        *b_leptEta;   //!
   TBranch        *b_neuRelIs;   //!
   TBranch        *b_chRelIso;   //!
   TBranch        *b_jetDR_in;   //!
   TBranch        *b_jetPtRatio_in;   //!
   TBranch        *b_jetBTagCSV_in;   //!
   TBranch        *b_sip3d;   //!
   TBranch        *b_mvaId;   //!
   TBranch        *b_innerHits;   //!
   TBranch        *b_logdxy;   //!
   TBranch        *b_logdz;   //!

   jetIdStudy(TTree *tree=0);
   virtual ~jetIdStudy();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef jetIdStudy_cxx
jetIdStudy::jetIdStudy(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/cmsrm/pc22/crovelli/data/topH/testJetV2/MC/tH125q_blvu_Yt1_H126toWW/merged/tH125q_blvu_Yt1_H126toWW_allDatasets.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/cmsrm/pc22/crovelli/data/topH/testJetV2/MC/tH125q_blvu_Yt1_H126toWW/merged/tH125q_blvu_Yt1_H126toWW_allDatasets.root");
      }
      f->GetObject("latino",tree);

   }
   Init(tree);
}

jetIdStudy::~jetIdStudy()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t jetIdStudy::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t jetIdStudy::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void jetIdStudy::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("hlt", &hlt, &b_hlt);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("pfMet", &pfMet, &b_pfMet);
   fChain->SetBranchAddress("caloMet", &caloMet, &b_caloMet);
   fChain->SetBranchAddress("projMet", &projMet, &b_projMet);
   fChain->SetBranchAddress("deltaPhi", &deltaPhi, &b_deltaPhi);
   fChain->SetBranchAddress("deltaR", &deltaR, &b_deltaR);
   fChain->SetBranchAddress("transvMass", &transvMass, &b_transvMass);
   fChain->SetBranchAddress("eleInvMass", &eleInvMass, &b_eleInvMass);
   fChain->SetBranchAddress("maxPtEle", &maxPtEle, &b_maxPtEle);
   fChain->SetBranchAddress("minPtEle", &minPtEle, &b_minPtEle);
   fChain->SetBranchAddress("maxEtaEle", &maxEtaEle, &b_maxEtaEle);
   fChain->SetBranchAddress("minEtaEle", &minEtaEle, &b_minEtaEle);
   fChain->SetBranchAddress("detaLeptons", &detaLeptons, &b_detaLeptons);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("finalLeptons", &finalLeptons, &b_finalLeptons);
   fChain->SetBranchAddress("jetVeto", &jetVeto, &b_jetVeto);
   fChain->SetBranchAddress("uncorrJetVeto", &uncorrJetVeto, &b_uncorrJetVeto);
   fChain->SetBranchAddress("preDeltaPhi", &preDeltaPhi, &b_preDeltaPhi);
   fChain->SetBranchAddress("finalSelection", &finalSelection, &b_finalSelection);
   fChain->SetBranchAddress("promptDecay", &promptDecay, &b_promptDecay);
   fChain->SetBranchAddress("genmll", &genmll, &b_genmll);
   fChain->SetBranchAddress("genptll", &genptll, &b_genptll);
   fChain->SetBranchAddress("genyll", &genyll, &b_genyll);
   fChain->SetBranchAddress("cteq66W", cteq66W, &b_cteq66W);
   fChain->SetBranchAddress("mstwW", mstwW, &b_mstwW);
   fChain->SetBranchAddress("nnpdfW", nnpdfW, &b_nnpdfW);
   fChain->SetBranchAddress("genHiggsPt", &genHiggsPt, &b_genHiggsPt);
   fChain->SetBranchAddress("genHiggsEta", &genHiggsEta, &b_genHiggsEta);
   fChain->SetBranchAddress("genTopPt", &genTopPt, &b_genTopPt);
   fChain->SetBranchAddress("genTopEta", &genTopEta, &b_genTopEta);
   fChain->SetBranchAddress("genWpfromH_Pt", &genWpfromH_Pt, &b_genWpfromH_Pt);
   fChain->SetBranchAddress("genWpfromH_Eta", &genWpfromH_Eta, &b_genWpfromH_Eta);
   fChain->SetBranchAddress("genWmfromH_Pt", &genWmfromH_Pt, &b_genWmfromH_Pt);
   fChain->SetBranchAddress("genWmfromH_Eta", &genWmfromH_Eta, &b_genWmfromH_Eta);
   fChain->SetBranchAddress("genWfromT_Pt", &genWfromT_Pt, &b_genWfromT_Pt);
   fChain->SetBranchAddress("genWfromT_Eta", &genWfromT_Eta, &b_genWfromT_Eta);
   fChain->SetBranchAddress("genLeptonPlusfromWfromH_Pt", &genLeptonPlusfromWfromH_Pt, &b_genLeptonPlusfromWfromH_Pt);
   fChain->SetBranchAddress("genLeptonPlusfromWfromH_Eta", &genLeptonPlusfromWfromH_Eta, &b_genLeptonPlusfromWfromH_Eta);
   fChain->SetBranchAddress("genLeptonMinusfromWfromH_Pt", &genLeptonMinusfromWfromH_Pt, &b_genLeptonMinusfromWfromH_Pt);
   fChain->SetBranchAddress("genLeptonMinusfromWfromH_Eta", &genLeptonMinusfromWfromH_Eta, &b_genLeptonMinusfromWfromH_Eta);
   fChain->SetBranchAddress("genLeptonfromWfromT_Pt", &genLeptonfromWfromT_Pt, &b_genLeptonfromWfromT_Pt);
   fChain->SetBranchAddress("genLeptonfromWfromT_Eta", &genLeptonfromWfromT_Eta, &b_genLeptonfromWfromT_Eta);
   fChain->SetBranchAddress("genNeutrinoPlusfromWfromH_Pt", &genNeutrinoPlusfromWfromH_Pt, &b_genNeutrinoPlusfromWfromH_Pt);
   fChain->SetBranchAddress("genNeutrinoPlusfromWfromH_Eta", &genNeutrinoPlusfromWfromH_Eta, &b_genNeutrinoPlusfromWfromH_Eta);
   fChain->SetBranchAddress("genNeutrinoMinusfromWfromH_Pt", &genNeutrinoMinusfromWfromH_Pt, &b_genNeutrinoMinusfromWfromH_Pt);
   fChain->SetBranchAddress("genNeutrinoMinusfromWfromH_Eta", &genNeutrinoMinusfromWfromH_Eta, &b_genNeutrinoMinusfromWfromH_Eta);
   fChain->SetBranchAddress("genNeutrinofromWfromT_Pt", &genNeutrinofromWfromT_Pt, &b_genNeutrinofromWfromT_Pt);
   fChain->SetBranchAddress("genNeutrinofromWfromT_Eta", &genNeutrinofromWfromT_Eta, &b_genNeutrinofromWfromT_Eta);
   fChain->SetBranchAddress("genForwardQuark_Pt", &genForwardQuark_Pt, &b_genForwardQuark_Pt);
   fChain->SetBranchAddress("genForwardQuark_Eta", &genForwardQuark_Eta, &b_genForwardQuark_Eta);
   fChain->SetBranchAddress("genbQuark_Pt", &genbQuark_Pt, &b_genbQuark_Pt);
   fChain->SetBranchAddress("genbQuark_Eta", &genbQuark_Eta, &b_genbQuark_Eta);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("npu", npu, &b_npu);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("ncbIDjets", &ncbIDjets, &b_ncbIDjets);
   fChain->SetBranchAddress("nuncorrjets", &nuncorrjets, &b_nuncorrjets);
   fChain->SetBranchAddress("dxyEVT", &dxyEVT, &b_dxyEVT);
   fChain->SetBranchAddress("dszEVT", &dszEVT, &b_dszEVT);
   fChain->SetBranchAddress("softtche", &softtche, &b_softtcche);
   fChain->SetBranchAddress("hardbjpb", &hardbjpb, &b_hardbjbp);
   fChain->SetBranchAddress("bTagSecVertex", &bTagSecVertex, &b_bTagSecVertex);
   fChain->SetBranchAddress("nSoftMu", &nSoftMu, &b_nSoftMu);
   fChain->SetBranchAddress("nSoftMuNoJets", &nSoftMuNoJets, &b_nSoftMuNoJets);
   fChain->SetBranchAddress("leadingJetBTagTrackCount", &leadingJetBTagTrackCount, &b_leadingJetBTagTrackCount);
   fChain->SetBranchAddress("subleadingJetBTagTrackCount", &subleadingJetBTagTrackCount, &b_subleadingJetBTagTrackCount);
   fChain->SetBranchAddress("subleadingJetsMaxBTagTrackCount", &subleadingJetsMaxBTagTrackCount, &b_subleadingJetsMaxBTagTrackCount);
   fChain->SetBranchAddress("leadingJetBTagJetBProb", &leadingJetBTagJetBProb, &b_leadingJetBTagJetBProb);
   fChain->SetBranchAddress("subleadingJetBTagJetBProb", &subleadingJetBTagJetBProb, &b_subleadingJetBTagJetBProb);
   fChain->SetBranchAddress("subleadingJetsMaxBTagJetBProb", &subleadingJetsMaxBTagJetBProb, &b_subleadingJetsMaxBTagJetBProb);
   fChain->SetBranchAddress("numExtraLep", &numExtraLep, &b_numExtraLep);
   fChain->SetBranchAddress("nsoftjet", &nsoftjet, &b_nsoftjet);
   fChain->SetBranchAddress("nsoftbjet", &nsoftbjet, &b_nsoftbjet);
   fChain->SetBranchAddress("numbtagCSVLcbIDaccepjets", &numbtagCSVLcbIDaccepjets, &b_numbtagCSVLcbIDaccepjets);
   fChain->SetBranchAddress("numbtagCSVMcbIDaccepjets", &numbtagCSVMcbIDaccepjets, &b_numbtagCSVMcbIDaccepjets);
   fChain->SetBranchAddress("numbtagCSVTcbIDaccepjets", &numbtagCSVTcbIDaccepjets, &b_numbtagCSVTcbIDaccepjets);
   fChain->SetBranchAddress("numcbIDcentralLjets", &numcbIDcentralLjets, &b_numcbIDcentralLjets);
   fChain->SetBranchAddress("numcbIDforwardLjets", &numcbIDforwardLjets, &b_numcbIDforwardLjets);
   fChain->SetBranchAddress("numcbIDcentralMjets", &numcbIDcentralMjets, &b_numcbIDcentralMjets);
   fChain->SetBranchAddress("numcbIDforwardMjets", &numcbIDforwardMjets, &b_numcbIDforwardMjets);
   fChain->SetBranchAddress("numbtagCSVLmvaIDaccepjets", &numbtagCSVLmvaIDaccepjets, &b_numbtagCSVLmvaIDaccepjets);
   fChain->SetBranchAddress("numbtagCSVMmvaIDaccepjets", &numbtagCSVMmvaIDaccepjets, &b_numbtagCSVMmvaIDaccepjets);
   fChain->SetBranchAddress("numbtagCSVTmvaIDaccepjets", &numbtagCSVTmvaIDaccepjets, &b_numbtagCSVTmvaIDaccepjets);
   fChain->SetBranchAddress("numbtagCSVLmvaIDcentraljets", &numbtagCSVLmvaIDcentraljets, &b_numbtagCSVLmvaIDcentraljets);
   fChain->SetBranchAddress("numbtagCSVMmvaIDcentraljets", &numbtagCSVMmvaIDcentraljets, &b_numbtagCSVMmvaIDcentraljets);
   fChain->SetBranchAddress("numbtagCSVTmvaIDcentraljets", &numbtagCSVTmvaIDcentraljets, &b_numbtagCSVTmvaIDcentraljets);
   fChain->SetBranchAddress("nummvaIDcentraljets", &nummvaIDcentraljets, &b_nummvaIDcentraljets);
   fChain->SetBranchAddress("nummvaIDforwardjets", &nummvaIDforwardjets, &b_nummvaIDforwardjets);
   fChain->SetBranchAddress("nummvaIDaccepINjets", &nummvaIDaccepINjets, &b_nummvaIDaccepINjets);
   fChain->SetBranchAddress("nummvaIDaccepOUTjets", &nummvaIDaccepOUTjets, &b_nummvaIDaccepOUTjets);
   fChain->SetBranchAddress("jetPuId_eta", jetPuId_eta, &b_jetPuId_eta);
   fChain->SetBranchAddress("jetPuId_phi", jetPuId_phi, &b_jetPuId_phi);
   fChain->SetBranchAddress("jetPuId_energy", jetPuId_energy, &b_jetPuId_energy);
   fChain->SetBranchAddress("jetPuId_pt", jetPuId_pt, &b_jetPuId_pt);
   fChain->SetBranchAddress("jetPuId_csv", jetPuId_csv, &b_jetPuId_csv);
   fChain->SetBranchAddress("jetPuId_betastar", jetPuId_betastar, &b_jetPuId_betastar);
   fChain->SetBranchAddress("jetPuId_rms", jetPuId_rms, &b_jetPuId_rms);
   fChain->SetBranchAddress("jetPuId_cutBased", jetPuId_cutBased, &b_jetPuId_cutBased);
   fChain->SetBranchAddress("jetPuId_mvaBased", jetPuId_mvaBased, &b_jetPuId_mvaBased);
   fChain->SetBranchAddress("jetPuId_associated", jetPuId_associated, &b_jetPuId_associated);
   fChain->SetBranchAddress("ptCVStaggedL", ptCVStaggedL, &b_ptCVStaggedL);
   fChain->SetBranchAddress("etaCVStaggedL", etaCVStaggedL, &b_etaCVStaggedL);
   fChain->SetBranchAddress("phiCVStaggedL", phiCVStaggedL, &b_phiCVStaggedL);
   fChain->SetBranchAddress("eneCVStaggedL", eneCVStaggedL, &b_eneCVStaggedL);
   fChain->SetBranchAddress("cvsCVStaggedL", cvsCVStaggedL, &b_cvsCVStaggedL);
   fChain->SetBranchAddress("ptForwardL", ptForwardL, &b_ptForwardL);
   fChain->SetBranchAddress("etaForwardL", etaForwardL, &b_etaForwardL);
   fChain->SetBranchAddress("phiForwardL", phiForwardL, &b_phiForwardL);
   fChain->SetBranchAddress("eneForwardL", eneForwardL, &b_eneForwardL);
   fChain->SetBranchAddress("cvsForwardL", cvsForwardL, &b_cvsForwardL);
   fChain->SetBranchAddress("ptCentralL", ptCentralL, &b_ptCentralL);
   fChain->SetBranchAddress("etaCentralL", etaCentralL, &b_etaCentralL);
   fChain->SetBranchAddress("phiCentralL", phiCentralL, &b_phiCentralL);
   fChain->SetBranchAddress("eneCentralL", eneCentralL, &b_eneCentralL);
   fChain->SetBranchAddress("cvsCentralL", cvsCentralL, &b_cvsCentralL);
   fChain->SetBranchAddress("ptCVStaggedM", ptCVStaggedM, &b_ptCVStaggedM);
   fChain->SetBranchAddress("etaCVStaggedM", etaCVStaggedM, &b_etaCVStaggedM);
   fChain->SetBranchAddress("phiCVStaggedM", phiCVStaggedM, &b_phiCVStaggedM);
   fChain->SetBranchAddress("eneCVStaggedM", eneCVStaggedM, &b_eneCVStaggedM);
   fChain->SetBranchAddress("cvsCVStaggedM", cvsCVStaggedM, &b_cvsCVStaggedM);
   fChain->SetBranchAddress("ptForwardM", ptForwardM, &b_ptForwardM);
   fChain->SetBranchAddress("etaForwardM", etaForwardM, &b_etaForwardM);
   fChain->SetBranchAddress("phiForwardM", phiForwardM, &b_phiForwardM);
   fChain->SetBranchAddress("eneForwardM", eneForwardM, &b_eneForwardM);
   fChain->SetBranchAddress("cvsForwardM", cvsForwardM, &b_cvsForwardM);
   fChain->SetBranchAddress("ptCentralM", ptCentralM, &b_ptCentralM);
   fChain->SetBranchAddress("etaCentralM", etaCentralM, &b_etaCentralM);
   fChain->SetBranchAddress("phiCentralM", phiCentralM, &b_phiCentralM);
   fChain->SetBranchAddress("eneCentralM", eneCentralM, &b_eneCentralM);
   fChain->SetBranchAddress("cvsCentralM", cvsCentralM, &b_cvsCentralM);
   fChain->SetBranchAddress("step", step, &b_step);
   fChain->SetBranchAddress("pxChMet", &pxChMet, &b_pxChMet);
   fChain->SetBranchAddress("pyChMet", &pyChMet, &b_pyChMet);
   fChain->SetBranchAddress("pzChMet", &pzChMet, &b_pzChMet);
   fChain->SetBranchAddress("pxLeadJet", pxLeadJet, &b_pxLeadJet);
   fChain->SetBranchAddress("pyLeadJet", pyLeadJet, &b_pyLeadJet);
   fChain->SetBranchAddress("pzLeadJet", pzLeadJet, &b_pzLeadJet);
   fChain->SetBranchAddress("pxSecondJet", pxSecondJet, &b_pxSecondJet);
   fChain->SetBranchAddress("pySecondJet", pySecondJet, &b_pySecondJet);
   fChain->SetBranchAddress("pzSecondJet", pzSecondJet, &b_pzSecondJet);
   fChain->SetBranchAddress("pxL1", &pxL1, &b_pxL1);
   fChain->SetBranchAddress("pyL1", &pyL1, &b_pyL1);
   fChain->SetBranchAddress("pzL1", &pzL1, &b_pzL1);
   fChain->SetBranchAddress("pxL2", &pxL2, &b_pxL2);
   fChain->SetBranchAddress("pyL2", &pyL2, &b_pyL2);
   fChain->SetBranchAddress("pzL2", &pzL2, &b_pzL2);
   fChain->SetBranchAddress("pxL3", &pxL3, &b_pxL3);
   fChain->SetBranchAddress("pyL3", &pyL3, &b_pyL3);
   fChain->SetBranchAddress("pzL3", &pzL3, &b_pzL3);
   fChain->SetBranchAddress("ch", ch, &b_ch);
   fChain->SetBranchAddress("bdt", bdt, &b_bdt);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("eta", eta, &b_eta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("flavour", flavour, &b_flavour);
   fChain->SetBranchAddress("lepid", lepid, &b_lepid);
   fChain->SetBranchAddress("lepiso", lepiso, &b_lepiso);
   fChain->SetBranchAddress("lepconv", lepconv, &b_lepconv);
   fChain->SetBranchAddress("scEnergy", scEnergy, &b_scEnergy);
   fChain->SetBranchAddress("R9", R9, &b_R9);
   fChain->SetBranchAddress("eneL1", &eneL1, &b_eneL1);
   fChain->SetBranchAddress("eneL2", &eneL2, &b_eneL2);
   fChain->SetBranchAddress("eneL3", &eneL3, &b_eneL3);
   fChain->SetBranchAddress("typeL1", &typeL1, &b_typeL1);
   fChain->SetBranchAddress("typeL2", &typeL2, &b_typeL2);
   fChain->SetBranchAddress("typeL3", &typeL3, &b_typeL3);
   fChain->SetBranchAddress("transvMassUp", &transvMassUp, &b_transvMassUp);
   fChain->SetBranchAddress("transvMassDown", &transvMassDown, &b_transvMassDown);
   fChain->SetBranchAddress("projPFMet", &projPFMet, &b_projPFMet);
   fChain->SetBranchAddress("projPFChargedMet", &projPFChargedMet, &b_projPFChargedMet);
   fChain->SetBranchAddress("signPFMet", &signPFMet, &b_signPFMet);
   fChain->SetBranchAddress("signPFChargedMet", &signPFChargedMet, &b_signPFChargedMet);
   fChain->SetBranchAddress("mtrchargedMet", &mtrchargedMet, &b_mtrchargedMet);
   fChain->SetBranchAddress("dymva1", &dymva1, &b_dymva1);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoJet", &rhoJet, &b_rhoJet);
   fChain->SetBranchAddress("leptPt", leptPt, &b_leptPt);
   fChain->SetBranchAddress("leptEta", leptEta, &b_leptEta);
   fChain->SetBranchAddress("neuRelIs", neuRelIs, &b_neuRelIs);
   fChain->SetBranchAddress("chRelIso", chRelIso, &b_chRelIso);
   fChain->SetBranchAddress("jetDR_in", jetDR_in, &b_jetDR_in);
   fChain->SetBranchAddress("jetPtRatio_in", jetPtRatio_in, &b_jetPtRatio_in);
   fChain->SetBranchAddress("jetBTagCSV_in", jetBTagCSV_in, &b_jetBTagCSV_in);
   fChain->SetBranchAddress("sip3d", sip3d, &b_sip3d);
   fChain->SetBranchAddress("mvaId", mvaId, &b_mvaId);
   fChain->SetBranchAddress("innerHits", innerHits, &b_innerHits);
   fChain->SetBranchAddress("logdxy", logdxy, &b_logdxy);
   fChain->SetBranchAddress("logdz", logdz, &b_logdz);
   Notify();
}

Bool_t jetIdStudy::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void jetIdStudy::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t jetIdStudy::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef jetIdStudy_cxx
