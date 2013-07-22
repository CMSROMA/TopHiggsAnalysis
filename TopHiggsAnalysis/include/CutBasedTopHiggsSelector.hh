
#ifndef CutBasedTopHiggsSelector_h
#define CutBasedTopHiggsSelector_h

#include "CommonTools/include/Selection.hh"
#include "CommonTools/include/Counters.hh"

class CutBasedTopHiggsSelector {

public:

  //! constructor
  CutBasedTopHiggsSelector();

  //! copy constructor
  CutBasedTopHiggsSelector( const CutBasedTopHiggsSelector& selector );

  //! destructor
  virtual ~CutBasedTopHiggsSelector();   

  //! configure from files
  void Configure(const char *fileCuts, const char* fileSwitches, const char *theTitle);

  //! configure pre-selection (it is not necessarily applied)
  void AppendPreselection(Selection *preselection) { _selection->append(preselection); }

  //! get the applied selection
  Selection* GetSelection() { return _selection; }

  //! set event by event observables
  void SetMcTruth(bool foundMcTree)   { m_foundMcTree = foundMcTree; }
  void SetHLT(bool passedHLT)         { m_passedHLT   = passedHLT;   }
  void SetIsChannel(bool channelsel)  { m_isThisChannel = channelsel; }
  void SetChannel(int channel)        { m_channel = channel; }

  void SetLepton1ID(int id1)          { m_lepton1id = id1; }
  void SetLepton2ID(int id2)          { m_lepton2id = id2; }
  void SetLepton3ID(int id3)          { m_lepton3id = id3; }

  void SetProcessID(int processID)       { m_processID     = processID; }
  void SetWeight(float weight)           { m_weight        = weight;    }
  void SetHighElePt(float highPt)        { m_highPt        = highPt;    }
  void SetLowElePt(float lowPt)          { m_lowPt         = lowPt;     }
  void SetInvMass(float mll)             { m_invMass       = mll;       }
  void SetElectronId(int isEleId)       { m_isElectronId  = isEleId; }
  void SetPositronId(int isPosId)       { m_isPositronId  = isPosId; }
  void SetElectronIsolation(int isEleIsol)   { m_isElectronIsol  = isEleIsol; }
  void SetPositronIsolation(int isPosIsol)   { m_isPositronIsol  = isPosIsol; }
  void SetElectronConvRejection(int isEleConvRej)   { m_isElectronConvRej  = isEleConvRej; }
  void SetPositronConvRejection(int isPosConvRej)   { m_isPositronConvRej  = isPosConvRej; }
  void SetElectronIp(int isEleIp)       { m_isElectronIp  = isEleIp; }
  void SetPositronIp(int isPosIp)       { m_isPositronIp  = isPosIp; }
  void SetJetVeto(bool passedCJV)        { m_passedJetVeto = passedCJV; }
  void SetUncorrJetVeto(bool passedUncorrCJV) { m_passedUncorrJetVeto = passedUncorrCJV; }
  void SetNJets(int njets)               { m_nJets           = njets; } // MVA ID for jets
  void SetNcbIDJets(int ncbjets)         { m_nCutBasedIDJets = ncbjets; } // Cut Based ID for jets
  void SetNUncorrJets(int nuncorrjets)   { m_nUncorrJets   = nuncorrjets; }
  void SetBTagJets(float btag)           { m_btagJets      = btag; }
  void SetNSoftMuons(int nsoftmu)        { m_nSoftMuons    = nsoftmu; }
  void SetNExtraLeptons(int nextralep)   { m_nExtraLeptons = nextralep; }
  void SetMet(float met)                 { m_met           = met;}
  void SetProjectedMet(float projmet)    { m_projectedMet  = projmet;}
  void SetMetOverPtLL(float metopt)      { m_metOverPtLL   = metopt;}
  void SetDYMVA(float dymva)             { m_dymva         = dymva; }
  void SetDeltaPhiLLJet(float deltaphiLLJet) { m_deltaPhiLLJet = deltaphiLLJet; }
  void SetDeltaPhi(float deltaPhi)       { m_deltaPhi      = deltaPhi;}
  void SetDetaLeptons(float deltaEta)    { m_detaLeptons   = deltaEta;}
  void SetPtll(float ptll) { m_ptll = ptll; }
  void SetExtraSlowLeptonPTCut(float ptMin) { m_extraSlowLeptonPTMin = ptMin; }
  void SetWWInvMass(float wwmass) { m_WWInvMass = wwmass; }
  void SetNvtx(int nvtx) { m_nVtx = nvtx; }

  //! get output of the selector
  bool output();

  //! latinos 
  bool outputStep0() { return m_step0; }
  bool outputStep1() { return m_step1; }

  //! display the electron efficiency
  void displayEfficiencies(std::string datasetName);

private:
  
  int m_channel; // the current channel {eee = 0, mmm = 1, eem = 2, mme = 3}
  int m_lepton1id; //id for lepton 1 {electron == 0, muon == 1}
  int m_lepton2id; //id for lepton 2
  int m_lepton3id; //id for lepton 3

  float m_weight;
  bool m_foundMcTree;
  bool m_passedHLT;
  bool m_isThisChannel;
  float m_highPt, m_lowPt;
  int m_isElectronId, m_isPositronId;
  int m_isElectronIsol, m_isPositronIsol;
  int m_isElectronConvRej, m_isPositronConvRej;
  int m_isElectronIp, m_isPositronIp;
  float m_invMass;
  bool m_passedJetVeto;
  bool m_passedUncorrJetVeto;
  int m_nJets, m_nCutBasedIDJets, m_nUncorrJets, m_nSoftMuons, m_nExtraLeptons;
  float m_btagJets;
  float m_met, m_projectedMet, m_deltaPhi, m_detaLeptons, m_metOverPtLL, m_deltaPhiLLJet, m_dymva;
  float m_maxPtElectron, m_minPtElectron, m_WWInvMass, m_ptll;
  float m_extraSlowLeptonPTMin;
  int m_nVtx;
  int m_processID;

  //! contains the preselection cuts
  Selection* _selection;

  //! counters for the efficiencies display, based on electron candidates
  Counters* globalCounter;

  //! Different cut steps
  bool m_step0, m_step1;  

  //! this is to do an efficiency for each process in the sample 
  //! (if more than one is present)
  //! to turn on it, use SetProcessID(int processID) with processID=!-1
  std::map<int, Counters*> multiProcessCounter;

};

#endif
