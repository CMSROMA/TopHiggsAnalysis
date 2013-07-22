#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

// Number of MC samples we want to check
#define NSAMPLES 19

using namespace std;

double weight(double ngen, double xsec, double filtereff, double lumi = 1);

void countEvents() {

  char nametree[500];
  sprintf(nametree,"FULL_SELECTION_EVENT_COUNTER_EEE");

  cout << "nametree = " << nametree << endl;

  // backgrounds
  TChain *chains[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    chains[isample] = new TChain(nametree);
  }

  TString pathname = "/cmsrm/pc24_2/jorda/data/results_V05/";

  // nominal sample first, then the systematics ones
  chains[0]->Add(pathname+"MC/tH125q_blvu_Yt1_H126toWW-Counters.root");
  chains[1]->Add(pathname+"MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-Counters.root");
  chains[2]->Add(pathname+"MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-Counters.root");
  chains[3]->Add(pathname+"MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-Counters.root");
  chains[4]->Add(pathname+"MC/TTWJets-Counters.root");
  chains[5]->Add(pathname+"MC/TTZJets-Counters.root");
  chains[6]->Add(pathname+"MC/WWWJets_8TeV-madgraph-Counters.root");
  chains[7]->Add(pathname+"MC/WWZNoGstarJets_8TeV-madgraph-Counters.root");
  chains[8]->Add(pathname+"MC/WZZNoGstarJets_8TeV-madgraph-Counters.root");
  chains[9]->Add(pathname+"MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-Counters.root");
  chains[10]->Add(pathname+"MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[11]->Add(pathname+"MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[12]->Add(pathname+"MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[13]->Add(pathname+"MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[14]->Add(pathname+"MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[15]->Add(pathname+"MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-Counters.root");
  chains[16]->Add(pathname+"MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-Counters.root");
  chains[17]->Add(pathname+"MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-Counters.root");
  chains[18]->Add(pathname+"MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-Counters.root");


  cout << "chains added. " << endl;

  std::vector<TString> sampleName;
  sampleName.push_back(pathname+"MC/tH125q_blvu_Yt1_H126toWW-datasetEEE.root"); // 0
  sampleName.push_back(pathname+"MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root"); // 1
  sampleName.push_back(pathname+"MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetEEE.root"); // 2
  sampleName.push_back(pathname+"MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root"); // 3
  sampleName.push_back(pathname+"MC/TTWJets-datasetEEE.root"); // 4
  sampleName.push_back(pathname+"MC/TTZJets-datasetEEE.root"); // 5
  sampleName.push_back(pathname+"MC/WWWJets_8TeV-madgraph-datasetEEE.root"); // 6
  sampleName.push_back(pathname+"MC/WWZNoGstarJets_8TeV-madgraph-datasetEEE.root"); // 7
  sampleName.push_back(pathname+"MC/WZZNoGstarJets_8TeV-madgraph-datasetEEE.root"); // 8
  sampleName.push_back(pathname+"MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetEEE.root"); // 9
  sampleName.push_back(pathname+"MC/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 10 
  sampleName.push_back(pathname+"MC/T_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 11
  sampleName.push_back(pathname+"MC/T_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 12
  sampleName.push_back(pathname+"MC/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 13
  sampleName.push_back(pathname+"MC/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 14
  sampleName.push_back(pathname+"MC/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola-datasetEEE.root"); // 15
  sampleName.push_back(pathname+"MC/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-datasetEEE.root"); // 16
  sampleName.push_back(pathname+"MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetEEE.root"); // 17
  sampleName.push_back(pathname+"MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetEEE.root"); // 18

  // xs from: https://docs.google.com/spreadsheet/ccc?key=0Aq5OAopf_dtsdDJqREJReGQyY21wbERldVFSZVJHbFE&hl=en#gid=56
  std::vector<float> sampleXsec;
  sampleXsec.push_back(0.0010449);      // 0: tH125q_blvu_Yt1_H126toWW xs = 15 fb * (BR W->ln)*(BR H->WW) = 15 fb*(3*0.108)*(0.215)
  sampleXsec.push_back(5.8123);         // 1: WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola
  sampleXsec.push_back(1.057512972672); // 2: WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola2
  sampleXsec.push_back(0.364868);       // 3: ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola  
  sampleXsec.push_back(0.232);          // 4: TTWJets
  sampleXsec.push_back(0.174);          // 5: TTZJets
  sampleXsec.push_back(0.0822);         // 6: WWWJets_8TeV-madgraph
  sampleXsec.push_back(0.0633);         // 7: WWZNoGstarJets_8TeV-madgraph
  sampleXsec.push_back(0.0192);         // 8: WZZNoGstarJets_8TeV-madgraph
  sampleXsec.push_back(234.0);          // 9: TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola
  sampleXsec.push_back(11.1773);        // 10: T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(59.5364608);     // 11: T_t-channel_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(1.0);            // 12: T_s-channel_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(11.1773);        // 13: Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(32.168436);      // 14: Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(1.0);            // 15: Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola
  sampleXsec.push_back(37509.0);        // 16: WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball
  sampleXsec.push_back(860.5);          // 17: DYJetsToLL_M-10To50filter_8TeV-madgraph
  sampleXsec.push_back(3532.8);         // 18: DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball

  std::vector<int> sampleProcessId; 
  // ids are taken from: https://docs.google.com/spreadsheet/ccc?key=0Ankm0DuoD0h0dHRhV1VNSlV1NEdhNFdyOXh3eFpSMHc&hl=en_US#gid=31
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 
  sampleProcessId.push_back(1); // 

  cout << "debug: start bkg" << endl;

  // backgrounds
  float nEv[NSAMPLES];
  for(int isample=0; isample<NSAMPLES; isample++) {
    nEv[isample] = 0.0;
  }

  for(int isample=0; isample<NSAMPLES; isample++) {

    cout << "processing sample " << isample << endl;
    Int_t           nCuts;
    Float_t         nSel[29];   //[nCuts]
    
    // List of branches
    TBranch        *b_nCuts;   //!
    TBranch        *b_nSel;   //!
    
    chains[isample]->SetBranchAddress("nCuts", &nCuts, &b_nCuts);
    chains[isample]->SetBranchAddress("nSel", nSel, &b_nSel);
    
    Long64_t nentries = chains[isample]->GetEntries();
    
    Long64_t nbytes = 0, nb = 0;
    // loop over files (>1 if VecBos in batch, splitted in many jobs)
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      
      nb = chains[isample]->GetEntry(jentry);   nbytes += nb;

      nEv[isample] += nSel[0];
    }
  }

  cout << "Processed all the chains." << endl; 

  if (sampleXsec.size() != sampleName.size() ) cout << "nasty error! check sizes..." << endl;

  std::ofstream weightsFile;
  weightsFile.open("weightTreestHV05bis.sh");
  weightsFile << "#! /bin/sh\n\n" << std::endl;
  weightsFile << "mkdir -p ../results/MC" << std::endl;
  weightsFile << "lumi=$1" << std::endl;
  
  // now write eee
  weightsFile << "echo \"Adding weights for eee datasets for \" $lumi \" pb-1...\"" << std::endl;
  weightsFile << "make LatinosAnalyzer" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    weightsFile << "./LatinosAnalyzer " << sampleName[isample].Data() 
		<< " "   << w << "*$lumi " << sampleProcessId[isample] << " 0 " << std::endl;
  }
  
  // now write mm
  weightsFile << "echo \"Adding weights for mmm datasets for \" $lumi \" pb-1...\"" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMMM = sampleName[isample].ReplaceAll("EEE","MMM");
    weightsFile << "./LatinosAnalyzer " << sampleNameMMM.Data()     
		<< " "   << w << "*$lumi " << sampleProcessId[isample] << " 1 " << std::endl;
  }

  // now write eem
  weightsFile << "echo \"Adding weights for eem datasets for \" $lumi \" pb-1...\"" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEEM = sampleName[isample].ReplaceAll("MMM","EEM");
    weightsFile << "./LatinosAnalyzer " << sampleNameEEM.Data()     
		<< " " << w << "*$lumi " << sampleProcessId[isample] << " 2 " << std::endl;
  }

  // now write me
  weightsFile << "echo \"Adding weights for mme datasets for \" $lumi \" pb-1...\"" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMME = sampleName[isample].ReplaceAll("EEM","MME");
    weightsFile << "./LatinosAnalyzer " << sampleNameMME.Data() 
		<< " " << w << "*$lumi " << sampleProcessId[isample] << " 3 " << std::endl;
  }
  weightsFile << "echo \"done weighting.\"" << std::endl;

  // now write ee
  weightsFile << "echo \"Adding weights for ee datasets for \" $lumi \" pb-1...\"" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameEE = sampleName[isample].ReplaceAll("MME","EE");
    weightsFile << "./LatinosAnalyzer " << sampleNameEE.Data() 
		<< " " << w << "*$lumi " << sampleProcessId[isample] << " 4 " << std::endl;
  }
  weightsFile << "echo \"done weighting.\"" << std::endl;

  // now write mm
  weightsFile << "echo \"Adding weights for mm datasets for \" $lumi \" pb-1...\"" << std::endl;

  for(int isample=0; isample<NSAMPLES; isample++) {
    //    cout << "Events processed for sample: " << sampleName[isample] << " = " << nEv[isample] << endl;
    float w = weight(nEv[isample], sampleXsec[isample], 1., 1.);
    TString sampleNameMM = sampleName[isample].ReplaceAll("EE","MM");
    weightsFile << "./LatinosAnalyzer " << sampleNameMM.Data() 
		<< " " << w << "*$lumi " << sampleProcessId[isample] << " 5 " << std::endl;
  }
  weightsFile << "echo \"done weighting.\"" << std::endl;

}

double weight(double ngen, double xsec, double filtereff, double lumi) {

  if(ngen==0) return 0;
  return xsec * filtereff * lumi / ngen;

}
