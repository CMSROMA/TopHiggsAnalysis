#define jetIdStudy_cxx
#include "jetIdStudy.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "../src/FiguresOfMeritEvaluator.cc"  

void jetIdStudy::Loop() { 

  // variables, with association                            
  TH1F *Hs_eta      = new TH1F("Hs_eta",      "Hs_eta",       30, -5.,5.);
  TH1F *Hs_phi      = new TH1F("Hs_phi",      "Hs_phi",       30, -3.2,3.2);
  TH1F *Hs_pt       = new TH1F("Hs_pt",       "Hs_pt",        30,  0.,200.);
  TH1F *Hs_csv      = new TH1F("Hs_csv",      "Hs_csv",       30,  0.,1.);
  TH1F *Hs_betastar = new TH1F("Hs_betastar", "Hs_betastar",  30,  0.,1.);
  TH1F *Hs_rms      = new TH1F("Hs_rms", "     Hs_rms",       30,  0.,0.17);
  TH1F *Hs_mvaBased = new TH1F("Hs_mvaBased", "Hs_mvaBased",  30, -1.,1.);
  Hs_eta      -> GetXaxis() -> SetTitle("jet #eta");
  Hs_phi      -> GetXaxis() -> SetTitle("jet #phi");
  Hs_pt       -> GetXaxis() -> SetTitle("jet p_{T} [GeV]");
  Hs_csv      -> GetXaxis() -> SetTitle("jet CSV output");
  Hs_betastar -> GetXaxis() -> SetTitle("jet #beta^{*}");
  Hs_rms      -> GetXaxis() -> SetTitle("jet RMS");
  Hs_mvaBased -> GetXaxis() -> SetTitle("jet Id MVA");

  // variables, without association                            
  TH1F *Hb_eta      = new TH1F("Hb_eta",      "Hb_eta",       30, -5.,5.);
  TH1F *Hb_phi      = new TH1F("Hb_phi",      "Hb_phi",       30, -3.2,3.2);
  TH1F *Hb_pt       = new TH1F("Hb_pt",       "Hb_pt",        30,  0.,200.);
  TH1F *Hb_csv      = new TH1F("Hb_csv",      "Hb_csv",       30,  0.,1.);
  TH1F *Hb_betastar = new TH1F("Hb_betastar", "Hb_betastar",  30,  0.,1.);
  TH1F *Hb_rms      = new TH1F("Hb_rms", "     Hb_rms",       30,  0.,0.17);
  TH1F *Hb_mvaBased = new TH1F("Hb_mvaBased", "Hb_mvaBased",  30, -1.,1.);
  Hb_eta      -> GetXaxis() -> SetTitle("jet #eta");
  Hb_phi      -> GetXaxis() -> SetTitle("jet #phi");
  Hb_pt       -> GetXaxis() -> SetTitle("jet p_{T} [GeV]");
  Hb_csv      -> GetXaxis() -> SetTitle("jet CSV output");
  Hb_betastar -> GetXaxis() -> SetTitle("jet #beta^{*}");
  Hb_rms      -> GetXaxis() -> SetTitle("jet RMS");
  Hb_mvaBased -> GetXaxis() -> SetTitle("jet Id MVA");

  // mva output, split by eta region, per ROCs - follow cut based regions
  TH1F *Hs_mvaBased_25 = new TH1F("Hs_mvaBased_25", "Hs_mvaBased_25", 30, -1.,1.);
  TH1F *Hs_mvaBased_30 = new TH1F("Hs_mvaBased_30", "Hs_mvaBased_30", 30, -1.,1.);
  TH1F *Hs_mvaBased_up = new TH1F("Hs_mvaBased_up", "Hs_mvaBased_up", 30, -1.,1.);
  TH1F *Hb_mvaBased_25 = new TH1F("Hb_mvaBased_25", "Hb_mvaBased_25", 30, -1.,1.);
  TH1F *Hb_mvaBased_30 = new TH1F("Hb_mvaBased_30", "Hb_mvaBased_30", 30, -1.,1.);
  TH1F *Hb_mvaBased_up = new TH1F("Hb_mvaBased_up", "Hb_mvaBased_up", 30, -1.,1.);

  // cut based output, split by eta region, per ROCs - follow cut based regions
  TH1F *Hs_cutBased_25 = new TH1F("Hs_cutBased_25", "Hs_cutBased_25", 30, -1.,1.);
  TH1F *Hs_cutBased_30 = new TH1F("Hs_cutBased_30", "Hs_cutBased_30", 30, -1.,1.);
  TH1F *Hs_cutBased_up = new TH1F("Hs_cutBased_up", "Hs_cutBased_up", 30, -1.,1.);
  TH1F *Hb_cutBased_25 = new TH1F("Hb_cutBased_25", "Hb_cutBased_25", 30, -1.,1.);
  TH1F *Hb_cutBased_30 = new TH1F("Hb_cutBased_30", "Hb_cutBased_30", 30, -1.,1.);
  TH1F *Hb_cutBased_up = new TH1F("Hb_cutBased_up", "Hb_cutBased_up", 30, -1.,1.);

  // mva output, split by eta region, per ROCs - follow tH definition
  TH1F *Hs_mvaBased_cenBtagM = new TH1F("Hs_mvaBased_cenBtagM", "Hs_mvaBased_cenBtagM", 30, -1.,1.);
  TH1F *Hs_mvaBased_cenBtagL = new TH1F("Hs_mvaBased_cenBtagL", "Hs_mvaBased_cenBtagL", 30, -1.,1.);
  TH1F *Hs_mvaBased_fwdM     = new TH1F("Hs_mvaBased_fwdM",     "Hs_mvaBased_fwdM",     30, -1.,1.);
  TH1F *Hs_mvaBased_fwdL     = new TH1F("Hs_mvaBased_fwdL",     "Hs_mvaBased_fwdL",     30, -1.,1.);
  TH1F *Hb_mvaBased_cenBtagM = new TH1F("Hb_mvaBased_cenBtagM", "Hb_mvaBased_cenBtagM", 30, -1.,1.);
  TH1F *Hb_mvaBased_cenBtagL = new TH1F("Hb_mvaBased_cenBtagL", "Hb_mvaBased_cenBtagL", 30, -1.,1.);
  TH1F *Hb_mvaBased_fwdM     = new TH1F("Hb_mvaBased_fwdM",     "Hb_mvaBased_fwdM",     30, -1.,1.);
  TH1F *Hb_mvaBased_fwdL     = new TH1F("Hb_mvaBased_fwdL",     "Hb_mvaBased_fwdL",     30, -1.,1.);

  // cut based output, split by eta region, per ROCs - follow tH definition
  TH1F *Hs_cutBased_cenBtagM = new TH1F("Hs_cutBased_cenBtagM", "Hs_cutBased_cenBtagM", 30, -1.,1.);
  TH1F *Hs_cutBased_cenBtagL = new TH1F("Hs_cutBased_cenBtagL", "Hs_cutBased_cenBtagL", 30, -1.,1.);
  TH1F *Hs_cutBased_fwdM     = new TH1F("Hs_cutBased_fwdM",     "Hs_cutBased_fwdM",     30, -1.,1.);
  TH1F *Hs_cutBased_fwdL     = new TH1F("Hs_cutBased_fwdL",     "Hs_cutBased_fwdL",     30, -1.,1.);
  TH1F *Hb_cutBased_cenBtagM = new TH1F("Hb_cutBased_cenBtagM", "Hb_cutBased_cenBtagM", 30, -1.,1.);
  TH1F *Hb_cutBased_cenBtagL = new TH1F("Hb_cutBased_cenBtagL", "Hb_cutBased_cenBtagL", 30, -1.,1.);
  TH1F *Hb_cutBased_fwdM     = new TH1F("Hb_cutBased_fwdM",     "Hb_cutBased_fwdM",     30, -1.,1.);
  TH1F *Hb_cutBased_fwdL     = new TH1F("Hb_cutBased_fwdL",     "Hb_cutBased_fwdL",     30, -1.,1.);
  
  // efficiencies                                  

  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // making simple variables plots
    for (int ii=0; ii<20; ii++) {

      if (jetPuId_pt[ii]<0) continue;

      bool matched = jetPuId_associated[ii];
      if (matched) {

	Hs_eta -> Fill(jetPuId_eta[ii]);
	Hs_phi -> Fill(jetPuId_phi[ii]);
	Hs_pt  -> Fill(jetPuId_pt[ii]);
	Hs_csv -> Fill(jetPuId_csv[ii]);   // chiara: occhio che spesso e' -1
	Hs_rms -> Fill(jetPuId_rms[ii]);
	Hs_betastar -> Fill(jetPuId_betastar[ii]);
	Hs_mvaBased -> Fill(jetPuId_mvaBased[ii]);

	int isCutBased = jetPuId_cutBased[ii];

	// for ROCs, MVA: same eta regions as in cut-based approach
	if (fabs(jetPuId_eta[ii])<2.5)      Hs_mvaBased_25 -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<3.0) Hs_mvaBased_30 -> Fill(jetPuId_mvaBased[ii]);
	else                                Hs_mvaBased_up -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in cut-based approach
	if (fabs(jetPuId_eta[ii])<2.5)      Hs_cutBased_25 -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<3.0) Hs_cutBased_30 -> Fill(isCutBased);
	else                                Hs_cutBased_up -> Fill(isCutBased);

	// for ROCs, MVA: same eta regions as in tH analysis, medium WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.679)                             Hs_mvaBased_cenBtagM -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.679) Hs_mvaBased_fwdM     -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hs_mvaBased_fwdM     -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in tH analysis, medium WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.679)                             Hs_cutBased_cenBtagM -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.679) Hs_cutBased_fwdM     -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hs_cutBased_fwdM     -> Fill(isCutBased);

	// for ROCs, MVA: same eta regions as in tH analysis, loose WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.244)                             Hs_mvaBased_cenBtagL -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.244) Hs_mvaBased_fwdL     -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hs_mvaBased_fwdL     -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in tH analysis, loose WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.244)                             Hs_cutBased_cenBtagL -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.244) Hs_cutBased_fwdL     -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hs_cutBased_fwdL     -> Fill(isCutBased);
	
      } else {

	Hb_eta -> Fill(jetPuId_eta[ii]);
	Hb_phi -> Fill(jetPuId_phi[ii]);
	Hb_pt  -> Fill(jetPuId_pt[ii]);
	Hb_csv -> Fill(jetPuId_csv[ii]);
	Hb_rms -> Fill(jetPuId_rms[ii]);
	Hb_betastar -> Fill(jetPuId_betastar[ii]);
	Hb_mvaBased -> Fill(jetPuId_mvaBased[ii]);

	int isCutBased = jetPuId_cutBased[ii];

	// for ROCs, MVA: same eta regions as in cut-based approach
	if (fabs(jetPuId_eta[ii])<2.5)      Hb_mvaBased_25 -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<3.0) Hb_mvaBased_30 -> Fill(jetPuId_mvaBased[ii]);
	else                                Hb_mvaBased_up -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in cut-based approach
	if (fabs(jetPuId_eta[ii])<2.5)      Hb_cutBased_25 -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<3.0) Hb_cutBased_30 -> Fill(isCutBased);
	else                                Hb_cutBased_up -> Fill(isCutBased);

	// for ROCs, MVA: same eta regions as in tH analysis, medium WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.679)                             Hb_mvaBased_cenBtagM -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.679) Hb_mvaBased_fwdM     -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hb_mvaBased_fwdM     -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in tH analysis, medium WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.679)                             Hb_cutBased_cenBtagM -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.679) Hb_cutBased_fwdM     -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hb_cutBased_fwdM     -> Fill(isCutBased);

	// for ROCs, MVA: same eta regions as in tH analysis, loose WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.244)                             Hb_mvaBased_cenBtagL -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.244) Hb_mvaBased_fwdL     -> Fill(jetPuId_mvaBased[ii]);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hb_mvaBased_fwdL     -> Fill(jetPuId_mvaBased[ii]);

	// for ROCs, cut based: same eta regions as in tH analysis, loose WP
	if      (fabs(jetPuId_eta[ii])<2.4 && jetPuId_csv[ii]>0.244)                             Hb_cutBased_cenBtagL -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])<2.4 && fabs(jetPuId_eta[ii])>1. && jetPuId_csv[ii]<0.244) Hb_cutBased_fwdL     -> Fill(isCutBased);
	else if (fabs(jetPuId_eta[ii])>2.4)                                                      Hb_cutBased_fwdL     -> Fill(isCutBased);
      }

    }  // loop over jets

  } // loop over entries


  cout << endl; 
  cout << "preparing FOMs" << endl; 

  FiguresOfMeritEvaluator rocComp25;
  rocComp25.setRange(0.6,1.0,0.3,1.0);
  rocComp25.addSignal("PU jetId, MVA", Hs_mvaBased_25); 
  rocComp25.addBackgrounds(Hb_mvaBased_25);   
  rocComp25.setCutDirection(">");
  rocComp25.addSignal("PU jetId, cut based", Hs_cutBased_25); 
  rocComp25.addBackgrounds(Hb_cutBased_25);   
  rocComp25.setCutDirection(">");
  rocComp25.drawResults("rocComp25");

  FiguresOfMeritEvaluator rocComp30;
  rocComp30.addSignal("PU jetId, MVA", Hs_mvaBased_30); 
  rocComp30.addBackgrounds(Hb_mvaBased_30);   
  rocComp30.setCutDirection(">");
  rocComp30.addSignal("PU jetId, cut based", Hs_cutBased_30); 
  rocComp30.addBackgrounds(Hb_cutBased_30);   
  rocComp30.setCutDirection(">");
  rocComp30.drawResults("rocComp25To3");

  FiguresOfMeritEvaluator rocCompup;
  rocCompup.addSignal("PU jetId, MVA", Hs_mvaBased_up); 
  rocCompup.addBackgrounds(Hb_mvaBased_up);   
  rocCompup.setCutDirection(">");
  rocCompup.addSignal("PU jetId, cut based", Hs_cutBased_up); 
  rocCompup.addBackgrounds(Hb_cutBased_up);   
  rocCompup.setCutDirection(">");
  rocCompup.drawResults("rocCompGt3");

  FiguresOfMeritEvaluator rocCompcBtagM;
  rocCompcBtagM.setRange(0.7,1.0,0.,0.5);
  rocCompcBtagM.addSignal("PU jetId, MVA", Hs_mvaBased_cenBtagM); 
  rocCompcBtagM.addBackgrounds(Hb_mvaBased_cenBtagM);   
  rocCompcBtagM.setCutDirection(">");
  rocCompcBtagM.addSignal("PU jetId, cut based", Hs_cutBased_cenBtagM); 
  rocCompcBtagM.addBackgrounds(Hb_cutBased_cenBtagM);   
  rocCompcBtagM.setCutDirection(">");
  rocCompcBtagM.drawResults("rocCompcBtagM");

  FiguresOfMeritEvaluator rocCompcBtagL;
  rocCompcBtagL.setRange(0.6,1.0,0.,0.8);
  rocCompcBtagL.addSignal("PU jetId, MVA", Hs_mvaBased_cenBtagL); 
  rocCompcBtagL.addBackgrounds(Hb_mvaBased_cenBtagL);   
  rocCompcBtagL.setCutDirection(">");
  rocCompcBtagL.addSignal("PU jetId, cut based", Hs_cutBased_cenBtagL); 
  rocCompcBtagL.addBackgrounds(Hb_cutBased_cenBtagL);   
  rocCompcBtagL.setCutDirection(">");
  rocCompcBtagL.drawResults("rocCompcBtagL");

  FiguresOfMeritEvaluator rocCompFwdM;
  rocCompFwdM.setRange(0.6,1.0,0.,0.8);
  rocCompFwdM.addSignal("PU jetId, MVA", Hs_mvaBased_fwdM); 
  rocCompFwdM.addBackgrounds(Hb_mvaBased_fwdM);   
  rocCompFwdM.setCutDirection(">");
  rocCompFwdM.addSignal("PU jetId, cut based", Hs_cutBased_fwdM); 
  rocCompFwdM.addBackgrounds(Hb_cutBased_fwdM);   
  rocCompFwdM.setCutDirection(">");
  rocCompFwdM.drawResults("rocCompFwdM");

  FiguresOfMeritEvaluator rocCompFwdL;
  rocCompFwdL.setRange(0.6,1.0,0.,0.8);
  rocCompFwdL.addSignal("PU jetId, MVA", Hs_mvaBased_fwdL); 
  rocCompFwdL.addBackgrounds(Hb_mvaBased_fwdL);   
  rocCompFwdL.setCutDirection(">");
  rocCompFwdL.addSignal("PU jetId, cut based", Hs_cutBased_fwdL); 
  rocCompFwdL.addBackgrounds(Hb_cutBased_fwdL);   
  rocCompFwdL.setCutDirection(">");
  rocCompFwdL.drawResults("rocCompFwdL");


  // cosmetics   
  Hs_eta      -> SetLineColor(2);
  Hs_phi      -> SetLineColor(2);
  Hs_pt       -> SetLineColor(2);
  Hs_csv      -> SetLineColor(2);
  Hs_rms      -> SetLineColor(2);
  Hs_betastar -> SetLineColor(2);
  Hs_mvaBased -> SetLineColor(2);

  Hb_eta      -> SetLineColor(4);
  Hb_phi      -> SetLineColor(4);
  Hb_pt       -> SetLineColor(4);
  Hb_csv      -> SetLineColor(4);
  Hb_rms      -> SetLineColor(4);
  Hb_betastar -> SetLineColor(4);
  Hb_mvaBased -> SetLineColor(4);


  // Plots 
  TCanvas c1("c1","c1",1);
  Hb_eta->DrawNormalized();
  Hs_eta->DrawNormalized("same"); 
  c1.SaveAs("jet_eta.png");

  TCanvas c2("c2","c2",1);
  Hb_phi->DrawNormalized();
  Hs_phi->DrawNormalized("same"); 
  c2.SaveAs("jet_phi.png");

  TCanvas c3("c3","c3",1);
  Hb_pt->DrawNormalized();
  Hs_pt->DrawNormalized("same"); 
  c3.SaveAs("jet_pt.png");

  TCanvas c4("c4","c4",1);
  Hb_csv->DrawNormalized();
  Hs_csv->DrawNormalized("same");
  c4.SaveAs("jet_csv.png");

  TCanvas c5("c5","c5",1);
  Hs_betastar->DrawNormalized();
  Hb_betastar->DrawNormalized("same");   
  c5.SaveAs("jet_betastar.png");

  TCanvas c6("c6","c6",1);
  Hs_rms->DrawNormalized();
  Hb_rms->DrawNormalized("same"); 
  c6.SaveAs("jet_rms.png");

  TCanvas c7("c7","c7",1);
  Hs_mvaBased->DrawNormalized();
  Hb_mvaBased->DrawNormalized("same"); 
  c7.SaveAs("jet_mvaBased.png");                                                                                                                               
}
