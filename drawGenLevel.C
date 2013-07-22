#include "TMath.h"
#include "TTree.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TPaveStats.h"

void drawGenLevel(){

  TFile *file = TFile::Open("genLevelInfo_tHYminus1.root");

  TTree *tree = (TTree*)file->Get("latino");
  
  TH1F* hLeptonMinusfromWfromH = new TH1F("hLeptonMinusfromWfromH","",40,0.,200.);
  TH1F* hLeptonPlusfromWfromH  = new TH1F("hLeptonPlusfromWfromH" ,"",40,0.,200.);
  TH1F* hLeptonfromWfromT      = new TH1F("hLeptonfromWfromT"     ,"",40,0.,200.);

  TH1F* hGenHiggsPt = new TH1F("hGenHiggsPt","",40,0.,400.);
  TH1F* hGenTopPt   = new TH1F("hGenTopPt"  ,"",40,0.,400.);

  TH1F* hGenHiggsEta = new TH1F("hGenHiggsEta","",40,-5.0,5.0);
  TH1F* hGenTopEta   = new TH1F("hGenTopEta"  ,"",40,-5.0,5.0);

  TH1F* hGenForwardQuarkEta = new TH1F("hGenForwardQuarkEta","",40,-5.0,5.0);
  TH1F* hGenbQuarkEta       = new TH1F("hGenbQuarkEta"      ,"",40,-5.0,5.0);

  TH1F* hGenForwardQuarkPt = new TH1F("hGenForwardQuarkPt","",40,0.,200.0);
  TH1F* hGenbQuarkPt       = new TH1F("hGenbQuarkPt"      ,"",40,0.,200.0);

  TCanvas *c1 = new TCanvas("canvas1","canvas1");
  c1->cd();
  hLeptonMinusfromWfromH->SetNormFactor(1.);
  hLeptonMinusfromWfromH->SetLineColor(1);
  hLeptonMinusfromWfromH->SetLineWidth(2);
  hLeptonMinusfromWfromH->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hLeptonMinusfromWfromH->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genLeptonMinusfromWfromH_Pt>>hLeptonMinusfromWfromH","genLeptonMinusfromWfromH_Pt>0.");

  TCanvas *c2 = new TCanvas("canvas2","canvas2");
  c2->cd();
  hLeptonPlusfromWfromH->SetNormFactor(1.);
  hLeptonPlusfromWfromH->SetLineColor(2);
  hLeptonPlusfromWfromH->SetLineWidth(2);
  hLeptonPlusfromWfromH->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hLeptonPlusfromWfromH->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genLeptonPlusfromWfromH_Pt>>hLeptonPlusfromWfromH","genLeptonPlusfromWfromH_Pt>0.");

  TCanvas *c3 = new TCanvas("canvas3","canvas3");
  c3->cd();
  hLeptonfromWfromT->SetNormFactor(1.);
  hLeptonfromWfromT->SetLineColor(4);
  hLeptonfromWfromT->SetLineWidth(3);
  hLeptonfromWfromT->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hLeptonfromWfromT->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genLeptonfromWfromT_Pt>>hLeptonfromWfromT","genLeptonfromWfromT_Pt>0.");

  TCanvas *c4 = new TCanvas("canvas4","canvas4");
  c4->cd();

  l2 = new TLegend(0.50,0.6,0.65,0.7,NULL,"brNDC");
  l2->AddEntry(hLeptonMinusfromWfromH,  " lepton(+) from W(+)","l");
  l2->AddEntry(hLeptonPlusfromWfromH ,  " lepton(-) from W(-)","l");
  l2->AddEntry(hLeptonfromWfromT     ,  " lepton(\\pm) from t(\\pm)","l");
  l2->SetTextSize(0.03);
  l2->SetFillColor(kWhite);
  l2->SetBorderSize(0);  

  hLeptonMinusfromWfromH->Draw();
  hLeptonPlusfromWfromH->Draw("same");
  hLeptonfromWfromT->Draw("same");
  l2->Draw("same");

  TCanvas *c5 = new TCanvas("canvas5","canvas5");
  c5->cd();
  hGenHiggsPt->SetNormFactor(1.);
  hGenHiggsPt->SetLineColor(1);
  hGenHiggsPt->SetLineWidth(2);
  hGenHiggsPt->SetTitle("Higgs pT");
  hGenHiggsPt->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hGenHiggsPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genHiggsPt>>hGenHiggsPt","genHiggsPt>0.");

  TCanvas *c6 = new TCanvas("canvas6","canvas6");
  c6->cd();
  hGenTopPt->SetNormFactor(1.);
  hGenTopPt->SetLineColor(2);
  hGenTopPt->SetLineWidth(2);
  hGenTopPt->SetTitle("Top pT");
  hGenTopPt->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hGenTopPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genTopPt>>hGenTopPt","genTopPt>0.");

  TCanvas *c7 = new TCanvas("canvas7","canvas7");
  c7->cd();
  hGenTopPt->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hGenTopPt->Draw();
  hGenHiggsPt->Draw("same");

  l3 = new TLegend(0.50,0.6,0.60,0.7,NULL,"brNDC");
  l3->AddEntry(hGenTopPt  , " Top","l");
  l3->AddEntry(hGenHiggsPt, " Higgs","l");
  l3->SetTextSize(0.03);
  l3->SetFillColor(kWhite);
  l3->SetBorderSize(0);  
  l3->Draw("same");
  
  TCanvas *c8 = new TCanvas("canvas8","canvas8");
  c8->cd();
  hGenHiggsEta->SetNormFactor(1.);
  hGenHiggsEta->SetLineColor(1);
  hGenHiggsEta->SetLineWidth(2);
  //hGenHiggsEta->SetTitle("Higgs \\eta");
  hGenHiggsEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenHiggsEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genHiggsEta>>hGenHiggsEta");

  TCanvas *c9 = new TCanvas("canvas9","canvas9");
  c9->cd();
  hGenTopEta->SetNormFactor(1.);
  hGenTopEta->SetLineColor(2);
  hGenTopEta->SetLineWidth(2);
  //hGenTopEta->SetTitle("Top \\eta");
  hGenTopEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenTopEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genTopEta>>hGenTopEta");

  TCanvas *c10 = new TCanvas("canvas10","canvas10");
  c10->cd();
  hGenTopEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenHiggsEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenHiggsEta->Draw("");
  hGenTopEta->Draw("same");

  l4 = new TLegend(0.5,0.3,0.6,0.4,NULL,"brNDC");
  l4->AddEntry(hGenTopEta  , " Top","l");
  l4->AddEntry(hGenHiggsEta, " Higgs","l");
  l4->SetTextSize(0.03);
  l4->SetFillColor(kWhite);
  l4->SetBorderSize(0);  
  l4->Draw("same");

  TCanvas *c11 = new TCanvas("canvas11","canvas11");
  c11->cd();
  hGenForwardQuarkEta->SetNormFactor(1.);
  hGenForwardQuarkEta->SetLineWidth(2);
  hGenForwardQuarkEta->SetLineColor(kOrange+4);
  //hGenForwardQuarkEta->SetTitle("quark q \\eta");
  hGenForwardQuarkEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenForwardQuarkEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genForwardQuark_Eta>>hGenForwardQuarkEta");
  
  TCanvas *c12 = new TCanvas("canvas12","canvas12");
  c12->cd();
  hGenbQuarkEta->SetNormFactor(1.);
  hGenbQuarkEta->SetLineWidth(2);
  //hGenbQuarkEta->SetTitle("quark b \\eta");
  hGenbQuarkEta->SetLineColor(kOrange+7);
  hGenbQuarkEta->GetXaxis()->SetTitle("\\eta^{q}");
  hGenbQuarkEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genbQuark_Eta>>hGenbQuarkEta");

  TCanvas *c13 = new TCanvas("canvas13","canvas13");
  c13->cd(); 
  hGenbQuarkEta->Draw();
  hGenForwardQuarkEta->Draw("same");

  l5 = new TLegend(0.2,0.7,0.3,0.8,NULL,"brNDC");
  l5->AddEntry(hGenForwardQuarkEta, " q quark","l");
  l5->AddEntry(hGenbQuarkEta      , " b quark","l");
  l5->SetTextSize(0.03);
  l5->SetFillColor(kWhite);
  l5->SetBorderSize(0);  
  l5->Draw("same");

  TCanvas *c14 = new TCanvas("canvas14","canvas14");
  c14->cd();
  hGenForwardQuarkPt->SetNormFactor(1.);
  hGenForwardQuarkPt->SetLineWidth(2);
  hGenForwardQuarkPt->SetLineColor(kOrange+4);
  //hGenForwardQuarkPt->SetTitle("quark q \\eta");
  hGenForwardQuarkPt->GetXaxis()->SetTitle("p_{T}^{q}");
  hGenForwardQuarkPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genForwardQuark_Pt>>hGenForwardQuarkPt");
  
  TCanvas *c15 = new TCanvas("canvas15","canvas15");
  c15->cd();
  hGenbQuarkPt->SetNormFactor(1.);
  hGenbQuarkPt->SetLineWidth(2);
  //hGenbQuarkPt->SetTitle("quark b \\eta");
  hGenbQuarkPt->SetLineColor(kOrange+7);
  hGenbQuarkPt->GetXaxis()->SetTitle("p_{T}^{q}");
  hGenbQuarkPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genbQuark_Pt>>hGenbQuarkPt");

  TCanvas *c16 = new TCanvas("canvas16","canvas16");
  c16->cd(); 
  hGenForwardQuarkPt->Draw();
  hGenbQuarkPt->Draw("same");


  l6 = new TLegend(0.2,0.7,0.3,0.8,NULL,"brNDC");
  l6->AddEntry(hGenForwardQuarkPt, " q quark","l");
  l6->AddEntry(hGenbQuarkPt      , " b quark","l");
  l6->SetTextSize(0.03);
  l6->SetFillColor(kWhite);
  l6->SetBorderSize(0);  
  l6->Draw("same");

} 


