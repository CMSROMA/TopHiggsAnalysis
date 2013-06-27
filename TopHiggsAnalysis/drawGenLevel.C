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

  TFile *file = TFile::Open("genLevelInfo.root");

  TTree *tree = (TTree*)file->Get("latino");
  
  TH1F* hLeptonMinusfromWfromH = new TH1F("hLeptonMinusfromWfromH","",40,0.,200.);
  TH1F* hLeptonPlusfromWfromH  = new TH1F("hLeptonPlusfromWfromH" ,"",40,0.,200.);
  TH1F* hLeptonfromWfromT      = new TH1F("hLeptonfromWfromT"     ,"",40,0.,200.);

  TH1F* hGenHiggsPt = new TH1F("hGenHiggsPt","",40,0.,400.);
  TH1F* hGenTopPt   = new TH1F("hGenTopPt"  ,"",40,0.,400.);

  TH1F* hGenHiggsEta = new TH1F("hGenHiggsEta","",40,-5.0,5.0);
  TH1F* hGenTopEta   = new TH1F("hGenTopEta"  ,"",40,-5.0,5.0);
  /*
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

  l2 = new TLegend(0.50,0.6,0.70,0.87,NULL,"brNDC");
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
  */

  /*
  TCanvas *c5 = new TCanvas("canvas5","canvas5");
  c5->cd();
  hGenHiggsPt->SetNormFactor(1.);
  hGenHiggsPt->SetLineColor(1);
  hGenHiggsPt->SetLineWidth(2);
  hGenHiggsPt->GetXaxis()->SetTitle("p^{H}_{T} (GeV)");
  hGenHiggsPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genHiggsPt>>hGenHiggsPt","genHiggsPt>0.");

  TCanvas *c6 = new TCanvas("canvas6","canvas6");
  c6->cd();
  hGenTopPt->SetNormFactor(1.);
  hGenTopPt->SetLineColor(2);
  hGenTopPt->SetLineWidth(2);
  hGenTopPt->GetXaxis()->SetTitle("p^{top}_{T} (GeV)");
  hGenTopPt->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genTopPt>>hGenTopPt","genTopPt>0.");

  TCanvas *c7 = new TCanvas("canvas7","canvas7");
  c7->cd();
  hGenTopPt->GetXaxis()->SetTitle("p^{gen}_{T} (GeV)");
  hGenTopPt->Draw();
  hGenHiggsPt->Draw("same");

  l3 = new TLegend(0.50,0.6,0.70,0.87,NULL,"brNDC");
  l3->AddEntry(hGenTopPt  , " Top","l");
  l3->AddEntry(hGenHiggsPt, " Higgs","l");
  l3->SetTextSize(0.03);
  l3->SetFillColor(kWhite);
  l3->SetBorderSize(0);  

  l3->Draw("same");
  */
  
  TCanvas *c8 = new TCanvas("canvas8","canvas8");
  c8->cd();
  hGenHiggsEta->SetNormFactor(1.);
  hGenHiggsEta->SetLineColor(1);
  hGenHiggsEta->SetLineWidth(2);
  hGenHiggsEta->GetXaxis()->SetTitle("\\eta^{H}");
  hGenHiggsEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genHiggsEta>>hGenHiggsEta");

  TCanvas *c9 = new TCanvas("canvas9","canvas9");
  c9->cd();
  hGenTopEta->SetNormFactor(1.);
  hGenTopEta->SetLineColor(2);
  hGenTopEta->SetLineWidth(2);
  hGenTopEta->GetXaxis()->SetTitle("\\eta^{top}");
  hGenTopEta->GetYaxis()->SetTitle("normalized to 1");
  tree->Draw("genTopEta>>hGenTopEta");

  TCanvas *c10 = new TCanvas("canvas10","canvas10");
  c10->cd();
  hGenTopEta->GetXaxis()->SetTitle("\\eta^{gen}");
  hGenTopEta->Draw();
  hGenHiggsEta->Draw("same");

  l4 = new TLegend(0.50,0.3,0.70,0.57,NULL,"brNDC");
  l4->AddEntry(hGenTopEta  , " Top","l");
  l4->AddEntry(hGenHiggsEta, " Higgs","l");
  l4->SetTextSize(0.03);
  l4->SetFillColor(kWhite);
  l4->SetBorderSize(0);  
  l4->Draw("same");

}


