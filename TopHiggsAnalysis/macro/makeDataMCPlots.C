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
#include "LatinoPlot.C"

#include <iostream>

#define NSPECIES 13
#define NVARIABLES 7
#define NCUTS 0
#define JETBINS 1

void makeDataMCPlots(int mH, const char *finalstate, float lumi, bool blindData=false, int signalFactor=1){

  gROOT->SetBatch(1);
  
  lumi = lumi * 1000.;

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine(".x LatinoStyle.C");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);  // Show overflow, underflow + SumOfWeights 
  gStyle->SetOptFit(111110); 
  gStyle->SetOptFile(1); 
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetMarkerColor(1);
  
  TString suffix="";
  
  TString species[NSPECIES];
  species[0]="Data";
  species[1]="tH";
  species[2]="WZ";
  species[3]="ZZ";
  species[4]="WW";
  species[5]="tt";
  species[6]="ttW";
  species[7]="ttZ";
  species[8]="WWW";
  species[9]="WWZ";
  species[10]="WZZ";
  species[11]="DYl";
  species[12]="DY";
  
  TString scalefactor_datadriven[NSPECIES][JETBINS];
  scalefactor_datadriven[0][0] = "1.";
  scalefactor_datadriven[1][0] = "1.";
  scalefactor_datadriven[2][0] = "1.";
  scalefactor_datadriven[3][0] = "1.";
  scalefactor_datadriven[4][0] = "1.";
  scalefactor_datadriven[5][0] = "1.";
  scalefactor_datadriven[6][0] = "1.";
  scalefactor_datadriven[7][0] = "1.";
  scalefactor_datadriven[8][0] = "1.";
  scalefactor_datadriven[9][0] = "1.";
  scalefactor_datadriven[10][0] = "1.";
  scalefactor_datadriven[11][0] = "1.";
  scalefactor_datadriven[12][0] = "1.";
  
  Color_t colors[NSPECIES];
  colors[0]=kBlack;
  colors[1]=kRed;       
  colors[2]=kAzure-9;
  colors[3]=kAzure-5;
  colors[4]=kAzure-1;
  colors[5]=kGray;
  colors[6]=kOrange+7;
  colors[7]=kOrange+1;
  colors[8]=kSpring+9;
  colors[9]=kSpring+3;
  colors[10]=kSpring-7;
  colors[11]=kViolet-1;
  colors[12]=kViolet-1;
  
  Color_t lineColors[NSPECIES];
  lineColors[0]=kBlack;
  lineColors[1]=kRed;      
  lineColors[2]=kAzure-9;
  lineColors[3]=kAzure-5;
  lineColors[4]=kAzure-1;
  lineColors[5]=kGray;
  lineColors[6]=kOrange+7;
  lineColors[7]=kOrange+1;
  lineColors[8]=kSpring+9;
  lineColors[9]=kSpring+3;
  lineColors[10]=kSpring-7;
  lineColors[11]=kViolet-1;
  lineColors[12]=kViolet-1;
  
  int legendOrder[NSPECIES];
  legendOrder[0]=0;
  legendOrder[1]=1;
  legendOrder[2]=2;
  legendOrder[3]=3;
  legendOrder[4]=4;
  legendOrder[5]=5;
  legendOrder[6]=6;
  legendOrder[7]=7;
  legendOrder[8]=8;
  legendOrder[9]=9;
  legendOrder[10]=10;
  legendOrder[11]=11;
  legendOrder[12]=12;
  
  TString files[NSPECIES];
  files[0]="/cmsrm/pc24_2/jorda/data/finalresults_V04//Data/DataAll.root";  
  files[1]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/tH125q_blvu_Yt1_H126toWW-datasetAll.root";  
  files[2]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola-datasetAll.root";
  files[3]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetAll.root";
  files[4]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola-datasetAll.root";
  files[5]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola-datasetAll.root";
  files[6]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/TTWJets-datasetAll.root";
  files[7]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/TTZJets-datasetAll.root";
  files[8]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/WWWJets_8TeV-madgraph-datasetAll.root";
  files[9]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/WWZNoGstarJets_8TeV-madgraph-datasetAll.root";
  files[10]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/WZZNoGstarJets_8TeV-madgraph-datasetAll.root";
  files[11]="/cmsrm/pc24_2/jorda/data/finalresults_V04//MC/DYJetsToLL_M-10To50filter_8TeV-madgraph-datasetAll.root";
  files[12]="/cmsrm/pc24_2/jorda/data/finalresults_V04/MC/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball-datasetAll.root";
  
  // directory for the plots
  TString plotsDir="./tHplots_V04/";

  // create the folder
  TString createFolder    = "mkdir " + plotsDir + " > /dev/null";
  system ((char*) createFolder);

  // create the subfolder for the corresponding channel
  TString createSubFolder = "mkdir " + plotsDir + "/" + TString(finalstate) + " > /dev/null";
  system ((char*) createSubFolder);

  TFile* fOut=new TFile("tH_histos_"+suffix+".root","RECREATE");
  
  char icut[NCUTS][100];
  //TH1F* histos[NSPECIES][NCUTS][NVARIABLES];    // 5 species, 2 cut levels, 8 variables
  TH1F* histos[NSPECIES][NVARIABLES];    // 5 species, 2 cut levels, 8 variables
  
  TString variables[NVARIABLES];
  variables[0]="pfmet";
  variables[1]="mll";
  variables[2]="pt1";
  variables[3]="pt2";
  variables[4]="pt3";
  variables[5]="njet";
  variables[6]="nvtx";
  
  TString units[NVARIABLES];
  units[0]="GeV";
  units[1]="GeV";
  units[2]="GeV/c^{2}";
  units[3]="GeV/c^{2}";
  units[4]="GeV/c";
  units[5]="";
  units[6]="";
  
  int nbins[NVARIABLES];
  nbins[0]=25;
  nbins[1]=30;
  nbins[2]=20;  
  nbins[3]=20;  
  nbins[4]=20;  
  nbins[5]=10;  
  nbins[6]=30;  

  float range[NVARIABLES][2]; // 8 variables, min, max
  // met
  range[0][0]=0.;
  range[0][1]=150.;
  // mll
  range[1][0]=0.;
  range[1][1]=300.;
  // pt 1
  range[2][0]=20.;
  range[2][1]=300.;   
  // pt 2
  range[3][0]=10.;
  range[3][1]=200.;   
  // pt 3
  range[4][0]=10.;
  range[4][1]=200.;   
  // njets
  range[5][0]=0.;
  range[5][1]=10.;
  // nvtx
  range[6][0]=1.;
  range[6][1]=60.;
  
  //                         0,1 2 3 4 5 6 
  int doplot[NVARIABLES] = { 1,1,1,1,1,1,1 };
  
  TString xaxisLabel[NVARIABLES];
  xaxisLabel[0]="PF E_{T}^{miss}";
  xaxisLabel[1]="m_{ll}";
  xaxisLabel[2]="p_{T}^{l,1}";
  xaxisLabel[3]="p_{T}^{l,2}";
  xaxisLabel[4]="p_{T}^{l,3}";
  xaxisLabel[5]="n jets";
  xaxisLabel[6]="n vtx";
  
  // back-up :P
  //xaxisLabel[1]="min(pr. PF E_{T}^{miss}, tk E_{T}^{miss})";
  //xaxisLabel[2]="m_{T}^{ll E_{T}^{miss}}";
  //xaxisLabel[6]="#Delta #phi_{ll}";
  //xaxisLabel[9]="M_{R}";
  //xaxisLabel[10]="DY MVA";
  
  
  TString binSize[NVARIABLES];
  
  //for (int z=0;z<2;++z){
  for (int z=0;z<NVARIABLES;++z){
    for (int i=0;i<NSPECIES;++i){
      histos[i][z]=new TH1F(variables[z]+"_W_"+species[i],variables[z]+"_W_"+species[i],nbins[z],range[z][0],range[z][1]);
      if(i==0)
	histos[i][z]->Sumw2();
      char binsiz[10];
      sprintf(binsiz,"%2.0f",(range[z][1]-range[z][0])/nbins[z]);
      binSize[z]=TString(binsiz);
    }
  }

  TString cut; // requirements
  //cut="( mll>12. && pt1>20.0 && mlll > 100. && pfmet > 30. && njet > 0 && njet < 3 && nextra == 0 && zveto == 1)*";
  cut="( mll>12. && pt1>20.0 && pt2>20.0 && nextra == 0 )*";

  TString channelcut("1*");
  if(TString(finalstate).Contains("all")) 
    channelcut=TString("(channel==0 || channel==1 || channel==2 || channel==3)*");

  if(TString(finalstate).Contains("eee")) channelcut=TString("(channel==0)*");
  if(TString(finalstate).Contains("mmm")) channelcut=TString("(channel==1)*");
  if(TString(finalstate).Contains("eem")) channelcut=TString("(channel==2)*");
  if(TString(finalstate).Contains("mme")) channelcut=TString("(channel==3)*");

  if(TString(finalstate).Contains("Zee"))  channelcut=TString("(channel==4)*");
  if(TString(finalstate).Contains("Zmm"))  channelcut=TString("(channel==5)*");
  if(TString(finalstate).Contains("Zall")) channelcut=TString("(channel==4 || channel==5)*");
  
  char lumistr[5];
  sprintf(lumistr,"%.1f",lumi);
  TString intLumi=TString(lumistr);     
  TFile *_file[NSPECIES];
  TTree *T1[NSPECIES];
  
  char lumiwgt[10];
  sprintf(lumiwgt,"%f*",lumi);

  if(!blindData) {
    _file[0]=TFile::Open(files[0]);
    T1[0] = (TTree*)_file[0]->Get("latino");
  } else T1[0] = 0;
  
  for (int i=1;i<NSPECIES;++i) {
    _file[i]=TFile::Open(files[i]);
    T1[i] = (TTree*)_file[i]->Get("latino");
  }
  
  int nspeciesToRun=NSPECIES;
  
  for (int z=0;z<NVARIABLES;++z){
    //for (int z=0;z<2;++z){

    if(doplot[z]==0) continue;
    int firstSpecie = 0;
    if(blindData) firstSpecie = 1;
    for (int i=firstSpecie;i<nspeciesToRun;++i){

      fOut->cd();
      TString histoName=variables[z]+"_W_"+species[i];
      std::cout << "Producing " << histoName << std::endl;

      if (T1[i]==0){
	std::cout << "Species " << i << " Tree not found" << std::endl;
	return;
      }
      
      if(i>0) {
	
	T1[i]->Project(histoName,variables[z],channelcut+cut+TString(lumiwgt)+TString("baseW*puW*effW*")+scalefactor_datadriven[i][0]);
	std::cout << "Done " << histoName << std::endl;

      }else{

	T1[i]->Project(histoName,variables[z],channelcut+cut+TString("1"));
	std::cout << "Done " << histoName << std::endl;

      }
      
    }

    LatinoPlot* myPlot = new LatinoPlot();
    lumi = lumi/1000;
    myPlot->setLumi(lumi);
    lumi = lumi*1000;
    myPlot->addLabel("");
    myPlot->setLabel((xaxisLabel[z]).Data());
    myPlot->setUnits((units[z]).Data());
    myPlot->setMass(mH);

    myPlot->setMCHist(itH,      histos[1][z]);
    myPlot->setMCHist(iWZ,      histos[2][z]);
    myPlot->setMCHist(iZZ,      histos[3][z]);
    myPlot->setMCHist(iWW,      histos[4][z]);
    myPlot->setMCHist(itt,      histos[5][z]);
    myPlot->setMCHist(ittw,     histos[6][z]);
    myPlot->setMCHist(ittz,     histos[7][z]);
    myPlot->setMCHist(iwww,     histos[8][z]);
    myPlot->setMCHist(iwwz,     histos[9][z]);
    myPlot->setMCHist(iwzz,     histos[10][z]);
    myPlot->setMCHist(idyl,     histos[11][z]);
    myPlot->setMCHist(idy,      histos[12][z]);

    if(!blindData) 
      myPlot->setDataHist(histos[0][z]);

    // Draw
    //--------------------------------------------------------------------
    TCanvas* c1 = new TCanvas(Form("test_%d_%d_lin", z, 0), Form("test_%d_%d_lin", z, 0));

    c1->SetLogy(0);

    //myPlot.setNoStack();
    myPlot->Draw(1);
    c1->GetFrame()->DrawClone();
    c1->SaveAs(plotsDir+"/"+TString(finalstate)+"/"+variables[z]+suffix+".lin.png");
    c1->SaveAs(plotsDir+"/"+TString(finalstate)+"/"+variables[z]+suffix+".root");
    c1->SaveAs(plotsDir+"/"+TString(finalstate)+"/"+variables[z]+suffix+".lin.eps");
    
    TCanvas* c2 = new TCanvas(Form("test_%d_%d_log", z, 0),Form("test_%d_%d_log", z, 0));
    
    c2->SetLogy(1);
    
    myPlot->Draw(1);
    
    c2->GetFrame()->DrawClone();
    
    c2->SaveAs(plotsDir+"/"+TString(finalstate)+"/"+variables[z]+suffix+".log.png");
    
    delete myPlot;

  }
  
  
  fOut->Write();
  fOut->Close();
  
}
