#include <TH2.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <fstream>

#include "TDirectory.h"
#include "TROOT.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
//#include "TH2D.h"
#include "TH3D.h"
//#include "TH2F.h"
#include "TH1I.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TChain.h"
#include "TMath.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include <vector>
#include "GlobalVariables.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "vector"
#include "map"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "tdrstyle_mod14.C"

using namespace std;
using namespace ROOT::Math;
int getBin(float x, float boundaries[],int b);
TH1F* convertbining(TH1F* h1,TH1F* h2);
//____________________________________________________________________________________
//-------------------------------------------------------------
//https://root.cern.ch/root/html/tutorials/tree/hvector.C.html  
//
//-------------------------------------------------------------
//_________________________________________________________________________________________
//


void  EnergyFlow_PlotsProducer3p8T()

{
  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  setTDRStyle();
    
  gDirectory->DeleteAll();
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");
  //gStyle->SetOptStat(111);
  bool save =true;
  bool ratioMC=true;
  string filenames[8] ={"../EFlowTreeRootFile/2EFlow_DetLevel_DataZeroBias1AllBx_tree_3p8Tesla_r251721.root",
			"../EFlowTreeRootFile/4EFlow_DetLevel_DataExpressPhysicsAllBx_tree_3p8Tesla.root",
			"../EFlowTreeRootFile/7EFlow_DetLevel_pythia8_Monash_tree_3p8Tesla.root",
            "../EFlowTreeRootFile/5EFlow_DetLevel_cuetp8m1_tree_3p8Tesla.root",
            "../EFlowTreeRootFile/6EFlow_DetLevel_epos_tree_3p8Tesla.root",
      "../EFlowTreeRootFile/9EFlow_GenLevel_pythia8_Monash_tree_3p8Tesla.root",
      "../EFlowTreeRootFile/10EFlow_GenLevel_TuneCUETP8M1_tree_3p8Tesla.root",
      "../EFlowTreeRootFile/11EFlow_GenLevel_epos_tree_3p8Tesla.root"
			 
  };
    
    
  string energyname[4]={"RecHit","CaloTower","PFClusters","Gen"};
  string fname[8]={"ZeroBias3p8","ExpPhy3p8","P8Monash3p8","Cuetp8m1_3p8","Epos3p8","GenP8Monash3p8","GenCuetp8m1_3p8","GenEpos3p8"};
  TFile *fOutFile = new TFile("EFlow_Histo.root","RECREATE");
    
    
    
  static const Int_t ftyp =8 ;// Data and MC
  static const Int_t etyp = 4 ;
    
  char title[999];
    
  TH1F * Energy_[ftyp+1][etyp+1];
  TH1F * Energy[ftyp+1][etyp+1];  
  TH1F * DataDetLOverMCDetL[ftyp+1];
  TH1F * MCDetLOverMCGenL[ftyp+1];
  TH1F * DataGenLOverMCGenL[ftyp+1];
  TH1F * PFClustersEnergy[ftyp+1][nEtaBins+1];
  TH1F * RatioAllvsFat[ftyp+1];
  for (int f=0; f<ftyp; f++){

    sprintf(title,"%sDetLOver%sDetL",fname[0].c_str(),fname[f].c_str());
    DataDetLOverMCDetL[f]=new TH1F(title,title,nEtaBins,EtaBins);
    DataDetLOverMCDetL[f]->Sumw2();

    sprintf(title,"%sDetLOver%sGenL",fname[f].c_str(),fname[f].c_str());
    MCDetLOverMCGenL[f]=new TH1F(title,title,nEtaBins,EtaBins);
    MCDetLOverMCGenL[f]->Sumw2();

    sprintf(title,"%sGenLOver%sGenL",fname[0].c_str(),fname[f].c_str());
    DataGenLOverMCGenL[f]=new TH1F(title,title,nEtaBins,EtaBins);
    DataGenLOverMCGenL[f]->Sumw2();

      sprintf(title,"%sAllBx%sFatBx",fname[0].c_str(),fname[f].c_str());
      RatioAllvsFat[f]=new TH1F(title,title,nEtaBins,EtaBins);
      RatioAllvsFat[f]->Sumw2();

    for (int e=0; e<etyp; e++){
	    
      sprintf(title,"EnergyFlow__%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy_[f][e]=new TH1F(title,title,nHBins,HBins);
      Energy_[f][e]->Sumw2();
	    
      sprintf(title,"EnergyFlow_%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy[f][e]=new TH1F(title,title,nEtaBins,EtaBins);
      Energy[f][e]->Sumw2();
	    
    }
	/*
    for (int eta=0; eta<nEtaBins; eta++){
	    
      sprintf(title,"EnergyDist_%s_EtaBin%d",fname[f].c_str(),eta);
      if ((eta>=0 && eta<=3) || (eta>=25 && eta<=28)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,150,0.,1500);
      }
      if ((eta>=4 && eta<=8) || (eta>=20 && eta<=24)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,100,0.,1000);
      }
      if (eta>=7 && eta<=19){
	PFClustersEnergy[f][eta]=new TH1F(title,title,50,0.,500);
      }
      PFClustersEnergy[f][eta]->Sumw2();
    }*/
  }
  int decade = 0;
  TFile *file[ftyp+1];
  TTree *fChain[ftyp+1];
    
  int TotNofEvent[ftyp+1];
  for (int f=0; f<ftyp; f++){
    	
    file[f]   = new TFile(filenames[f].c_str(),"READ");

    cout<<"file : "<<filenames[f].c_str()<<endl;
    
    fChain[f] = (TTree*)file[f]->Get("EFlow");
    if( f>4) fChain[f]->SetBranchAddress("Gen", &Gen, &b_Gen);
    else {
      fChain[f]->SetBranchAddress("RecHit", &RecHit, &b_RecHit);
      fChain[f]->SetBranchAddress("CaloTower", &CaloTower, &b_CaloTower);
      fChain[f]->SetBranchAddress("PFClusters", &PFClusters, &b_PFClusters);
    }
   
    // ----------------------- Event -----------------------//
    Long64_t nentries = fChain[f]->GetEntriesFast();
    cout<<"Entries "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;

    TotNofEvent[f]= nentries;

    int maxevent=10000;
    for (Long64_t ev=0; ev<nentries;ev++) {
      
        //for (Long64_t ev=0; ev<maxevent;ev++) {
	
      Long64_t iev = fChain[f]->LoadTree(ev);
      if (iev < 0)break;
      nb = fChain[f]->GetEntry(ev);   nbytes += nb;
      // cout<<"EVENT "<<ev<<endl;
      double progress = 10.0*ev/(1.0*nentries);
      int k = TMath::FloorNint(progress); 
      if (k > decade) 
	cout<<10*k<<" %"<<endl;
      decade = k; 
      for (int i =0; i<nEtaBins;i++) {
	  
	RecHitEtaSums[i]=0.;
	CaloEtaSums[i]=0.;
	PFClustersEtaSums[i]=0.;
	GenEtaSums[i]=0.;
      }
      if (f>4){
	// -----------------------  Gen  -----------------//
	for(unsigned long gen=0; gen<Gen->size(); gen++) {
	  GenEtaSums[gen]= GenEtaSums[gen] + Gen->at(gen);
	}
	for (int k=0;k<nEtaBins;k++){
	Energy_[f][3]->Fill(k,GenEtaSums[k]);
	} 
      }
      
      else{
	  
	// ----------------------- RecHit --------------- //
	for(unsigned long rec=0; rec<RecHit->size(); rec++) {
	  RecHitEtaSums[rec]= RecHitECALEtaSums[rec] + RecHit->at(rec);
	}  
	  
	//----------------------- CaloTower  -------------//
	for(unsigned long cal=0; cal<CaloTower->size(); cal++) {
	  CaloEtaSums[cal]= CaloEtaSums[cal] + CaloTower->at(cal);
	}  
	//----------------------- PF Cluster --------------- //
	for(unsigned long pfc=0; pfc<PFClusters->size(); pfc++) {
	  PFClustersEtaSums[pfc]= PFClustersEtaSums[pfc] + PFClusters->at(pfc);
	}  
	//{"RecHit","CaloTower","PFClusters"};
	for (int k=0;k<nEtaBins;k++){
	  Energy_[f][0]->Fill(k,RecHitEtaSums[k]);
	  Energy_[f][1]->Fill(k,CaloEtaSums[k]);
	  Energy_[f][2]->Fill(k,PFClustersEtaSums[k]);
	  //PFClustersEnergy[f][k]->Fill(PFClustersEtaSums[k]);
	} 
	  
      }//else
    }//event 
  }//File
    
  // Plots style
    
  //"ZeroBias3p8","ExpPhy3p8","P8Monash3p8","Cuetp8m1_3p8","Epos3p8","GenP8Monash3p8","GenCuetp8m1_3p8","GenEpos3p8"
    
    //%%%%%%%%%EnergyDistributions///%%%%%%%%%%%%

 /* for (int f=0; f<4; f++){
    
    for (int eta=0; eta<nEtaBins; eta++){
      PFClustersEnergy[f][eta]->GetXaxis()->SetTitle("Energy");
      PFClustersEnergy[f][eta]->GetYaxis()->SetTitle("1/N");
      // PFClustersEnergy[f][eta]->SetMinimum(1e-5);
      // PFClustersEnergy[f][eta]->SetMaximum(10);
      if (f == 0){
	PFClustersEnergy[f][eta]->SetLineColor(1);
	PFClustersEnergy[f][eta]->SetMarkerColor(1);
	PFClustersEnergy[f][eta]->SetMarkerStyle(20);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==1){
	PFClustersEnergy[f][eta]->SetLineColor(kBlue);
	PFClustersEnergy[f][eta]->SetMarkerColor(kBlue);
	PFClustersEnergy[f][eta]->SetMarkerStyle(24);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==2 ){
	PFClustersEnergy[f][eta]->SetLineColor(kRed);
	PFClustersEnergy[f][eta]->SetMarkerColor(kRed);
	PFClustersEnergy[f][eta]->SetMarkerStyle(25);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==3) {
	PFClustersEnergy[f][eta]->SetLineColor(kYellow);
	PFClustersEnergy[f][eta]->SetMarkerColor(kYellow);
	PFClustersEnergy[f][eta]->SetMarkerStyle(26);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==4){
	PFClustersEnergy[f][eta]->SetLineColor(kMagenta);
	PFClustersEnergy[f][eta]->SetMarkerColor(kMagenta);
	PFClustersEnergy[f][eta]->SetMarkerStyle(27);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==5){
	PFClustersEnergy[f][eta]->SetLineColor(kGreen);
	PFClustersEnergy[f][eta]->SetMarkerColor(kGreen);
	PFClustersEnergy[f][eta]->SetMarkerStyle(28);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }


      PFClustersEnergy[f][eta]->Scale(1./TotNofEvent[f]);
      PFClustersEnergy[f][eta]->SetMinimum(1e-5);
      PFClustersEnergy[f][eta]->SetMaximum(10);
    }
	
  } */
    
    
  //string energyname[4]={"RecHit","CaloTower","PFClusters","Gen"};
  //string fname[11]={"ZeroBias3p8","ExpPhy3p8","P8Monash3p8","Cuetp8m1_3p8","Epos3p8","GenP8Monash3p8","GenCuetp8m1_3p8","GenEpos3p8"};
  for (int f=0; f<ftyp; f++){
      
    for (int e=0; e<etyp; e++){
      convertbining(Energy_[f][e],Energy[f][e]);
      Energy[f][e]->Scale(1./TotNofEvent[f]);
	  
      Energy[f][e]->GetXaxis()->SetTitle("#eta");
      Energy[f][e]->GetYaxis()->SetTitle("(1 / N) dE/d#eta (GeV)");
      Energy[f][e]->SetMinimum(0.4);
      Energy[f][e]->SetMaximum(2000.);
      //DataZeroBias3p8 RecHit
      if (f == 0 && (e ==0)){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(21);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //DataZeroBias3p8 CaloTower
      if (f == 0 && e == 1){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(22);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //DataZeroBias3p8 PFClusters
      if (f == 0 && e ==2){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(20);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //ExpressPhy RecHit
      if (f == 1 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //ExpressPhy CaloTower
      if (f == 1 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //ExpressPhy PFClusters
      if (f == 1 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(27);
	Energy[f][e]->SetMarkerSize(2.);
      }
      //P8Monash3p8 RecHit
      if (f == 2 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //P8Monash3p8 CaloTower
      if (f == 2 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //P8Monash3p8 PFClusters
      if (f == 2 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(2);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(28);
	Energy[f][e]->SetMarkerSize(1.5);
      }
        //Cuetp8m1_3p8 RecHit
        if (f == 3 && e == 0){
            Energy[f][e]->SetLineStyle(1);
            Energy[f][e]->SetLineColor(kRed);
            Energy[f][e]->SetMarkerColor(kRed);
            Energy[f][e]->SetMarkerStyle(25);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        //Cuetp8m1_3p8 CaloTower
        if (f == 3 && e==1 ){
            Energy[f][e]->SetLineStyle(2);
            Energy[f][e]->SetLineColor(kRed);
            Energy[f][e]->SetMarkerColor(kRed);
            Energy[f][e]->SetMarkerStyle(26);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        //Cuetp8m1_3p8 PFClusters
        if (f == 3 && e==2 ){
            Energy[f][e]->SetLineStyle(3);
            Energy[f][e]->SetLineWidth(1);
            Energy[f][e]->SetLineColor(kRed);
            Energy[f][e]->SetMarkerColor(kRed);
            Energy[f][e]->SetMarkerStyle(30);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        //Epos3p8 RecHit
        if (f == 4 && e == 0){
            Energy[f][e]->SetLineStyle(1);
            Energy[f][e]->SetLineColor(kBlue);
            Energy[f][e]->SetMarkerColor(kBlue);
            Energy[f][e]->SetMarkerStyle(25);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        //Epos3p8 CaloTower
        if (f == 4 && e==1 ){
            Energy[f][e]->SetLineStyle(2);
            Energy[f][e]->SetLineColor(kBlue);
            Energy[f][e]->SetMarkerColor(kBlue);
            Energy[f][e]->SetMarkerStyle(26);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        //Epos3p8 PFClusters
        if (f == 4 && e==2 ){
            Energy[f][e]->SetLineStyle(3);
            Energy[f][e]->SetLineColor(kBlue);
            Energy[f][e]->SetMarkerColor(kBlue);
            Energy[f][e]->SetMarkerStyle(32);
            Energy[f][e]->SetMarkerSize(1.5);
        }
        
        
        // Gen Plots
        if (f == 5 && e==3 ){
            Energy[f][e]->SetLineStyle(3);
            Energy[f][e]->SetLineWidth(5);
            Energy[f][e]->SetLineColor(kTeal+3);
            Energy[f][e]->SetMarkerColor(kTeal+3);
            //Energy[f][e]->SetMarkerStyle(27);
            //Energy[f][e]->SetMarkerSize(1.5);
        }
        if (f == 6 && e==3 ){
            
            Energy[f][e]->SetLineStyle(3);
            Energy[f][e]->SetLineWidth(5);
            Energy[f][e]->SetLineColor(kMagenta+2);
            Energy[f][e]->SetMarkerColor(kMagenta+2);
            //Energy[f][e]->SetMarkerStyle(28);
            //Energy[f][e]->SetMarkerSize(1.5);
        }
        if (f == 7 && e==3 ){
            Energy[f][e]->SetLineStyle(3);
            Energy[f][e]->SetLineWidth(5);
            Energy[f][e]->SetLineColor(kRed);
            Energy[f][e]->SetMarkerColor(kRed);
            //Energy[f][e]->SetMarkerStyle(30);
            //Energy[f][e]->SetMarkerSize(1.5);
        }
        
        
    
    /*
   

     
      if (f == 9 && e==3 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	//Energy[f][e]->SetMarkerStyle(32);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
	   
      if (f == 10 && e==3 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	//Energy[f][e]->SetMarkerStyle(3);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
    
    //Data Fat Bunch PFClusters
    if (f == 11 && e ==2){
       Energy[f][e]->SetLineColor(kMagenta);
       Energy[f][e]->SetLineStyle(2);
       Energy[f][e]->SetMarkerColor(kMagenta);
       Energy[f][e]->SetMarkerStyle(24);
       Energy[f][e]->SetMarkerSize(1.5);
    }
        
    //Data BPTX PlusOnly
    if (f == 12 && e ==2){
       Energy[f][e]->SetLineColor(kRed);
       Energy[f][e]->SetLineStyle(2);
       Energy[f][e]->SetMarkerColor(kRed);
       Energy[f][e]->SetMarkerStyle(20);
       Energy[f][e]->SetMarkerSize(1.5);
    }
    //Data BPTX MinusOnly
    if (f == 13 && e ==2){
       Energy[f][e]->SetLineColor(kBlue);
       Energy[f][e]->SetLineStyle(2);
       Energy[f][e]->SetMarkerColor(kBlue);
       Energy[f][e]->SetMarkerStyle(20);
       Energy[f][e]->SetMarkerSize(1.5);
    }
    //Data NoBPTX
    if (f == 14 && e ==2){
       Energy[f][e]->SetLineColor(kGreen);
       Energy[f][e]->SetLineStyle(2);
       Energy[f][e]->SetMarkerColor(kGreen);
       Energy[f][e]->SetMarkerStyle(20);
       Energy[f][e]->SetMarkerSize(1.5);
    }
        
        //Data ZBPlusL1
        if (f == 15 && e ==2){
            Energy[f][e]->SetLineColor(28);
            Energy[f][e]->SetLineStyle(2);
            Energy[f][e]->SetMarkerColor(28);
            Energy[f][e]->SetMarkerStyle(25);
            Energy[f][e]->SetMarkerSize(1.5);
        }
	   */
    }//Energy type
      //"ZeroBias3p8","ExpPhy3p8","P8Monash3p8","Cuetp8m1_3p8","Epos3p8","GenP8Monash3p8","GenCuetp8m1_3p8","GenEpos3p8"
      
    if(f>1 && f<5)DataDetLOverMCDetL[f-2]->Divide(Energy[0][2],Energy[f][2], 1.,1.,"B");
    if(f>4 && f<8)MCDetLOverMCGenL[f-5]->Divide(Energy[f-3][2],Energy[f][3], 1.,1.,"B");
    RatioAllvsFat[0]->Divide(Energy[0][2],Energy[1][2], 1.,1.,"B");
  }//File Type

  DataDetLOverMCDetL[0]->GetXaxis()->SetTitle("#eta");
  DataDetLOverMCDetL[0]->GetYaxis()->SetTitle("#frac{1}{N}#frac{dE}{d#eta}[DATA] / #frac{1}{N}#frac{dE}{d#eta}[MC]");
  DataDetLOverMCDetL[0]->SetMinimum(0.45);
  DataDetLOverMCDetL[0]->SetMaximum(2.);
  DataDetLOverMCDetL[0]->SetLineColor(kTeal+3);
  DataDetLOverMCDetL[0]->SetLineWidth(4);
  DataDetLOverMCDetL[0]->SetLineStyle(1);
  DataDetLOverMCDetL[0]->SetMarkerColor(kTeal+3);
    
  DataDetLOverMCDetL[1]->SetLineColor(kMagenta+2);
  DataDetLOverMCDetL[1]->SetLineWidth(4);
  DataDetLOverMCDetL[1]->SetLineStyle(2);
  DataDetLOverMCDetL[1]->SetMarkerColor(kMagenta+2);
    
  DataDetLOverMCDetL[2]->SetLineColor(kRed);
  DataDetLOverMCDetL[2]->SetLineWidth(4);
  DataDetLOverMCDetL[2]->SetLineStyle(3);
  DataDetLOverMCDetL[2]->SetMarkerColor(kRed);
    
 // DataDetLOverMCDetL[3]->SetLineColor(kBlue);
  //DataDetLOverMCDetL[3]->SetLineWidth(4);
  //DataDetLOverMCDetL[3]->SetLineStyle(4);
  //DataDetLOverMCDetL[3]->SetMarkerColor(kBlue);
    
  //DataDetLOverMCDetL[4]->SetLineColor(kOrange);
  //DataDetLOverMCDetL[4]->SetLineWidth(4);
  //DataDetLOverMCDetL[4]->SetLineStyle(5);
  //DataDetLOverMCDetL[4]->SetMarkerColor(kOrange);

  MCDetLOverMCGenL[0]->GetXaxis()->SetTitle("#eta");
  MCDetLOverMCGenL[0]->GetYaxis()->SetTitle("#frac{1}{N}#frac{dE}{d#eta}[PFClusters] / #frac{1}{N}#frac{dE}{d#eta}[GenLevel]");
  MCDetLOverMCGenL[0]->SetMinimum(0.2);
  MCDetLOverMCGenL[0]->SetMaximum(1.6);
  MCDetLOverMCGenL[0]->SetLineColor(kTeal+3);
  MCDetLOverMCGenL[0]->SetLineWidth(4);
  MCDetLOverMCGenL[0]->SetLineStyle(1);
  MCDetLOverMCGenL[0]->SetMarkerColor(kTeal+3);
    
  MCDetLOverMCGenL[1]->SetLineColor(kMagenta+2);
  MCDetLOverMCGenL[1]->SetLineWidth(4);
  MCDetLOverMCGenL[1]->SetLineStyle(2);
  MCDetLOverMCGenL[1]->SetMarkerColor(kMagenta+2);
    
  MCDetLOverMCGenL[2]->SetLineColor(kRed);
  MCDetLOverMCGenL[2]->SetLineWidth(4);
  MCDetLOverMCGenL[2]->SetLineStyle(3);
  MCDetLOverMCGenL[2]->SetMarkerColor(kRed);
    
  //MCDetLOverMCGenL[3]->SetLineColor(kBlue);
  //MCDetLOverMCGenL[3]->SetLineWidth(4);
  //MCDetLOverMCGenL[3]->SetLineStyle(4);
  //MCDetLOverMCGenL[3]->SetMarkerColor(kBlue);
    
  //MCDetLOverMCGenL[4]->SetLineColor(kOrange);
  //MCDetLOverMCGenL[4]->SetLineWidth(4);
  //MCDetLOverMCGenL[4]->SetLineStyle(5);
  //MCDetLOverMCGenL[4]->SetMarkerColor(kOrange);
    
    RatioAllvsFat[0]->GetXaxis()->SetTitle("#eta");
    RatioAllvsFat[0]->GetYaxis()->SetTitle("#frac{1}{N}#frac{dE}{d#eta}[AllBx] / #frac{1}{N}#frac{dE}{d#eta}[FatBx]");
    RatioAllvsFat[0]->SetMinimum(0.2);
    RatioAllvsFat[0]->SetMaximum(1.6);
    RatioAllvsFat[0]->SetLineColor(kTeal+3);
    RatioAllvsFat[0]->SetLineWidth(4);
    RatioAllvsFat[0]->SetLineStyle(1);
    RatioAllvsFat[0]->SetMarkerColor(kTeal+3);


  //////************************************************************************************////
  /////                           NOW PLOTS                                                 ///
  ////*************************************************************************************////
  TLine *line=new TLine(-5.191,1.,5.191,1.);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
      
  TLegend *leg[999];
  leg[0]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[0]->AddEntry(Energy[0][2],"ZB3p8","PL");
  leg[0]->AddEntry(Energy[1][2],"ExpPhy3p8","PL");
    
    
    leg[3]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
    leg[3]->AddEntry(DataDetLOverMCDetL[0],"pythia8Monash","L");
    leg[3]->AddEntry(DataDetLOverMCDetL[1],"CuetP8M1","L");
    leg[3]->AddEntry(DataDetLOverMCDetL[3],"epos","L");
    
    leg[4]= tdrLeg(0.35,0.62,0.57,0.90);
    leg[4]->AddEntry(Energy[1][2],"pythia8Monash","p");
    leg[4]->AddEntry(Energy[2][2],"pythia8MBR(#times 10)","p");
    leg[4]->AddEntry(Energy[3][2],"herwigpp(#times 10^{2})","p");
    leg[4]->AddEntry(Energy[4][2],"epos(#times 10^{3})","p");
    leg[4]->AddEntry(Energy[5][2],"qgsjetII(#times 10^{4})","p");

  
    
/*
  leg[1]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[1]->AddEntry(Energy[1][0],"RecHits","PL");
  leg[1]->AddEntry(Energy[1][1],"CaloTowers","PL");
  leg[1]->AddEntry(Energy[1][2],"PFClusters","PL");

  leg[2]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[2]->AddEntry(Energy[0][2],"ZB3p8","PL");
  leg[2]->AddEntry(Energy[1][2],"ZBlhcf","PL");
  leg[2]->AddEntry(Energy[2][2],"pythia8MBR","PL");
  leg[2]->AddEntry(Energy[3][2],"herwigpp","PL");
  leg[2]->AddEntry(Energy[4][2],"epos","PL");
  leg[2]->AddEntry(Energy[5][2],"qgsjetII","PL");

 

  leg[4]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[4]->AddEntry(Energy[1][2],"pythia8Monash","p");
  leg[4]->AddEntry(Energy[2][2],"pythia8MBR(#times 10)","p");
  leg[4]->AddEntry(Energy[3][2],"herwigpp(#times 10^{2})","p");
  leg[4]->AddEntry(Energy[4][2],"epos(#times 10^{3})","p");
  leg[4]->AddEntry(Energy[5][2],"qgsjetII(#times 10^{4})","p");

  leg[5]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[5]->AddEntry(Energy[0][2],"Data All Bx","p");
  leg[5]->AddEntry(Energy[11][2],"Data Fat Bx","p");
    
  leg[6]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[6]->AddEntry(Energy[0][2],"Data","p");
  leg[6]->AddEntry(Energy[12][2],"Bptx PlusOnly","p");
  leg[6]->AddEntry(Energy[13][2],"Bptx MinusOnly","p");
  leg[6]->AddEntry(Energy[14][2],"No Bptx","p");
    
    leg[7]= tdrLeg(0.35,0.62,0.57,0.90);
    leg[7]->AddEntry(RatioAllvsFat[0],"Data AllBx/FatBx","p");
    
    leg[8]= tdrLeg(0.35,0.62,0.57,0.90);
    leg[8]->AddEntry(Energy[0][2],"DataExpr ZB","p");
    leg[8]->AddEntry(Energy[15][2],"DataExpr ZBpL1","p");
    
 */
  TCanvas *c[999];

  
  //------------------ DATA --------------------//
  //Energy[0][0]->GetXaxis()->SetNdivisions(511);
 
  c[0] = tdrCanvas("EnergyFlow_0T_Data",Energy[0][0],4,11,kRectangular);
  c[0]->cd();
  c[0]->SetLogy();
  for (int f=0; f<2; f++){
    Energy[f][2]->Draw("same");
  }
  tex->DrawLatex(0.75,0.85,"DATA3p8");
  leg[0]->Draw();
  sprintf(title,"Plots/EnergyFlow_3p8_Data.pdf");
  if(save)c[0]->SaveAs(title);
    
    // //------------------ Ratio of Det Level DATA and MC --------------------//
    
    DataDetLOverMCDetL[0]->GetXaxis()->SetNdivisions(511);
    c[3] = tdrCanvas("DataDetLOverMCDetL",DataDetLOverMCDetL[0],4,11,kRectangular);
    c[3]->cd();
    for (int f=0; f<3; f++)DataDetLOverMCDetL[f]->Draw("same hist");
    tex->DrawLatex(0.75,0.20,"PFClusters");
    leg[3]->Draw(); line->Draw();
    sprintf(title,"Plots/DataDetLOverMCDetL.pdf");
    if(save)c[3]->SaveAs(title);
    
    //------------------ MC Detlevel and Genlevel Compare EFlow was Fixed-------------------//
    if(ratioMC){
        Energy[5][3]->SetMinimum(1.);
        Energy[5][3]->SetMaximum(1e8);
        Energy[5][3]->GetXaxis()->SetNdivisions(511);
        c[4] = tdrCanvas("EnergyFlow_0T_MC_Det_GenLevel_Comp",Energy[5][3],4,11,kRectangular);
        c[4]->cd();c[4]->SetLogy();
        
        for (int f=2; f<8; f++){
            
            if(f<5){
                Energy[f][2]->Scale( (pow(10,(f-2))) ) ;
                Energy[f][2]->Draw("same E0");
                //Energy[f][2]->Scale( 1./(pow(10,(f-1))) ) ;
                
            }
            if(f>4){
                Energy[f][3]->Scale( (pow(10,(f-5))) ) ;
                Energy[f][3]->Draw("same hist");
                //Energy[f][3]->Scale( 1./ (pow(10,(f-6))) ) ;
            }
        }
        leg[4]->Draw(); tex->DrawLatex(0.66,0.87,"Lines : Gen Level");
        leg[4]->Draw(); tex->DrawLatex(0.66,0.82,"Points: PFClusters");
        sprintf(title,"Plots/EnergyFlow_0T_MC_Det_GenLevel_Comp.pdf"); if(save)c[4]->SaveAs(title);
        
        
        //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
        
        MCDetLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
        c[5] = tdrCanvas("MCDetLOverMCGenL",MCDetLOverMCGenL[0],4,11,kRectangular);
        c[5]->cd();
        for (int f=0; f<3; f++)MCDetLOverMCGenL[f]->Draw("same hist");
        //tex->DrawLatex(0.75,0.20,"Detec/GenLevel");
        leg[3]->Draw(); line->Draw();
        sprintf(title,"Plots/MCDetLOverMCGenL.pdf");
        if(save)c[5]->SaveAs(title);
    }
    

    
   /*
    //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
    
    MCDetLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
    c[5] = tdrCanvas("MCDetLOverMCGenL",MCDetLOverMCGenL[0],4,11,kRectangular);
    c[5]->cd();
    for (int f=0; f<3; f++)MCDetLOverMCGenL[f]->Draw("same hist");
    //tex->DrawLatex(0.75,0.20,"Detec/GenLevel");
    leg[3]->Draw(); line->Draw();
    sprintf(title,"Plots/MCDetLOverMCGenL.pdf");
    if(save)c[5]->SaveAs(title);
*/
/*
 //------------------ MC --------------------//
  Energy[1][0]->GetXaxis()->SetNdivisions(511); 
   
  c[1] = tdrCanvas("EnergyFlow_0T_MC",Energy[1][0],4,11,kRectangular);
  c[1]->cd();
  c[1]->SetLogy();
  for (int e=0; e<3; e++){  
    Energy[1][e]->Draw("same");
  }
  leg[1]->Draw();
  tex->DrawLatex(0.73,0.85,"Pythia8Monash");
  sprintf(title,"Plots/EnergyFlow_0T_Pythia8Monash.pdf");
  if(save)c[1]->SaveAs(title);
  //------------------ DATA MC Compare EFlow--------------------//
  c[2] = tdrCanvas("EnergyFlow_Data_MC_PFClusters",Energy[0][0],4,11,kRectangular);
  c[2]->cd();c[2]->SetLogy();
  for (int f=0; f<6; f++){  
    Energy[f][2]->Draw("same");
  }
  leg[2]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_0T_Data_MC_PFClusters.pdf"); if(save)c[2]->SaveAs(title);
  
 
  //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
  if(ratioMC){
  Energy[6][3]->SetMinimum(1.);
  Energy[6][3]->SetMaximum(1e8);
  Energy[6][3]->GetXaxis()->SetNdivisions(511);
  c[4] = tdrCanvas("EnergyFlow_0T_MC_Det_GenLevel_Comp",Energy[6][3],4,11,kRectangular);
  c[4]->cd();c[4]->SetLogy();
  
  for (int f=1; f<11; f++){
    
    if(f<6){
      Energy[f][2]->Scale( (pow(10,(f-1))) ) ;
      Energy[f][2]->Draw("same E0");
      //Energy[f][2]->Scale( 1./(pow(10,(f-1))) ) ;

    }
    if(f>5 && f<11){
      Energy[f][3]->Scale( (pow(10,(f-6))) ) ;
      Energy[f][3]->Draw("same hist");
      //Energy[f][3]->Scale( 1./ (pow(10,(f-6))) ) ;
    }
  }
  leg[4]->Draw(); tex->DrawLatex(0.66,0.87,"Lines : Gen Level");
  leg[4]->Draw(); tex->DrawLatex(0.66,0.82,"Points: PFClusters");
  sprintf(title,"Plots/EnergyFlow_0T_MC_Det_GenLevel_Comp.pdf"); if(save)c[4]->SaveAs(title);
 
 
 
 //------------------ MC Detlevel and Genlevel Compare EFlow was Fixed-------------------//
 if(ratioMC){
 Energy[5][3]->SetMinimum(1.);
 Energy[5][3]->SetMaximum(1e8);
 Energy[5][3]->GetXaxis()->SetNdivisions(511);
 c[4] = tdrCanvas("EnergyFlow_0T_MC_Det_GenLevel_Comp",Energy[5][3],4,11,kRectangular);
 c[4]->cd();c[4]->SetLogy();
 
 for (int f=2; f<8; f++){
 
 if(f<5){
 Energy[f][2]->Scale( (pow(10,(f-1))) ) ;
 Energy[f][2]->Draw("same E0");
 //Energy[f][2]->Scale( 1./(pow(10,(f-1))) ) ;
 
 }
 if(f>4 && f<8){
 Energy[f][3]->Scale( (pow(10,(f-6))) ) ;
 Energy[f][3]->Draw("same hist");
 //Energy[f][3]->Scale( 1./ (pow(10,(f-6))) ) ;
 }
 }
 leg[4]->Draw(); tex->DrawLatex(0.66,0.87,"Lines : Gen Level");
 leg[4]->Draw(); tex->DrawLatex(0.66,0.82,"Points: PFClusters");
 sprintf(title,"Plots/EnergyFlow_0T_MC_Det_GenLevel_Comp.pdf"); if(save)c[4]->SaveAs(title);
 
 
 //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
 
 MCDetLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
 c[5] = tdrCanvas("MCDetLOverMCGenL",MCDetLOverMCGenL[0],4,11,kRectangular);
 c[5]->cd();
 for (int f=0; f<3; f++)MCDetLOverMCGenL[f]->Draw("same hist");
 //tex->DrawLatex(0.75,0.20,"Detec/GenLevel");
 leg[3]->Draw(); line->Draw();
 sprintf(title,"Plots/MCDetLOverMCGenL.pdf");
 if(save)c[5]->SaveAs(title);
 }
 

 

  //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
 
  MCDetLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
  c[5] = tdrCanvas("MCDetLOverMCGenL",MCDetLOverMCGenL[0],4,11,kRectangular);
  c[5]->cd();
  for (int f=0; f<5; f++)MCDetLOverMCGenL[f]->Draw("same hist");
  //tex->DrawLatex(0.75,0.20,"Detec/GenLevel");
  leg[3]->Draw(); line->Draw();
  sprintf(title,"Plots/MCDetLOverMCGenL.pdf");
  if(save)c[5]->SaveAs(title);
  }
    
  //------------------ All Bunch vs Fat Bunch--------------------//
  c[6] = tdrCanvas("EnergyFlow_Data_allBx_vs_fatBx_PFClusters",Energy[0][0],4,11,kRectangular);
  c[6]->cd();c[6]->SetLogy();
  Energy[0][2]->Draw("same");
  Energy[11][2]->Draw("same");
  leg[5]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_Data_allBx_vs_fatBx_PFClusters.pdf"); if(save)c[6]->SaveAs(title);
    
    */
    
    //------------------ Ratio of All Bunch vs Fat Bunch--------------------//
/*  TH1D * RatioAllvsFat = new TH1D("RatioAllvsFat","RatioAllvsFat",nEtaBins,EtaBins);
  RatioAllvsFat->GetXaxis()->SetTitle("#eta");
  RatioAllvsFat->GetYaxis()->SetTitle("AllBx/FatBx");
  RatioAllvsFat->Divide(Energy[0][2],Energy[11][2],1,1,"B");
    
  RatioAllvsFat->SetLineColor(kBlack);
  RatioAllvsFat->SetMarkerColor(kBlack);
  //RatioAllvsFat->Divide(Energy[11][2]);
  //RatioAllvsFat->SetMinimum(0.5);
  //RatioAllvsFat->SetMaximum(1.5);*/
  /*
    RatioAllvsFat[0]->GetXaxis()->SetNdivisions(511);
  c[7] = tdrCanvas("RatioAllvsFat",RatioAllvsFat[0],4,11,kRectangular);
  c[7]->cd();//c[7]->SetLogy();
  RatioAllvsFat[0]->Draw("same hist");
  leg[7]->Draw(); line->Draw();
  sprintf(title,"Plots/Ratio_of_EnergyFlow_Data_allBx_vs_fatBx_PFClusters.pdf"); if(save)c[7]->SaveAs(title);
*/
  //------------------ Single Bunch--------------------//
    /*
  c[8] = tdrCanvas("EnergyFlow_Data_BPTX_Trigger_PFClusters",Energy[0][0],4,11,kRectangular);
  c[8]->cd();c[8]->SetLogy();
  Energy[12][2]->Draw("same");
  Energy[12][2]->SetMinimum(0.1);
  Energy[12][2]->SetMaximum(1e3);
  Energy[0][2]->Draw("same");
  Energy[13][2]->Draw("same");
  Energy[14][2]->Draw("same");
  leg[6]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_Data_BPTX_Trigger_PFClusters.pdf"); if(save)c[8]->SaveAs(title);
*/
   
    /*
    //------------------ DATA ZBvsL1PlusZB Compare EFlow--------------------//
    c[9] = tdrCanvas("EnergyFlow_DataExpr_DataZBpL1_PFClusters",Energy[0][0],4,11,kRectangular);
    c[9]->cd();c[9]->SetLogy();

        Energy[0][2]->Draw("same");
        Energy[15][2]->Draw("same");
    
    leg[8]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
    sprintf(title,"Plots/EnergyFlow_0T_DataExpr_DataZBpL1_PFClusters.pdf"); if(save)c[9]->SaveAs(title);
   
    
  */
  // //------------------ DATA MC Compare PFClusterEnergy Dist--------------------//
  // int count=0;

  // leg[3]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  // leg[3]->AddEntry(PFClustersEnergy[0][0],"DATA","PL");
  // leg[3]->AddEntry(PFClustersEnergy[1][0],"pythia8Monash","PL");
  // leg[3]->AddEntry(PFClustersEnergy[2][0],"pythia8MBR","PL");
  // leg[3]->AddEntry(PFClustersEnergy[3][0],"herwigpp","PL");
  // leg[3]->AddEntry(PFClustersEnergy[4][0],"epos","PL");
  // leg[3]->AddEntry(PFClustersEnergy[5][0],"qgsjetII","PL");
    
  // for (int k=0; k<nEtaBins; k++){
  // 	sprintf(title,"cPFClusterEnergy_Comp_Pythia8Monash_EtaBin%d",k);
  // 	c[3+count] = tdrCanvas(title,PFClustersEnergy[0][0],4,11,kRectangular);
  // 	c[3+count]->cd();
  // 	c[3+count]->SetLogy();
  // 	for (int f=0; f<6; f++)PFClustersEnergy[f][k]->Draw("same");
  // 	sprintf(title,"%3.3f<#eta<%3.3f",EtaBins[k],EtaBins[k+1]);
  // 	tex->DrawLatex(0.70,0.86,title);
  // 	leg[3]->Draw();
  // 	sprintf(title,"Plots/cPFClusterEnergy_Comp_Data_MC_EtaBin%d.pdf",k);
  // 	if(save)c[3+count]->SaveAs(title);
  // 	count++;
  // }
  
  fOutFile->Write();
  //fOutFile->Close();
  
}
int getBin(float x, float boundaries[],int b)
{
  int i;
  int n = b;
  if (x<boundaries[0] || x>=boundaries[n])
    return -1;
  for(i=0;i<n;i++)
    {
      if (x>=boundaries[i] && x<boundaries[i+1])
	return i;
    }
  return 0;
}
TH1F* convertbining(TH1F* h1,TH1F* h2){                                                                                                                                                      
  double content;
  double error;
  double Width;
  const size_t N =h1->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=N;i++)                                                                                            
    {
      Width   = h2->GetBinWidth(i);
      content = h1->GetBinContent(i);
      error   = h1->GetBinError(i);
      h2->SetBinContent(i,content*(1./Width));
      h2->SetBinError(i, error*(1./Width));
           
    }
  return h2;
  
}
// for (int f=0; f<ftyp; f++){
//   sprintf(title,"PFClustersEnergyvsEta__onlyHCAL_%s",fname[f].c_str());
//   PFClustersEnergyvsEta_[f][0]=new TH1F(title,title,nHBins,HBins);
//   PFClustersEnergyvsEta_[f][0]->Sumw2();

//   sprintf(title,"PFClustersEnergyvsEta__HCALandECAL_%s",fname[f].c_str());
//   PFClustersEnergyvsEta_[f][1]=new TH1F(title,title,nHBins,HBins);
//   PFClustersEnergyvsEta_[f][1]->Sumw2();

//     ///Only ECAL
//     sprintf(title,"PFClustersEnergyvsEta__ECAL_%s",fname[f].c_str());
//     PFClustersEnergyvsEta_[f][2]=new TH1F(title,title,nHBins,HBins);
//     PFClustersEnergyvsEta_[f][2]->Sumw2();

//     ///Only HF
//     sprintf(title,"PFClustersEnergyvsEta__HF_%s",fname[f].c_str());
//     PFClustersEnergyvsEta_[f][3]=new TH1F(title,title,nHBins,HBins);
//     PFClustersEnergyvsEta_[f][3]->Sumw2();

      
//   sprintf(title,"RecHitEnergyvsEta__%s",fname[f].c_str());
//   RecHitEnergyvsEta_[f][0]=new TH1F(title,title,nHBins,HBins);
//   RecHitEnergyvsEta_[f][0]->Sumw2();

//     ///ECAL+HCAL
//     sprintf(title,"RecHitEnergyvsEtaEcalHcal__%s",fname[f].c_str());
//     RecHitEnergyvsEtaEcalHcal_[f][0]=new TH1F(title,title,nHBins,HBins);
//     RecHitEnergyvsEtaEcalHcal_[f][0]->Sumw2();
      
      
//   sprintf(title,"EnergyFlow__Gen_SelectEventforDetLevel_%s",fname[f].c_str());
//   Energy_GenDetLevel_[f][0]=new TH1F(title,title,nHBins,HBins);
//   Energy_GenDetLevel_[f][0]->Sumw2();

//     //HF
//     sprintf(title,"EnergyFlow__Gen_SelectEventforDetLevel_HF_%s",fname[f].c_str());
//     Energy_GenDetLevel_HF_[f][0]=new TH1F(title,title,nHBins,HBins);
//     Energy_GenDetLevel_HF_[f][0]->Sumw2();
      
      
//   for (int e=0; e<etyp; e++){
      
//     sprintf(title,"EnergyFlow__%s_%s",fname[f].c_str(),energyname[e].c_str());
//     Energy_[f][e]=new TH1F(title,title,nHBins,HBins);
//     Energy_[f][e]->Sumw2();
    
//   }
// }
// for (int f=0; f<ftyp; f++){
   
//   sprintf(title,"PFClustersEnergyvsEta_onlyHCAL_%s",fname[f].c_str());
//   PFClustersEnergyvsEta[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//   PFClustersEnergyvsEta[f][0]->Sumw2();

//   sprintf(title,"PFClustersEnergyvsEta_HCALandECAL_%s",fname[f].c_str());
//   PFClustersEnergyvsEta[f][1]=new TH1F(title,title,nEtaBins,EtaBins);
//   PFClustersEnergyvsEta[f][1]->Sumw2();

//     ///Only ECAL
//     sprintf(title,"PFClustersEnergyvsEta_ECAL_%s",fname[f].c_str());
//     PFClustersEnergyvsEta[f][2]=new TH1F(title,title,nEtaBins,EtaBins);
//     PFClustersEnergyvsEta[f][2]->Sumw2();

//     ///Only HF
//     sprintf(title,"PFClustersEnergyvsEta_HF_%s",fname[f].c_str());
//     PFClustersEnergyvsEta[f][3]=new TH1F(title,title,nEtaBins,EtaBins);
//     PFClustersEnergyvsEta[f][3]->Sumw2();
      
      
//   sprintf(title,"RecHitEnergyvsEta_%s",fname[f].c_str());
//   RecHitEnergyvsEta[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//   RecHitEnergyvsEta[f][0]->Sumw2();

//     ///ECAL+HCAL
//     sprintf(title,"RecHitEnergyvsEtaEcalHcal_%s",fname[f].c_str());
//     RecHitEnergyvsEtaEcalHcal[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//     RecHitEnergyvsEtaEcalHcal[f][0]->Sumw2();

    
//   sprintf(title,"EnergyFlow_Gen_SelectEventforDetLevel_%s",fname[f].c_str());
//   Energy_GenDetLevel[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//   Energy_GenDetLevel[f][0]->Sumw2();
    
//     //HF
//     sprintf(title,"EnergyFlow_Gen_SelectEventforDetLevel_HF_%s",fname[f].c_str());
//     Energy_GenDetLevel_HF[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//     Energy_GenDetLevel_HF[f][0]->Sumw2();
      
//   sprintf(title,"Energy_CorFactor_PFCluster_onlyHCAL_GendEdEtaDetlevelEvnSelectOverRecdEdEta_%s",fname[f].c_str());
//   Energy_CorFactor[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//   Energy_CorFactor[f][0]->Sumw2();
    
//   sprintf(title,"Energy_CorFactor_PFCluster_HCALandECAL_GendEdEtaDetlevelEvnSelectOverRecdEdEta_%s",fname[f].c_str());
//   Energy_CorFactor[f][1]=new TH1F(title,title,nEtaBins,EtaBins);
//   Energy_CorFactor[f][1]->Sumw2();
    
//   sprintf(title,"Energy_CorFactor_RecHit_HCAL_GendEdEtaDetlevelEvnSelectOverRecdEdEta_%s",fname[f].c_str());
//   Energy_CorFactor[f][2]=new TH1F(title,title,nEtaBins,EtaBins);
//   Energy_CorFactor[f][2]->Sumw2();

//     sprintf(title,"Energy_CorFactor_RecHit_HCALandECAL_GendEdEtaDetlevelEvnSelectOverRecdEdEta_%s",fname[f].c_str());
//     Energy_CorFactor[f][3]=new TH1F(title,title,nEtaBins,EtaBins);
//     Energy_CorFactor[f][3]->Sumw2();


//     //HF
//     sprintf(title,"Energy_CorFactor_PFCluster_HF_GendEdEtaDetlevelEvnSelectOverRecdEdEta_%s",fname[f].c_str());
//     Energy_CorFactor_HF[f][0]=new TH1F(title,title,nEtaBins,EtaBins);
//     Energy_CorFactor_HF[f][0]->Sumw2();
      

//   for (int e=0; e<etyp; e++){
      
//     sprintf(title,"EnergyFlow_%s_%s",fname[f].c_str(),energyname[e].c_str());
      
//     Energy[f][e]=new TH1F(title,title,nEtaBins,EtaBins);
//     Energy[f][e]->Sumw2();
        
//   }
// }
  

  
//  ///Sadece ECAL
//   for (int f=0; f<ftyp; f++){
        
//       for (int eta=0; eta<nEtaBins; eta++){
            
//           sprintf(title,"PFClusterEnergy_ECAL_%s_EtaBin%d",fname[f].c_str(),eta);
//           if ((eta>=0 && eta<=3) || (eta>=24 && eta<=28)){
//               PFClustersEnergy_ECAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           if ((eta>=4 && eta<=7) || (eta>=20 && eta<=23)){
//               PFClustersEnergy_ECAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           if (eta>=8 && eta<=19){
//               PFClustersEnergy_ECAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           PFClustersEnergy_ECAL[f][eta]->Sumw2();
//       }
//   }
//   //////

//   ///Sadece HCAL
//   for (int f=0; f<ftyp; f++){
        
//       for (int eta=0; eta<nEtaBins; eta++){
            
//           sprintf(title,"PFClusterEnergy_HCAL_%s_EtaBin%d",fname[f].c_str(),eta);
//           if ((eta>=0 && eta<=3) || (eta>=24 && eta<=28)){
//               PFClustersEnergy_HCAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           if ((eta>=4 && eta<=7) || (eta>=20 && eta<=23)){
//               PFClustersEnergy_HCAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           if (eta>=8 && eta<=19){
//               PFClustersEnergy_HCAL[f][eta]=new TH1F(title,title,50,0.,5);
//           }
//           PFClustersEnergy_HCAL[f][eta]->Sumw2();
//       }
//   }
//   //////

// RecHitHCALEtaSums[i]=0.;
// RecHitECALEtaSums[i]=0.;
// RecHitHFEtaSums[i]=0.;
	
// PFCandEtaSums[i]=0.;
	
// PFClustersECALEtaSums[i]=0.;
// PFClustersHCALEtaSums[i]=0.;
// PFClustersHFEtaSums[i]=0.;
	
// PFClustersEtaSumsUp[i]=0.;
// PFClustersECALEtaSumsUp[i]=0.;
// PFClustersHCALEtaSumsUp[i]=0.;
// PFClustersHFEtaSumsUp[i]=0.;
	
// PFClustersEtaSumsLo[i]=0.;
// PFClustersECALEtaSumsLo[i]=0.;
// PFClustersHCALEtaSumsLo[i]=0.;
// PFClustersHFEtaSumsLo[i]=0.;
	
// PFClustersEtaSumsRaw[i]=0.;
// PFClustersECALEtaSumsRaw[i]=0.;
// PFClustersHCALEtaSumsRaw[i]=0.;
// PFClustersHFEtaSumsRaw[i]=0.;
	
// GenEtaSums_1GeV[i]=0.;
// GenEtaSums_2GeV[i]=0.;
// GenEtaSums_HadronElectronCut[i]=0.;
// GenEtaSums_EM[i]=0.;
// GenEtaSums_EM_withCuts[i]=0.;
// GenEtaSums_Had[i]=0.;
// GenEtaSums_Had_withCuts[i]=0.;
// GenEtaSumsDetEvntSelct[i]=0.;
// GenEtaSumsDetEvntSelctHF[i]=0.;
////Sadece ECAL
/*
  for (int f=0; f<ftyp; f++){
        
  for (int eta=0; eta<nEtaBins; eta++){
  PFClustersEnergy_ECAL[f][eta]->GetXaxis()->SetTitle("Energy");
  PFClustersEnergy_ECAL[f][eta]->GetYaxis()->SetTitle("1/N");
  // PFClustersEnergy_ECAL[f][eta]->SetMinimum(1e-5);
  // PFClustersEnergy_ECAL[f][eta]->SetMaximum(10);
  if (f == 0){
  PFClustersEnergy_ECAL[f][eta]->SetLineColor(kGreen);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerColor(kGreen);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerStyle(21);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerSize(0.5);
  }
  if (f ==1){
  PFClustersEnergy_ECAL[f][eta]->SetLineColor(kRed);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerColor(kRed);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerStyle(25);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerSize(0.5);
  }
  if (f ==2){
  PFClustersEnergy_ECAL[f][eta]->SetLineColor(kRed);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerColor(kRed);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerStyle(26);
  PFClustersEnergy_ECAL[f][eta]->SetMarkerSize(0.5);
  }
            
  PFClustersEnergy_ECAL[f][eta]->Scale(1./TotNofEvent[f]);
  PFClustersEnergy_ECAL[f][eta]->SetMinimum(1e-5);
  PFClustersEnergy_ECAL[f][eta]->SetMaximum(10);
  }
  }
  ////
  
    
  ////Sadece HCAL
  for (int f=0; f<ftyp; f++){
        
  for (int eta=0; eta<nEtaBins; eta++){
  PFClustersEnergy_HCAL[f][eta]->GetXaxis()->SetTitle("Energy");
  PFClustersEnergy_HCAL[f][eta]->GetYaxis()->SetTitle("1/N");
  // PFClustersEnergy_HCAL[f][eta]->SetMinimum(1e-5);
  // PFClustersEnergy_HCAL[f][eta]->SetMaximum(10);
  if (f == 0){
  PFClustersEnergy_HCAL[f][eta]->SetLineColor(kMagenta);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerColor(kMagenta);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerStyle(22);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerSize(0.5);
  }
  if (f ==1){
  PFClustersEnergy_HCAL[f][eta]->SetLineColor(30);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerColor(30);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerStyle(25);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerSize(0.5);
  }
  if (f ==2){
  PFClustersEnergy_HCAL[f][eta]->SetLineColor(kRed);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerColor(kRed);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerStyle(26);
  PFClustersEnergy_HCAL[f][eta]->SetMarkerSize(0.5);
  }
            
  PFClustersEnergy_HCAL[f][eta]->Scale(1./TotNofEvent[f]);
  PFClustersEnergy_HCAL[f][eta]->SetMinimum(1e-5);
  PFClustersEnergy_HCAL[f][eta]->SetMaximum(10);
  }
  }
*/
