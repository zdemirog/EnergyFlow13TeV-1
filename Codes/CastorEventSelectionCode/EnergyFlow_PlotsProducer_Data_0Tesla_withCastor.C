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


void EnergyFlow_PlotsProducer_Data_0Tesla_withCastor(){
  gROOT->ProcessLine(".L tdrstyle_mod14.C");
  setTDRStyle();
    
  gDirectory->DeleteAll();
  gROOT->ForceStyle();
  gStyle->SetOptFit(0);
 // gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");
  gStyle->SetOptStat(111111);
  bool save =true;
  bool ratioMC=false;
  string filenames[4] ={"EFlowTreeRootFile/EFlow_DetLevel_DataZeroBiasAllBx_tree_0Tesla.root",
			"EFlowTreeRootFile/EFlow_DetLevel_herwigpp_tree_0Tesla.root",
			"root2/EFlow_DetLevel_DataZeroBiasAllBx_tree_0Tesla_withCastor.root",
			"root2/EFlow_DetLevel_herwigpp_tree_0Tesla_withCastor.root"	 
  };
    
    
  string energyname[4]={"RecHit","CaloTower","PFClusters","Gen"};
  string fname[4]={"ZeroBias1_0T","pythia8Monash_0T","ZeroBias1_0T_withCastor","pythia8Monash_0T_withCastor"};
  TFile *fOutFile = new TFile("EFlow_Histo.root","RECREATE");
    
    
    
  static const Int_t ftyp =4 ;// Data and MC
  static const Int_t etyp = 4 ;
    
  char title[999];
    
  TH1F * Energy_[ftyp+1][etyp+1];
  TH1F * Energy[ftyp+1][etyp+1];  
  TH1F * DataDetLOverMCDetL[ftyp+1];
  TH1F * MCDetLOverMCGenL[ftyp+1];
  TH1F * DataGenLOverMCGenL[ftyp+1];
  TH1F * PFClustersEnergy[ftyp+1][nEtaBins+1];
  TH1F * RatioZBiasvsMonash[ftyp+1];
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

    sprintf(title,"%sExpress%sZeroBias",fname[0].c_str(),fname[f].c_str());
    RatioZBiasvsMonash[f]=new TH1F(title,title,nEtaBins,EtaBins);
    RatioZBiasvsMonash[f]->Sumw2();


    for (int e=0; e<etyp; e++){
	    
      sprintf(title,"EnergyFlow__%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy_[f][e]=new TH1F(title,title,nHBins,HBins);
      Energy_[f][e]->Sumw2();
	    
      sprintf(title,"EnergyFlow_%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy[f][e]=new TH1F(title,title,nEtaBins,EtaBins);
      Energy[f][e]->Sumw2();
	    
    }
	
   /* for (int eta=0; eta<nEtaBins; eta++){
	    
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
   // if( f>5 && f<11) fChain[f]->SetBranchAddress("Gen", &Gen, &b_Gen);
   // else {
      fChain[f]->SetBranchAddress("RecHit", &RecHit, &b_RecHit);
      fChain[f]->SetBranchAddress("CaloTower", &CaloTower, &b_CaloTower);
      fChain[f]->SetBranchAddress("PFClusters", &PFClusters, &b_PFClusters);
     // fChain[f]->SetBranchAddress("CastorTower", &CastorTower, &CastorTower);
   // }
   
    // ----------------------- Event -----------------------//
    Long64_t nentries = fChain[f]->GetEntriesFast();
    cout<<"Entries "<<nentries<<endl;
    Long64_t nbytes = 0, nb = 0;

    TotNofEvent[f]= nentries;

    int maxevent=10000;
    for (Long64_t ev=0; ev<nentries;ev++) {
      
      //  for (Long64_t ev=0; ev<maxevent;ev++) {
	
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
     /* if (f>5 && f<11){
	// -----------------------  Gen  -----------------//
	for(unsigned long gen=0; gen<Gen->size(); gen++) {
	  GenEtaSums[gen]= GenEtaSums[gen] + Gen->at(gen);
	}
	for (int k=0;k<nEtaBins;k++){
	Energy_[f][3]->Fill(k,GenEtaSums[k]);
	} 
      } */ 
      
    //  else{
	


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
	///---------------Castor Tower---------------//
	//for(unsigned long cstr=0; cstr<CastorTower->size(); cstr++) {
	 // CastorTowerEtaSums[cstr]= CastorTowerEtaSums[cstr] + CastorTower->at(cstr);
//	} 
	//{"RecHit","CaloTower","PFClusters"};
	for (int k=0;k<nEtaBins;k++){
	  Energy_[f][0]->Fill(k,RecHitEtaSums[k]);
	  Energy_[f][1]->Fill(k,CaloEtaSums[k]);
	  Energy_[f][2]->Fill(k,PFClustersEtaSums[k]);
  //        Energy_[f][3]->Fill(k,CastorTowerEtaSums[k]);
	  ///PFClustersEnergy[f][k]->Fill(PFClustersEtaSums[k]);
	} 
	  
      //}//else
    }//event 
  }//File
    
  // Plots style
    
  //"Data","pythia8_Monash","pythia8_TuneMBR","herwigpp","epos","qgsjetII","Genpythia8","Genherwigpp","Genepos","GenqgsjetII"
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
      if (f == 1){
	PFClustersEnergy[f][eta]->SetLineColor(2);
	PFClustersEnergy[f][eta]->SetMarkerColor(2);
	PFClustersEnergy[f][eta]->SetMarkerStyle(20);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
    if (f == 2){
	PFClustersEnergy[f][eta]->SetLineColor(1);
	PFClustersEnergy[f][eta]->SetMarkerColor(1);
	PFClustersEnergy[f][eta]->SetMarkerStyle(20);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f == 3){
	PFClustersEnergy[f][eta]->SetLineColor(2);
	PFClustersEnergy[f][eta]->SetMarkerColor(2);
	PFClustersEnergy[f][eta]->SetMarkerStyle(20);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }

      PFClustersEnergy[f][eta]->Scale(1./TotNofEvent[f]);
      PFClustersEnergy[f][eta]->SetMinimum(1e-5);
      PFClustersEnergy[f][eta]->SetMaximum(10);
    }
	
  }*/
    
    
  //string energyname[4]={"RecHit","CaloTower","PFClusters","Gen"};
  //string fname[11]={"Data","pythia8_Monash","pythia8_TuneMBR","herwigpp","epos","qgsjetII","Genpythia8Monash","Genpythia8MBR","Genherwigpp","Genepos","GenqgsjetII"};
  for (int f=0; f<ftyp; f++){
      
    for (int e=0; e<etyp; e++){
      convertbining(Energy_[f][e],Energy[f][e]);
      Energy[f][e]->Scale(1./TotNofEvent[f]);
	  
      Energy[f][e]->GetXaxis()->SetTitle("#eta");
      Energy[f][e]->GetYaxis()->SetTitle("(1 / N) dE/d#eta (GeV)");
      Energy[f][e]->SetMinimum(0.4);
      Energy[f][e]->SetMaximum(2000.);
      //Data ZeroBias1 RecHit 3.8T
      if (f == 0 && (e ==0)){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(21);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Data ZeroBias1 CaloTower 3.8T
      if (f == 0 && e == 1){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(22);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Data ZeroBias1 PFClusters 3.8T
      if (f == 0 && e ==2){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(20);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash RecHit 3.8T
      if (f == 1 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash CaloTower 3.8T
      if (f == 1 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash PFClusters 3.8T
	  
      if (f == 1 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(27);
	Energy[f][e]->SetMarkerSize(2.);
      }
      //Data ZeroBias1 RecHit 0T
      if (f == 2 && (e ==0)){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(21);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Data ZeroBias1 CaloTower 0T
      if (f == 2 && e == 1){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(22);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Data ZeroBias1 PFClusters 0T
      if (f == 2 && e ==2){
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(20);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash RecHit 0T
      if (f == 3 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash  CaloTower 0T
      if (f == 3 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia8Monash PFClusters 0T
	  
      if (f == 3 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(27);
	Energy[f][e]->SetMarkerSize(2.);
      }
    }//Energy type

    RatioZBiasvsMonash[0]->Divide(Energy[0][2],Energy[1][2], 1.,1.,"B");
    RatioZBiasvsMonash[1]->Divide(Energy[2][2],Energy[3][2], 1.,1.,"B");
  }//File Type

 
    RatioZBiasvsMonash[0]->GetXaxis()->SetTitle("#eta");
    RatioZBiasvsMonash[0]->GetYaxis()->SetTitle("#frac{1}{N}#frac{dE}{d#eta}[ZeroBias1] / #frac{1}{N}#frac{dE}{d#eta}[Pythia8Monash]");
    RatioZBiasvsMonash[0]->SetMinimum(0.2);
    RatioZBiasvsMonash[0]->SetMaximum(3.1);
    RatioZBiasvsMonash[0]->SetLineColor(kBlue+3);
    RatioZBiasvsMonash[0]->SetLineWidth(4);
    RatioZBiasvsMonash[0]->SetLineStyle(1);
    RatioZBiasvsMonash[0]->SetMarkerColor(kBlue+3);
    RatioZBiasvsMonash[0]->GetYaxis()->SetTitleSize(0.052);
   
    RatioZBiasvsMonash[1]->GetXaxis()->SetTitle("#eta");
    RatioZBiasvsMonash[1]->GetYaxis()->SetTitle("#frac{1}{N}#frac{dE}{d#eta}[ZeroBias1] / #frac{1}{N}#frac{dE}{d#eta}[Pythia8Monash]");
    RatioZBiasvsMonash[1]->SetMinimum(0.2);
    RatioZBiasvsMonash[1]->SetMaximum(3.1);
    RatioZBiasvsMonash[1]->SetLineColor(kRed+2);
    RatioZBiasvsMonash[1]->SetLineWidth(4);
    RatioZBiasvsMonash[1]->SetLineStyle(1);
    RatioZBiasvsMonash[1]->SetMarkerColor(kRed+2);
    RatioZBiasvsMonash[1]->GetYaxis()->SetTitleSize(0.052);

  //////************************************************************************************////
  /////                           NOW PLOTS                                                 ///
  ////*************************************************************************************////
  TLine *line=new TLine(-5.191,1.,5.191,1.);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
      
  TLegend *leg[999];
  leg[0]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[0]->AddEntry(Energy[0][0],"RecHits","PL");
  leg[0]->AddEntry(Energy[0][1],"CaloTowers","PL");
  leg[0]->AddEntry(Energy[0][2],"PFClusters","PL");

  leg[1]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[1]->AddEntry(Energy[1][0],"RecHits","PL");
  leg[1]->AddEntry(Energy[1][1],"CaloTowers","PL");
  leg[1]->AddEntry(Energy[1][2],"PFClusters","PL");

  leg[2]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[2]->AddEntry(Energy[0][2],"Data ZeroBias1 0T","PL");
  leg[2]->AddEntry(Energy[1][2],"Pythia8Monash 0T","PL");


  leg[3]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[3]->AddEntry(RatioZBiasvsMonash[0],"Data ZeroBias1/Pythia8Monash","p");
    
  leg[4]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[4]->AddEntry(Energy[2][0],"RecHits","PL");
  leg[4]->AddEntry(Energy[2][1],"CaloTowers","PL");
  leg[4]->AddEntry(Energy[2][2],"PFClusters","PL");

  leg[5]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[5]->AddEntry(Energy[3][0],"RecHits","PL");
  leg[5]->AddEntry(Energy[3][1],"CaloTowers","PL");
  leg[5]->AddEntry(Energy[3][2],"PFClusters","PL");

  leg[6]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[6]->AddEntry(RatioZBiasvsMonash[0],"ZeroBias1/Pythia8Monash0T","p");
  leg[6]->AddEntry(RatioZBiasvsMonash[1],"ZeroBias1/Pythia8Monash0T with Castor","p");

  leg[7]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[7]->AddEntry(Energy[2][2],"Data ZeroBias1 0T with Castor","PL");
  leg[7]->AddEntry(Energy[3][2],"Pythia8Monash 0T with Castor","PL");

  leg[8]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[8]->AddEntry(Energy[0][2],"Data ZeroBias1 0T","PL");
  leg[8]->AddEntry(Energy[1][2],"Pythia8Monash 0T","PL");
  leg[8]->AddEntry(Energy[2][2],"Data ZeroBias1 0T with Castor","PL");
  leg[8]->AddEntry(Energy[3][2],"Pythia8Monash 0T with Castor","PL");


  TCanvas *c[999];

  
  //------------------ DATA 3p8T --------------------//
  Energy[0][0]->GetXaxis()->SetNdivisions(511); 
 
  c[0] = tdrCanvas("EnergyFlow_0T_DataZeroBias",Energy[0][0],4,11,kRectangular);
  c[0]->cd();
  c[0]->SetLogy();
  for (int e=0; e<3; e++){  
    Energy[0][e]->Draw("same");
  }
  tex->DrawLatex(0.75,0.85,"ZeroBias1 0T");
  leg[0]->Draw();
  sprintf(title,"Plots/EnergyFlow_0T_DataZeroBias.pdf");
  if(save)c[0]->SaveAs(title);
  //------------------ Zerobias 3p8T --------------------//
  Energy[1][0]->GetXaxis()->SetNdivisions(511); 
   
  c[1] = tdrCanvas("EnergyFlow_0T_Monash",Energy[1][0],4,11,kRectangular);
  c[1]->cd();
  c[1]->SetLogy();
  for (int e=0; e<3; e++){  
    Energy[1][e]->Draw("same");
  }
  leg[1]->Draw();
  tex->DrawLatex(0.73,0.85,"Pythia8Monash 0T");
  sprintf(title,"Plots/EnergyFlow_0T_Monash.pdf");
  if(save)c[1]->SaveAs(title);
  //------------------ DATA Express Zerobias  Compare EFlow--------------------//
  c[2] = tdrCanvas("EnergyFlow_Data_PFClusters_0T",Energy[0][0],4,11,kRectangular);
  c[2]->cd();c[2]->SetLogy();
  for (int f=0; f<2; f++){  
    Energy[f][2]->Draw("same");
  }
  leg[2]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_Data_PFClusters_0T.pdf"); if(save)c[2]->SaveAs(title);
 
  
  RatioZBiasvsMonash[0]->GetXaxis()->SetNdivisions(511);
  c[3] = tdrCanvas("RatioZeroBiasvsMonash_0T",RatioZBiasvsMonash[0],4,11,kRectangular);
  c[3]->cd();//c[7]->SetLogy();
  RatioZBiasvsMonash[0]->Draw("same hist");
  leg[3]->Draw(); line->Draw();
  sprintf(title,"Plots/RatioZeroBiasvsMonash_0T.pdf"); if(save)c[3]->SaveAs(title);


  c[4] = tdrCanvas("EnergyFlow_0T_ZeroBias_withCastor",Energy[0][0],4,11,kRectangular);
  c[4]->cd();
  c[4]->SetLogy();
  for (int e=0; e<3; e++){  
    Energy[2][e]->Draw("same");
  }
  tex->DrawLatex(0.75,0.85,"ZeroBias1 0T_withCastor");
  leg[4]->Draw();
  sprintf(title,"Plots/EnergyFlow_0T_ZeroBias_withCastor.pdf");
  if(save)c[4]->SaveAs(title);

  //------------------ Zero Bias 0T --------------------//
  Energy[1][0]->GetXaxis()->SetNdivisions(511); 
   
  c[5] = tdrCanvas("EnergyFlow_0T_Monash_withCastor",Energy[1][0],4,11,kRectangular);
  c[5]->cd();
  c[5]->SetLogy();
  for (int e=0; e<3; e++){  
    Energy[3][e]->Draw("same");
  }
  leg[5]->Draw();
  tex->DrawLatex(0.73,0.85,"Pythia8Monash 0T_withCastor");
  sprintf(title,"Plots/EnergyFlow_0T_Monash_withCastor.pdf");
  if(save)c[5]->SaveAs(title);

 
  
  RatioZBiasvsMonash[1]->GetXaxis()->SetNdivisions(511);
  c[6] = tdrCanvas("RatioZeroBiasvsMonash_0T_withCastor",RatioZBiasvsMonash[1],4,11,kRectangular);
  c[6]->cd();//c[7]->SetLogy();
  RatioZBiasvsMonash[0]->Draw("same hist");
  RatioZBiasvsMonash[1]->Draw("same hist");
  leg[6]->Draw(); line->Draw();
  sprintf(title,"Plots/RatioZeroBiasvsMonash_0T_withCastor.pdf"); if(save)c[6]->SaveAs(title);


  //------------------ DATA Express Zerobias  Compare EFlow--------------------//
  c[8] = tdrCanvas("EnergyFlow_Data_MC_PFClusters_0T_withCastor",Energy[0][0],4,11,kRectangular);
  c[8]->cd();c[8]->SetLogy();
    Energy[2][2]->Draw("same");
    Energy[3][2]->Draw("same");
  leg[7]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_Data_MC_PFClusters_0T_withCastor.pdf"); if(save)c[8]->SaveAs(title);



 //------------------ DATA Express Zerobias  Compare EFlow--------------------//
  c[9] = tdrCanvas("EnergyFlow_Data_MC_PFClusters_0T_compare",Energy[0][0],4,11,kRectangular);
  c[9]->cd();c[9]->SetLogy();
    Energy[0][2]->Draw("same");
    Energy[1][2]->Draw("same");
    Energy[2][2]->Draw("same");
    Energy[3][2]->Draw("same");
  leg[8]->Draw(); tex->DrawLatex(0.75,0.20,"PFClusters");
  sprintf(title,"Plots/EnergyFlow_Data_MC_PFClusters_0T_compare.pdf"); if(save)c[9]->SaveAs(title);




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
