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
#include "TGraphAsymmErrors.h"
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
TH1F* convertbiningMCGen(TH1F* h1,TH1F* h2);
TH1F* convertbiningMCDet(TH1F* h1,TH1F* h2);
TH1F* convertbiningData(TH1F* h1,TH1F* h2);
TH1F* unfoldbinbybin(TH1F* h1,TH1F* h2,TH1F* h3,TH1F* h4,TH1F* h5,TH1F* h6);
TH1F* modelunct(TH1F* h1,TH1F* h2,TH1F* h3,TH1F* h4);
//____________________________________________________________________________________
//-------------------------------------------------------------
//https://root.cern.ch/root/html/tutorials/tree/hvector.C.html  
//
//-------------------------------------------------------------
//_________________________________________________________________________________________
//


void  EnergyFlow_PlotsProducerHFSDEnhancedVeto()

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
  bool DataMCCompare=false;
  bool ratioMC=false;
  bool energydist=true;
  bool unfold =false;
  bool unfoldratio=false;
  string filenames[15] ={
      "/Volumes/BACKUP2/SDEnhancedNew/Data/EFlow_DetLevel_MinBias_DataZeroBias1_13TeVAllBx_0T_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_TuneMonash_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_TuneMBR_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_Herwig_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_Epos0T_13TeV-pythia8_AllCollect_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_qgsjetII_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",
      
      "/Volumes/BACKUP2/SDEnhancedNew/GenLevelMC/EFlow_GenLevel_MinBias_TuneMonash_13TeV-pythia8_0Tesla_SDEnhanced_FromBenoit_Full.root",
      "/Volumes/BACKUP2/SDEnhancedNew/GenLevelMC/EFlow_GenLevel_MinBias_TuneMBR_13TeV-pythia8_0Tesla_SDEnhanced_FromBenoit_Full.root",
      "/Volumes/BACKUP2/SDEnhancedNew/GenLevelMC/EFlow_GenLevel_MinBias_Herwig_13TeV-pythia8_0Tesla_SDEnhanced_FromBenoit_Full.root",
      "/Volumes/BACKUP2/SDEnhancedNew/GenLevelMC/EFlow_GenLevel_MinBias_Epos0T_13TeV-pythia8_0Tesla_SDEnhanced_FromBenoit_Full.root",
      "/Volumes/BACKUP2/SDEnhancedNew/GenLevelMC/EFlow_GenLevel_MinBias_qgsjetII_13TeV-pythia8_0Tesla_SDEnhanced_FromBenoit_Full.root",
      
      "/Volumes/BACKUP2/SDEnhancedNew/Data/EFlow_DetLevel_MinBias_DataZeroBias1_13TeVAllBx_0T_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",//5 GeV cut Data
      "/Volumes/BACKUP2/SDEnhancedNew/Data/EFlow_DetLevel_MinBias_DataZeroBias1_13TeVAllBx_0T_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5p5GeV.root",//5.5 GeV cut Data
      
      "/Volumes/BACKUP2/SDEnhancedNew/MCDetLevel/EFlow_DetLevel_MinBias_TuneMonash_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5GeV.root",// 5 GeV cut MC
      "/Volumes/BACKUP2/SDEnhancedNew/Monash5p5/EFlow_DetLevel_MinBias_TuneMonash_13TeV-pythia8_JINSTBinning_SDEnhanced_NoiseCut4GeV_Emin5p5GeV.root"// 5.5 GeV cut MC
      //"EFlow_Histo_MagnetON.root"
      
      
  };
    
    
  string energyname[6]={"RecHit","CaloTower","PFClusters","Gen","PFCandidate","DetGenGen"};
  string fname[15]={"Data","pythia8_Monash","pythia8_TuneMBR","herwigpp","epos","qgsjetII","Genpythia8Monash","Genpythia8MBR","GenTuneCUETP8M1withoutMBR","Genepos","GenqgsjetII",
		    "DataCenter","DataUp","MCCenter", "MCUp"};
 TFile *fOutFile = new TFile("EFlow_VetoHisto.root","RECREATE");
 char e1[999]; sprintf(e1,"EnergyFlow_Data_CaloTower");
 char e2[999]; sprintf(e2,"EnergyFlow_epos_CaloTower");
  static const Int_t ftyp =15 ;// Data and MC
  static const Int_t etyp = 6 ;
  float etanorm = 1.;
  //https://indico.cern.ch/event/452328/contribution/1/attachments/1165639/1680640/Inelastic_Cross_Section_-_151006.pdf
  
  

  char title[999];
  char hname[999];
    
  TH1F * Energy_[ftyp+1][etyp+1];
  TH1F * Energy[ftyp+1][etyp+1];
  TH1F * DataDetLOverMCDetL[ftyp+1];
  TH1F * MCDetLOverMCGenL[ftyp+1];
  TH1F * DataGenLOverMCGenL[ftyp+1];
  TH1F * PFClustersEnergy[ftyp+1][nEtaBins+1];
  TF1  *fPFClustersEnergy[ftyp+1][nEtaBins+1];
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
      
    for (int e=0; e<etyp; e++){
	    
      sprintf(title,"EnergyFlow__%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy_[f][e]=new TH1F(title,title,nHBins,HBins);
      Energy_[f][e]->Sumw2();
	    
      sprintf(title,"EnergyFlow_%s_%s",fname[f].c_str(),energyname[e].c_str());
      Energy[f][e]=new TH1F(title,title,nEtaBins,EtaBins);
      Energy[f][e]->Sumw2();
	    
    }
	
    for (int eta=0; eta<nEtaBins; eta++){
	    
      sprintf(title,"EnergyDist_%s_EtaBin%d",fname[f].c_str(),eta);
      if ((eta>=0 && eta<=3) || (eta>=24 && eta<=30)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,800,0.,800.);
	fPFClustersEnergy[f][eta]=new TF1(title,"expo",5,500);
      }
      if ((eta>=4 && eta<=5) || (eta>=22 && eta<=23)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,600,0.,600.);
	fPFClustersEnergy[f][eta]=new TF1(title,"expo",5,400);
      }
      if ((eta>=6 && eta<=7) || (eta>=20 && eta<=21)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,300,0.,300.);
	fPFClustersEnergy[f][eta]=new TF1(title,"expo",5,200);
      }
      if ((eta>=8 && eta<=9) || (eta>=18 && eta<=19)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,200,0.,200.);
	fPFClustersEnergy[f][eta]=new TF1(title,"expo",5,150);
      }
      if ((eta>=10 && eta<=17)){
	PFClustersEnergy[f][eta]=new TH1F(title,title,50,0.,50.);
	fPFClustersEnergy[f][eta]=new TF1(title,"expo",5,40);
      }
      PFClustersEnergy[f][eta]->Sumw2();
    }
  }
  //-------------------------------------------------
  // Add CUETM1 + MBR from Sercan http://ssen.web.cern.ch/ssen/CMS-FSQ-15-006-RivetPlots/plotsMBRcuetm1/index.html 
  TH1F * EnergyCUETM1[1];
  sprintf(title,"EnergyFlow_PythiaCUETM1_MCGen");
  EnergyCUETM1[0]=new TH1F(title,title,nEtaBins,EtaBins);
  EnergyCUETM1[0]->Sumw2();
  EnergyCUETM1[0]->SetLineStyle(6);
  EnergyCUETM1[0]->SetLineWidth(5);
  EnergyCUETM1[0]->SetLineColor(7);
  EnergyCUETM1[0]->SetMarkerColor(7);

  double Valuecuetm1[]={8.415735e+01,1.141035e+02,1.525285e+02,2.023164e+02,2.659359e+02,3.388043e+02,5.986713e+02};
  double Errorcuetm1[]={1.897476e-01,2.609960e-01,3.538606e-01,4.864738e-01,6.441896e-01,9.039275e-01,8.508348e-01};
  
  const size_t NN =EnergyCUETM1[0]->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=NN;i++)                                                                                            
    {
      
      EnergyCUETM1[0]->SetBinContent(i,Valuecuetm1[i-1]);
      EnergyCUETM1[0]->SetBinError(i,Errorcuetm1[i-1]);
     

     // cout <<NN<<"  "<< Errorcuetm1[i-1]<<" "<<Valuecuetm1[i-1]<<" "<<i<<endl;
           
    }
  TH1F * DataGenLOverMCGenLCUETM1[1];
  sprintf(title,"DataGenLOverPythiaCUETM1GenL");
  DataGenLOverMCGenLCUETM1[0]=new TH1F(title,title,nEtaBins,EtaBins);
  DataGenLOverMCGenLCUETM1[0]->Sumw2();
  
  DataGenLOverMCGenLCUETM1[0]->SetLineColor(7);
  DataGenLOverMCGenLCUETM1[0]->SetLineWidth(5);
  DataGenLOverMCGenLCUETM1[0]->SetLineStyle(6);
  DataGenLOverMCGenLCUETM1[0]->SetMarkerColor(7);

  // Add CUETS1 from Paulo 
  TH1F * EnergyCUETS1[3];
  sprintf(title,"EnergyFlow_PythiaCUETS1_MCGen");
  EnergyCUETS1[0]=new TH1F(title,title,nEtaBins,EtaBins);
  sprintf(title,"EnergyFlow_PythiaCUETS1_MCGenUP");
  EnergyCUETS1[1]=new TH1F(title,title,nEtaBins,EtaBins);
  sprintf(title,"EnergyFlow_PythiaCUETS1_MCGenDown");
  EnergyCUETS1[2]=new TH1F(title,title,nEtaBins,EtaBins);
  EnergyCUETS1[0]->Sumw2();
  EnergyCUETS1[1]->Sumw2();
  EnergyCUETS1[2]->Sumw2();
  EnergyCUETS1[0]->SetLineStyle(7);
  EnergyCUETS1[0]->SetLineWidth(5);
  EnergyCUETS1[0]->SetLineColor(9);
  EnergyCUETS1[0]->SetMarkerColor(9);

  double Valuecuets1[]={88.01872,116.9152,155.7445,209.6226,279.0925,352.8993,618.8548};
  double Errorcuets1[]={4.597351,6.591072,8.087782,8.926368,15.00443,14.84799,17.22231};
  
  const size_t NNN =EnergyCUETS1[0]->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=NNN;i++)                                                                                            
    {
      
      EnergyCUETS1[0]->SetBinContent(i,Valuecuets1[i-1]);
      EnergyCUETS1[0]->SetBinError(i,Errorcuets1[i-1]);

      EnergyCUETS1[1]->SetBinContent(i,(Valuecuets1[i-1]+Errorcuets1[i-1]));
      EnergyCUETS1[1]->SetBinError(i,0);

      EnergyCUETS1[2]->SetBinContent(i,(Valuecuets1[i-1]-Errorcuets1[i-1]));
      EnergyCUETS1[2]->SetBinError(i,0);

     

      //cout <<NN<<"  "<< Errorcuets1[i-1]<<" "<<Valuecuets1[i-1]<<" "<<i<<endl;
           
    }
  TH1F * DataGenLOverMCGenLCUETS1[3];
  sprintf(title,"DataGenLOverPythiaCUETS1GenL");
  DataGenLOverMCGenLCUETS1[0]=new TH1F(title,title,nEtaBins,EtaBins);
  sprintf(title,"DataGenLOverPythiaCUETS1GenLUP");
  DataGenLOverMCGenLCUETS1[1]=new TH1F(title,title,nEtaBins,EtaBins);
  sprintf(title,"DataGenLOverPythiaCUETS1GenLDwn");
  DataGenLOverMCGenLCUETS1[2]=new TH1F(title,title,nEtaBins,EtaBins);
  DataGenLOverMCGenLCUETS1[0]->Sumw2();
  DataGenLOverMCGenLCUETS1[1]->Sumw2();
  DataGenLOverMCGenLCUETS1[2]->Sumw2();
  DataGenLOverMCGenLCUETS1[0]->SetLineColor(9);
  DataGenLOverMCGenLCUETS1[0]->SetLineWidth(5);
  DataGenLOverMCGenLCUETS1[0]->SetLineStyle(7);
  DataGenLOverMCGenLCUETS1[0]->SetMarkerColor(9);
  
  
  
  //--------------------------------------------------
  //--------------------------------------------------
  
  int decade = 0;
  TFile *file[ftyp+1];
  TTree *fChain[ftyp+1];
  
  float PileUp = 0.06; // Run Number 247324  Magnet Of
  float punorm  = 1. / (PileUp / (1-(TMath::Exp((-1)*PileUp))));
    
  //float punorm =1;
  int TotNofEvent[ftyp+1];
  for (int f=0; f<ftyp; f++){

    

    cout << " PU norm "<<punorm<<endl;
    // --- open file in EOS
    file[f]   = TFile::Open(filenames[f].c_str(),"READ");
   // file[f]   = new TFile(filenames[f].c_str(),"READ");

    cout<<"file : "<<filenames[f].c_str()<<endl;
  
    fChain[f] = (TTree*)file[f]->Get("EFlow");
    
    if (f<6 || f>10) {
	fChain[f]->SetBranchAddress("RecHit", &RecHit, &b_RecHit);
	fChain[f]->SetBranchAddress("CaloTower", &CaloTower, &b_CaloTower);
	fChain[f]->SetBranchAddress("SDHFSignalTower",&SDHFSignalTower);//New Part
	fChain[f]->SetBranchAddress("SDHFVetoTower",&SDHFVetoTower);//New Part
	fChain[f]->SetBranchAddress("PFClusters", &PFClusters, &b_PFClusters);
	fChain[f]->SetBranchAddress("PFCandidate", &PFCandidate, &b_PFCandidate);
	fChain[f]->SetBranchAddress("GenDetEventSelect", &GenDetEventSelect, &b_GenDetEventSelect);
      }
      if( f>5 && f<11 ) {
          fChain[f]->SetBranchAddress("Gen", &Gen, &b_Gen);
          fChain[f]->SetBranchAddress("SignalGen",&SignalGen);
          fChain[f]->SetBranchAddress("SignalVetoGen",&SignalVetoGen);
          
      }
     //else cout <<"No File for tree"<<endl;
   
      // ----------------------- Event -----------------------//
      Long64_t nentries = fChain[f]->GetEntriesFast();
      cout<<"Entries "<<nentries<<endl;
      Long64_t nbytes = 0, nb = 0;

      TotNofEvent[f]= nentries;
      float HFPlusHFMinusRecHit =0;
      float HFPlusHFMinusPFCl =0;
      float HFPlusHFMinusCalo =0;
      float HFPlusHFMinusGen =0;
      float CastorGen =0;
      float HFPlusHFMinusCandidate =0;
      float HFPlusHFMinusGenDetEtaSums=0;
      int maxevent=10000;
   
      for (Long64_t ev=0; ev<nentries;ev++) {
      
	//	for (Long64_t ev=0; ev<maxevent;ev++) {
	
	Long64_t iev = fChain[f]->LoadTree(ev);
	if (iev < 0)break;
	nb = fChain[f]->GetEntry(ev);   nbytes += nb;
	// cout<<"EVENT "<<ev<<endl;
	double progress = 10.0*ev/(1.0*nentries);
	int k = TMath::FloorNint(progress); 
	if (k > decade) 
	  cout<<10*k<<" %"<<endl;
	decade = k; 
	for (int i =0; i<BIN;i++) {
	  
	  RecHitEtaSums[i]=0.;
	  CaloEtaSums[i]=0.;
	  PFClustersEtaSums[i]=0.;
	  GenEtaSums[i]=0.;
	  PFCandEtaSums[i]=0.;
	  GenEtaSumsDetEvntSelct[i]=0.;
	}
      
	//HF Bin
	//0.,1.,2.,3.,4.,5.,6.,7.,23.,24.,25.,26.,27.,28.,29.,30.
      
	HFPlusHFMinusRecHit =0;
	HFPlusHFMinusPFCl =0;
	HFPlusHFMinusCalo =0;
	HFPlusHFMinusGen =0;
	HFPlusHFMinusCandidate =0;
	HFPlusHFMinusGenDetEtaSums=0;
	//cout <<RecHit->size()<<" "<<CaloTower->size()<<endl;
      
      
	if(f<6 || f>10){
	  
	  // ----------------------- RecHit --------------- //
	  for(unsigned long rec=0; rec<RecHit->size(); rec++) {
	    //for(unsigned long rec=1; (rec<RecHit->size())-1; rec++) {
	    RecHitEtaSums[rec]= RecHit->at(rec);
	    //RecHitEtaSums[rec-1]= RecHitECALEtaSums[rec-1] + RecHit->at(rec-1);
	  }  
	  //----------------------- SDHFVetoTower -------------//
	  for(unsigned long cal=0; cal<SDHFVetoTower->size(); cal++) {
	    //  for(unsigned long cal=1; cal<(CaloTower->size())-1; cal++) {
	    CaloEtaSums[cal]=SDHFVetoTower->at(cal);
	    //CaloEtaSums[cal-1]= CaloEtaSums[cal-1] + CaloTower->at(cal-1);
	  }  
	  //----------------------- PF Cluster --------------- //
	  for(unsigned long pfc=0; pfc<PFClusters->size(); pfc++) {
	    //  for(unsigned long pfc=1; pfc<(PFClusters->size())-1; pfc++) {
	    PFClustersEtaSums[pfc]=PFClusters->at(pfc);
	    //PFClustersEtaSums[pfc-1]= PFClustersEtaSums[pfc-1] + PFClusters->at(pfc-1);
	  }
	  //----------------------- PFCandidate --------------- //
	  for(unsigned long pfcan=0; pfcan<PFCandidate->size(); pfcan++) {
	  
	    PFCandEtaSums[pfcan]=PFCandidate->at(pfcan);
	  
	  }
	  //----------------------- DetGen --------------- //
	  for(unsigned long det=0; det<GenDetEventSelect->size(); det++) {
	  
	    GenEtaSumsDetEvntSelct[det]= GenDetEventSelect->at(det);
	  
	  } 
	  //{"RecHit","CaloTower","PFClusters"};
	
	  for (int k=0;k<nEtaBins;k++){
	    

	    HFPlusHFMinusRecHit = etanorm*(RecHitEtaSums[k] + RecHitEtaSums[29-k]);
	    HFPlusHFMinusPFCl   = etanorm*(PFClustersEtaSums[k] + PFClustersEtaSums[29-k]);
	    HFPlusHFMinusCalo   = etanorm*(CaloEtaSums[k] + CaloEtaSums[29-k]);
	    HFPlusHFMinusCandidate= etanorm*(PFCandEtaSums[k] + PFCandEtaSums[29-k]);
	    HFPlusHFMinusGenDetEtaSums= etanorm*(GenEtaSumsDetEvntSelct[k] + GenEtaSumsDetEvntSelct[29-k]);
	    Energy_[f][0]->Fill(k,HFPlusHFMinusRecHit);
	    Energy_[f][1]->Fill(k,HFPlusHFMinusCalo);
	    Energy_[f][2]->Fill(k,HFPlusHFMinusPFCl);
	    Energy_[f][4]->Fill(k,HFPlusHFMinusCandidate);
	    Energy_[f][5]->Fill(k,HFPlusHFMinusGenDetEtaSums);
	  
	  
	    // Energy_[f][0]->Fill(k,(etanorm * (RecHitEtaSums[k]+RecHitEtaSums[30-k])));
	    // Energy_[f][1]->Fill(k,(etanorm * (CaloEtaSums[k]+CaloEtaSums[30-k])));
	    // Energy_[f][2]->Fill(k,(etanorm * (PFClustersEtaSums[k] + PFClustersEtaSums[30-k])));
	  
	  
	    if(HFPlusHFMinusCalo>0){
	      PFClustersEnergy[f][k]->Fill(HFPlusHFMinusCalo);
	      }
	  
	  } 
	  
	}
	if (f>5 && f<11 ){
	  // -----------------------  Gen  -----------------//
	
        
        for(unsigned long gen=0; gen<SignalGen->size() ; gen++) {
            GenEtaSums[gen]= SignalGen->at(gen);
        }
        for(unsigned long gen=0; gen<Gen->size() ; gen++) {
            GenEtaSums_[gen]= Gen->at(gen);
        }
        for (int k=0;k<nEtaBins;k++){
            if (k==0){
                CastorGen=1./2. * (GenEtaSums_[k] + GenEtaSums_[29-k]);
                Energy_[f][3]->Fill(k,CastorGen);
            }
            if(k>0) {
                HFPlusHFMinusGen =etanorm * (GenEtaSums[k] + GenEtaSums[29-k]);
                Energy_[f][3]->Fill(k,HFPlusHFMinusGen);// HFPlus or HFMinus .
            }
        }
        
    }		
      }//event
  }//File
    
    // Plots style

 
  



  //"Data","pythia8_Monash","pythia8_TuneMBR","herwigpp","epos","qgsjetII","Genpythia8","Genherwigpp","Genepos","GenqgsjetII"
    
  for (int f=0; f<ftyp; f++){
    
    for (int eta=0; eta<nEtaBins; eta++){
      PFClustersEnergy[f][eta]->GetXaxis()->SetTitle("Energy");
      PFClustersEnergy[f][eta]->GetYaxis()->SetTitle("1/N");
      // PFClustersEnergy[f][eta]->SetMinimum(1e-5);
      // PFClustersEnergy[f][eta]->SetMaximum(10);
      if (f == 0){
	
    PFClustersEnergy[f][eta]->Scale(punorm);
	PFClustersEnergy[f][eta]->SetLineColor(1);
	PFClustersEnergy[f][eta]->SetMarkerColor(1);
	PFClustersEnergy[f][eta]->SetMarkerStyle(20);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==1 || f==10){
	PFClustersEnergy[f][eta]->SetLineColor(kBlue);
	PFClustersEnergy[f][eta]->SetMarkerColor(kBlue);
	PFClustersEnergy[f][eta]->SetMarkerStyle(24);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==2 || f==9){
	PFClustersEnergy[f][eta]->SetLineColor(kRed);
	PFClustersEnergy[f][eta]->SetMarkerColor(kRed);
	PFClustersEnergy[f][eta]->SetMarkerStyle(25);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==3 || f==8){
	PFClustersEnergy[f][eta]->SetLineColor(kYellow);
	PFClustersEnergy[f][eta]->SetMarkerColor(kYellow);
	PFClustersEnergy[f][eta]->SetMarkerStyle(26);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==4 || f==7){
	PFClustersEnergy[f][eta]->SetLineColor(kMagenta);
	PFClustersEnergy[f][eta]->SetMarkerColor(kMagenta);
	PFClustersEnergy[f][eta]->SetMarkerStyle(27);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }
      if (f ==5 || f==6){
	PFClustersEnergy[f][eta]->SetLineColor(kGreen);
	PFClustersEnergy[f][eta]->SetMarkerColor(kGreen);
	PFClustersEnergy[f][eta]->SetMarkerStyle(28);
	PFClustersEnergy[f][eta]->SetMarkerSize(0.5);
      }


      PFClustersEnergy[f][eta]->Scale(1./(TotNofEvent[f]));
      
      PFClustersEnergy[f][eta]->SetMinimum(1e-5);
      PFClustersEnergy[f][eta]->SetMaximum(10);
    }
	
  }
    
    
  //string energyname[4]={"RecHit","CaloTower","PFClusters","Gen"};
  //string fname[11]={"Data","pythia8_Monash","pythia8_TuneMBR","herwigpp","epos","qgsjetII","Genpythia8Monash","Genpythia8MBR","Genherwigpp","Genepos","GenqgsjetII"};
  for (int f=0; f<ftyp; f++){
      
    for (int e=0; e<etyp; e++){
      Energy_[f][e]->Scale(1./(TotNofEvent[f]));
     
      if ( f ==0 || f ==11 || f==12)convertbiningData(Energy_[f][e],Energy[f][e]);
      if ( f>0 && f<6)convertbiningMCDet(Energy_[f][e],Energy[f][e]);
      if ( f>12)convertbiningMCDet(Energy_[f][e],Energy[f][e]); //new MC samples
      if ( f>5 && f<11)convertbiningMCGen(Energy_[f][e],Energy[f][e]);
     
      Energy[f][e]->GetXaxis()->SetTitle("|#eta|");
      Energy[f][e]->GetYaxis()->SetTitle("dE/d|#eta| (GeV)");
      Energy[f][e]->SetMinimum(1.); // was changed from 10. to 1.
      Energy[f][e]->SetMaximum(1e9); // was changed from 10000. to 1e9
      //Data RecHit
      if (f == 0 && (e ==0) ){
	Energy[f][e]->Scale(punorm);
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(21);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      
      //Data CaloTower Magnet off
      if (f == 0 && e == 1){
	Energy[f][e]->Scale(punorm);
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineWidth(2);
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(20);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 11 && (e ==1) ){
	Energy[f][e]->Scale(punorm);

      }
      if (f == 12 && (e ==1) ){
	Energy[f][e]->Scale(punorm);

      }
 
      //Data PFClusters
      if ((f == 0 )&& e ==2){
	Energy[f][e]->Scale(punorm);
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(22);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Data PFCandidate
      if ((f == 0 )&& e ==4){
	Energy[f][e]->Scale(punorm);
	Energy[f][e]->SetLineColor(1);
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetMarkerColor(1);
	Energy[f][e]->SetMarkerStyle(23);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 Monash RecHit
      if (f == 1 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 Monash CaloTower
      if (f == 1 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 Monash PFClusters	  
      if (f == 1 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(27);
	Energy[f][e]->SetMarkerSize(2.);
      }
      //Pythia 8 Monash PFCandidate	  
      if (f == 1 && e==4 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	Energy[f][e]->SetMarkerStyle(34);
	Energy[f][e]->SetMarkerSize(2.);
      }
	  
      //Pythia 8 MBR RecHit
      if (f == 2 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 MBR CaloTower
      if (f == 2 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 MBR PFClusters
	  
      if (f == 2 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(2);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(28);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Pythia 8 MBR PFCandidate
	  
      if (f == 2 && e==4 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(2);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	Energy[f][e]->SetMarkerStyle(34);
	Energy[f][e]->SetMarkerSize(1.5);
      }	  
      //Herwigpp RecHit
      if (f == 3 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Herwigpp CaloTower
      if (f == 3 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Herwigpp PFClusters
	  
      if (f == 3 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	Energy[f][e]->SetMarkerStyle(30);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Herwigpp PFCandidate
	  
      if (f == 3 && e==4 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(1);
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	Energy[f][e]->SetMarkerStyle(34);
	Energy[f][e]->SetMarkerSize(1.5);
      }

      //Epos RecHit
      if (f == 4 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Epos CaloTower
      if (f == 4 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Epos PFClusters
	  
      if (f == 4 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	Energy[f][e]->SetMarkerStyle(32);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //Epos PFCandidate
	  
      if (f == 4 && e==4 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	Energy[f][e]->SetMarkerStyle(34);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //qgsjetII RecHit
      if (f == 5 && e == 0){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	Energy[f][e]->SetMarkerStyle(25);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //qgsjetII CaloTower
      if (f == 5 && e==1 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	Energy[f][e]->SetMarkerStyle(26);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //qgsjetII PFClusters
      if (f == 5 && e==2 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	Energy[f][e]->SetMarkerStyle(3);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      //qgsjetII PFCandidate
      if (f == 5 && e==4 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	Energy[f][e]->SetMarkerStyle(34);
	Energy[f][e]->SetMarkerSize(1.5);
      }
      // Gen-Det Plots
     
      if (f == 1 && e==5 ){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	//Energy[f][e]->SetMarkerStyle(27);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 2 && e==5 ){
	
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	//Energy[f][e]->SetMarkerStyle(28);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 3 && e==5 ){
	Energy[f][e]->SetLineStyle(3);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kOrange);
	Energy[f][e]->SetMarkerColor(kOrange);
	//Energy[f][e]->SetMarkerStyle(30);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 4 && e==5 ){
	Energy[f][e]->SetLineStyle(4);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	//Energy[f][e]->SetMarkerStyle(32);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
	   
      if (f == 5 && e==5 ){
	Energy[f][e]->SetLineStyle(5);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	//Energy[f][e]->SetMarkerStyle(3);
	//Energy[f][e]->SetMarkerSize(1.5);
      }

      
      // Gen Plots
      if (f == 6 && e==3 ){
	Energy[f][e]->SetLineStyle(1);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kTeal+3);
	Energy[f][e]->SetMarkerColor(kTeal+3);
	//Energy[f][e]->SetMarkerStyle(27);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 7 && e==3 ){
	
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	//Energy[f][e]->SetMarkerStyle(28);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 8 && e==3 ){
	Energy[f][e]->SetLineStyle(2);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kMagenta+2);
	Energy[f][e]->SetMarkerColor(kMagenta+2);
	//Energy[f][e]->SetMarkerStyle(30);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
      if (f == 9 && e==3 ){
	Energy[f][e]->SetLineStyle(4);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kBlue);
	Energy[f][e]->SetMarkerColor(kBlue);
	//Energy[f][e]->SetMarkerStyle(32);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
	   
      if (f == 10 && e==3 ){
	Energy[f][e]->SetLineStyle(5);
	Energy[f][e]->SetLineWidth(5);
	Energy[f][e]->SetLineColor(kRed);
	Energy[f][e]->SetMarkerColor(kRed);
	
	//Energy[f][e]->SetMarkerStyle(3);
	//Energy[f][e]->SetMarkerSize(1.5);
      }
	   
    }//Energy type
    // Gen
    if(f>0 && f<6)DataDetLOverMCDetL[f-1]->Divide(Energy[0][1],Energy[f][1], 1.,1.,"B");
    // Ratio to MCDetLevel-MCHadLevel
    //if(f>5)MCDetLOverMCGenL[f-6]->Divide(Energy[f-5][1],Energy[f][3], 1.,1.,"B");
    // Ratio to MCHadLevel-MCDetLevel
    if(f>5 && f<11)MCDetLOverMCGenL[f-6]->Divide(Energy[f][3],Energy[f-5][1], 1.,1.,"B");

      
      if(f>5 && f<11){
          int bining1 = Energy[f][3]->GetNbinsX();
          for (int i=1;i<=bining1;i++){
              cout<<"divide: " << (Energy[f][3]->GetBinContent(i))/(Energy[f-5][1]->GetBinContent(i)) << endl;
              
              
              // }
              //int bining2 = Energy[f-5][1]->GetNbinsX();
              //for (int i=1;i<=bining2;i++){
              // cout<<"monash det content: " << Energy[f-5][1]->GetBinContent(i) << endl;
              
              
          }
      }
   
    // Det select Gen
    // if(f>0){
    //   DataDetLOverMCDetL[f-1]->Divide(Energy[0][1],Energy[f][1], 1.,1.,"B");
    //   MCDetLOverMCGenL[f-1]->Divide(Energy[f][1],Energy[f][5], 1.,1.,"B");
    // }
  }//File Type
  for (int j=0 ; j<5 ;j++){
  DataDetLOverMCDetL[j]->GetXaxis()->SetTitle("|#eta|");
  DataDetLOverMCDetL[j]->GetYaxis()->SetTitle("#frac{dE}{d|#eta|}[DATA] / #frac{dE}{d|#eta|}[MC]");
  DataDetLOverMCDetL[j]->SetMinimum(0.1);
  DataDetLOverMCDetL[j]->SetMaximum(3.);
  }
  DataDetLOverMCDetL[0]->SetLineColor(kTeal+3);
  DataDetLOverMCDetL[0]->SetLineWidth(5);
  DataDetLOverMCDetL[0]->SetLineStyle(1);
  DataDetLOverMCDetL[0]->SetMarkerColor(kTeal+3);

  
  DataDetLOverMCDetL[1]->SetLineColor(kMagenta+2);
  DataDetLOverMCDetL[1]->SetLineWidth(5);
  DataDetLOverMCDetL[1]->SetLineStyle(2);
  DataDetLOverMCDetL[1]->SetMarkerColor(kMagenta+2);
    
  DataDetLOverMCDetL[2]->SetLineColor(kOrange);
  DataDetLOverMCDetL[2]->SetLineWidth(5);
  DataDetLOverMCDetL[2]->SetLineStyle(3);
  DataDetLOverMCDetL[2]->SetMarkerColor(kOrange);
    
  DataDetLOverMCDetL[3]->SetLineColor(kBlue);
  DataDetLOverMCDetL[3]->SetLineWidth(5);
  DataDetLOverMCDetL[3]->SetLineStyle(4);
  DataDetLOverMCDetL[3]->SetMarkerColor(kBlue);
    
  DataDetLOverMCDetL[4]->SetLineColor(kRed);
  DataDetLOverMCDetL[4]->SetLineWidth(5);
  DataDetLOverMCDetL[4]->SetLineStyle(5);
  DataDetLOverMCDetL[4]->SetMarkerColor(kRed);
  for (int j=0 ; j<5 ;j++){
  DataGenLOverMCGenL[j]->GetXaxis()->SetTitle("|#eta|");
  DataGenLOverMCGenL[j]->GetYaxis()->SetTitle("#frac{dE}{d|#eta|}[DATA] / #frac{dE}{d|#eta|}[MC]");
  DataGenLOverMCGenL[j]->SetMinimum(0.1);
  DataGenLOverMCGenL[j]->SetMaximum(3.);
  }
  DataGenLOverMCGenL[0]->SetLineColor(kTeal+3);
  DataGenLOverMCGenL[0]->SetLineWidth(5);
  DataGenLOverMCGenL[0]->SetLineStyle(1);
  DataGenLOverMCGenL[0]->SetMarkerColor(kTeal+3);
    
  DataGenLOverMCGenL[1]->SetLineColor(kMagenta+2); //from kCyan+2 to kMagenta+2
  DataGenLOverMCGenL[1]->SetLineWidth(5);
  DataGenLOverMCGenL[1]->SetLineStyle(2);
  DataGenLOverMCGenL[1]->SetMarkerColor(kMagenta+2); //from kCyan+2 to kMagenta+2
    
  DataGenLOverMCGenL[2]->SetLineColor(kMagenta+2);
  DataGenLOverMCGenL[2]->SetLineWidth(5);
  DataGenLOverMCGenL[2]->SetLineStyle(2);
  DataGenLOverMCGenL[2]->SetMarkerColor(kMagenta+2);
    
  DataGenLOverMCGenL[3]->SetLineColor(kBlue);
  DataGenLOverMCGenL[3]->SetLineWidth(5);
  DataGenLOverMCGenL[3]->SetLineStyle(4);
  DataGenLOverMCGenL[3]->SetMarkerColor(kBlue);
    
  DataGenLOverMCGenL[4]->SetLineColor(kRed);
  DataGenLOverMCGenL[4]->SetLineWidth(5);
  DataGenLOverMCGenL[4]->SetLineStyle(5);
  DataGenLOverMCGenL[4]->SetMarkerColor(kRed);


  
  
  

  for (int j=0 ; j<5 ;j++){
  MCDetLOverMCGenL[j]->GetXaxis()->SetTitle("|#eta|");
  MCDetLOverMCGenL[j]->GetYaxis()->SetTitle("#frac{dE}{d|#eta|}[had] / #frac{dE}{d|#eta|}[det]");
  MCDetLOverMCGenL[j]->SetMinimum(0.1);
  MCDetLOverMCGenL[j]->SetMaximum(3.);
  }
  MCDetLOverMCGenL[0]->SetLineColor(kTeal+3);
  MCDetLOverMCGenL[0]->SetLineWidth(5);
  MCDetLOverMCGenL[0]->SetLineStyle(1);
  MCDetLOverMCGenL[0]->SetMarkerColor(kTeal+3);
    
  MCDetLOverMCGenL[1]->SetLineColor(kMagenta+2);
  MCDetLOverMCGenL[1]->SetLineWidth(5);
  MCDetLOverMCGenL[1]->SetLineStyle(2);
  MCDetLOverMCGenL[1]->SetMarkerColor(kMagenta+2);
    
  MCDetLOverMCGenL[2]->SetLineColor(kOrange);
  MCDetLOverMCGenL[2]->SetLineWidth(5);
  MCDetLOverMCGenL[2]->SetLineStyle(3);
  MCDetLOverMCGenL[2]->SetMarkerColor(kOrange);
    
  MCDetLOverMCGenL[3]->SetLineColor(kBlue);
  MCDetLOverMCGenL[3]->SetLineWidth(5);
  MCDetLOverMCGenL[3]->SetLineStyle(4);
  MCDetLOverMCGenL[3]->SetMarkerColor(kBlue);
    
  MCDetLOverMCGenL[4]->SetLineColor(kRed);
  MCDetLOverMCGenL[4]->SetLineWidth(5);
  MCDetLOverMCGenL[4]->SetLineStyle(5);
  MCDetLOverMCGenL[4]->SetMarkerColor(kRed);



  
  //////************************************************************************************////
  /////                           NOW PLOTS                                                 ///
  ////*************************************************************************************////
  TLine *line=new TLine(3.139,1.,6.6,1.);
  TLatex *tex = new TLatex();
  tex->SetNDC();
  tex->SetTextSize(0.045);
      
  TLegend *leg[999];
  leg[0]= tdrLeg(0.3822055,0.626087,0.6015038,0.8991304);
  leg[0]->AddEntry(Energy[0][0],"RecHits","PL");
  leg[0]->AddEntry(Energy[0][1],"CaloTowers","PL");
  leg[0]->AddEntry(Energy[0][2],"PFClusters","PL");
  leg[0]->AddEntry(Energy[0][4],"PFCandidate","PL");

  leg[1]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[1]->AddEntry(Energy[1][0],"RecHits","PL");
  leg[1]->AddEntry(Energy[1][1],"CaloTowers","PL");
  leg[1]->AddEntry(Energy[1][2],"PFClusters","PL");
  leg[1]->AddEntry(Energy[1][4],"PFCandidate","PL");
  
  leg[2]= tdrLeg(0.3270677,0.6278261,0.5463659,0.9008696);
  leg[2]->AddEntry(Energy[0][1],"DATA","PL");
  leg[2]->AddEntry(Energy[1][1],"Pythia8 Monash","PL");
  leg[2]->AddEntry(Energy[2][1],"Pythia8 MBR","PL");
  //leg[2]->AddEntry(Energy[3][1],"herwigpp","PL");
  leg[2]->AddEntry(Energy[4][1],"EPOS-LHC","PL");
  leg[2]->AddEntry(Energy[5][1],"QGSJET II","PL");
  //leg[2]->AddEntry(EnCalErrBar[0],"EvSelErrBar","F");
  //leg[2]->AddEntry(EnCalErrBar[0],"EnCalErrBar","F");
  
  leg[3]= tdrLeg(0.6444724,0.6311189,0.8643216,0.9055944);
  leg[3]->AddEntry(DataDetLOverMCDetL[0],"Pythia8 Monash","L");
  leg[3]->AddEntry(DataDetLOverMCDetL[1],"Pythia8 MBR","L");
  //leg[3]->AddEntry(DataDetLOverMCDetL[2],"herwigpp","L");
  leg[3]->AddEntry(DataDetLOverMCDetL[3],"EPOS-LHC","L");
  leg[3]->AddEntry(DataDetLOverMCDetL[4],"QGSJET II","L");

  leg[4]= tdrLeg(0.35,0.62,0.57,0.90);
  leg[4]->AddEntry(Energy[1][1],"Pythia8 Monash","p");
  leg[4]->AddEntry(Energy[2][1],"Pythia8 MBR (#times 10)","p");
  //leg[4]->AddEntry(Energy[3][1],"herwigpp(#times 10^{2})","p");
  leg[4]->AddEntry(Energy[4][1],"EPOS-LHC (#times 10^{2})","p");
  leg[4]->AddEntry(Energy[5][1],"QGSJET II (#times 10^{3})","p");

  leg[5]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[5]->AddEntry(PFClustersEnergy[0][1],"DATA","PL");
  leg[5]->AddEntry(PFClustersEnergy[1][1],"Pythia8 Monash","PL");
  leg[5]->AddEntry(PFClustersEnergy[2][1],"Pythia8 MBR","PL");
  //leg[5]->AddEntry(PFClustersEnergy[3][0],"herwigpp","PL");
  leg[5]->AddEntry(PFClustersEnergy[4][1],"EPOS-LHC","PL");
  leg[5]->AddEntry(PFClustersEnergy[5][1],"QGSJET II","PL");

  leg[8]= tdrLeg(0.40,0.6195462,0.62,0.8935428);
  leg[8]->AddEntry(Energy[0][1],"Data 0T","PL");
  //leg[8]->AddEntry(EnergyFlow_Data_CaloTower3p8T,"Data 3.8T","PL");
  leg[8]->AddEntry(Energy[4][1],"pythia8 MBR 0T","PL");
  //leg[8]->AddEntry(EnergyFlow_pythia8_TuneMBR_CaloTower3p8T,"Pythia8 MBR 3.8T","PL");
  
  
  TCanvas *c[999];

 ///{"RecHit","CaloTower","PFClusters","Gen","PFCandidate","DetGenGen"}
  //------------------ DATA RecHit,CaloTower, PFClstr,PFcandidate Status--------------------//
  Energy[0][0]->GetXaxis()->SetNdivisions(511);
  Energy[0][1]->GetXaxis()->SetNdivisions(511); 
  Energy[1][0]->GetXaxis()->SetNdivisions(511);
  Energy[1][1]->GetXaxis()->SetNdivisions(511); 
  if (DataMCCompare){
    c[0] = tdrCanvas("EnergyFlow_0T_Data",Energy[0][1],4,11,kRectangular);
    c[0]->cd();
    c[0]->SetLogy();
    for (int e=0; e<3; e++){  
      Energy[0][e]->Draw("same");
      
    }
    Energy[0][4]->Draw("same");
    tex->DrawLatex(0.75,0.85,"DATA");
    tex->DrawLatex(0.75,0.15,"B = 0 T");
    leg[0]->Draw();
    sprintf(title,"Plots_SDEnhancedVeto_v2/EnergyFlow_0T_Data.pdf");
    if(save)c[0]->SaveAs(title);
    //---------------- Magnet On Compare with Magnet off
   /*
    c[16] = tdrCanvas("EnergyFlow_Data_MC_Magnet_on_off",Energy[0][1],4,11,kRectangular);
    c[16]->cd();
    c[16]->SetLogy();
    Energy[0][1]->Draw("same");// Magnet On
   // EnergyFlow_Data_CaloTower3p8T->Draw("same");//Magnet Of
    Energy[4][1]->Draw("same");// Magnet On
   // EnergyFlow_pythia8_TuneMBR_CaloTower3p8T->Draw("same");//Magnet Of
    //tex->DrawLatex(0.75,0.85,"DATA");
    leg[8]->Draw();
    tex->DrawLatex(0.75,0.15,"B = 0 T");
    sprintf(title,"Plots_SDEnhancedVeto_v2/EnergyFlow_Data_MC_Magnet_on_off.pdf");
    if(save)c[16]->SaveAs(title);

    */
    //---------- MC for Monash RecHit,CaloTower, PFClstr,PFcandidate Status -------------//
    
    
    c[1] = tdrCanvas("EnergyFlow_0T_MC",Energy[1][1],4,11,kRectangular);
    c[1]->cd();
    c[1]->SetLogy();
    for (int e=0; e<3; e++){  
      Energy[1][e]->Draw("same");
    }
    Energy[1][4]->Draw("same");
    leg[1]->Draw();
    tex->DrawLatex(0.72,0.85,"Pythia8 Monash");
    tex->DrawLatex(0.72,0.15,"B = 0 T");
    sprintf(title,"Plots_SDEnhancedVeto_v2/EnergyFlow_0T_Pythia8Monash.pdf");
    if(save)c[1]->SaveAs(title);
    //------------- DATA & MC Detlevel Compare EFlow for Calotower ------------------//
    c[2] = tdrCanvas("EnergyFlow_Data_MC_CaloTower",Energy[0][1],4,11,kRectangular);
    c[2]->cd();c[2]->SetLogy();
    //EvSelErrBar[0]->Draw("SAME2");
    //EnCalErrBar[0]->Draw("SAME2");
    for (int f=0; f<6; f++){  
      if (f!=3)Energy[f][1]->Draw("same");
    }
    
    leg[2]->Draw(); tex->DrawLatex(0.75,0.15,"B = 0 T"); tex->DrawLatex(0.75,0.20,"CaloTower");
    sprintf(title,"Plots_SDEnhancedVeto_v2/EnergyFlow_0T_Data_MC_CaloTower.pdf"); if(save)c[2]->SaveAs(title);
    
    // //------------------ Ratio of Det Level DATA and MC --------------------//
    
    DataDetLOverMCDetL[0]->GetXaxis()->SetNdivisions(511);
    DataDetLOverMCDetL[1]->GetXaxis()->SetNdivisions(511); 
    c[3] = tdrCanvas("DataDetLOverMCDetL",DataDetLOverMCDetL[1],4,11,kRectangular);
    c[3]->cd();
    for (int f=0; f<5; f++)   //for da bir döngü yok
      if (f!=2)DataDetLOverMCDetL[f]->Draw("same hist");
    tex->DrawLatex(0.75,0.20,"CaloTower");
    tex->DrawLatex(0.75,0.15,"B = 0 T");
    leg[3]->Draw(); line->Draw();
    sprintf(title,"Plots_SDEnhancedVeto_v2/DataDetLOverMCDetL.pdf");
    if(save)c[3]->SaveAs(title);
  }
  //------------------ MC Detlevel and Genlevel Compare EFlow--------------------//
  if(ratioMC){
      modelunct(MCDetLOverMCGenL[0],MCDetLOverMCGenL[1],MCDetLOverMCGenL[3],MCDetLOverMCGenL[4]);
    //Energy[6][3]->SetMinimum(1.);
    //Energy[6][3]->SetMaximum(1e9);
    Energy[6][3]->GetXaxis()->SetNdivisions(511);
    c[4] = tdrCanvas("EnergyFlow_0T_MC_Det_GenLevel_Comp",Energy[6][3],4,11,kRectangular);
    c[4]->cd();c[4]->SetLogy();

    // for (int f=1; f<6; f++){
    
    //   if (f!=3){
    //      Energy[f][1]->Scale( (pow(10,(f-1))) ) ;
    //      Energy[f][1]->Draw("same E0");
   
    //      Energy[f][5]->Scale( (pow(10,(f-1))) ) ;
    //      Energy[f][5]->Draw("same hist");
   
    //   }
    //  }
    int count =1;
    int count1 =1;
    for (int f=1; f<11; f++){
      if (f!=3 && f!=8){
      
	if(f<6){
	  Energy[f][1]->Scale( (pow(10,(count-1))) ) ;
	  Energy[f][1]->Draw("same E0");
	  //Energy[f][2]->Scale( 1./(pow(10,(f-1))) ) ;
	  count ++;
	}
	if(f>5){
	  Energy[f][3]->Scale( (pow(10,(count1-1))) ) ;
	  Energy[f][3]->Draw("same hist");
	  count1++;
	  //Energy[f][3]->Scale( 1./ (pow(10,(f-6))) ) ;
	}
     
      }
    }
    leg[4]->Draw(); tex->DrawLatex(0.66,0.87,"Lines : Gen Level");
    leg[4]->Draw(); tex->DrawLatex(0.66,0.82,"Points: CaloTower");
    tex->DrawLatex(0.75,0.15,"B = 0 T");
    sprintf(title,"Plots_SDEnhancedVeto_v2/EnergyFlow_0T_MC_Det_GenLevel_Comp.pdf"); if(save)c[4]->SaveAs(title);
    

    //-------- MC Detlevel and Genlevel Compare EFlow CorrectionFactor ---------//
      
    MCDetLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
    //MCDetLOverMCGenL[1]->GetXaxis()->SetNdivisions(511);
    c[5] = tdrCanvas("MCDetLOverMCGenL",MCDetLOverMCGenL[1],4,11,kRectangular);
    c[5]->cd();
      for (int f=0; f<5; f++){
      if (f!=2)MCDetLOverMCGenL[f]->Draw("same hist");
      }
    //tex->DrawLatex(0.75,0.20,"Detec/GenLevel");
    leg[3]->Draw(); line->Draw();tex->DrawLatex(0.75,0.15,"B = 0 T");
    sprintf(title,"Plots_SDEnhancedVeto_v2/MCDetLOverMCGenL.pdf");
    if(save)c[5]->SaveAs(title);
      
  }
  
  // //------------------ DATA MC Compare CaloTower Dist--------------------//
    
  if(energydist){
    int count=0;
    for (int k=0; k<nEtaBins; k++){
      sprintf(title,"cCaloTowerEnergy_Comp_MC_EtaBin%d",k);
      if((k>=0 && k<=3) ||  (k>=24 && k<=30))c[6+count] = tdrCanvas(title,PFClustersEnergy[0][0],4,11,kRectangular);
      if ((k>=4 && k<=5) || (k>=22 && k<=23))c[6+count] = tdrCanvas(title,PFClustersEnergy[0][4],4,11,kRectangular);
      if ((k>=6 && k<=7) || (k>=20 && k<=21))c[6+count] = tdrCanvas(title,PFClustersEnergy[0][6],4,11,kRectangular);
      if ((k>=8 && k<=9) || (k>=18 && k<=19))c[6+count] = tdrCanvas(title,PFClustersEnergy[0][8],4,11,kRectangular);
      if ((k>=10 && k<=17))c[6+count] = tdrCanvas(title,PFClustersEnergy[0][10],4,11,kRectangular);
      c[6+count]->cd();
      c[6+count]->SetLogy();
      for (int f=0; f<6; f++){
	//PFClustersEnergy[f][k]->Fit(fPFClustersEnergy[f][k],"RQN");
	if( f!=3)PFClustersEnergy[f][k]->Draw("same");
	// else PFClustersEnergy[f][k]->Draw("histo same");
	//fPFClustersEnergy[f][k]->Draw("same");
	
      }
      sprintf(title,"%3.3f<|#eta|<%3.3f",EtaBins[(7-k)-1],EtaBins[(7-k)]);
      tex->DrawLatex(0.70,0.86,title);
      leg[5]->Draw();tex->DrawLatex(0.75,0.15,"B = 0 T");
      sprintf(title,"cCaloTowerEnergy_Comp_Data_MC_EtaBin%d.pdf",k);
      if(save)c[6+count]->SaveAs(title);
      count++;
    }
  }
  if(unfold){
    unfoldbinbybin(MCDetLOverMCGenL[0],MCDetLOverMCGenL[1],MCDetLOverMCGenL[3],MCDetLOverMCGenL[4],Energy[0][1],Energy[0][1]);
    unfoldbinbybin(MCDetLOverMCGenL[0],MCDetLOverMCGenL[1],MCDetLOverMCGenL[3],MCDetLOverMCGenL[4],Energy[11][1],Energy[11][1]);
    unfoldbinbybin(MCDetLOverMCGenL[0],MCDetLOverMCGenL[1],MCDetLOverMCGenL[3],MCDetLOverMCGenL[4],Energy[12][1],Energy[12][1]);
    //for (int f=1; f<ftyp; f++){DataGenLOverMCGenL[f-1]->Divide(Energy[0][1],Energy[f][3], 1.,1.,"B");}
    for (int f=6; f<11; f++){DataGenLOverMCGenL[f-6]->Divide(Energy[0][1],Energy[f][3], 1.,1.,"B");}
    // Error bar
    DataGenLOverMCGenLCUETM1[0]->Divide(Energy[0][1],EnergyCUETM1[0], 1.,1.,"B");
    DataGenLOverMCGenLCUETS1[0]->Divide(Energy[0][1],EnergyCUETS1[0], 1.,1.,"B");
    DataGenLOverMCGenLCUETS1[1]->Divide(Energy[0][1],EnergyCUETS1[1], 1.,1.,"B");
    DataGenLOverMCGenLCUETS1[2]->Divide(Energy[0][1],EnergyCUETS1[2], 1.,1.,"B");
  
    TGraphAsymmErrors *EvSelErrBar[2];
    TGraphAsymmErrors *EnCalErrBar[2];
    TGraphAsymmErrors *Totalsysunc[2];
    TGraphAsymmErrors *MagnetOnOffUnc[1];
    TGraphAsymmErrors *gModelUnc[2];
    /// ???  how many units x has  ???///
    int bining = Energy[0][0]->GetNbinsX();
    Float_t x[bining],y[bining],yr[bining];
    Float_t xl[bining],xh[bining];
    
    Float_t ylEvSel[bining],ylEnCal[bining],ylMOnOffUnc[bining],ylTotal[bining],ylS1[bining];
    Float_t ylrEvSel[bining],ylrEnCal[bining],ylrMOnOffUnc[bining],ylrmodelUnc[bining],ylrTotal[bining],ylrS1[bining];
    
    Float_t yhEvSel[bining],yhEnCal[bining],yhMOnOffUnc[bining],yhTotal[bining],yhS1[bining];
    Float_t yhrEvSel[bining],yhrEnCal[bining],yhrMOnOffUnc[bining],yhrmodelUnc[bining],yhrTotal[bining],yhrS1[bining];
    
    // Calculation UNC unc for Energy Thresholds
    // for (int i=1;i<=bining;i++) 
    //   {
    // 	float eta= Energy[0][1]->GetBinCenter(i);
    // 	x[i-1]= eta;
    // 	xl[i-1]=(Energy[0][1]->GetBinWidth(i))/2.;
    // 	xh[i-1]=(Energy[0][1]->GetBinWidth(i))/2.;
    // 	y[i-1] =Energy[0][1]->GetBinContent(i);
    // 	yr[i-1] =1.;
	
    // 	yhEvSel[i-1]=Energy[11][1]->GetBinContent(i) - Energy[0][1]->GetBinContent(i);
    // 	ylEvSel[i-1]=Energy[0][1]->GetBinContent(i) - Energy[12][1]->GetBinContent(i);	

    // 	yhrEvSel[i-1]=(1.-(Energy[11][1]->GetBinContent(i)/ Energy[0][1]->GetBinContent(i)));
    // 	ylrEvSel[i-1]=(1.-(Energy[0][1]->GetBinContent(i) / Energy[12][1]->GetBinContent(i)));
    //   }
    // EvSelErrBar[0]= new TGraphAsymmErrors(bining,x,y,xl,xh,ylEvSel,yhEvSel);
    // EvSelErrBar[0]->SetFillColor(kOrange+1);
    // EvSelErrBar[0]->SetFillStyle(1001);
    // EvSelErrBar[0]->SetName("EnergyThreshold");
    // EvSelErrBar[1]= new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrEvSel,yhrEvSel);
    // EvSelErrBar[1]->SetFillColor(kOrange+1);
    // EvSelErrBar[1]->SetFillStyle(1001);
    // EvSelErrBar[0]->SetName("EnergyThreshold_ratio");


    // Calculation UNC unc for Energy Thresholds
    for (int i=1;i<=bining;i++) 
      {
	float eta= Energy[0][1]->GetBinCenter(i); //i think that it must be  GetBinCenter(i-1)
	x[i-1]= eta;
	xl[i-1]=(Energy[0][1]->GetBinWidth(i))/2.;
	xh[i-1]=(Energy[0][1]->GetBinWidth(i))/2.;
	y[i-1] =Energy[0][1]->GetBinContent(i);
	yr[i-1] =1.;
	
          
          //Energy[12][1] has the energy of emin = 5.5 GeV
          //Energy[14][1] MC (5.5GeV) Monash  and  Energy[13][1] Monash 5 GeV
	yhEvSel[i-1]=Energy[0][1]->GetBinContent(i) + (abs(Energy[0][1]->GetBinContent(i) - Energy[12][1]->GetBinContent(i)));
	ylEvSel[i-1]=Energy[0][1]->GetBinContent(i) - (abs(Energy[0][1]->GetBinContent(i) - Energy[12][1]->GetBinContent(i)));	

	yhrEvSel[i-1]=1-((Energy[12][1]->GetBinContent(i) / Energy[0][1]->GetBinContent(i)) / (Energy[14][1]->GetBinContent(i) / Energy[13][1]->GetBinContent(i)));// (5.5Gev / 5Gev )_data / (5.5Gev / 5Gev )_MC
	ylrEvSel[i-1]=1-((Energy[12][1]->GetBinContent(i) / Energy[0][1]->GetBinContent(i)) / (Energy[14][1]->GetBinContent(i) / Energy[13][1]->GetBinContent(i)));// (5.5Gev / 5Gev )_data / (5.5Gev / 5Gev )_MC(1.-(Energy[0][1]->GetBinContent(i) / Energy[12][1]->GetBinContent(i)));
          cout << " yhr "<<yhrEvSel[i-1]<<" ylr "<< ylrEvSel[i-1]<<endl;
      }
    EvSelErrBar[0]= new TGraphAsymmErrors(bining,x,y,xl,xh,ylEvSel,yhEvSel);
    EvSelErrBar[0]->SetFillColor(kOrange+1);
    EvSelErrBar[0]->SetFillStyle(1001);
    EvSelErrBar[0]->SetName("EnergyThreshold");
    EvSelErrBar[1]= new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrEvSel,yhrEvSel);
    EvSelErrBar[1]->SetFillColor(kOrange+1);
    EvSelErrBar[1]->SetFillStyle(1001);
    EvSelErrBar[1]->SetName("EnergyThreshold_ratio"); //  EvSelErrBar[0] -->must be   EvSelErrBar[1]
    
   
    // Calculation UNC unc for Energy Scale
    float unc  =0.10;
    for (int i=1;i<=bining;i++) 
      {
	if (i==bining)unc = 0.17;//PointOfCastor
	else unc = unc;
	yhEnCal[i-1]=(Energy[0][1]->GetBinContent(i))*unc ;
	ylEnCal[i-1]=(Energy[0][1]->GetBinContent(i))*unc ;
	yhrEnCal[i-1]= unc;
	ylrEnCal[i-1]= unc;
	
      }
    
    EnCalErrBar[0]= new TGraphAsymmErrors(bining,x,y,xl,xh,ylEnCal,yhEnCal);
    EnCalErrBar[0]->SetFillColor(kRed-6);
    EnCalErrBar[0]->SetFillStyle(1001);
    EnCalErrBar[0]->SetName("EnergyScale");
    EnCalErrBar[1]= new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrEnCal,yhrEnCal);
    EnCalErrBar[1]->SetFillColor(kRed-6);
    EnCalErrBar[1]->SetFillStyle(1001);
    EnCalErrBar[1]->SetName("EnergyScale_ratio");
    
    
    Float_t MagOnOffUnc[]={0.01,0.005,0.011,0.005,0.007,0.002,0};// From CMSNote119_v4 fig 7 Calotower
    // Calculation UNC unc for MagnetOff and ON
    for (int i=1;i<=bining;i++) 
      {
	yhrMOnOffUnc[i-1]=MagOnOffUnc[i-1];
	ylrMOnOffUnc[i-1]=MagOnOffUnc[i-1];
      }

    MagnetOnOffUnc[0]=new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrMOnOffUnc,yhrMOnOffUnc);
    MagnetOnOffUnc[0]->SetFillColor(kCyan+1);
    MagnetOnOffUnc[0]->SetFillStyle(1001);
    MagnetOnOffUnc[0]->SetName("MagnetOnOffUnc");

      /////////////////////// Model unc calculations From my side////////////////////
   /*
       float A,B,C,D,tpl,xx,yy,zz,tt;
       Float_t enb[7];
       const size_t Nunc = MCDetLOverMCGenL[0]->GetNbinsX();
      
      for (int k=0; k<1; k++) {
       for (size_t i=1; i<=Nunc; i++) {
           for (int j=0; j<7; j++) {
               A =MCDetLOverMCGenL[0]->GetBinContent(i);
               B =MCDetLOverMCGenL[1]->GetBinContent(i);
               C =MCDetLOverMCGenL[2]->GetBinContent(i);
               D =MCDetLOverMCGenL[3]->GetBinContent(i);
               tpl = (A+B+C+D)/4;
       
               xx = abs(tpl-A);
               yy = abs(tpl-B);
               zz = abs(tpl-C);
               tt = abs(tpl-D);
       
               if (xx>yy && xx>zz && xx>tt) enb[j]=xx;
               else if(yy>xx && yy>zz && yy>tt) enb[j]=yy;
               else if(zz>xx && zz>yy && zz>tt) enb[j]=zz;
               else if(tt>xx && tt>yy && tt>zz) enb[j]=tt;
               //cout << "ModelUnc" << enb[j]<<endl;
               }
           cout << "ModelUnc = "<< "Binning0-7 "<< enb[k]<<endl;
       }
      
      }

    */
      
      ///////////////////////////////////////////////////////////////////////


     //Float_t ModelUnc[]={0.027,0.017,0.016,0.016,0.016,0.016,0.021};//  old calc
     //Float_t ModelUnc[]={0.665,0.520,0.470,0.446,0.568,0.233,0.160};// my calc
      //Float_t ModelUnc[]={0.2969135,0.2576705,0.2578085,0.261371,0.37896775,0.2928915,0.169176};
      Float_t ModelUnc[]={31.5194,46.1851,67.8845,100.953,177.889,186.876,0.169722};

    for (int i=1;i<=bining;i++)
      {
	yhrmodelUnc[i-1]=ModelUnc[i-1];
	ylrmodelUnc[i-1]=ModelUnc[i-1];
      }
    
    gModelUnc[0]=new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrmodelUnc,yhrmodelUnc);
    gModelUnc[0]->SetFillColor(kCyan+1);
    gModelUnc[0]->SetFillStyle(1001);
    gModelUnc[0]->SetName("gModelUnc");

    
    // Calculation TOATAL UNC
    for (int i=1;i<=bining;i++) 
      {	
	// yhrTotal[i-1]=TMath::Sqrt(pow((yhrEvSel[i-1]),2) + pow(yhrEnCal[i-1],2) + pow(yhrMOnOffUnc[i-1],2) +pow(0.01,2) );// 1% from MC Model 
	// ylrTotal[i-1]=TMath::Sqrt(pow((ylrEvSel[i-1]),2) + pow(ylrEnCal[i-1],2) + pow(ylrMOnOffUnc[i-1],2)+pow(0.01,2));// 1% from MC Model

	yhrTotal[i-1]=TMath::Sqrt(pow((yhrEvSel[i-1]),2) + pow(yhrEnCal[i-1],2) +pow(ModelUnc[i-1],2) );// 1% from MC Model
	ylrTotal[i-1]=TMath::Sqrt(pow((ylrEvSel[i-1]),2) + pow(ylrEnCal[i-1],2) +pow(ModelUnc[i-1],2));// 1% from MC Model
	
	yhTotal[i-1]=(Energy[0][1]->GetBinContent(i))*yhrTotal[i-1];
	ylTotal[i-1]=(Energy[0][1]->GetBinContent(i))*ylrTotal[i-1];

      }

    Totalsysunc[0]=new TGraphAsymmErrors(bining,x,y,xl,xh,ylTotal,yhTotal);
    Totalsysunc[0]->SetFillColor(kGray+1);
    Totalsysunc[0]->SetFillStyle(1001);
    Totalsysunc[0]->SetName("TotalUnc");

    Totalsysunc[1]=new TGraphAsymmErrors(bining,x,yr,xl,xh,ylrTotal,yhrTotal);
    Totalsysunc[1]->SetFillColor(kGray+1);
    Totalsysunc[1]->SetFillStyle(1001);
    Totalsysunc[1]->SetName("TotalUnc_Ratio");


   //CUETS1 UNC 
    TGraphAsymmErrors * CUETS1Unc[2];
    Float_t cutes1R[bining];
    Float_t cutes1Y[bining];
    for (int i=1;i<=bining;i++) 
      {
	cutes1Y[i-1]= EnergyCUETS1[0]->GetBinContent(i);
	cutes1R[i-1]= DataGenLOverMCGenLCUETS1[0]->GetBinContent(i);
	yhS1[i-1]= Errorcuets1[i-1];
	ylS1[i-1]= Errorcuets1[i-1];
	yhrS1[i-1]= DataGenLOverMCGenLCUETS1[1]->GetBinContent(i) - DataGenLOverMCGenLCUETS1[0]->GetBinContent(i);
	ylrS1[i-1]= DataGenLOverMCGenLCUETS1[0]->GetBinContent(i) - DataGenLOverMCGenLCUETS1[2]->GetBinContent(i);
	
      }
    
    CUETS1Unc[0]= new TGraphAsymmErrors(bining,x,cutes1Y,xl,xh,ylS1,yhS1);
    CUETS1Unc[0]->SetFillColor(kRed-6);
    CUETS1Unc[0]->SetFillStyle(1001);
    CUETS1Unc[0]->SetName("CUETS1");
    CUETS1Unc[1]= new TGraphAsymmErrors(bining,x,cutes1R,xl,xh,ylrS1,yhrS1);
    CUETS1Unc[1]->SetFillColor(kRed-6);
    CUETS1Unc[1]->SetFillStyle(1001);
    CUETS1Unc[1]->SetName("CUETS1R");
    
    
    //-----------------------------------------------------------//
    Energy[0][1]->SetFillColor(kGray+1);
    Energy[0][1]->SetFillStyle(1001);
    
   //------LEGENDS--------------------//
    leg[6]= tdrLeg(0.6328321,0.1826087,0.8521303,0.4747826);
    leg[6]->AddEntry(Energy[0][1],"Data","PLF");
    leg[6]->AddEntry(Energy[6][3],"Pythia8 Monash","L");
    leg[6]->AddEntry(Energy[9][3],"EPOS-LHC","L");
    leg[6]->AddEntry(Energy[10][3],"QGSJET II","L");
      
    //-----------------------------------------//
    leg[7]= tdrLeg(0.3157895,0.6730435,0.5162907,0.9043478);
    leg[7]->AddEntry(DataGenLOverMCGenL[0],"Pythia8 Monash","L");
    leg[7]->AddEntry(DataGenLOverMCGenL[3],"EPOS-LHC","L");
    leg[7]->AddEntry(DataGenLOverMCGenL[4],"QGSJET II","L");
    leg[7]->AddEntry(Totalsysunc[1],"Total Exp. Unc.","F");
      
    //-------------------------------------------//
    leg[8]= tdrLeg(0.6328321,0.1826087,0.8521303,0.4747826);
    leg[8]->AddEntry(Energy[0][1],"Data","PLF");
    leg[8]->AddEntry(Energy[8][3],"Pythia8 CUETP8M1","L");
    //leg[8]->AddEntry(EnergyCUETM1[0],"Pythia8 CUETP8M1+MBR","L");
    //leg[8]->AddEntry(CUETS1Unc[0],"Pythia8 CUETP8S1","F");
    //------------------------------------------//
    leg[9]= tdrLeg(0.3157895,0.6730435,0.5162907,0.9043478);
    leg[9]->AddEntry(DataGenLOverMCGenL[2],"Pythia8 CUETP8M1","L");
    //leg[9]->AddEntry(DataGenLOverMCGenLCUETM1[0],"Pythia8 CUETP8M1+MBR","L");
    //leg[9]->AddEntry(CUETS1Unc[1],"Pythia8 CUETP8S1","F");
    leg[9]->AddEntry(Totalsysunc[1],"Total Exp. Unc.","F");
      
      //-----------------------------------------------------------------//
      //-------------------- Now Plots for Unfold Status ----------------//
      //-----------------------------------------------------------------//
    
      // ---------------------  Data & MCs Gen ---------------------------//
      
    c[14] = tdrCanvas("Unfolded_Data",Energy[0][1],4,11,kRectangular);
    c[14]->cd();c[14]->SetLogy();
    
    // for (int f=6; f<11; f++){  
    //   if (f!=7)Energy[f][3]->Draw("same hist");
    // }
      Totalsysunc[0]->Draw("SAME2");
      Energy[0][1]->Draw("same");      // Data Calotower
      Energy[6][3]->Draw("same hist"); //Monash Gen
      Energy[9][3]->Draw("same hist"); //Epos Gen
      Energy[10][3]->Draw("same hist");// QgsJetII Gen
    
    leg[6]->Draw(); tex->DrawLatex(0.75,0.15,"B = 0 T");tex->DrawLatex(0.23,0.15,"SD-Enhanced");
    sprintf(title,"Plots_SDEnhancedVeto_v2/Unfolded_Data.pdf"); if(save)c[14]->SaveAs(title);
    
      
      //------- Ratio of Monash Epos and QgsJet is divided the Unfold data ---//
    DataGenLOverMCGenL[1]->SetMinimum(0.1);
    DataGenLOverMCGenL[1]->SetMaximum(3.);
    DataGenLOverMCGenL[0]->GetXaxis()->SetNdivisions(511);
    c[15] = tdrCanvas("DataGenLOverMCGenL",DataGenLOverMCGenL[1],4,11,kRectangular);
    c[15]->cd();
    Totalsysunc[1]->Draw("SAME2");
     //for (int f=0; f<5; f++){
      //if (f!=1 && f!=2)DataGenLOverMCGenL[f]->Draw("same");}
    DataGenLOverMCGenL[0]->Draw("same hist"); //Ratio for Monash (Monash / UnfoldData)
    DataGenLOverMCGenL[3]->Draw("same hist"); //Ratio for Epos (Epos / UnfoldData)
    DataGenLOverMCGenL[4]->Draw("same hist"); //Ratio for qgsjetII (qgsjetII / UnfoldData)
    //if (f!=1 && f!=2)
    //DataGenLOverMCGenLCUETM1[0]->Draw("same hist");
    tex->DrawLatex(0.75,0.15,"B = 0 T");tex->DrawLatex(0.20,0.15,"SD-Enhanced");
    //tex->DrawLatex(0.75,0.20,"CaloTower");tex->DrawLatex(0.25,0.85,"B = 0 T");
    leg[7]->Draw();line->Draw();
    
    sprintf(title,"Plots_SDEnhancedVeto_v2/DataGenLOverMCGenL1.pdf");
    if(save)c[15]->SaveAs(title);
      

    
      //- Data MC 2 From Sercan and Paolo  our Cuetp8 Sercan's Cuet and UnfoldData + Unc ---//
      //------------------------------------------------------------------------------------//
    c[16] = tdrCanvas("Unfolded_Data2",Energy[0][1],4,11,kRectangular);
    c[16]->cd();c[16]->SetLogy();
      Totalsysunc[0]->Draw("SAME2");
      Energy[0][1]->Draw("same");
      Energy[8][3]->Draw("same hist");
      //EnergyCUETM1[0]->Draw("same hist");
      //CUETS1Unc[0]->Draw("SAME2");

    // for (int f=6; f<11; f++){  
    //   if (f!=7)Energy[f][3]->Draw("same hist");
    // }
    //CUETS1Unc[0]->Draw("SAME2");
    //Energy[8][3]->Draw("same hist");
    //EnergyCUETM1[0]->Draw("same hist");
    //EnergyCUETS1[0]->Draw("same hist");
    //Energy[0][1]->Draw("same");
    leg[8]->Draw(); tex->DrawLatex(0.75,0.15,"B = 0 T");tex->DrawLatex(0.20,0.15,"SD-Enhanced");
    sprintf(title,"Plots_SDEnhancedVeto_v2/Unfolded_Data2.pdf"); if(save)c[16]->SaveAs(title);
    
    
      //--------------------- Ratio 2 For Cuet's ----------------------------//
    DataGenLOverMCGenL[2]->SetMinimum(0.1);
    DataGenLOverMCGenL[2]->SetMaximum(3.);
    DataGenLOverMCGenL[2]->GetXaxis()->SetNdivisions(511); 
    c[17] = tdrCanvas("DataGenLOverMCGenL2",DataGenLOverMCGenL[2],4,11,kRectangular);
    c[17]->cd();
    Totalsysunc[1]->Draw("SAME2");
    
    // for (int f=0; f<5; f++){
    //   if (f!=1)DataGenLOverMCGenL[f]->Draw("same hist");}
    //CUETS1Unc[1]->Draw("SAME2");
    DataGenLOverMCGenL[2]->Draw("same hist");
    //DataGenLOverMCGenLCUETM1[0]->Draw("same hist");
    //DataGenLOverMCGenLCUETS1[0]->Draw("same hist");
    tex->DrawLatex(0.75,0.15,"B = 0 T");tex->DrawLatex(0.20,0.15,"SD-Enhanced");
    //tex->DrawLatex(0.75,0.20,"CaloTower");tex->DrawLatex(0.25,0.85,"B = 0 T");
    leg[9]->Draw();line->Draw();
    
    sprintf(title,"Plots_SDEnhancedVeto_v2/DataGenLOverMCGenL.pdf");
    if(save)c[17]->SaveAs(title);
    
  }
  


  
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
TH1F* convertbiningData(TH1F* h1,TH1F* h2){                                                                                                                                                      
  double content;
  double error;
  double Width;
  double CASTORCorrection =1.;
  const size_t N =h1->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=N;i++)                                                                                            
    {
      if (i==1) CASTORCorrection = 1.039;
      else CASTORCorrection = 1.;
      
      Width   = h2->GetBinWidth((N+1)-i);
      content = h1->GetBinContent(i)* CASTORCorrection;
      error   = h1->GetBinError(i);
      h2->SetBinContent(((N+1)-i),content*(1./Width));
      h2->SetBinError(((N+1)-i), error*(1./Width));

      cout <<N<<"  "<<content<<" "<<(N+1)-i<<" "<<i<<endl;
           
    }
  return h2;
  
}
TH1F* convertbiningMCGen(TH1F* h1,TH1F* h2){                                                                                                                                                      
  double content;
  double error;
  double Width;
  double CASTORCorrection =1.;
  const size_t N =h1->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=N;i++)                                                                                            
    {
      if (i>1) CASTORCorrection = 1.;
      if (i==1) CASTORCorrection = 1.;
      
      Width   = h2->GetBinWidth((N+1)-i);
      content = h1->GetBinContent(i);
      error   = h1->GetBinError(i);
      h2->SetBinContent(((N+1)-i),content*(1./Width));
      h2->SetBinError(((N+1)-i), error*(1./Width));

      cout <<N<<"  "<<content<<" "<<(N+1)-i<<" "<<i<<endl;
           
    }
  return h2;
  
}

TH1F* convertbiningMCDet(TH1F* h1,TH1F* h2){                                                                                                                                                      
  double content;
  double error;
  double Width;
  double CASTORCorrection =1.;
  const size_t N =h1->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=N;i++)                                                                                            
    {
      
      if (i==1) CASTORCorrection = 0.96;
      else CASTORCorrection = 1.;
      
      Width   = h2->GetBinWidth((N+1)-i);
      content = h1->GetBinContent(i) * CASTORCorrection;
      error   = h1->GetBinError(i);
      h2->SetBinContent(((N+1)-i),content*(1./Width));
      h2->SetBinError(((N+1)-i), error*(1./Width));

      cout <<N<<"  "<<content<<" "<<(N+1)-i<<" "<<i<<endl;
           
    }
  return h2;
  
}
TH1F* unfoldbinbybin(TH1F* h1,TH1F* h2,TH1F* h3,TH1F* h4,TH1F* h5,TH1F* h6){                                                                                                                                                      
  double content1,content2,content3,content4,content5,content6;
  double error;
  double Width;
  const size_t N =h1->GetNbinsX();                                                                                                                                          
  for (size_t i=1;i<=N;i++)                                                                                            
    {
      
      content1 = h1->GetBinContent(i);
      content2 = h2->GetBinContent(i);
      content3 = h3->GetBinContent(i);
      content4 = h4->GetBinContent(i);
      error    = h5->GetBinError(i) * ((content1 + content2 + content3 +content4) / 4.);
      content5 = h5->GetBinContent(i);
      content6 = content5 * ((content1 + content2 + content3 +content4) / 4.);
      
      h6->SetBinContent(i,content6);
      h6->SetBinError(i, error);

     
           
    }
  return h6;
  
}
 /// the part added by Zuhal  for model unc //

TH1F* modelunct(TH1F* h1,TH1F* h2,TH1F* h3,TH1F* h4){
    float A,B,C,D,tpl,xx,yy,zz,tt;
    //float A,B,C,D,tpl;
    //Float_t enb[7];
    const size_t Nunc = h1->GetNbinsX();
    ///*
    //for (int k=0; k<1; k++) {
    for (size_t i=1; i<=Nunc; i++) {
        /// for (int j=0; j<7; j++) {
        A =h1->GetBinContent(i);
        B =h2->GetBinContent(i);
        C =h3->GetBinContent(i);
        D =h4->GetBinContent(i);
        tpl = (A+B+C+D)/4;
        
        xx = abs(tpl-A);
        yy = abs(tpl-B);
        zz = abs(tpl-C);
        tt = abs(tpl-D);
        float biggest[4]={xx,yy,zz,tt};
        float temp=0.;
        for(int i=0;i<4;i++){
            if(biggest[i]>temp) temp=biggest[i];
        }
        cout << "The biggest value is: " << temp << endl;
        // if (xx>yy && xx>zz && xx>tt) enb[j]=xx;
        // else if(yy>xx && yy>zz && yy>tt) enb[j]=yy;
        // else if(zz>xx && zz>yy && zz>tt) enb[j]=zz;
        // else if(tt>xx && tt>yy && tt>zz) enb[j]=tt;
        //cout << "ModelUnc" << enb[j]<<endl;
        //}
        // cout << "ModelUnc = "<< "Binning0-7 "<< enb[k]<<endl;
        cout << "A: "<< A << " B: "<< B << " C: " << C << " D: " << D <<" tpl: "<< tpl << " xx: " << xx << " yy:" << yy << " zz: " << zz << " tt: "<< tt <<endl;
    }
    return 0;
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




/*svn co -N svn+ssh://svn.cern.ch/reps/tdr2 AN-15-119
  cd AN-15-119
  svn update utils
  svn update -N notes
  svn update notes/AN-15-119
  eval `./notes/tdr runtime -csh`
  cd notes/AN-15-119/trunk
  tdr --style=an b AN-15-119*/
