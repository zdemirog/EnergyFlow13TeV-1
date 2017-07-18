#include <TH2.h>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
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
#include "GlobalVariables.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "vector"
#include "map"
#include "TLorentzVector.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "tdrstyle_mod14.C"
#include "Run2015A_L1TechBPTXQuiet_0T.h"
using namespace std;
using namespace ROOT::Math;
int getBin(float x, float boundaries[],int b);
//____________________________________________________________________________________
//-------------------------------------------------------------
//https://root.cern.ch/root/html/tutorials/tree/hvector.C.html  
//
//-------------------------------------------------------------
//_________________________________________________________________________________________
//


//void  EnergyFlow_DetLevel_TreeProducer(int FileNumber=0)
void  EnergyFlow_DetLevel_L1TechBPTXQuiet_AllBx_E4N4_HFANDSelection_TreeProducer_Run247324()

{
    gROOT->ProcessLine("#include <vector>"); 
    //gROOT->ProcessLine(".L MC.C+");
    //setTDRStyle();
    
    gDirectory->DeleteAll();
    gROOT->ForceStyle();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0); 
    
    cout<<"file dir :   "<<readfilesdir.c_str()<<"     file number : "<<readfilesnumber_in_dir<<endl;
    static const Int_t ftyp =1;
    static const Int_t fnumber =readfilesnumber_in_dir ;

    //string fname[6]={"Data","herwigpp","pythia8_Monash","pythia8_MBR","epos","qgsjetII"};
    for (int i=0; i<ftyp;i++){
	cout << " file :"<<i<<" ==> "<<filetype.c_str()<<endl;
    };
    //gStyle->SetOptStat(111);
    bool save =true;
    string filenames ="";
    string energyname[15]={"RecHitHCAL","RecHitHF","RecHitECAL","CaloTower",
			   "hcalPFClustersCorrEnergy","ecalPFClustersCorrEnergy","HFPFClustersCorrEnergy",
			   "hcalPFClustersRawEnergy","ecalPFClustersRawEnergy","HFPFClustersRawEnergy",
			   "Gen","Events","CastorTower","CastorRecHit","PFCandidate"};
    
    
   
    
    
    static const Int_t etyp = 15 ;
   
    char title[999];
  
    vector<float> RecHitHCAL; 
    vector<float> RecHitHF;
    vector<float> RecHitECAL; 
    vector<float> RecHit;
    vector<float> CaloTower;
    vector<float> HCALTower;
    vector<float> hcalPFClustersCorrEnergy; 
    vector<float> ecalPFClustersCorrEnergy;
    vector<float> HFPFClustersCorrEnergy;
    vector<float> hcalPFClustersRawEnergy;
    vector<float> ecalPFClustersRawEnergy;
    vector<float> HFPFClustersRawEnergy;
    vector<float> PFClusters;
    vector<float> GenDetEventSelect;
    vector<float> CastorTower;
    vector<float> RecHitCastor; 
    vector<float> PFCandidate;

    TFile *fOutFile[ftyp+1];
    TTree *EnergyFlow_Det[ftyp+1];
    int TotNofEvent[ftyp+1];
    int NofEvtPassTrg[ftyp+1];
    int NofEvtPassTrg2[ftyp+1];
    
    TTree *fChain[fnumber+1];
    TFile *file[fnumber+1];

    int decade = 0;

    float Norm=1.0;
    if (filetype=="MC") Norm=1.0/1.117;
    else Norm=Norm;

    for (int f=0; f<ftyp; f++){
	//for (int f=FileNumber; f<FileNumber+1; f++){
	//----------------------Creating tree for output--------------//
	sprintf(title,"EFlow_DetLevel_%s_JINSTBinning_L1TechBPTXQuiet_AllBx_HFAnd_TiggerBit7_Emin4GeV_NoiseCut4GeV_Run247324.root",readfilesname.c_str());
	fOutFile[f]= new TFile(title,"RECREATE");
	//sprintf(title,"%s",fname.c_str());
	sprintf(title,"EFlow");
	EnergyFlow_Det[f]= new TTree(title,title);
	EnergyFlow_Det[f]->Branch("RecHitHCAL",&RecHitHCAL); 
	EnergyFlow_Det[f]->Branch("RecHitHF",&RecHitHF);
	EnergyFlow_Det[f]->Branch("RecHitECAL",&RecHitECAL); 
	EnergyFlow_Det[f]->Branch("CaloTower",&CaloTower);
	EnergyFlow_Det[f]->Branch("HCALTower",&HCALTower);
	EnergyFlow_Det[f]->Branch("RecHit",&RecHit);
	EnergyFlow_Det[f]->Branch("hcalPFClustersCorrEnergy",&hcalPFClustersCorrEnergy); 
	EnergyFlow_Det[f]->Branch("ecalPFClustersCorrEnergy",&ecalPFClustersCorrEnergy);
	EnergyFlow_Det[f]->Branch("HFPFClustersCorrEnergy",&HFPFClustersCorrEnergy);
	EnergyFlow_Det[f]->Branch("hcalPFClustersRawEnergy",&hcalPFClustersRawEnergy);
	EnergyFlow_Det[f]->Branch("ecalPFClustersRawEnergy",&ecalPFClustersRawEnergy);
	EnergyFlow_Det[f]->Branch("HFPFClustersRawEnergy",&HFPFClustersRawEnergy);
	EnergyFlow_Det[f]->Branch("PFClusters",&PFClusters);
	EnergyFlow_Det[f]->Branch("GenDetEventSelect",&GenDetEventSelect);
        EnergyFlow_Det[f]->Branch("CastorTower",&CastorTower);
        EnergyFlow_Det[f]->Branch("RecHitCastor",&RecHitCastor);
	EnergyFlow_Det[f]->Branch("PFCandidate",&PFCandidate);

	TotNofEvent[f]= 0;
    NofEvtPassTrg[f]=0;
    NofEvtPassTrg2[f]=0;
	//------- File Number ---------//                                                                                                                                                
        for (int fn=0; fn < fnumber ; fn++){
            // --- open file in EOS                                                                                                                                                      
            char Fname[999];
            sprintf(Fname,"trees_%d.root",(fn+1));
            string sFname(Fname);
            string filenames=readfilesdir+sFname;
            cout<<"file : "<<filenames.c_str()<<endl;
            file[fn]   = TFile::Open(filenames.c_str(),"READ");
	    
	    // --- read tree
	    fChain[fn] = (TTree*)file[fn]->Get("EflowTree/data");
	    //Event
	    fChain[fn]->SetBranchAddress("run", &run, &b_run);
	    fChain[fn]->SetBranchAddress("lumi", &lumi, &b_lumi);
	    fChain[fn]->SetBranchAddress("event", &event, &b_event);
	    fChain[fn]->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
	    fChain[fn]->SetBranchAddress("processID", &processID, &b_processID);
	    fChain[fn]->SetBranchAddress("cmenergy", &cmenergy, &b_cmenergy);
	    fChain[fn]->SetBranchAddress("puTrueNumInteractions", &puTrueNumInteractions, &b_puTrueNumInteractions);
	    fChain[fn]->SetBranchAddress("PUNumInteractions", &PUNumInteractions, &b_PUNumInteractions);
	    
	    std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * genParticlesp4_ = 0;
	    if(filetype =="MC"){
		fChain[fn]->SetBranchAddress("genParticlesp4", &genParticlesp4_);
		fChain[fn]->SetBranchAddress("genParticlescharge", &genParticlescharge, &b_genParticlescharge);
		fChain[fn]->SetBranchAddress("genParticlespdg", &genParticlespdg, &b_genParticlespdg);
		fChain[fn]->SetBranchAddress("genParticlesstatus", &genParticlesstatus, &b_genParticlesstatus);
		fChain[fn]->SetBranchAddress("Xix", &Xix, &b_Xix);
		fChain[fn]->SetBranchAddress("Xiy", &Xiy, &b_Xiy);
		fChain[fn]->SetBranchAddress("XiSD", &XiSD, &b_XiSD);
		fChain[fn]->SetBranchAddress("XiDD", &XiDD, &b_XiDD);
	    }
	    fChain[fn]->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
	    fChain[fn]->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
	    fChain[fn]->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
	    fChain[fn]->SetBranchAddress("vtxxErr", &vtxxErr, &b_vtxxErr);
	    fChain[fn]->SetBranchAddress("vtxyErr", &vtxyErr, &b_vtxyErr);
	    fChain[fn]->SetBranchAddress("vtxzErr", &vtxzErr, &b_vtxzErr);
	    fChain[fn]->SetBranchAddress("vtxisValid", &vtxisValid, &b_vtxisValid);
	    fChain[fn]->SetBranchAddress("vtxisFake", &vtxisFake, &b_vtxisFake);
	    fChain[fn]->SetBranchAddress("vtxchi2", &vtxchi2, &b_vtxchi2);
	    fChain[fn]->SetBranchAddress("vtxndof", &vtxndof, &b_vtxndof);
	    fChain[fn]->SetBranchAddress("vtxnTracks", &vtxnTracks, &b_vtxnTracks);
	    fChain[fn]->SetBranchAddress("simvtxx", &simvtxx, &b_simvtxx);
	    fChain[fn]->SetBranchAddress("simvtxy", &simvtxy, &b_simvtxy);
	    fChain[fn]->SetBranchAddress("simvtxz", &simvtxz, &b_simvtxz);
	    
	    //Triger
	    if (filetype=="DATA"){
	    fChain[fn]->SetBranchAddress("trgl1L1GTAlgo", &trgl1L1GTAlgo, &b_trgl1L1GTAlgo);
	    fChain[fn]->SetBranchAddress("trgl1L1GTTech", &trgl1L1GTTech, &b_trgl1L1GTTech);
	    fChain[fn]->SetBranchAddress("trgZeroBias", &trgZeroBias, &b_trgZeroBias);
	    }
	    
	    // // PF Candidate
	    std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * PFCandidatesp4_ = 0;
	    fChain[fn]->SetBranchAddress("PFCandidatesparticleId", &PFCandidatesparticleId, &b_PFCandidatesparticleId);
	    fChain[fn]->SetBranchAddress("PFCandidatesp4", &PFCandidatesp4_);
	    fChain[fn]->SetBranchAddress("PFCandidatesrawEcalEnergy", &PFCandidatesrawEcalEnergy, &b_PFCandidatesrawEcalEnergy);
	    fChain[fn]->SetBranchAddress("PFCandidatesrawHcalEnergy", &PFCandidatesrawHcalEnergy, &b_PFCandidatesrawHcalEnergy);
	    //Calo Tower
	    std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * CaloTowersp4_ = 0;
	    fChain[fn]->SetBranchAddress("CaloTowersp4", &CaloTowersp4_);
	    fChain[fn]->SetBranchAddress("CaloTowersemEnergy", &CaloTowersemEnergy, &b_CaloTowersemEnergy);
	    fChain[fn]->SetBranchAddress("CaloTowershadEnergy", &CaloTowershadEnergy, &b_CaloTowershadEnergy);
	    fChain[fn]->SetBranchAddress("CaloTowershasEB", &CaloTowershasEB, &b_CaloTowershasEB);
	    fChain[fn]->SetBranchAddress("CaloTowershasEE", &CaloTowershasEE, &b_CaloTowershasEE);
	    fChain[fn]->SetBranchAddress("CaloTowershasHB", &CaloTowershasHB, &b_CaloTowershasHB);
	    fChain[fn]->SetBranchAddress("CaloTowershasHE", &CaloTowershasHE, &b_CaloTowershasHE);
	    fChain[fn]->SetBranchAddress("CaloTowershasHF", &CaloTowershasHF, &b_CaloTowershasHF);
	    // Rec Hit
	    fChain[fn]->SetBranchAddress("EcalRecHitenergy", &EcalRecHitenergy, &b_EcalRecHitenergy);
	    fChain[fn]->SetBranchAddress("EcalRecHitEt", &EcalRecHitEt, &b_EcalRecHitEt);
	    fChain[fn]->SetBranchAddress("EcalRecHittime", &EcalRecHittime, &b_EcalRecHittime);
	    fChain[fn]->SetBranchAddress("EcalRecHitieta", &EcalRecHitieta, &b_EcalRecHitieta);
	    fChain[fn]->SetBranchAddress("EcalRecHitiphi", &EcalRecHitiphi, &b_EcalRecHitiphi);
	    fChain[fn]->SetBranchAddress("EcalRecHiteta", &EcalRecHiteta, &b_EcalRecHiteta);
	    fChain[fn]->SetBranchAddress("EcalRecHitphi", &EcalRecHitphi, &b_EcalRecHitphi);
	    
	    fChain[fn]->SetBranchAddress("HBHERecHitenergy", &HBHERecHitenergy, &b_HBHERecHitenergy);
	    fChain[fn]->SetBranchAddress("HBHERecHitEt", &HBHERecHitEt, &b_HBHERecHitEt);
	    fChain[fn]->SetBranchAddress("HBHERecHittime", &HBHERecHittime, &b_HBHERecHittime);
	    fChain[fn]->SetBranchAddress("HBHERecHitieta", &HBHERecHitieta, &b_HBHERecHitieta);
	    fChain[fn]->SetBranchAddress("HBHERecHitiphi", &HBHERecHitiphi, &b_HBHERecHitiphi);
	    fChain[fn]->SetBranchAddress("HBHERecHitdepth", &HBHERecHitdepth, &b_HBHERecHitdepth);
	    //fChain[fn]->SetBranchAddress("HBHERecHitHBHENumRecHits", &HBHERecHitHBHENumRecHits, &b_HBHERecHitHBHENumRecHits);
	    fChain[fn]->SetBranchAddress("HBHERecHiteta", &HBHERecHiteta, &b_HBHERecHiteta);
	    fChain[fn]->SetBranchAddress("HBHERecHitphi", &HBHERecHitphi, &b_HBHERecHitphi);
	    fChain[fn]->SetBranchAddress("HFRecHitenergy", &HFRecHitenergy, &b_HFRecHitenergy);
	    fChain[fn]->SetBranchAddress("HFRecHitEt", &HFRecHitEt, &b_HFRecHitEt);
	    fChain[fn]->SetBranchAddress("HFRecHittime", &HFRecHittime, &b_HFRecHittime);
	    fChain[fn]->SetBranchAddress("HFRecHitieta", &HFRecHitieta, &b_HFRecHitieta);
	    fChain[fn]->SetBranchAddress("HFRecHitiphi", &HFRecHitiphi, &b_HFRecHitiphi);
	    fChain[fn]->SetBranchAddress("HFRecHitdepth", &HFRecHitdepth, &b_HFRecHitdepth);
	    //fChain[fn]->SetBranchAddress("HFRecHitHFNumRecHits", &HFRecHitHFNumRecHits, &b_HFRecHitHFNumRecHits);
	    fChain[fn]->SetBranchAddress("HFRecHiteta", &HFRecHiteta, &b_HFRecHiteta);
	    fChain[fn]->SetBranchAddress("HFRecHitphi", &HFRecHitphi, &b_HFRecHitphi);
	    
	    
	    fChain[fn]->SetBranchAddress("CastorRecHitEnergy", &CastorRecHitEnergy, &b_CastorRecHitEnergy);
	    fChain[fn]->SetBranchAddress("CastorRecHitSector", &CastorRecHitSector, &b_CastorRecHitSector);
	    fChain[fn]->SetBranchAddress("CastorRecHitModule", &CastorRecHitModule, &b_CastorRecHitModule);
	    fChain[fn]->SetBranchAddress("CastorRecHitisBad", &CastorRecHitisBad, &b_CastorRecHitisBad);
	    fChain[fn]->SetBranchAddress("CastorRecHitisSaturated", &CastorRecHitisSaturated, &b_CastorRecHitisSaturated);
	    fChain[fn]->SetBranchAddress("CastorRecHitisDesaturated", &CastorRecHitisDesaturated, &b_CastorRecHitisDesaturated);
	    
	    std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * CastorTowerp4_ = 0;
	    fChain[fn]->SetBranchAddress("CastorTowerp4", &CastorTowerp4_);
	    fChain[fn]->SetBranchAddress("CastorToweremEnergy", &CastorToweremEnergy, &b_CastorToweremEnergy);
	    fChain[fn]->SetBranchAddress("CastorTowerhadEnergy", &CastorTowerhadEnergy, &b_CastorTowerhadEnergy);
	    fChain[fn]->SetBranchAddress("CastorTowerNrechits", &CastorTowerNrechits, &b_CastorTowerNrechits);
	    
	    //PF Clusters
	    fChain[fn]->SetBranchAddress("hcalPFClustersenergy", &hcalPFClustersenergy, &b_hcalPFClustersenergy);
	    fChain[fn]->SetBranchAddress("hcalPFClusterscorrectedEnergy", &hcalPFClusterscorrectedEnergy, &b_hcalPFClusterscorrectedEnergy);
	    fChain[fn]->SetBranchAddress("hcalPFClusterscorrectedEnergyUncertainty", &hcalPFClusterscorrectedEnergyUncertainty, &b_hcalPFClusterscorrectedEnergyUncertainty);
	    fChain[fn]->SetBranchAddress("hcalPFClusterstime", &hcalPFClusterstime, &b_hcalPFClusterstime);
	    fChain[fn]->SetBranchAddress("hcalPFClustersdepth", &hcalPFClustersdepth, &b_hcalPFClustersdepth);
	    fChain[fn]->SetBranchAddress("hcalPFClusterspt", &hcalPFClusterspt, &b_hcalPFClusterspt);
	    fChain[fn]->SetBranchAddress("hcalPFClustersEt", &hcalPFClustersEt, &b_hcalPFClustersEt);
	    fChain[fn]->SetBranchAddress("hcalPFClusterseta", &hcalPFClusterseta, &b_hcalPFClusterseta);
	    fChain[fn]->SetBranchAddress("hcalPFClustersphi", &hcalPFClustersphi, &b_hcalPFClustersphi);
	    fChain[fn]->SetBranchAddress("hcalPFClusterssize", &hcalPFClusterssize, &b_hcalPFClusterssize);
	    fChain[fn]->SetBranchAddress("hcalPFClustersisInClean", &hcalPFClustersisInClean, &b_hcalPFClustersisInClean);
	    fChain[fn]->SetBranchAddress("hcalPFClustersisInUnClean", &hcalPFClustersisInUnClean, &b_hcalPFClustersisInUnClean);
	    fChain[fn]->SetBranchAddress("hfPFClustersenergy", &hfPFClustersenergy, &b_hfPFClustersenergy);
	    fChain[fn]->SetBranchAddress("hfPFClusterscorrectedEnergy", &hfPFClusterscorrectedEnergy, &b_hfPFClusterscorrectedEnergy);
	    fChain[fn]->SetBranchAddress("hfPFClusterscorrectedEnergyUncertainty", &hfPFClusterscorrectedEnergyUncertainty, &b_hfPFClusterscorrectedEnergyUncertainty);
	    fChain[fn]->SetBranchAddress("hfPFClusterstime", &hfPFClusterstime, &b_hfPFClusterstime);
	    fChain[fn]->SetBranchAddress("hfPFClustersdepth", &hfPFClustersdepth, &b_hfPFClustersdepth);
	    fChain[fn]->SetBranchAddress("hfPFClusterspt", &hfPFClusterspt, &b_hfPFClusterspt);
	    fChain[fn]->SetBranchAddress("hfPFClustersEt", &hfPFClustersEt, &b_hfPFClustersEt);
	    fChain[fn]->SetBranchAddress("hfPFClusterseta", &hfPFClusterseta, &b_hfPFClusterseta);
	    fChain[fn]->SetBranchAddress("hfPFClustersphi", &hfPFClustersphi, &b_hfPFClustersphi);
	    fChain[fn]->SetBranchAddress("hfPFClusterssize", &hfPFClusterssize, &b_hfPFClusterssize);
	    fChain[fn]->SetBranchAddress("hfPFClustersisInClean", &hfPFClustersisInClean, &b_hfPFClustersisInClean);
	    fChain[fn]->SetBranchAddress("hfPFClustersisInUnClean", &hfPFClustersisInUnClean, &b_hfPFClustersisInUnClean);
	    fChain[fn]->SetBranchAddress("ecalPFClustersenergy", &ecalPFClustersenergy, &b_ecalPFClustersenergy);
	    fChain[fn]->SetBranchAddress("ecalPFClusterscorrectedEnergy", &ecalPFClusterscorrectedEnergy, &b_ecalPFClusterscorrectedEnergy);
	    fChain[fn]->SetBranchAddress("ecalPFClusterscorrectedEnergyUncertainty", &ecalPFClusterscorrectedEnergyUncertainty, &b_ecalPFClusterscorrectedEnergyUncertainty);
	    fChain[fn]->SetBranchAddress("ecalPFClusterstime", &ecalPFClusterstime, &b_ecalPFClusterstime);
	    fChain[fn]->SetBranchAddress("ecalPFClustersdepth", &ecalPFClustersdepth, &b_ecalPFClustersdepth);
	    fChain[fn]->SetBranchAddress("ecalPFClusterspt", &ecalPFClusterspt, &b_ecalPFClusterspt);
	    fChain[fn]->SetBranchAddress("ecalPFClustersEt", &ecalPFClustersEt, &b_ecalPFClustersEt);
	    fChain[fn]->SetBranchAddress("ecalPFClusterseta", &ecalPFClusterseta, &b_ecalPFClusterseta);
	    fChain[fn]->SetBranchAddress("ecalPFClustersphi", &ecalPFClustersphi, &b_ecalPFClustersphi);
	    fChain[fn]->SetBranchAddress("ecalPFClusterssize", &ecalPFClusterssize, &b_ecalPFClusterssize);
	    fChain[fn]->SetBranchAddress("ecalPFClustersisInClean", &ecalPFClustersisInClean, &b_ecalPFClustersisInClean);
	    fChain[fn]->SetBranchAddress("ecalPFClustersisInUnClean", &ecalPFClustersisInUnClean, &b_ecalPFClustersisInUnClean);
	    fChain[fn]->SetBranchAddress("bx", &bx, &b_bx);
	    
            // ----------------------- Cut-------------------------//
            lumimin =97;
            lumimax =311;
            emin = 4.;
            NoiseCut = 4.;
            Etabnd = 5.205;//Max Eta HF for CMS
            cmseta = 6.6;//Max Eta for CMS
            etamin =3.152;
            //etamax =4.889;//event select
            etamax =5.205;//event select
            etamaxMinus =-3.152;
            etaminMinus =-5.205;//event select
            minXiSD=1e-2;
	    
	    
	    // event 
	    Long64_t nentries = fChain[fn]->GetEntriesFast();
	    cout<<"Entries "<<nentries<<endl;
	    Long64_t nbytes = 0, nb = 0;
	    int maxevent=1000;
	    for (Long64_t ev=0; ev<nentries;ev++) {
		
		//for (Long64_t ev=0; ev<maxevent;ev++) {
		
		Long64_t iev = fChain[fn]->LoadTree(ev);
		if (iev < 0)break;
		nb = fChain[fn]->GetEntry(ev);   nbytes += nb;
		
		double progress = 10.0*ev/(1.0*nentries);
		int k = TMath::FloorNint(progress); 
		if (k > decade) 
		    cout<<10*k<<" %"<<endl;
		decade = k; 
		RecHit.clear(); 
		RecHitHCAL.clear(); 
		RecHitHF.clear();
		RecHitECAL.clear(); 
		CaloTower.clear();
		HCALTower.clear();
		hcalPFClustersCorrEnergy.clear(); 
		ecalPFClustersCorrEnergy.clear();
		HFPFClustersCorrEnergy.clear();
		hcalPFClustersRawEnergy.clear();
		ecalPFClustersRawEnergy.clear();
		HFPFClustersRawEnergy.clear();
		PFClusters.clear();
		GenDetEventSelect.clear();
                CastorTower.clear();
		RecHitCastor.clear();
		PFCandidate.clear();
                
		
		for (int i =0; i<nEtaBins;i++) {
		    RecHitHCALEtaSums[i]=0.;
		    
		    RecHitECALEtaSums[i]=0.;
		    RecHitHFEtaSums[i]=0.;
		    HCALTowerEtaSums[i]=0.;
		    CaloEtaSums[i]=0.;
		    PFCandEtaSums[i]=0.;
		    
		    PFClustersEtaSums[i]=0.;
		    PFClustersECALEtaSums[i]=0.;
		    PFClustersHCALEtaSums[i]=0.;
		    PFClustersHFEtaSums[i]=0.;
		    
		    PFClustersEtaSumsUp[i]=0.;
		    PFClustersECALEtaSumsUp[i]=0.;
		    PFClustersHCALEtaSumsUp[i]=0.;
		    PFClustersHFEtaSumsUp[i]=0.;
		    
		    PFClustersEtaSumsLo[i]=0.;
		    PFClustersECALEtaSumsLo[i]=0.;
		    PFClustersHCALEtaSumsLo[i]=0.;
		    PFClustersHFEtaSumsLo[i]=0.;
		    
		    PFClustersEtaSumsRaw[i]=0.;
		    PFClustersECALEtaSumsRaw[i]=0.;
		    PFClustersHCALEtaSumsRaw[i]=0.;
		    PFClustersHFEtaSumsRaw[i]=0.;
		    
		    GenEtaSumsDetEvntSelct[i]=0.;
                    CastorTowerEtaSums[i]=0.;
		    RecHitCastorEtaSums[i]=0.;
		    
		}
		
		// Ini. for each event
		EnergyCutRecHit=false;
		EnergyCutRecHitMinus=false;
		EnergyCutTowerMinus=false;
		EnergyCutTowerPlus=false;
		EnergyCutTower=false;
		EnergyCutRecHitPlus=false;
		EnergyCutPFCluster=false;
		EnergyCutRecHitHFAnd=false;
		RunCut=false;
		LumiCut=false;
		bxCut=false;
			    
		if (filetype =="DATA" && run==247324) RunCut=true;
		if(filetype =="MC") RunCut =true;// for MC
		
		// ----------------------- fat bunch ------------------//
		if (filetype =="DATA" && bx==208) bxCut=true;
		if(filetype =="MC") bxCut =true;// for MC
		
		// ----------------------- lumi cut ------------------//
		
		if (filetype =="DATA")LumiCut=(lumi>=lumimin &&lumi<lumimax);
		else LumiCut=true;// for MC
		
		//if(LumiCut && RunCut ){
		if(RunCut && LumiCut){
		    
		    // ----------------------- triger cut --------------//
		    if (filetype =="DATA")TrigerPass =  (trgl1L1GTTech->at(7)== 1);
		    else TrigerPass=true;
		    
		    if(TrigerPass){
			//-- Event Selection for Det Level --//
			//--------- HF OR-------//
/*			for(unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
			    XYZTVector caltwr = (*CaloTowersp4_)[cal];
			    
			    if((abs(caltwr.Eta())>=etamin && abs(caltwr.Eta())<=etamax) && ((caltwr.E())*Norm) > emin) {
				EnergyCutTower=true;
			    }
			}
 */
                //------------ HF_AND ------------//
                //----------HF Minus-------//
                for (unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
                    
                    XYZTVector caltwr = (*CaloTowersp4_)[cal];
                    
                    if((caltwr.Eta()>=etaminMinus && caltwr.Eta()<=etamaxMinus) && ((caltwr.E())*Norm) > emin) {
                        
                        EnergyCutTowerMinus=true;
                    }
                }
                
                
                //----------HF Plus-------//
                for (unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
                    
                    XYZTVector caltwr = (*CaloTowersp4_)[cal];
                    
                    if((caltwr.Eta()>=etamin && caltwr.Eta()<=etamax) && ((caltwr.E())*Norm) > emin) {
                        
                        EnergyCutTowerPlus=true;
                    }
                }
                
                
			if (EnergyCutTowerMinus && EnergyCutTowerPlus){
			    // Gen Level Event selection for Detlevel
			    if (filetype =="MC"){
				for (unsigned long gd=0; gd<genParticlesp4_->size(); gd++){
				    XYZTVector genDet = (*genParticlesp4_)[gd];
				    
				    if((abs(genDet.Eta())<= cmseta)  ) {
					selcetetabin=getBin(genDet.Eta(),EtaBins, nEtaBins);
					GenEtaSumsDetEvntSelct[selcetetabin] = GenEtaSumsDetEvntSelct[selcetetabin] + genDet.E();
				    }
				}//Gen
			    }
			    // ----------------------- ECALRecHit --------------- //
			    for(unsigned long j=0; j<EcalRecHitenergy->size(); j++) {
				selcetetabin=getBin(EcalRecHiteta->at(j),EtaBins,nEtaBins);
				RecHitECALEtaSums[selcetetabin]= RecHitECALEtaSums[selcetetabin] + EcalRecHitenergy->at(j);
			    }  
			    
			    // ----------------------- HCALRecHit --------------- //
			    for(unsigned long j=0; j<HBHERecHitenergy->size(); j++) {
				selcetetabin=getBin(HBHERecHiteta->at(j),EtaBins,nEtaBins);
				RecHitHCALEtaSums[selcetetabin]= RecHitHCALEtaSums[selcetetabin] + HBHERecHitenergy->at(j);
			    }
			    
			    // ----------------------- HFRecHit --------------- //
			    for(unsigned long j=0; j<HFRecHitenergy->size(); j++) {
				if((abs(HFRecHiteta->at(j))>=etamin && abs(HFRecHiteta->at(j))<=Etabnd)) {
				    selcetetabin=getBin(HFRecHiteta->at(j),EtaBins,nEtaBins);
				    RecHitHFEtaSums[selcetetabin]= RecHitHFEtaSums[selcetetabin] + ((HFRecHitenergy->at(j))*Norm);
				}
			    }  
                //----------------------- CaloTower  -------------//
                for (unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
                    
                    XYZTVector caltwr = (*CaloTowersp4_)[cal];
                    selcetetabin=getBin(caltwr.Eta(),EtaBins,nEtaBins);
                    if(abs(caltwr.Eta()) >=etamin && abs(caltwr.Eta()) <=Etabnd && caltwr.E()*Norm > NoiseCut) HCALTowerEtaSums[selcetetabin]= HCALTowerEtaSums[selcetetabin] + (caltwr.E())*Norm;//HF MC Norm
                    else HCALTowerEtaSums[selcetetabin]= HCALTowerEtaSums[selcetetabin] + (caltwr.E());
                }
                
			    // // -----------------------  PF cand  -----------------------//
			    // //PFCand
			    for (unsigned long pf=0; pf<PFCandidatesp4_->size(); pf++){
			    XYZTVector pfcand = (*PFCandidatesp4_)[pf];
			    //cout << " pf energy "<<pfcand.E()<< " pf eta "<< pfcand.Eta()<<endl;
			    
			    // if(pfcand.E() > 10.) {
			    //cout <<"event "<< ev <<" PFCandidates "<<PFCandidatesp4_->size()<<endl; 
			    selcetetabin=getBin(pfcand.Eta(),EtaBins, nEtaBins);
			    
			    PFCandEtaSums[selcetetabin]= PFCandEtaSums[selcetetabin] + pfcand.E();
			    
			    //cout<<" eta "<<pfcand.Eta()<<" bin "<<selcetetabin<<" energy " << pfcand.E() <<" total PF energy "<<PFCandEtaSums[selcetetabin]<<endl;
			    // }
			    }//PFCand
			    
			    //----------------------- PF Cluster --------------- //
			    //HCALcorrected energy
			    for(unsigned long pfc=0; pfc<hcalPFClusterscorrectedEnergy->size(); pfc++) {
				selcetetabin=getBin(hcalPFClusterseta->at(pfc),EtaBins, nEtaBins);
				PFClustersHCALEtaSums[selcetetabin]= PFClustersHCALEtaSums[selcetetabin]+hcalPFClusterscorrectedEnergy->at(pfc);
			    }
			    //Ecal corrected energy
			    for(unsigned long pfc=0; pfc<ecalPFClusterscorrectedEnergy->size(); pfc++) {
				selcetetabin=getBin(ecalPFClusterseta->at(pfc),EtaBins, nEtaBins);
				PFClustersECALEtaSums[selcetetabin]= PFClustersECALEtaSums[selcetetabin]+ ecalPFClusterscorrectedEnergy->at(pfc);
			    }
			    //Hf corrected energy	   
			    for(unsigned long pfc=0; pfc<hfPFClusterscorrectedEnergy->size(); pfc++) {
				if((abs(hfPFClusterseta->at(pfc))>=etamin && abs(hfPFClusterseta->at(pfc))<=Etabnd)) {
				    selcetetabin=getBin(hfPFClusterseta->at(pfc),EtaBins, nEtaBins);
				    PFClustersHFEtaSums[selcetetabin]= PFClustersHFEtaSums[selcetetabin]+((hfPFClusterscorrectedEnergy->at(pfc))*Norm);
				}
			    }
			    
			    //HCAL raw energy
			    for(unsigned long pfc=0; pfc<hcalPFClustersenergy->size(); pfc++) {
				selcetetabin=getBin(hcalPFClusterseta->at(pfc),EtaBins, nEtaBins);
				PFClustersHCALEtaSumsRaw[selcetetabin]= PFClustersHCALEtaSumsRaw[selcetetabin]+hcalPFClustersenergy->at(pfc);
				
			    }
			    //ECAL raw energy
			    for(unsigned long pfc=0; pfc< ecalPFClustersenergy->size(); pfc++) {
				selcetetabin=getBin(ecalPFClusterseta->at(pfc),EtaBins, nEtaBins);
				PFClustersECALEtaSumsRaw[selcetetabin]= PFClustersECALEtaSumsRaw[selcetetabin]+ecalPFClustersenergy->at(pfc);
			    }
			    //HF raw energy
			    for(unsigned long pfc=0; pfc<hfPFClustersenergy->size(); pfc++) {
				
				if((abs(hfPFClusterseta->at(pfc))>=etamin && abs(hfPFClusterseta->at(pfc))<=Etabnd) ) { 
				    selcetetabin=getBin(hfPFClusterseta->at(pfc),EtaBins, nEtaBins);
				    PFClustersHFEtaSumsRaw[selcetetabin]= PFClustersHFEtaSumsRaw[selcetetabin]+((hfPFClustersenergy->at(pfc))*Norm);
				}
			    }
			    // --------------------Castor Tower-----------------//
			    for (unsigned long cas=0; cas<CastorTowerp4_->size(); cas++){
				
				XYZTVector castwr = (*CastorTowerp4_)[cas];
				selcetetabin=getBin(castwr.Eta(),EtaBins,nEtaBins);
				CastorTowerEtaSums[selcetetabin]= CastorTowerEtaSums[selcetetabin] + castwr.E();
			    }
                
			    // ----------------------- Castor RecHit --------------- //                                                                                                               
                            for(unsigned long j=0; j<CastorRecHitEnergy->size(); j++) {
				//selcetetabin=getBin(HFRecHiteta->at(j),EtaBins,nEtaBins);
				RecHitCastorEtaSums[0]= RecHitCastorEtaSums[0] + (CastorRecHitEnergy->at(j));
                               
                            }
                
			    
			    // Filling tree
			    float totenergycalotower,totenergyrechit,totenergypfclusters,totenergypfcandidate;
			    for (int k=0;k<nEtaBins;k++){
				// (HFPFClustersRawEnergy[k]).push_back(PFClustersHFEtaSumsRaw[k]);
				// (Gen_DetEvntSelct[k]).push_back(GenEtaSumsDetEvntSelct[k]);
				totenergyrechit =0.;
				totenergypfclusters=0.;
				totenergycalotower=0.;
				totenergypfcandidate=0.;

				RecHitHCAL.push_back(RecHitHCALEtaSums[k]);
				RecHitHF.push_back(RecHitHFEtaSums[k]);
				RecHitECAL.push_back(RecHitECALEtaSums[k]);
				RecHitCastor.push_back(RecHitCastorEtaSums[0]);
				if (k==0)totenergyrechit = RecHitHCALEtaSums[k]+RecHitECALEtaSums[k]+RecHitHFEtaSums[k]+RecHitCastorEtaSums[0];
				else totenergyrechit = RecHitHCALEtaSums[k]+RecHitECALEtaSums[k]+RecHitHFEtaSums[k];
				RecHit.push_back(totenergyrechit);

				HCALTower.push_back(HCALTowerEtaSums[k]);
			        CastorTower.push_back(CastorTowerEtaSums[k]);
				totenergycalotower = HCALTowerEtaSums[k]+CastorTowerEtaSums[k];
				CaloTower.push_back(totenergycalotower);

				hcalPFClustersCorrEnergy.push_back(PFClustersHCALEtaSums[k]); 
				ecalPFClustersCorrEnergy.push_back(PFClustersECALEtaSums[k]);
				HFPFClustersCorrEnergy.push_back(PFClustersHFEtaSums[k]);
				hcalPFClustersRawEnergy.push_back(PFClustersHCALEtaSumsRaw[k]);
				ecalPFClustersRawEnergy.push_back(PFClustersECALEtaSumsRaw[k]);
				HFPFClustersRawEnergy.push_back(PFClustersHFEtaSumsRaw[k]);
				//totenergypfclusters = PFClustersECALEtaSums[k]+PFClustersHCALEtaSumsRaw[k]+PFClustersHFEtaSumsRaw[k]+CastorTowerEtaSums[k];
				totenergypfclusters = PFClustersECALEtaSums[k]+PFClustersHCALEtaSumsRaw[k]+PFClustersHFEtaSumsRaw[k];
				PFClusters.push_back(totenergypfclusters);

				GenDetEventSelect.push_back(GenEtaSumsDetEvntSelct[k]);

				PFCandidate.push_back(PFCandEtaSums[k]);
			    }
			    
			    TotNofEvent[f]++;
			    EnergyFlow_Det[f]->Fill(); 
			}//HF OR Energy Cut
                NofEvtPassTrg[f]++;
		    }//Trigerpass
                NofEvtPassTrg2[f]++;
		}//Lumi && Run && bnx Cut	
	    }//Event
	    cout <<"  Total event : "<<TotNofEvent[f]<<endl;
        cout <<"  Trg Passed Evt : "<<NofEvtPassTrg[f]<<endl;
        cout <<"  Trg Passed Evt2 : "<<NofEvtPassTrg2[f]<<endl;
	}//file number
	fOutFile[f]->Write();
	delete fOutFile[f];
	cout<<"Filling tree and write in root file with "<<TotNofEvent[f]<<" events"<<" for file : "<<readfilesname.c_str()<<endl;
    }//File
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
