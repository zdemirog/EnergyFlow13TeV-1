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


void  EnergyFlow_DataZeroBias_TreeProducer_withCastor_ZuhalHFeventSelection.C(int FileNumber=0)

{
    
    //gROOT->ProcessLine(".L tdrstyle_mod14.C");
    //setTDRStyle();
    
    gDirectory->DeleteAll();
    gROOT->ForceStyle();
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0); 
    gROOT->ProcessLine("#include <vector>"); 

    static const Int_t ftyp =6 ;// Data and MC
    string fname[6]={"DataZeroBiasAllBx","herwigpp","pythia8_Monash","pythia8_MBR","epos","qgsjetII"};
    for (int i=0; i<ftyp;i++){
	cout << " file :"<<i<<" ==> "<<fname[i].c_str()<<endl;
    };
    //gStyle->SetOptStat(111);
    bool save =true;
    
    string filenames[6] ={"root://eoscms.cern.ch//eos/cms/store/user/higgs313/CFFtrees/ZeroBias_byHans/data_RunIILowPU_0T_01072015_data_ZeroBias1_Run2015A_Run247324_byHans.root",
			  "root://eoscms.cern.ch//eos/cms/store/user/zdemirog/13TeV/MC/MinBias_TuneEE5C_13TeV-herwigpp_MagnetOff_trees.root",
			  "root://eoscms.cern.ch//eos/cms/store/user/zdemirog/13TeV/MC/TuneMonash_MagnetOff_TuneMonashMagnetOff.root",
			  "root://eoscms.cern.ch//eos/cms/store/user/zdemirog/13TeV/MC/MinBias_TuneMBR_13TeV-pythia8_MagnetOff_trees.root",
			  "root://eoscms.cern.ch//eos/cms/store/user/zdemirog/13TeV/MC/ReggeGribovPartonMC_13TeV-EPOS_MagnetOff_trees.root",
			  "root://eoscms.cern.ch//eos/cms/store/user/zdemirog/13TeV/MC/ReggeGribovPartonMC_13TeV-QGSJetII_MagnetOff_trees.root"
			  
			  
    };
    
    
    string energyname[13]={"RecHitHCAL","RecHitHF","RecHitECAL","CaloTower",
			   "hcalPFClustersCorrEnergy","ecalPFClustersCorrEnergy","HFPFClustersCorrEnergy",
			   "hcalPFClustersRawEnergy","ecalPFClustersRawEnergy","HFPFClustersRawEnergy",
			   "Gen","Events","CastorTower"};
    
    
   
    
    
    static const Int_t etyp = 13 ;
   
    char title[999];
  
    vector<float> RecHitHCAL; 
    vector<float> RecHitHF;
    vector<float> RecHitECAL; 
    vector<float> RecHit;
    vector<float> CaloTower;
    vector<float> hcalPFClustersCorrEnergy; 
    vector<float> ecalPFClustersCorrEnergy;
    vector<float> HFPFClustersCorrEnergy;
    vector<float> hcalPFClustersRawEnergy;
    vector<float> ecalPFClustersRawEnergy;
    vector<float> HFPFClustersRawEnergy;
    vector<float> PFClusters;
    vector<float> GenDetEventSelect;
    vector<float> CastorTower;

    TFile *file[ftyp+1];
    TFile *fOutFile[ftyp+1];
    TTree *fChain[ftyp+1];

    TTree *EnergyFlow_Det[ftyp+1];
  
    int TotNofEvent[ftyp+1];
  


    int decade = 0;
    //for (int f=0; f<ftyp; f++){
    for (int f=FileNumber; f<FileNumber+1; f++){
	//----------------------Creating tree for output--------------//
	sprintf(title,"EFlow_DetLevel_%s_tree_0Tesla_withCastor.root",fname[f].c_str());
	fOutFile[f]= new TFile(title,"RECREATE");
	//sprintf(title,"%s",fname[f].c_str());
	sprintf(title,"EFlow");
	EnergyFlow_Det[f]= new TTree(title,title);
	EnergyFlow_Det[f]->Branch("RecHitHCAL",&RecHitHCAL); 
	EnergyFlow_Det[f]->Branch("RecHitHF",&RecHitHF);
	EnergyFlow_Det[f]->Branch("RecHitECAL",&RecHitECAL); 
	EnergyFlow_Det[f]->Branch("CaloTower",&CaloTower);
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
	// --- open file in EOS
	file[f]   = TFile::Open(filenames[f].c_str(),"READ");
      
	cout<<"file : "<<filenames[f].c_str()<<endl;
	// --- read tree
	fChain[f] = (TTree*)file[f]->Get("EflowTree/data");
	//Event
	fChain[f]->SetBranchAddress("run", &run, &b_run);
	fChain[f]->SetBranchAddress("lumi", &lumi, &b_lumi);
	fChain[f]->SetBranchAddress("event", &event, &b_event);
	fChain[f]->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
	fChain[f]->SetBranchAddress("processID", &processID, &b_processID);
	fChain[f]->SetBranchAddress("cmenergy", &cmenergy, &b_cmenergy);
	fChain[f]->SetBranchAddress("puTrueNumInteractions", &puTrueNumInteractions, &b_puTrueNumInteractions);
	fChain[f]->SetBranchAddress("PUNumInteractions", &PUNumInteractions, &b_PUNumInteractions);
	
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * genParticlesp4_ = 0;
	if(f>0){
	fChain[f]->SetBranchAddress("genParticlesp4", &genParticlesp4_);
	fChain[f]->SetBranchAddress("genParticlescharge", &genParticlescharge, &b_genParticlescharge);
	fChain[f]->SetBranchAddress("genParticlespdg", &genParticlespdg, &b_genParticlespdg);
	fChain[f]->SetBranchAddress("genParticlesstatus", &genParticlesstatus, &b_genParticlesstatus);
	fChain[f]->SetBranchAddress("Xix", &Xix, &b_Xix);
	fChain[f]->SetBranchAddress("Xiy", &Xiy, &b_Xiy);
	fChain[f]->SetBranchAddress("XiSD", &XiSD, &b_XiSD);
	fChain[f]->SetBranchAddress("XiDD", &XiDD, &b_XiDD);
	}
	fChain[f]->SetBranchAddress("vtxx", &vtxx, &b_vtxx);
	fChain[f]->SetBranchAddress("vtxy", &vtxy, &b_vtxy);
	fChain[f]->SetBranchAddress("vtxz", &vtxz, &b_vtxz);
	fChain[f]->SetBranchAddress("vtxxErr", &vtxxErr, &b_vtxxErr);
	fChain[f]->SetBranchAddress("vtxyErr", &vtxyErr, &b_vtxyErr);
	fChain[f]->SetBranchAddress("vtxzErr", &vtxzErr, &b_vtxzErr);
	fChain[f]->SetBranchAddress("vtxisValid", &vtxisValid, &b_vtxisValid);
	fChain[f]->SetBranchAddress("vtxisFake", &vtxisFake, &b_vtxisFake);
	fChain[f]->SetBranchAddress("vtxchi2", &vtxchi2, &b_vtxchi2);
	fChain[f]->SetBranchAddress("vtxndof", &vtxndof, &b_vtxndof);
	fChain[f]->SetBranchAddress("vtxnTracks", &vtxnTracks, &b_vtxnTracks);
	fChain[f]->SetBranchAddress("simvtxx", &simvtxx, &b_simvtxx);
	fChain[f]->SetBranchAddress("simvtxy", &simvtxy, &b_simvtxy);
	fChain[f]->SetBranchAddress("simvtxz", &simvtxz, &b_simvtxz);
      
	//Triger
	//if (f==0){
	fChain[f]->SetBranchAddress("trgl1L1GTAlgo", &trgl1L1GTAlgo, &b_trgl1L1GTAlgo);
	fChain[f]->SetBranchAddress("trgl1L1GTTech", &trgl1L1GTTech, &b_trgl1L1GTTech);
	fChain[f]->SetBranchAddress("trgZeroBias", &trgZeroBias, &b_trgZeroBias);
	//}
      
	// // PF Candidate
	// std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * PFCandidatesp4_ = 0;
	// fChain[f]->SetBranchAddress("PFCandidatesparticleId", &PFCandidatesparticleId, &b_PFCandidatesparticleId);
	// fChain[f]->SetBranchAddress("PFCandidatesp4", &PFCandidatesp4_);
	// fChain[f]->SetBranchAddress("PFCandidatesrawEcalEnergy", &PFCandidatesrawEcalEnergy, &b_PFCandidatesrawEcalEnergy);
	// fChain[f]->SetBranchAddress("PFCandidatesrawHcalEnergy", &PFCandidatesrawHcalEnergy, &b_PFCandidatesrawHcalEnergy);
	//Calo Tower
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * CaloTowersp4_ = 0;
	fChain[f]->SetBranchAddress("CaloTowersp4", &CaloTowersp4_);
	fChain[f]->SetBranchAddress("CaloTowersemEnergy", &CaloTowersemEnergy, &b_CaloTowersemEnergy);
	fChain[f]->SetBranchAddress("CaloTowershadEnergy", &CaloTowershadEnergy, &b_CaloTowershadEnergy);
	fChain[f]->SetBranchAddress("CaloTowershasEB", &CaloTowershasEB, &b_CaloTowershasEB);
	fChain[f]->SetBranchAddress("CaloTowershasEE", &CaloTowershasEE, &b_CaloTowershasEE);
	fChain[f]->SetBranchAddress("CaloTowershasHB", &CaloTowershasHB, &b_CaloTowershasHB);
	fChain[f]->SetBranchAddress("CaloTowershasHE", &CaloTowershasHE, &b_CaloTowershasHE);
	fChain[f]->SetBranchAddress("CaloTowershasHF", &CaloTowershasHF, &b_CaloTowershasHF);
	// Rec Hit
	fChain[f]->SetBranchAddress("EcalRecHitenergy", &EcalRecHitenergy, &b_EcalRecHitenergy);
	fChain[f]->SetBranchAddress("EcalRecHitEt", &EcalRecHitEt, &b_EcalRecHitEt);
	fChain[f]->SetBranchAddress("EcalRecHittime", &EcalRecHittime, &b_EcalRecHittime);
	fChain[f]->SetBranchAddress("EcalRecHitieta", &EcalRecHitieta, &b_EcalRecHitieta);
	fChain[f]->SetBranchAddress("EcalRecHitiphi", &EcalRecHitiphi, &b_EcalRecHitiphi);
	fChain[f]->SetBranchAddress("EcalRecHiteta", &EcalRecHiteta, &b_EcalRecHiteta);
	fChain[f]->SetBranchAddress("EcalRecHitphi", &EcalRecHitphi, &b_EcalRecHitphi);
      
	fChain[f]->SetBranchAddress("HBHERecHitenergy", &HBHERecHitenergy, &b_HBHERecHitenergy);
	fChain[f]->SetBranchAddress("HBHERecHitEt", &HBHERecHitEt, &b_HBHERecHitEt);
	fChain[f]->SetBranchAddress("HBHERecHittime", &HBHERecHittime, &b_HBHERecHittime);
	fChain[f]->SetBranchAddress("HBHERecHitieta", &HBHERecHitieta, &b_HBHERecHitieta);
	fChain[f]->SetBranchAddress("HBHERecHitiphi", &HBHERecHitiphi, &b_HBHERecHitiphi);
	fChain[f]->SetBranchAddress("HBHERecHitdepth", &HBHERecHitdepth, &b_HBHERecHitdepth);
	//fChain[f]->SetBranchAddress("HBHERecHitHBHENumRecHits", &HBHERecHitHBHENumRecHits, &b_HBHERecHitHBHENumRecHits);
	fChain[f]->SetBranchAddress("HBHERecHiteta", &HBHERecHiteta, &b_HBHERecHiteta);
	fChain[f]->SetBranchAddress("HBHERecHitphi", &HBHERecHitphi, &b_HBHERecHitphi);
	fChain[f]->SetBranchAddress("HFRecHitenergy", &HFRecHitenergy, &b_HFRecHitenergy);
	fChain[f]->SetBranchAddress("HFRecHitEt", &HFRecHitEt, &b_HFRecHitEt);
	fChain[f]->SetBranchAddress("HFRecHittime", &HFRecHittime, &b_HFRecHittime);
	fChain[f]->SetBranchAddress("HFRecHitieta", &HFRecHitieta, &b_HFRecHitieta);
	fChain[f]->SetBranchAddress("HFRecHitiphi", &HFRecHitiphi, &b_HFRecHitiphi);
	fChain[f]->SetBranchAddress("HFRecHitdepth", &HFRecHitdepth, &b_HFRecHitdepth);
	//fChain[f]->SetBranchAddress("HFRecHitHFNumRecHits", &HFRecHitHFNumRecHits, &b_HFRecHitHFNumRecHits);
	fChain[f]->SetBranchAddress("HFRecHiteta", &HFRecHiteta, &b_HFRecHiteta);
	fChain[f]->SetBranchAddress("HFRecHitphi", &HFRecHitphi, &b_HFRecHitphi);

	//Castor RecHit
	fChain[f]->SetBranchAddress("CastorRecHitEnergy", &CastorRecHitEnergy, &b_CastorRecHitEnergy);
	fChain[f]->SetBranchAddress("CastorRecHitSector", &CastorRecHitSector, &b_CastorRecHitSector);
	fChain[f]->SetBranchAddress("CastorRecHitModule", &CastorRecHitModule, &b_CastorRecHitModule);
        fChain[f]->SetBranchAddress("CastorRecHitisBad", &CastorRecHitisBad, &b_CastorRecHitisBad);
	fChain[f]->SetBranchAddress("CastorRecHitisSaturated", &CastorRecHitisSaturated, &b_CastorRecHitisSaturated);
	fChain[f]->SetBranchAddress("CastorRecHitisDesaturated", &CastorRecHitisDesaturated, &b_CastorRecHitisDesaturated);
	//Castor Tower
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * CastorTowerp4_ = 0;
	fChain[f]->SetBranchAddress("CastorTowerp4", &CastorTowerp4_);
	fChain[f]->SetBranchAddress("CastorToweremEnergy", &CastorToweremEnergy, &b_CastorToweremEnergy);
	fChain[f]->SetBranchAddress("CastorTowerhadEnergy", &CastorTowerhadEnergy, &b_CastorTowerhadEnergy);
	fChain[f]->SetBranchAddress("CastorTowerNrechits", &CastorTowerNrechits, &b_CastorTowerNrechits);

	//PF Clusters
	fChain[f]->SetBranchAddress("hcalPFClustersenergy", &hcalPFClustersenergy, &b_hcalPFClustersenergy);
	fChain[f]->SetBranchAddress("hcalPFClusterscorrectedEnergy", &hcalPFClusterscorrectedEnergy, &b_hcalPFClusterscorrectedEnergy);
	fChain[f]->SetBranchAddress("hcalPFClusterscorrectedEnergyUncertainty", &hcalPFClusterscorrectedEnergyUncertainty, &b_hcalPFClusterscorrectedEnergyUncertainty);
	fChain[f]->SetBranchAddress("hcalPFClusterstime", &hcalPFClusterstime, &b_hcalPFClusterstime);
	fChain[f]->SetBranchAddress("hcalPFClustersdepth", &hcalPFClustersdepth, &b_hcalPFClustersdepth);
	fChain[f]->SetBranchAddress("hcalPFClusterspt", &hcalPFClusterspt, &b_hcalPFClusterspt);
	fChain[f]->SetBranchAddress("hcalPFClustersEt", &hcalPFClustersEt, &b_hcalPFClustersEt);
	fChain[f]->SetBranchAddress("hcalPFClusterseta", &hcalPFClusterseta, &b_hcalPFClusterseta);
	fChain[f]->SetBranchAddress("hcalPFClustersphi", &hcalPFClustersphi, &b_hcalPFClustersphi);
	fChain[f]->SetBranchAddress("hcalPFClusterssize", &hcalPFClusterssize, &b_hcalPFClusterssize);
	fChain[f]->SetBranchAddress("hcalPFClustersisInClean", &hcalPFClustersisInClean, &b_hcalPFClustersisInClean);
	fChain[f]->SetBranchAddress("hcalPFClustersisInUnClean", &hcalPFClustersisInUnClean, &b_hcalPFClustersisInUnClean);
	fChain[f]->SetBranchAddress("hfPFClustersenergy", &hfPFClustersenergy, &b_hfPFClustersenergy);
	fChain[f]->SetBranchAddress("hfPFClusterscorrectedEnergy", &hfPFClusterscorrectedEnergy, &b_hfPFClusterscorrectedEnergy);
	fChain[f]->SetBranchAddress("hfPFClusterscorrectedEnergyUncertainty", &hfPFClusterscorrectedEnergyUncertainty, &b_hfPFClusterscorrectedEnergyUncertainty);
	fChain[f]->SetBranchAddress("hfPFClusterstime", &hfPFClusterstime, &b_hfPFClusterstime);
	fChain[f]->SetBranchAddress("hfPFClustersdepth", &hfPFClustersdepth, &b_hfPFClustersdepth);
	fChain[f]->SetBranchAddress("hfPFClusterspt", &hfPFClusterspt, &b_hfPFClusterspt);
	fChain[f]->SetBranchAddress("hfPFClustersEt", &hfPFClustersEt, &b_hfPFClustersEt);
	fChain[f]->SetBranchAddress("hfPFClusterseta", &hfPFClusterseta, &b_hfPFClusterseta);
	fChain[f]->SetBranchAddress("hfPFClustersphi", &hfPFClustersphi, &b_hfPFClustersphi);
	fChain[f]->SetBranchAddress("hfPFClusterssize", &hfPFClusterssize, &b_hfPFClusterssize);
	fChain[f]->SetBranchAddress("hfPFClustersisInClean", &hfPFClustersisInClean, &b_hfPFClustersisInClean);
	fChain[f]->SetBranchAddress("hfPFClustersisInUnClean", &hfPFClustersisInUnClean, &b_hfPFClustersisInUnClean);
	fChain[f]->SetBranchAddress("ecalPFClustersenergy", &ecalPFClustersenergy, &b_ecalPFClustersenergy);
	fChain[f]->SetBranchAddress("ecalPFClusterscorrectedEnergy", &ecalPFClusterscorrectedEnergy, &b_ecalPFClusterscorrectedEnergy);
	fChain[f]->SetBranchAddress("ecalPFClusterscorrectedEnergyUncertainty", &ecalPFClusterscorrectedEnergyUncertainty, &b_ecalPFClusterscorrectedEnergyUncertainty);
	fChain[f]->SetBranchAddress("ecalPFClusterstime", &ecalPFClusterstime, &b_ecalPFClusterstime);
	fChain[f]->SetBranchAddress("ecalPFClustersdepth", &ecalPFClustersdepth, &b_ecalPFClustersdepth);
	fChain[f]->SetBranchAddress("ecalPFClusterspt", &ecalPFClusterspt, &b_ecalPFClusterspt);
	fChain[f]->SetBranchAddress("ecalPFClustersEt", &ecalPFClustersEt, &b_ecalPFClustersEt);
	fChain[f]->SetBranchAddress("ecalPFClusterseta", &ecalPFClusterseta, &b_ecalPFClusterseta);
	fChain[f]->SetBranchAddress("ecalPFClustersphi", &ecalPFClustersphi, &b_ecalPFClustersphi);
	fChain[f]->SetBranchAddress("ecalPFClusterssize", &ecalPFClusterssize, &b_ecalPFClusterssize);
	fChain[f]->SetBranchAddress("ecalPFClustersisInClean", &ecalPFClustersisInClean, &b_ecalPFClustersisInClean);
	fChain[f]->SetBranchAddress("ecalPFClustersisInUnClean", &ecalPFClustersisInUnClean, &b_ecalPFClustersisInUnClean);
	fChain[f]->SetBranchAddress("bx", &bx, &b_bx);

	// ----------------------- Cut-------------------------//
	lumimin =88;
	lumimax =1000;
	emin = 5.;
	CastorNoiseThreshold = 2./sqrt(5);
	Etabnd = 5.191;//Max Eta for CMS
	//Etabnd = 6.6;//Max Eta for CMS
	etamin =3.139;
	etamax =4.889;//event select
	etamaxMinus =-3.139;
	etaminMinus =-4.889;//event select
	//etamax =5.191;//event select 
	minXiSD=1e-6;
	// ----------------------- Event -----------------------//
  
	Long64_t nentries = fChain[f]->GetEntriesFast();
	cout<<"Entries "<<nentries<<endl;
  
  
	Long64_t nbytes = 0, nb = 0;

	TotNofEvent[f]= 0;
	

	int maxevent=10000;
	//for (Long64_t ev=0; ev<nentries;ev++) {
      
	    for (Long64_t ev=0; ev<maxevent;ev++) {
	
	    Long64_t iev = fChain[f]->LoadTree(ev);
	    if (iev < 0)break;
	    nb = fChain[f]->GetEntry(ev);   nbytes += nb;

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
	    hcalPFClustersCorrEnergy.clear(); 
	    ecalPFClustersCorrEnergy.clear();
	    HFPFClustersCorrEnergy.clear();
	    hcalPFClustersRawEnergy.clear();
	    ecalPFClustersRawEnergy.clear();
	    HFPFClustersRawEnergy.clear();
	    PFClusters.clear();
	    GenDetEventSelect.clear();
	    CastorTower.clear();	   
	    for (int i =0; i<nEtaBins;i++) {
		RecHitHCALEtaSums[i]=0.;
	   
		RecHitECALEtaSums[i]=0.;
		RecHitHFEtaSums[i]=0.;
		CaloEtaSums[i]=0.;
		CastorTowerEtaSums[i]=0.;
		//PFCandEtaSums[i]=0.;
 
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
      
      
	    // if (f==0 && run==251721){
	    if (f==0 && run==247324) RunCut=true;
	    if(f>0) RunCut =true;// for MC
            
	    // ----------------------- fat bunch ------------------//
	    // if (f==0 && bx==208) bxCut=true;
	    bxCut =true;// for MC
            
	    // ----------------------- lumi cut ------------------//
      
	    if (f==0)LumiCut=(lumi>=lumimin &&lumi<lumimax);
	    else LumiCut=true;// for MC
      
	    //if(LumiCut && RunCut ){
	    if(RunCut && LumiCut && bxCut){
	  
		// ----------------------- triger cut --------------//
		if (f==0)TrigerPass =  (trgZeroBias== 1);
		else TrigerPass=true;
	  
		if(TrigerPass){
		    //-- Event Selection for Det Level --//
			
			/////-------from Sebastian Castor event selection ------/////	
			//castorEventSelection = False
			//for i in xrange(0, self.fChain.CastorTowerp4.size()):
			//if self.fChain.CastorTowerp4.at(i).energy() > self.CastorNoiseThreshold*math.sqrt(self.fChain.CastorTowerNrechits.at(i)):
       			//castorEventSelection = True
       			///break

		    /*    for(unsigned long cstr=0; cstr<CastorTowerp4_->size(); cstr++){
			XYZTVector cstrtwr = (*CastorTowerp4_)[cstr];
   			if(cstrtwr.E() > CastorNoiseThreshold * sqrt(CastorTowerNrechits->at(cstr))) {
			    EnergyCutTower=true;
			}
		    }
*/
		    //--------- HF OR-------//
		    for(unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
		    XYZTVector caltwr = (*CaloTowersp4_)[cal];
		    if((abs(caltwr.Eta())>=etamin && abs(caltwr.Eta())<=etamax) && caltwr.E() > emin) {
		      EnergyCutTower=true;
		    	}
		      }
		    if (EnergyCutTower){
			// Gen Level Event selection for Detlevel
			if (f>0){
			    for (unsigned long gd=0; gd<genParticlesp4_->size(); gd++){
				XYZTVector genDet = (*genParticlesp4_)[gd];
			  
				if((abs(genDet.Eta())<= Etabnd)  ) {
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
				RecHitHFEtaSums[selcetetabin]= RecHitHFEtaSums[selcetetabin] + HFRecHitenergy->at(j);
			    }
			}  
			//----------------------- CaloTower  -------------//
			for (unsigned long cal=0; cal<CaloTowersp4_->size(); cal++){
		      
			    XYZTVector caltwr = (*CaloTowersp4_)[cal];
			    selcetetabin=getBin(caltwr.Eta(),EtaBins,nEtaBins);
			    CaloEtaSums[selcetetabin]= CaloEtaSums[selcetetabin] + caltwr.E();
			}
		  
			// // -----------------------  PF cand  -----------------------//
			// //PFCand
			// for (unsigned long pf=0; pf<PFCandidatesp4_->size(); pf++){
			// 	XYZTVector pfcand = (*PFCandidatesp4_)[pf];
			// 	//cout << " pf energy "<<pfcand.E()<< " pf eta "<< pfcand.Eta()<<endl;
		  
			// 	// if(pfcand.E() > 10.) {
			// 	//cout <<"event "<< ev <<" PFCandidates "<<PFCandidatesp4_->size()<<endl; 
			// 	selcetetabin=getBin(pfcand.Eta(),EtaBins, nEtaBins);
		  
			// 	PFCandEtaSums[selcetetabin]= PFCandEtaSums[selcetetabin] + pfcand.E();
		  
			// 	//cout<<" eta "<<pfcand.Eta()<<" bin "<<selcetetabin<<" energy " << pfcand.E() <<" total PF energy "<<PFCandEtaSums[selcetetabin]<<endl;
			// 	// }
			// }//PFCand
		  
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
				PFClustersHFEtaSums[selcetetabin]= PFClustersHFEtaSums[selcetetabin]+hfPFClusterscorrectedEnergy->at(pfc);
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
				PFClustersHFEtaSumsRaw[selcetetabin]= PFClustersHFEtaSumsRaw[selcetetabin]+hfPFClustersenergy->at(pfc);
			    }
			}

			//----------------------- Castor Tower  -------------//

			//totalCastorEnergy = 0.
			//for i in xrange(0, self.fChain.CastorTowerp4.size()):
  			//if self.fChain.CastorTowerp4.at(i).energy() > self.CastorNoiseThreshold*math.sqrt(self.fChain.CastorTowerNrechits.at(i)):
  			//totalCastorEnergy += self.fChain.CastorTowerp4.at(i).energy()


			for (unsigned long cas=0; cas<CastorTowerp4_->size(); cas++){
		      
			    XYZTVector castwr = (*CastorTowerp4_)[cas];
			    selcetetabin=getBin(castwr.Eta(),EtaBins,nEtaBins);
			    CastorTowerEtaSums[selcetetabin]= CastorTowerEtaSums[selcetetabin] + castwr.E();
			}

		  
			// Filling tree
			float totenergyrechit,totenergypfclusters;
			for (int k=0;k<nEtaBins;k++){
			   			    // (HFPFClustersRawEnergy[k]).push_back(PFClustersHFEtaSumsRaw[k]);
			    // (Gen_DetEvntSelct[k]).push_back(GenEtaSumsDetEvntSelct[k]);
			    totenergyrechit =0.;
			    totenergypfclusters=0.;
			    RecHitHCAL.push_back(RecHitHCALEtaSums[k]);
			    RecHitHF.push_back(RecHitHFEtaSums[k]);
			    RecHitECAL.push_back(RecHitECALEtaSums[k]); 
			    totenergyrechit = RecHitHCALEtaSums[k]+RecHitECALEtaSums[k]+RecHitHFEtaSums[k];
			    RecHit.push_back(totenergyrechit);
			    CaloTower.push_back(CaloEtaSums[k]); 
			    CastorTower.push_back(CastorTowerEtaSums[k]); 
			    hcalPFClustersCorrEnergy.push_back(PFClustersHCALEtaSums[k]); 
			    ecalPFClustersCorrEnergy.push_back(PFClustersECALEtaSums[k]);
			    HFPFClustersCorrEnergy.push_back(PFClustersHFEtaSums[k]);
			    hcalPFClustersRawEnergy.push_back(PFClustersHCALEtaSumsRaw[k]);
			    ecalPFClustersRawEnergy.push_back(PFClustersECALEtaSumsRaw[k]);
			    HFPFClustersRawEnergy.push_back(PFClustersHFEtaSumsRaw[k]);
			    totenergypfclusters = PFClustersECALEtaSums[k]+PFClustersHCALEtaSumsRaw[k]+PFClustersHFEtaSumsRaw[k]+CastorTowerEtaSums[k];
			    PFClusters.push_back(totenergypfclusters);
			    GenDetEventSelect.push_back(GenEtaSumsDetEvntSelct[k]);
			    
			}
		  
			TotNofEvent[f]++;
			EnergyFlow_Det[f]->Fill(); 
		    }//HF OR Energy Cut
		}//Trigerpass
	    }//Lumi && Run && bnx Cut
     

	   
	}//Event
	//cout <<filenames[f].c_str()<<"  Total event : "<<TotNofEvent[f]<<endl;
	fOutFile[f]->Write();
	delete fOutFile[f];
	cout<<"Filling tree and write in root file with "<<TotNofEvent[f]<<" events"<<" for file : "<<fname[f].c_str()<<endl;
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
