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
#include "MinBias_TuneMonash_13TeV-pythia8.h"
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
void  EnergyFlow_GenLevel_TuneMonash0T_TreeProducer()

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
    string energyname[13]={"RecHitHCAL","RecHitHF","RecHitECAL","CaloTower",
			   "hcalPFClustersCorrEnergy","ecalPFClustersCorrEnergy","HFPFClustersCorrEnergy",
			   "hcalPFClustersRawEnergy","ecalPFClustersRawEnergy","HFPFClustersRawEnergy",
			   "Gen","Events","CastorTower"};
    
    
   
    
    
    static const Int_t etyp = 13 ;
   
    char title[999];
    
    vector<float> Gen;
    vector<float> Gen_1GeV;
    vector<float> Gen_2GeV;
    vector<float> Gen_HadronElectronCut;
    vector<float> Gen_EM;
    vector<float> Gen_EM_withEnergyCuts;
    vector<float> Gen_Had;
    vector<float> Gen_Had_withEnergyCuts;
    
    TFile *fOutFile[ftyp+1];
    TTree *EnergyFlow_Gen[ftyp+1];
    
    TTree *fChain[fnumber+1];
    TFile *file[fnumber+1];
    int TotNofEventGen[ftyp+1];

    int decade = 0;

    //float Norm=1.0;
    //if (filetype=="MC") Norm=1.0/1.117;
    //else Norm=Norm;

    for (int f=0; f<ftyp; f++){
	//for (int f=FileNumber; f<FileNumber+1; f++){
	//----------------------Creating tree for output--------------//
	sprintf(title,"EFlow_GenLevel_%s_0Tesla_WithNewEta.root",readfilesname.c_str());
	fOutFile[f]= new TFile(title,"RECREATE");
	//sprintf(title,"%s",fname.c_str());
	sprintf(title,"EFlow");
        EnergyFlow_Gen[f]= new TTree(title,title);
        EnergyFlow_Gen[f]->Branch("Gen",&Gen);
        EnergyFlow_Gen[f]->Branch("Gen_1GeV", &Gen_1GeV);
        EnergyFlow_Gen[f]->Branch("Gen_2GeV", &Gen_2GeV);
        EnergyFlow_Gen[f]->Branch("Gen_HadronElectronCut",&Gen_HadronElectronCut);
        EnergyFlow_Gen[f]->Branch("Gen_EM",&Gen_EM);
        EnergyFlow_Gen[f]->Branch("Gen_EM_withEnergyCuts",&Gen_EM_withEnergyCuts);
        EnergyFlow_Gen[f]->Branch("Gen_Had",&Gen_Had);
        EnergyFlow_Gen[f]->Branch("Gen_Had_withEnergyCuts",&Gen_Had_withEnergyCuts);
	TotNofEventGen[f]= 0;
	//------- File Number ---------//                                                                                                                                                
        for (int fn=0; fn < fnumber ; fn++){
            // --- open file in EOS                                                                                                                                                      
            char Fname[999];
            sprintf(Fname,"MinBias_TuneMonash%d3_13TeV-pythia8_MagnetOff_trees.root",(fn+1));
            string sFname(Fname);
            string filenames=readfilesdir+sFname;
            cout<<"file : "<<filenames.c_str()<<endl;
            file[fn]   = TFile::Open(filenames.c_str(),"READ");
	    
	    // --- read tree
	    fChain[fn] = (TTree*)file[fn]->Get("EflowTree/data");
	    //Event
	    //if(filetype =="MC"){
            fChain[f]->SetBranchAddress("run", &run, &b_run);
            fChain[f]->SetBranchAddress("lumi", &lumi, &b_lumi);
            fChain[f]->SetBranchAddress("event", &event, &b_event);
            fChain[f]->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
            fChain[f]->SetBranchAddress("processID", &processID, &b_processID);
            fChain[f]->SetBranchAddress("cmenergy", &cmenergy, &b_cmenergy);
            fChain[f]->SetBranchAddress("puTrueNumInteractions", &puTrueNumInteractions, &b_puTrueNumInteractions);
            fChain[f]->SetBranchAddress("PUNumInteractions", &PUNumInteractions, &b_PUNumInteractions);
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
            std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > * genParticlesp4_ = 0;
            
            fChain[f]->SetBranchAddress("genParticlesp4", &genParticlesp4_);
            fChain[f]->SetBranchAddress("genParticlescharge", &genParticlescharge, &b_genParticlescharge);
            fChain[f]->SetBranchAddress("genParticlespdg", &genParticlespdg, &b_genParticlespdg);
            fChain[f]->SetBranchAddress("genParticlesstatus", &genParticlesstatus, &b_genParticlesstatus);
            fChain[f]->SetBranchAddress("Xix", &Xix, &b_Xix);
            fChain[f]->SetBranchAddress("Xiy", &Xiy, &b_Xiy);
            fChain[f]->SetBranchAddress("XiSD", &XiSD, &b_XiSD);
            fChain[f]->SetBranchAddress("XiDD", &XiDD, &b_XiDD);
            
       // }
	    // ----------------------- Cut-------------------------//
            lumimin =97;
            lumimax =311;
            emin = 5.;
            NoiseCut = 4.;
            Etabnd = 5.205;//Max Eta for CMS
            cmseta = 6.6;//Max Eta for CMS
            etamin =3.152;
            //etamax =4.889;//event select
            etamax =5.205;//event select
            etamaxMinus =-3.152;
            etaminMinus =-5.205;//event select
            minXiSD=1e-6;
            
	    
	    // event 
	    Long64_t nentries = fChain[fn]->GetEntriesFast();
	    cout<<"Entries "<<nentries<<endl;
	    Long64_t nbytes = 0, nb = 0;
	    int maxevent=10000;
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
            
            Gen.clear();
            Gen_1GeV.clear();
            Gen_2GeV.clear();
            Gen_HadronElectronCut.clear();
            Gen_EM.clear();
            Gen_EM_withEnergyCuts.clear();
            Gen_Had.clear();
            Gen_Had_withEnergyCuts.clear();
            
		for (int i =0; i<nEtaBins;i++) {
            
            GenEtaSums[i]=0.;
            GenEtaSums_1GeV[i]=0.;
            GenEtaSums_2GeV[i]=0.;
            GenEtaSums_HadronElectronCut[i]=0.;
            GenEtaSums_EM[i]=0.;
            GenEtaSums_EM_withEnergyCuts[i]=0.;
            GenEtaSums_Had[i]=0.;
            GenEtaSums_Had_withEnergyCuts[i]=0.;
            
		}
		
		// Ini. for each event
            XiCutGen=false;
            EnergyCutGen1GeV=false;
            EnergyCutGen2GeV=false;
        //--------------------XiCut-----------------------//
		if (filetype =="MC" && ( (Xix > minXiSD && Xiy < minXiSD) || (Xix < minXiSD && Xiy > minXiSD) )) XiCutGen=true;

		
        if(XiCutGen){
            
            for (unsigned long g=0; g<genParticlesp4_->size(); g++){
                XYZTVector gen = (*genParticlesp4_)[g];
                ElectronPhotonCut=false;
                EBEnergyCut=false;
                EEEnergyCut=false;
                HBEnergyCut=false;
                HEEnergyCut=false;
                HFEnergyCut=false;
                HF_EPEnergyCut=false;
                HadronCut=false;
                //Electron-Photon Cut begins
                if(genParticlespdg->at(g) ==11 || genParticlespdg->at(g) ==22) {
                    ElectronPhotonCut=true;
                    //EB
                    if((abs(gen.Eta())<= 1.479)  ) {
                        if(gen.E()> 0.23 ) {
                            EBEnergyCut=true;
                        }
                    }
                    //EE
                    if((abs(gen.Eta())>=1.479) && (abs(gen.Eta())<= 3.0 )  ) {
                        if(gen.E()> 0.6 && gen.Et()> 0.15) {
                            EEEnergyCut=true;
                        }
                    }
                    //HF
                    if((abs(gen.Eta()) >= 2.8) && (abs(gen.Eta()) <= 5.191)   ) {
                        if(gen.E()> 1.4 ) {
                            HF_EPEnergyCut=true;
                        }
                    }
                    
                } //Electron-Photon Cut ends
                //Hadron Cut begins
                if(!ElectronPhotonCut ) {
                    HadronCut=true;
                    //HB
                    if((abs(gen.Eta())<= 1.4)  ) {
                        if(gen.E()> 1. ) {
                            HBEnergyCut=true;
                        }
                    }
                    //HE
                    if((abs(gen.Eta())>=1.3) && (abs(gen.Eta())<= 3.0 ) ) {
                        if(gen.E()> 1.1 ) {
                            HEEnergyCut=true;
                        }
                    }
                    //HF
                    if((abs(gen.Eta()) >= 2.8) && (abs(gen.Eta()) <= 5.191)   ) {
                        if(gen.E()> 1.4 ) {
                            HFEnergyCut=true;
                        }
                    }
                    
                } ////Hadron Cut ends
                // -- Gen level apply only XiSD cut --//
                if((abs(gen.Eta())<= cmseta)  ) {
                    
                    selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                    GenEtaSums[selcetetabin]= GenEtaSums[selcetetabin] + gen.E();
                    
                    //-- Gen level apply Energy cut for each sub det both  EM and HAD particles--//
                    if (HFEnergyCut || HBEnergyCut || HEEnergyCut || EEEnergyCut || EBEnergyCut || HF_EPEnergyCut){
                        selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                        GenEtaSums_HadronElectronCut[selcetetabin]= GenEtaSums_HadronElectronCut[selcetetabin] + gen.E();
                    }
                    
                    //--Gen level apply only EM particle cut --//
                    if (ElectronPhotonCut){
                        
                        selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                        GenEtaSums_EM[selcetetabin]= GenEtaSums_EM[selcetetabin] + gen.E();
                        //--Gen level apply EM particle & energy cut --//
                        if (EBEnergyCut || EEEnergyCut || HF_EPEnergyCut){
                            selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                            GenEtaSums_EM_withEnergyCuts[selcetetabin]= GenEtaSums_EM_withEnergyCuts[selcetetabin] + gen.E();
                        }
                    }
                    
                    //--Gen level apply only HAD particle cut --//
                    if (HadronCut){
                        
                        selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                        GenEtaSums_Had[selcetetabin]= GenEtaSums_Had[selcetetabin] + gen.E();
                        //--Gen level apply HAD particle & energy cut --//
                        if (HBEnergyCut || HEEnergyCut || HFEnergyCut){
                            
                            selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                            GenEtaSums_Had_withEnergyCuts[selcetetabin]= GenEtaSums_Had_withEnergyCuts[selcetetabin] + gen.E();
                            
                        }
                    }
                    
                    
                    //--Gen level 1 GeV cut --// 
                    if(gen.E()> 1. ) {
                        EnergyCutGen1GeV=true;
                        selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                        
                        GenEtaSums_1GeV[selcetetabin]= GenEtaSums_1GeV[selcetetabin] + gen.E();
                        
                    }
                    
                    //--Gen level 2 GeV cut --//
                    if(gen.E()> 2. ) {
                        EnergyCutGen2GeV=true;
                        selcetetabin=getBin(gen.Eta(),EtaBins, nEtaBins);
                        
                        GenEtaSums_2GeV[selcetetabin]= GenEtaSums_2GeV[selcetetabin] + gen.E();
                        
                    }
                    
                    //cout<<" eta "<<gen.Eta()<<" bin "<<selcetetabin<<" energy " << gen.E() <<" total Gen energy "<<GenEtaSums[selcetetabin]<<endl;
                }// Gen particle eta cut
            }//Gen particle loop 
            TotNofEventGen[f]++;
            
            for (int k=0;k<nEtaBins;k++){
                
                Gen.push_back(GenEtaSums[k]);
                Gen_1GeV.push_back(GenEtaSums_1GeV[k]);
                Gen_2GeV.push_back(GenEtaSums_1GeV[k]);
                Gen_HadronElectronCut.push_back(GenEtaSums_HadronElectronCut[k]);
                Gen_EM.push_back(GenEtaSums_EM[k]);
                Gen_EM_withEnergyCuts.push_back(GenEtaSums_EM_withEnergyCuts[k]);
                Gen_Had.push_back(GenEtaSums_Had[k]);
                Gen_Had_withEnergyCuts.push_back(GenEtaSums_Had_withEnergyCuts[k]);
                
                
            }
            EnergyFlow_Gen[f]->Fill();
            
        }//xi cut
            
	    }//Event
	    cout <<"  Total event : "<<TotNofEventGen[f]<<endl;
	}//file number
	fOutFile[f]->Write();
	delete fOutFile[f];
	cout<<"Filling tree and write in root file with "<<TotNofEventGen[f]<<" events"<<" for file : "<<readfilesname.c_str()<<endl;
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
