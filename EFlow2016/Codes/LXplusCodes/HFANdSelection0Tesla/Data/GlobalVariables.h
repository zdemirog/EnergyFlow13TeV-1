vector<float>   *CaloTower=0;
vector<float>   *HCALTower=0;
vector<float>   *RecHit=0;
vector<float>   *PFClusters=0;
vector<float>   *Gen=0;
vector<float>   *CastorTower=0;
vector<float>   *RecHitCastor=0;
vector<float>   *PFCandidate=0;

vector<float>   *CaloTowersemEnergy=0;
vector<float>   *CaloTowershadEnergy=0;
vector<int>     *CaloTowershasEB=0;
vector<int>     *CaloTowershasEE=0;
vector<int>     *CaloTowershasHB=0;
vector<int>     *CaloTowershasHE=0;
vector<int>     *CaloTowershasHF=0;

vector<float>   *CastorRecHitEnergy=0;
vector<int>     *CastorRecHitSector=0;
vector<int>     *CastorRecHitModule=0;
vector<int>     *CastorRecHitisBad=0;
vector<int>     *CastorRecHitisSaturated=0;
vector<int>     *CastorRecHitisDesaturated=0;
vector<float>   *CastorToweremEnergy=0;
vector<float>   *CastorTowerhadEnergy=0;
vector<int>     *CastorTowerNrechits=0;
vector<float>   *EcalRecHitenergy=0;
vector<float>   *EcalRecHitEt=0;
vector<float>   *EcalRecHittime=0;
vector<int>     *EcalRecHitieta=0;
vector<int>     *EcalRecHitiphi=0;
vector<float>   *EcalRecHiteta=0;
vector<float>   *EcalRecHitphi=0;

vector<float>   *HBHERecHitenergy=0;
vector<float>   *HBHERecHitEt=0;
vector<float>   *HBHERecHittime=0;
vector<int>     *HBHERecHitieta=0;
vector<int>     *HBHERecHitiphi=0;
vector<int>     *HBHERecHitdepth=0;
vector<float>   *HBHERecHiteta=0;
vector<float>   *HBHERecHitphi=0;
vector<float>   *HFRecHitenergy=0;
vector<float>   *HFRecHitEt=0;
vector<float>   *HFRecHittime=0;
vector<int>     *HFRecHitieta=0;
vector<int>     *HFRecHitiphi=0;
vector<int>     *HFRecHitdepth=0;
vector<float>   *HFRecHiteta=0;
vector<float>   *HFRecHitphi=0;
vector<int>     *PFCandidatesparticleId=0;
vector<float>   *PFCandidatesrawEcalEnergy=0;
vector<float>   *PFCandidatesrawHcalEnergy=0;

vector<int>     *trgl1L1GTAlgo=0;
vector<int>     *trgl1L1GTTech=0;
vector<float>   *vtxx=0;
vector<float>   *vtxy=0;
vector<float>   *vtxz=0;
vector<float>   *vtxxErr=0;
vector<float>   *vtxyErr=0;
vector<float>   *vtxzErr=0;
vector<int>     *vtxisValid=0;
vector<int>     *vtxisFake=0;
vector<float>   *vtxchi2=0;
vector<int>     *vtxndof=0;
vector<int>     *vtxnTracks=0;
Int_t           trgZeroBias;

vector<int>     *genParticlescharge=0;
vector<int>     *genParticlespdg=0;
vector<int>     *genParticlesstatus=0;
vector<float>   *ecalPFClustersenergy=0;
vector<float>   *ecalPFClusterscorrectedEnergy=0;
vector<float>   *ecalPFClusterscorrectedEnergyUncertainty=0;
vector<float>   *ecalPFClusterstime=0;
vector<float>   *ecalPFClustersdepth=0;
vector<float>   *ecalPFClusterspt=0;
vector<float>   *ecalPFClustersEt=0;
vector<float>   *ecalPFClusterseta=0;
vector<float>   *ecalPFClustersphi=0;
vector<int>     *ecalPFClusterssize=0;
vector<int>     *ecalPFClustersisInClean=0;
vector<int>     *ecalPFClustersisInUnClean=0;
vector<float>   *hcalPFClustersenergy=0;
vector<float>   *hcalPFClusterscorrectedEnergy=0;
vector<float>   *hcalPFClusterscorrectedEnergyUncertainty=0;
vector<float>   *hcalPFClusterstime=0;
vector<float>   *hcalPFClustersdepth=0;
vector<float>   *hcalPFClusterspt=0;
vector<float>   *hcalPFClustersEt=0;
vector<float>   *hcalPFClusterseta=0;
vector<float>   *hcalPFClustersphi=0;
vector<int>     *hcalPFClusterssize=0;
vector<int>     *hcalPFClustersisInClean=0;
vector<int>     *hcalPFClustersisInUnClean=0;
vector<float>   *hfPFClustersenergy=0;
vector<float>   *hfPFClusterscorrectedEnergy=0;
vector<float>   *hfPFClusterscorrectedEnergyUncertainty=0;
vector<float>   *hfPFClusterstime=0;
vector<float>   *hfPFClustersdepth=0;
vector<float>   *hfPFClusterspt=0;
vector<float>   *hfPFClustersEt=0;
vector<float>   *hfPFClusterseta=0;
vector<float>   *hfPFClustersphi=0;
vector<int>     *hfPFClusterssize=0;
vector<int>     *hfPFClustersisInClean=0;
vector<int>     *hfPFClustersisInUnClean=0;

Int_t           run;
Int_t           lumi;
Int_t           event;
Int_t           bx;
Float_t         genWeight;
Float_t         alphaQCD;
Float_t         qScale;
Int_t           processID;
Float_t         Xix;
Float_t         Xiy;
Float_t         XiSD;
Float_t         XiDD;
Float_t         cmenergy;
Float_t         instLumiPerBX;
Float_t         puTrueNumInteractions;
Float_t         PUNumInteractions;

Float_t         simvtxx;
Float_t         simvtxy;
Float_t         simvtxz;


// List of branches
TBranch        *b_CaloTower;
TBranch        *b_HCALTower;
TBranch        *b_RecHit;
TBranch        *b_PFClusters;
TBranch        *b_Gen;
TBranch        *b_CastorTower;
TBranch        *b_RecHitCastor;
TBranch        *b_PFCandidate;

TBranch        *b_CaloTowersemEnergy;   //!
TBranch        *b_CaloTowershadEnergy;   //!
TBranch        *b_CaloTowershasEB;   //!
TBranch        *b_CaloTowershasEE;   //!
TBranch        *b_CaloTowershasHB;   //!
TBranch        *b_CaloTowershasHE;   //!
TBranch        *b_CaloTowershasHF;   //!

TBranch        *b_CastorRecHitEnergy;   //!
TBranch        *b_CastorRecHitSector;   //!
TBranch        *b_CastorRecHitModule;   //!
TBranch        *b_CastorRecHitisBad;   //!
TBranch        *b_CastorRecHitisSaturated;   //!
TBranch        *b_CastorRecHitisDesaturated;   //!
TBranch        *b_CastorToweremEnergy;   //!
TBranch        *b_CastorTowerhadEnergy;   //!
TBranch        *b_CastorTowerNrechits;   //!
TBranch        *b_EcalRecHitenergy;   //!
TBranch        *b_EcalRecHitEt;   //!
TBranch        *b_EcalRecHittime;   //!
TBranch        *b_EcalRecHitieta;   //!
TBranch        *b_EcalRecHitiphi;   //!
TBranch        *b_EcalRecHiteta;   //!
TBranch        *b_EcalRecHitphi;   //!


TBranch        *b_genParticlescharge;   //!
TBranch        *b_genParticlespdg;   //!
TBranch        *b_genParticlesstatus;   //!
TBranch        *b_HBHERecHitenergy;   //!
TBranch        *b_HBHERecHitEt;   //!
TBranch        *b_HBHERecHittime;   //!
TBranch        *b_HBHERecHitieta;   //!
TBranch        *b_HBHERecHitiphi;   //!
TBranch        *b_HBHERecHitdepth;   //!
TBranch        *b_HBHERecHiteta;   //!
TBranch        *b_HBHERecHitphi;   //!
TBranch        *b_HFRecHitenergy;   //!
TBranch        *b_HFRecHitEt;   //!
TBranch        *b_HFRecHittime;   //!
TBranch        *b_HFRecHitieta;   //!
TBranch        *b_HFRecHitiphi;   //!
TBranch        *b_HFRecHitdepth;   //!
TBranch        *b_HFRecHiteta;   //!
TBranch        *b_HFRecHitphi;   //!
TBranch        *b_HFRecHitLong_energy;   //!
TBranch        *b_HFRecHitShort_energy;   //!
TBranch        *b_PFCandidatesrawEcalEnergy;   //!
TBranch        *b_PFCandidatesrawHcalEnergy;   //!
TBranch        *b_PFCandidatesparticleId;   //!

TBranch        *b_trgl1L1GTAlgo;   //!
TBranch        *b_trgl1L1GTTech;   //!
TBranch        *b_vtxx;   //!
TBranch        *b_vtxy;   //!
TBranch        *b_vtxz;   //!
TBranch        *b_vtxxErr;   //!
TBranch        *b_vtxyErr;   //!
TBranch        *b_vtxzErr;   //!
TBranch        *b_vtxisValid;   //!
TBranch        *b_vtxisFake;   //!
TBranch        *b_vtxchi2;   //!
TBranch        *b_vtxndof;   //!
TBranch        *b_vtxnTracks;   //!
TBranch        *b_trgZeroBias;   //!


TBranch        *b_ecalPFClustersenergy;   //!
TBranch        *b_ecalPFClusterscorrectedEnergy;   //!
TBranch        *b_ecalPFClusterscorrectedEnergyUncertainty;   //!
TBranch        *b_ecalPFClusterstime;   //!
TBranch        *b_ecalPFClustersdepth;   //!
TBranch        *b_ecalPFClusterspt;   //!
TBranch        *b_ecalPFClustersEt;   //!
TBranch        *b_ecalPFClusterseta;   //!
TBranch        *b_ecalPFClustersphi;   //!
TBranch        *b_ecalPFClusterssize;   //!
TBranch        *b_ecalPFClustersisInClean;   //!
TBranch        *b_ecalPFClustersisInUnClean;   //!
TBranch        *b_hcalPFClustersenergy;   //!
TBranch        *b_hcalPFClusterscorrectedEnergy;   //!
TBranch        *b_hcalPFClusterscorrectedEnergyUncertainty;   //!
TBranch        *b_hcalPFClusterstime;   //!
TBranch        *b_hcalPFClustersdepth;   //!
TBranch        *b_hcalPFClusterspt;   //!
TBranch        *b_hcalPFClustersEt;   //!
TBranch        *b_hcalPFClusterseta;   //!
TBranch        *b_hcalPFClustersphi;   //!
TBranch        *b_hcalPFClusterssize;   //!
TBranch        *b_hcalPFClustersisInClean;   //!
TBranch        *b_hcalPFClustersisInUnClean;   //!
TBranch        *b_hfPFClustersenergy;   //!
TBranch        *b_hfPFClusterscorrectedEnergy;   //!
TBranch        *b_hfPFClusterscorrectedEnergyUncertainty;   //!
TBranch        *b_hfPFClusterstime;   //!
TBranch        *b_hfPFClustersdepth;   //!
TBranch        *b_hfPFClusterspt;   //!
TBranch        *b_hfPFClustersEt;   //!
TBranch        *b_hfPFClusterseta;   //!
TBranch        *b_hfPFClustersphi;   //!
TBranch        *b_hfPFClusterssize;   //!
TBranch        *b_hfPFClustersisInClean;   //!
TBranch        *b_hfPFClustersisInUnClean;   //!


TBranch        *b_run;   //!
TBranch        *b_lumi;   //!
TBranch        *b_event;   //!
TBranch        *b_bx;   //!
TBranch        *b_genWeight;   //!
TBranch        *b_alphaQCD;   //!
TBranch        *b_qScale;   //!
TBranch        *b_processID;   //!
TBranch        *b_Xix;   //!
TBranch        *b_Xiy;   //!
TBranch        *b_XiSD;   //!
TBranch        *b_XiDD;   //!
TBranch        *b_cmenergy;   //!
TBranch        *b_instLumiPerBX;   //!
TBranch        *b_puTrueNumInteractions;   //!
TBranch        *b_PUNumInteractions;   //!
TBranch        *b_simvtxx;   //!
TBranch        *b_simvtxy;   //!
TBranch        *b_simvtxz;   //!


/* int nEtaBins =82; */
/* float EtaBins[83]={-5.191,-4.889,-4.716,-4.538,-4.363,-4.191,-4.013 */
/* 		    ,-3.839,-3.664,-3.489,-3.314,-3.139,-2.964,-2.853 */
/* 		    ,-2.650,-2.500,-2.322,-2.172,-2.043,-1.930,-1.830 */
/* 		    ,-1.740,-1.653,-1.566,-1.479,-1.392,-1.305,-1.218 */
/* 		    ,-1.131,-1.044,-0.957,-0.879,-0.783,-0.696,-0.609 */
/* 		    ,-0.522,-0.435,-0.348,-0.261,-0.174,-0.087,0.000 */
/* 		    ,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696 */
/* 		    ,0.783,0.879,0.957,1.044,1.131,1.218,1.305,1.392 */
/* 		    ,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172 */
/* 		    ,2.322,2.500,2.650,2.853,2.964,3.139,3.314,3.489 */
/* 		    ,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889 */
/* 		    ,5.191}; */
/* double abseta_bin[42]={0.000 */
/* 		       ,0.087,0.174,0.261,0.348,0.435,0.522,0.609,0.696 */
/* 		       ,0.783,0.879,0.957,1.044,1.131,1.218,1.305,1.392 */
/* 		       ,1.479,1.566,1.653,1.740,1.830,1.930,2.043,2.172 */
/* 		       ,2.322,2.500,2.650,2.853,2.964,3.139,3.314,3.489 */
/* 		       ,3.664,3.839,4.013,4.191,4.363,4.538,4.716,4.889 */
/* 		       ,5.191}; */
//vector<double> tmp;
//map<int, vector<double> > vCaloEtaSums;
int nEtaBins =30;
float EtaBins[31]={-6.6,-5.191,-4.889,-4.538,-4.191,-3.839,-3.489,-3.139,-2.650,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,
		   0,0.348,0.696,1.044,1.392,1.74,2.172,2.650,3.139,3.489,3.839,4.191,4.538,4.889,5.191, 6.6};

int nHBins =30;
float HBins[31]={0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.};

static const Int_t BIN= 30;
Double_t RecHitEtaSums[BIN];

Double_t RecHitHCALEtaSums[BIN];
//ECAL
Double_t RecHitECALEtaSums[BIN];
//
Double_t RecHitHFEtaSums[BIN];

Double_t CaloEtaSums[BIN];
Double_t HCALTowerEtaSums[BIN];
Double_t CastorTowerEtaSums[BIN];
Double_t RecHitCastorEtaSums[BIN];
Double_t PFCandEtaSums[BIN];

Double_t PFClustersECALEtaSums[BIN];
Double_t PFClustersHCALEtaSums[BIN];
Double_t PFClustersHFEtaSums[BIN];
Double_t PFClustersEtaSums[BIN];

Double_t PFClustersECALEtaSumsUp[BIN];
Double_t PFClustersHCALEtaSumsUp[BIN];
Double_t PFClustersHFEtaSumsUp[BIN];
Double_t PFClustersEtaSumsUp[BIN];

Double_t PFClustersECALEtaSumsLo[BIN];
Double_t PFClustersHCALEtaSumsLo[BIN];
Double_t PFClustersHFEtaSumsLo[BIN];
Double_t PFClustersEtaSumsLo[BIN];


Double_t PFClustersECALEtaSumsRaw[BIN];
Double_t PFClustersHCALEtaSumsRaw[BIN];
Double_t PFClustersHFEtaSumsRaw[BIN];
Double_t PFClustersEtaSumsRaw[BIN];

Double_t GenEtaSums[BIN];
Double_t GenEtaSums_1GeV[BIN];
Double_t GenEtaSums_2GeV[BIN];
Double_t GenEtaSums_HadronElectronCut[BIN];
Double_t GenEtaSums_EM[BIN];
Double_t GenEtaSums_EM_withEnergyCuts[BIN];
Double_t GenEtaSums_Had[BIN];
Double_t GenEtaSums_Had_withEnergyCuts[BIN];
Double_t GenEtaSumsDetEvntSelct[BIN];
Double_t GenEtaSumsDetEvntSelctHF[BIN];
         
int selcetetabin;
float Etabnd;
float cmseta;
bool EnergyCutRecHit(false);
bool EnergyCutRecHitMinus(false);
bool EnergyCutTowerMinus(false);
bool EnergyCutTowerPlus(false);
bool EnergyCutTower(false);
bool EnergyCutRecHitPlus(false);
bool EnergyCutRecHitHFAnd(false);
bool EnergyCutPFCluster(false);
bool XiCutGen(false);
bool EnergyCutGen1GeV(false);
bool EnergyCutGen2GeV(false);

bool ElectronPhotonCut(false);
bool EBEnergyCut(false);
bool EEEnergyCut(false);
bool HBEnergyCut(false);
bool HEEnergyCut(false);
bool HFEnergyCut(false);
bool HF_EPEnergyCut(false);
bool HadronCut(false);


bool EtaCut(false);
bool RunCut(false);
bool bxCut(false);
//bool VertexCut(false);
bool LumiCut(false);
bool TrigerPass(false);
int lumimin,lumimax;
float etamin,etamax,emin,minXiSD,etaminMinus,etamaxMinus, CastorNoiseThreshold;
