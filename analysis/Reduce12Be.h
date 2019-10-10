#ifndef Reduce_h
#define Reduce_h

#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TPolyLine3D.h>
#include <TRandom3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>

#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>

#include "TypeDef.h"

// Header file for the classes stored in the TTree if any.

class Reduce {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t           YdMul;
    Int_t           YdChannel[128];   //[YdMul]
    Int_t           YdEnergyRaw[128];   //[YdMul]
    Float_t         YdEnergy[128];   //[YdMul]
    Int_t           YdRing[128];   //[YdMul]
    Int_t           YdSector[128];   //[YdMul]
    Int_t           YuMul;
    Int_t           YuChannel[128];   //[YuMul]
    Int_t           YuEnergyRaw[128];   //[YuMul]
    Float_t         YuEnergy[128];   //[YuMul]
    Int_t           YuRing[128];   //[YuMul]
    Int_t           YuSector[128];   //[YuMul]
    Int_t           CsI1Mul;
    Int_t           CsI1Channel[128];   //[CsI1Mul]
    Float_t         CsI1EnergyRaw[128];   //[CsI1Mul]
    Int_t           CsI2Mul;
    Int_t           CsI2Channel[128];   //[CsI2Mul]
    Float_t         CsI2EnergyRaw[128];   //[CsI2Mul]
    Int_t           ICChannel;
    Float_t         ICEnergyRaw;
    Int_t           Sd1rMul;
    Int_t           Sd1rChannel[128];   //[Sd1rMul]
    Int_t           Sd1rEnergyRaw[128];   //[Sd1rMul]
    Float_t         Sd1rEnergy[128];   //[Sd1rMul]
    Int_t           Sd1sMul;
    Int_t           Sd1sChannel[128];   //[Sd1sMul]
    Int_t           Sd1sEnergyRaw[128];   //[Sd1sMul]
    Float_t         Sd1sEnergy[128];   //[Sd1sMul]
    Int_t           Sd2rMul;
    Int_t           Sd2rChannel[128];   //[Sd2rMul]
    Int_t           Sd2rEnergyRaw[128];   //[Sd2rMul]
    Float_t         Sd2rEnergy[128];   //[Sd2rMul]
    Int_t           Sd2sMul;
    Int_t           Sd2sChannel[128];   //[Sd2sMul]
    Int_t           Sd2sEnergyRaw[128];   //[Sd2sMul]
    Float_t         Sd2sEnergy[128];   //[Sd2sMul]
    Int_t           SurMul;
    Int_t           SurChannel[128];   //[SurMul]
    Int_t           SurEnergyRaw[128];   //[SurMul]
    Int_t           SusMul;
    Int_t           SusChannel[128];   //[SusMul]
    Int_t           SusEnergyRaw[128];   //[SusMul]

    // List of branches
    TBranch        *b_YdMul;   //!
    TBranch        *b_YdChannel;   //!
    TBranch        *b_YdEnergyRaw;   //!
    TBranch        *b_YdEnergy;   //!
    TBranch        *b_YdRing;   //!
    TBranch        *b_YdSector;   //!
    TBranch        *b_YuMul;   //!
    TBranch        *b_YuChannel;   //!
    TBranch        *b_YuEnergyRaw;   //!
    TBranch        *b_YuEnergy;   //!
    TBranch        *b_YuRing;   //!
    TBranch        *b_YuSector;   //!
    TBranch        *b_CsI1Mul;   //!
    TBranch        *b_CsI1Channel;   //!
    TBranch        *b_CsI1EnergyRaw;   //!
    TBranch        *b_CsI2Mul;   //!
    TBranch        *b_CsI2Channel;   //!
    TBranch        *b_CsI2EnergyRaw;   //!
    TBranch        *b_ICChannel;   //!
    TBranch        *b_ICEnergyRaw;   //!
    TBranch        *b_Sd1rMul;   //!
    TBranch        *b_Sd1rChannel;   //!
    TBranch        *b_Sd1rEnergyRaw;   //!
    TBranch        *b_Sd1rEnergy;   //!
    TBranch        *b_Sd1sMul;   //!
    TBranch        *b_Sd1sChannel;   //!
    TBranch        *b_Sd1sEnergyRaw;   //!
    TBranch        *b_Sd1sEnergy;   //!
    TBranch        *b_Sd2rMul;   //!
    TBranch        *b_Sd2rChannel;   //!
    TBranch        *b_Sd2rEnergyRaw;   //!
    TBranch        *b_Sd2rEnergy;   //!
    TBranch        *b_Sd2sMul;   //!
    TBranch        *b_Sd2sChannel;   //!
    TBranch        *b_Sd2sEnergyRaw;   //!
    TBranch        *b_Sd2sEnergy;   //!
    TBranch        *b_SurMul;   //!
    TBranch        *b_SurChannel;   //!
    TBranch        *b_SurEnergyRaw;   //!
    TBranch        *b_SusMul;   //!
    TBranch        *b_SusChannel;   //!
    TBranch        *b_SusEnergyRaw;   //!

    Reduce(TTree *tree=0);
    virtual ~Reduce();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);

    // Calibration files
    void ReadCalibrations();
    std::pair<double, double> sd1rCalibrations[24];
    std::pair<double, double> sd1sCalibrations[32];
    std::pair<double, double> sd2rCalibrations[24];
    std::pair<double, double> sd2sCalibrations[32];
    std::pair<double, double> ydCalibrations[128];
    std::pair<double, double> yuCalibrations[128];

    // Tree variables
    void InitTree();
    TTree *outTree;
    int fYuMul;
    int fYuChannel[128];
    int fYuRing[128];
    int fYuSector[128];
    int fYuEnergyRaw[128];
    float fYuEnergy[128];
    int fSd1Mul;
    int fSd1ChannelRing[128];
    int fSd1ChannelSector[128];
    int fSd1Ring[128];
    int fSd1Sector[128];
    int fSd1EnergyRawRing[128];
    int fSd1EnergyRawSector[128];
    float fSd1EnergyRing[128];
    float fSd1EnergySector[128];
    int fSd2Mul;
    int fSd2ChannelRing[128];
    int fSd2ChannelSector[128];
    int fSd2Ring[128];
    int fSd2Sector[128];
    int fSd2EnergyRawRing[128];
    int fSd2EnergyRawSector[128];
    float fSd2EnergyRing[128];
    float fSd2EnergySector[128];
    int fICChannel;
    float fICEnergyRaw;
};

#endif

#ifdef Reduce_cxx
Reduce::Reduce(TTree *tree) : fChain(0) {
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if(tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/joshhooker/Desktop/data/12Be_dp/NewDecode_5021.root");
        if(!f || !f->IsOpen()) {
            f = new TFile("/Users/joshhooker/Desktop/data/12Be_dp/NewDecode_5021.root");
        }
        f->GetObject("AutoTree",tree);

    }
    Init(tree);
}

Reduce::~Reduce() {
    if(!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t Reduce::GetEntry(Long64_t entry) {
    // Read contents of entry.
    if(!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t Reduce::LoadTree(Long64_t entry) {
    // Set the environment to read one entry
    if(!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if(centry < 0) return centry;
    if(fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
        Notify();
    }
    return centry;
}

void Reduce::Init(TTree *tree) {
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if(!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("YdMul", &YdMul, &b_YdMul);
    fChain->SetBranchAddress("YdChannel", YdChannel, &b_YdChannel);
    fChain->SetBranchAddress("YdEnergyRaw", YdEnergyRaw, &b_YdEnergyRaw);
    fChain->SetBranchAddress("YdEnergy", YdEnergy, &b_YdEnergy);
    fChain->SetBranchAddress("YdRing", YdRing, &b_YdRing);
    fChain->SetBranchAddress("YdSector", YdSector, &b_YdSector);
    fChain->SetBranchAddress("YuMul", &YuMul, &b_YuMul);
    fChain->SetBranchAddress("YuChannel", YuChannel, &b_YuChannel);
    fChain->SetBranchAddress("YuEnergyRaw", YuEnergyRaw, &b_YuEnergyRaw);
    fChain->SetBranchAddress("YuEnergy", YuEnergy, &b_YuEnergy);
    fChain->SetBranchAddress("YuRing", YuRing, &b_YuRing);
    fChain->SetBranchAddress("YuSector", YuSector, &b_YuSector);
    fChain->SetBranchAddress("CsI1Mul", &CsI1Mul, &b_CsI1Mul);
    fChain->SetBranchAddress("CsI1Channel", CsI1Channel, &b_CsI1Channel);
    fChain->SetBranchAddress("CsI1EnergyRaw", CsI1EnergyRaw, &b_CsI1EnergyRaw);
    fChain->SetBranchAddress("CsI2Mul", &CsI2Mul, &b_CsI2Mul);
    fChain->SetBranchAddress("CsI2Channel", CsI2Channel, &b_CsI2Channel);
    fChain->SetBranchAddress("CsI2EnergyRaw", CsI2EnergyRaw, &b_CsI2EnergyRaw);
    fChain->SetBranchAddress("ICChannel", &ICChannel, &b_ICChannel);
    fChain->SetBranchAddress("ICEnergyRaw", &ICEnergyRaw, &b_ICEnergyRaw);
    fChain->SetBranchAddress("Sd1rMul", &Sd1rMul, &b_Sd1rMul);
    fChain->SetBranchAddress("Sd1rChannel", Sd1rChannel, &b_Sd1rChannel);
    fChain->SetBranchAddress("Sd1rEnergyRaw", Sd1rEnergyRaw, &b_Sd1rEnergyRaw);
    fChain->SetBranchAddress("Sd1rEnergy", Sd1rEnergy, &b_Sd1rEnergy);
    fChain->SetBranchAddress("Sd1sMul", &Sd1sMul, &b_Sd1sMul);
    fChain->SetBranchAddress("Sd1sChannel", Sd1sChannel, &b_Sd1sChannel);
    fChain->SetBranchAddress("Sd1sEnergyRaw", Sd1sEnergyRaw, &b_Sd1sEnergyRaw);
    fChain->SetBranchAddress("Sd1sEnergy", Sd1sEnergy, &b_Sd1sEnergy);
    fChain->SetBranchAddress("Sd2rMul", &Sd2rMul, &b_Sd2rMul);
    fChain->SetBranchAddress("Sd2rChannel", Sd2rChannel, &b_Sd2rChannel);
    fChain->SetBranchAddress("Sd2rEnergyRaw", Sd2rEnergyRaw, &b_Sd2rEnergyRaw);
    fChain->SetBranchAddress("Sd2rEnergy", Sd2rEnergy, &b_Sd2rEnergy);
    fChain->SetBranchAddress("Sd2sMul", &Sd2sMul, &b_Sd2sMul);
    fChain->SetBranchAddress("Sd2sChannel", Sd2sChannel, &b_Sd2sChannel);
    fChain->SetBranchAddress("Sd2sEnergyRaw", Sd2sEnergyRaw, &b_Sd2sEnergyRaw);
    fChain->SetBranchAddress("Sd2sEnergy", Sd2sEnergy, &b_Sd2sEnergy);
    fChain->SetBranchAddress("SurMul", &SurMul, &b_SurMul);
    fChain->SetBranchAddress("SurChannel", SurChannel, &b_SurChannel);
    fChain->SetBranchAddress("SurEnergyRaw", SurEnergyRaw, &b_SurEnergyRaw);
    fChain->SetBranchAddress("SusMul", &SusMul, &b_SusMul);
    fChain->SetBranchAddress("SusChannel", SusChannel, &b_SusChannel);
    fChain->SetBranchAddress("SusEnergyRaw", SusEnergyRaw, &b_SusEnergyRaw);
    Notify();
}

Bool_t Reduce::Notify() {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void Reduce::Show(Long64_t entry) {
    // Print contents of entry.
    // If entry is not specified, print current entry
    if(!fChain) return;
    fChain->Show(entry);
}

Int_t Reduce::Cut(Long64_t entry) {
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

void Reduce::ReadCalibrations() {
    printf("Reading Sd1r Calibration File\n");
    std::ifstream inSd1r("Sd1r_AlphaCalibration.txt");
    std::string dummyLine;
    getline(inSd1r, dummyLine);
    int var1;
    double var2, var3;
    while(inSd1r >> var1 >> var2 >> var3) {
        sd1rCalibrations[var1] = std::make_pair(var2, var3);
    }
    inSd1r.close();

    printf("Reading Sd1s Calibration File\n");
    std::ifstream inSd1s("Sd1s_Calibration.txt");
    getline(inSd1s, dummyLine);
    while(inSd1s >> var1 >> var2 >> var3) {
        sd1sCalibrations[var1] = std::make_pair(var2, var3);
    }
    inSd1s.close();

    printf("Reading Sd2r Calibration File\n");
    std::ifstream inSd2r("Sd2r_Calibration.txt");
    getline(inSd2r, dummyLine);
    while(inSd2r >> var1 >> var2 >> var3) {
        sd2rCalibrations[var1] = std::make_pair(var2, var3);
    }
    inSd2r.close();

    printf("Reading Sd2s Calibration File\n");
    std::ifstream inSd2s("Sd2s_Calibration.txt");
    getline(inSd2s, dummyLine);
    while(inSd2s >> var1 >> var2 >> var3) {
        sd2sCalibrations[var1] = std::make_pair(var2, var3);
    }
    inSd2s.close();

    printf("Reading Yd Calibration File\n");
    std::ifstream inYd("Yd_Calibration.txt");
    getline(inYd, dummyLine);
    while(inYd >> var1 >> var2 >> var3) {
        ydCalibrations[var1] = std::make_pair(var2, var3);
    }
    inYd.close();

    printf("Reading Yu Calibration File\n");
    std::ifstream inYu("Yu_Calibration.txt");
    getline(inYu, dummyLine);
    while(inYu >> var1 >> var2 >> var3) {
        yuCalibrations[var1] = std::make_pair(var2, var3);
    }
    inYu.close();
}

void Reduce::InitTree() {
    outTree = new TTree("AutoTree", "12Be(d, p) events after IC cut and PID in Sd detectors");

    outTree->Branch("YuMul", &fYuMul, "YuMul/I");
    outTree->Branch("YuChannel", &fYuChannel, "YuChannel[YuMul]/I");
    outTree->Branch("YuRing", &fYuRing, "YuRing[YuMul]/I");
    outTree->Branch("YuSector", &fYuSector, "YuSector[YuMul]/I");
    outTree->Branch("YuEnergyRaw", &fYuEnergyRaw, "YuEnergyRaw[YuMul]/I");
    outTree->Branch("YuEnergy", &fYuEnergy, "YuEnergy[YuMul]/F");

    outTree->Branch("Sd1Mul", &fSd1Mul, "Sd1Mul/I");
    outTree->Branch("Sd1Ring", &fSd1Ring, "Sd1Ring[Sd1Mul]/I");
    outTree->Branch("Sd1Sector", &fSd1Sector, "Sd1Sector[Sd1Mul]/I");
    outTree->Branch("Sd1EnergyRawRing", &fSd1EnergyRawRing, "Sd1EnergyRawRing[Sd1Mul]/I");
    outTree->Branch("Sd1EnergyRawSector", &fSd1EnergyRawSector, "Sd1EnergyRawSector[Sd1Mul]/I");
    outTree->Branch("Sd1EnergyRing", &fSd1EnergyRing, "Sd1EnergyRing[Sd1Mul]/F");
    outTree->Branch("Sd1EnergySector", &fSd1EnergySector, "Sd1EnergySector[Sd1Mul]/F");

    outTree->Branch("Sd2Mul", &fSd2Mul, "Sd2Mul/I");
    outTree->Branch("Sd2Ring", &fSd2Ring, "Sd2Ring[Sd2Mul]/I");
    outTree->Branch("Sd2Sector", &fSd2Sector, "Sd2Sector[Sd2Mul]/I");
    outTree->Branch("Sd2EnergyRawRing", &fSd2EnergyRawRing, "Sd2EnergyRawRing[Sd2Mul]/I");
    outTree->Branch("Sd2EnergyRawSector", &fSd2EnergyRawSector, "Sd2EnergyRawSector[Sd2Mul]/I");
    outTree->Branch("Sd2EnergyRing", &fSd2EnergyRing, "Sd2EnergyRing[Sd2Mul]/F");
    outTree->Branch("Sd2EnergySector", &fSd2EnergySector, "Sd2EnergySector[Sd2Mul]/F");

    outTree->Branch("ICEnergyRaw", &ICEnergyRaw);
}

#endif // #ifdef Reduce_cxx
