//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  9 14:14:15 2019 by ROOT version 6.19/01
// from TTree AutoTree/12Be(d, p) events after IC cut and PID in Sd detectors
// found on file: reducedFile.root
//////////////////////////////////////////////////////////

#ifndef Analysis12Be_h
#define Analysis12Be_h

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

class Analysis12Be {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t           YuMul;
    Int_t           YuChannel[128];   //[YuMul]
    Int_t           YuRing[128];   //[YuMul]
    Int_t           YuSector[128];   //[YuMul]
    Int_t           YuEnergyRaw[128];   //[YuMul]
    Float_t         YuEnergy[128];   //[YuMul]
    Int_t           Sd1Mul;
    Int_t           Sd1Ring[128];   //[Sd1Mul]
    Int_t           Sd1Sector[128];   //[Sd1Mul]
    Int_t           Sd1EnergyRawRing[128];   //[Sd1Mul]
    Int_t           Sd1EnergyRawSector[128];   //[Sd1Mul]
    Float_t         Sd1EnergyRing[128];   //[Sd1Mul]
    Float_t         Sd1EnergySector[128];   //[Sd1Mul]
    Int_t           Sd2Mul;
    Int_t           Sd2Ring[128];   //[Sd2Mul]
    Int_t           Sd2Sector[128];   //[Sd2Mul]
    Int_t           Sd2EnergyRawRing[128];   //[Sd2Mul]
    Int_t           Sd2EnergyRawSector[128];   //[Sd2Mul]
    Float_t         Sd2EnergyRing[128];   //[Sd2Mul]
    Float_t         Sd2EnergySector[128];   //[Sd2Mul]
    Float_t         ICEnergyRaw;

    // List of branches
    TBranch        *b_YuMul;   //!
    TBranch        *b_YuChannel;   //!
    TBranch        *b_YuRing;   //!
    TBranch        *b_YuSector;   //!
    TBranch        *b_YuEnergyRaw;   //!
    TBranch        *b_YuEnergy;   //!
    TBranch        *b_Sd1Mul;   //!
    TBranch        *b_Sd1Ring;   //!
    TBranch        *b_Sd1Sector;   //!
    TBranch        *b_Sd1EnergyRawRing;   //!
    TBranch        *b_Sd1EnergyRawSector;   //!
    TBranch        *b_Sd1EnergyRing;   //!
    TBranch        *b_Sd1EnergySector;   //!
    TBranch        *b_Sd2Mul;   //!
    TBranch        *b_Sd2Ring;   //!
    TBranch        *b_Sd2Sector;   //!
    TBranch        *b_Sd2EnergyRawRing;   //!
    TBranch        *b_Sd2EnergyRawSector;   //!
    TBranch        *b_Sd2EnergyRing;   //!
    TBranch        *b_Sd2EnergySector;   //!
    TBranch        *b_ICEnergyRaw;   //!

    Analysis12Be(TTree *tree=0);
    virtual ~Analysis12Be();
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
};

#endif

#ifdef Analysis12Be_cxx
Analysis12Be::Analysis12Be(TTree *tree) : fChain(0) {
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if(tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("reducedFile.root");
        if(!f || !f->IsOpen()) {
            f = new TFile("reducedFile.root");
        }
        f->GetObject("AutoTree",tree);
    }
    Init(tree);
}

Analysis12Be::~Analysis12Be() {
   if(!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis12Be::GetEntry(Long64_t entry) {
    // Read contents of entry.
    if(!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t Analysis12Be::LoadTree(Long64_t entry) {
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

void Analysis12Be::Init(TTree *tree) {
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

    fChain->SetBranchAddress("YuMul", &YuMul, &b_YuMul);
    fChain->SetBranchAddress("YuChannel", YuChannel, &b_YuChannel);
    fChain->SetBranchAddress("YuRing", YuRing, &b_YuRing);
    fChain->SetBranchAddress("YuSector", YuSector, &b_YuSector);
    fChain->SetBranchAddress("YuEnergyRaw", YuEnergyRaw, &b_YuEnergyRaw);
    fChain->SetBranchAddress("YuEnergy", YuEnergy, &b_YuEnergy);
    fChain->SetBranchAddress("Sd1Mul", &Sd1Mul, &b_Sd1Mul);
    fChain->SetBranchAddress("Sd1Ring", Sd1Ring, &b_Sd1Ring);
    fChain->SetBranchAddress("Sd1Sector", Sd1Sector, &b_Sd1Sector);
    fChain->SetBranchAddress("Sd1EnergyRawRing", Sd1EnergyRawRing, &b_Sd1EnergyRawRing);
    fChain->SetBranchAddress("Sd1EnergyRawSector", Sd1EnergyRawSector, &b_Sd1EnergyRawSector);
    fChain->SetBranchAddress("Sd1EnergyRing", Sd1EnergyRing, &b_Sd1EnergyRing);
    fChain->SetBranchAddress("Sd1EnergySector", Sd1EnergySector, &b_Sd1EnergySector);
    fChain->SetBranchAddress("Sd2Mul", &Sd2Mul, &b_Sd2Mul);
    fChain->SetBranchAddress("Sd2Ring", Sd2Ring, &b_Sd2Ring);
    fChain->SetBranchAddress("Sd2Sector", Sd2Sector, &b_Sd2Sector);
    fChain->SetBranchAddress("Sd2EnergyRawRing", Sd2EnergyRawRing, &b_Sd2EnergyRawRing);
    fChain->SetBranchAddress("Sd2EnergyRawSector", Sd2EnergyRawSector, &b_Sd2EnergyRawSector);
    fChain->SetBranchAddress("Sd2EnergyRing", Sd2EnergyRing, &b_Sd2EnergyRing);
    fChain->SetBranchAddress("Sd2EnergySector", Sd2EnergySector, &b_Sd2EnergySector);
    fChain->SetBranchAddress("ICEnergyRaw", &ICEnergyRaw, &b_ICEnergyRaw);
    Notify();
}

Bool_t Analysis12Be::Notify() {
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void Analysis12Be::Show(Long64_t entry) {
    // Print contents of entry.
    // If entry is not specified, print current entry
    if(!fChain) return;
    fChain->Show(entry);
}

Int_t Analysis12Be::Cut(Long64_t entry) {
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}

void Analysis12Be::ReadCalibrations() {
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

#endif // #ifdef Analysis12Be_cxx
