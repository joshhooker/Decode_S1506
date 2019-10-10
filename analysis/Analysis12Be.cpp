#define Analysis12Be_cxx
#include "Analysis12Be.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();
    Analysis12Be t(chain);
    t.Loop();
    return 0;
}

TChain* MakeChain() {
    auto *chain = new TChain("AutoTree");

    TString PathToFiles = "/Users/joshhooker/Desktop/data/12Be_dp/";
    chain->Add(PathToFiles + "reducedFile12Be.root");

    return chain;
}

void Analysis12Be::Loop() {
    if (fChain == 0) return;

    TFile *file = new TFile("output.root", "recreate");

    Long64_t nentries = fChain->GetEntries();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        if(jentry != 0 && jentry % 100000 == 0) printf("Processed %lld/%lld events\n", jentry, nentries);

        std::vector<yyDetect> yuDetect_;
        for(int i = 0; i < YuMul; i++) {
            float energy = yuCalibrations[YuChannel[i]].first*YuEnergyRaw[i] + yuCalibrations[YuChannel[i]].second;
            yyDetect hit = {YuChannel[i], YuRing[i], YuSector[i], YuEnergyRaw[i], static_cast<float>(energy/1000.)};
            yuDetect_.push_back(hit);
        }

        std::vector<s3Detect> sd1Detect_;
        for(int i = 0; i < Sd1Mul; i++) { // Just use ring multiplicity and the first sector
            float energyRing = sd1rCalibrations[Sd1Ring[i]].first*Sd1EnergyRawRing[i] + sd1rCalibrations[Sd1Ring[i]].second;
            float energySector = sd1sCalibrations[Sd1Sector[i]].first*Sd1EnergyRawSector[i] + sd1sCalibrations[Sd1Sector[i]].second;
            s3Detect hit = {Sd1Ring[i], Sd1Sector[i], Sd1EnergyRawRing[i], Sd1EnergyRawSector[i], static_cast<float>(energyRing/1000.), static_cast<float>(energySector/1000.)};
            sd1Detect_.push_back(hit);
        }

        std::vector<s3Detect> sd2Detect_;
        for(int i = 0; i < Sd2Mul; i++) { // Just use ring multiplicity and the first sector
            float energyRing = sd2rCalibrations[Sd2Ring[i]].first*Sd2EnergyRawRing[i] + sd2rCalibrations[Sd2Ring[i]].second;
            float energySector = sd2sCalibrations[Sd2Sector[i]].first*Sd2EnergyRawSector[i] + sd2sCalibrations[Sd2Sector[i]].second;
            s3Detect hit = {Sd2Ring[i], Sd2Sector[i], Sd2EnergyRawRing[i], Sd2EnergyRawSector[i], static_cast<float>(energyRing/1000.), static_cast<float>(energySector/1000.)};
            sd2Detect_.push_back(hit);
        }

    }


    file->Close();
}
