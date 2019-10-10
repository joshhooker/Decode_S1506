#define Reduce_cxx
#include "Reduce12Be.h"

TChain* MakeChain();

int main() {
    TChain *chain = MakeChain();
    Reduce t(chain);
    t.Loop();
    return 0;
}

TChain* MakeChain() {
    auto *chain = new TChain("AutoTree");

    TString PathToFiles = "/Users/joshhooker/Desktop/data/12Be_dp/NewDecode_";

    std::vector<int> run_{5021, 5022, 5023, 5024, 5025, 5027, 5028, 5029, 5030,
                          5031, 5032, 5033, 5034, 5035, 5036, 5037, 5038, 5039,
                          5041, 5042, 5044, 5045, 5048, 5049, 5050, 5051, 5052,
                          5053, 5054, 5055, 5056, 5057, 5058, 5060, 5061, 5064,
                          5065, 5066, 5067, 5068, 5069, 5070, 5071, 5072, 5073,
                          5074, 5075, 5076, 5077, 5078, 5079, 5080, 5081, 5082,
                          5083, 5084, 5085, 5086, 5087, 5088, 5089, 5090, 5091,
                          5092, 5094, 5095, 5096, 5097, 5098, 5100, 5102, 5103,
                          5104, 5105, 5106, 5107, 5108, 5109, 5110, 5111, 5112,
                          5115, 5116, 5117, 5118, 5119, 5120, 5121, 5122, 5123,
                          5124, 5125, 5126, 5127, 5128, 5129, 5130, 5131, 5132,
                          5180, 5181, 5182, 5183, 5184, 5185, 5186, 5188, 5189,
                          5190, 5191, 5192, 5193, 5194, 5195, 5196, 5197, 5198,
                          5199, 5201, 5202, 5203, 5204, 5205, 5206, 5207, 5208,
                          5209, 5210, 5211, 5212, 5213, 5214, 5215, 5216, 5217,
                          5218, 5219, 5220, 5221};

    for(size_t i = 0; i < run_.size(); i++) {
        chain->Add(PathToFiles + Form("%d.root", run_[i]));
    }

    return chain;
}

void Reduce::Loop() {
    if(fChain == 0) return;

    TFile *f_cut = TFile::Open("BeCutIC.root"); //PID cut for Be
    TCutG *pidcut = (TCutG*) f_cut->Get("BeCut_NewCalSd1rSd2r");
    f_cut->Close();

    TFile *file = new TFile("reducedFile12Be.root", "recreate");

    ReadCalibrations();
    InitTree();

    TH2F* hPID = new TH2F("sdPID", "siPID; sd2 Energy [MeV]; sd1 Energy [MeV]", 500, 0, 130, 500, 0, 40);
    TH2F* hPIDPostCut = new TH2F("sdPIDPostCut", "siPIDPostCut; sd2 Energy [MeV]; sd1 Energy [MeV]", 500, 0, 130, 500, 0, 40);

    Long64_t nentries = fChain->GetEntries();

    Long64_t nbytes = 0, nb = 0;
    for(Long64_t jentry = 0; jentry < nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if(ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        if(jentry != 0 && jentry % 100000 == 0) printf("Processed %lld/%lld events\n", jentry, nentries);

        if(Sd1rMul == 0 || Sd1sMul == 0 || Sd2rMul == 0 || Sd2sMul == 0) continue;
        if(YuMul == 0) continue;
        if(ICEnergyRaw < 620 || ICEnergyRaw > 1100) continue;

        // Add detectors to structures and recalibrate so we don't have to run the decode again
        std::vector<yyDetect> ydDetect_;
        for(int i = 0; i < YdMul; i++) {
            float energy = ydCalibrations[YdChannel[i]].first*YdEnergyRaw[i] + ydCalibrations[YdChannel[i]].second;
            yyDetect hit = {YdChannel[i], YdRing[i], YdSector[i], YdEnergyRaw[i], static_cast<float>(energy)};
            ydDetect_.push_back(hit);
        }

        std::vector<yyDetect> yuDetect_;
        for(int i = 0; i < YuMul; i++) {
            float energy = yuCalibrations[YuChannel[i]].first*YuEnergyRaw[i] + yuCalibrations[YuChannel[i]].second;
            yyDetect hit = {YuChannel[i], YuRing[i], YuSector[i], YuEnergyRaw[i], static_cast<float>(energy/1000.)};
            yuDetect_.push_back(hit);
        }

        std::vector<s3Detect> sd1Detect_;
        for(int i = 0; i < Sd1rMul; i++) { // Just use ring multiplicity and the first sector
            float energyRing = sd1rCalibrations[Sd1rChannel[i]].first*Sd1rEnergyRaw[i] + sd1rCalibrations[Sd1rChannel[i]].second;
            float energySector = sd1sCalibrations[Sd1sChannel[0]].first*Sd1sEnergyRaw[0] + sd1sCalibrations[Sd1sChannel[0]].second;
            s3Detect hit = {Sd1rChannel[i], Sd1sChannel[0], Sd1rEnergyRaw[i], Sd1sEnergyRaw[0], static_cast<float>(energyRing/1000.), static_cast<float>(energySector/1000.)};
            sd1Detect_.push_back(hit);
        }

        std::vector<s3Detect> sd2Detect_;
        for(int i = 0; i < Sd2rMul; i++) { // Just use ring multiplicity and the first sector
            float energyRing = sd2rCalibrations[Sd2rChannel[i]].first*Sd2rEnergyRaw[i] + sd2rCalibrations[Sd2rChannel[i]].second;
            float energySector = sd2sCalibrations[Sd2sChannel[0]].first*Sd2sEnergyRaw[0] + sd2sCalibrations[Sd2sChannel[0]].second;
            s3Detect hit = {Sd2rChannel[i], Sd2sChannel[0], Sd2rEnergyRaw[i], Sd2sEnergyRaw[0], static_cast<float>(energyRing/1000.), static_cast<float>(energySector/1000.)};
            sd2Detect_.push_back(hit);
        }


        hPID->Fill(sd2Detect_[0].energyRing, sd1Detect_[0].energyRing);
        if(!pidcut->IsInside(sd2Detect_[0].energyRing, sd1Detect_[0].energyRing)) continue;
        hPIDPostCut->Fill(sd2Detect_[0].energyRing, sd1Detect_[0].energyRing);

        fYuMul = 0;
        for(auto yy: yuDetect_) {
            fYuChannel[fYuMul] = yy.channel;
            fYuRing[fYuMul] = yy.ring;
            fYuSector[fYuMul] = yy.sector;
            fYuEnergyRaw[fYuMul] = yy.energyRaw;
            fYuEnergy[fYuMul] = yy.energy;
            fYuMul++;
        }

        fSd1Mul = 0;
        for(auto s3: sd1Detect_) {
            fSd1Ring[fSd1Mul] = s3.ring;
            fSd1Sector[fSd1Mul] = s3.sector;
            fSd1EnergyRawRing[fSd1Mul] = s3.energyRawRing;
            fSd1EnergyRawSector[fSd1Mul] = s3.energyRawSector;
            fSd1EnergyRing[fSd1Mul] = s3.energyRing;
            fSd1EnergySector[fSd1Mul] = s3.energySector;
            fSd1Mul++;
        }

        fSd2Mul = 0;
        for(auto s3: sd2Detect_) {
            fSd2Ring[fSd2Mul] = s3.ring;
            fSd2Sector[fSd2Mul] = s3.sector;
            fSd2EnergyRawRing[fSd2Mul] = s3.energyRawRing;
            fSd2EnergyRawSector[fSd2Mul] = s3.energyRawSector;
            fSd2EnergyRing[fSd2Mul] = s3.energyRing;
            fSd2EnergySector[fSd2Mul] = s3.energySector;
            fSd2Mul++;
        }

        outTree->Fill();

    }

    ///////////////////////////////////////
    // Done with event by event analysis //
    ///////////////////////////////////////

    outTree->Write();
    hPID->Write();
    hPIDPostCut->Write();

    file->Close();
}