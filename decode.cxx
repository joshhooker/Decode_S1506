using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "/home/reactions/treeIris/include/IDet.h" //the treeIris library dependancy - You need to change the path to IDet.h depending on your computer!
#include "iostream"
#include "fstream"

struct yyDetect {
    int channel;
    int energyRaw;
    float energy;
    int sector;
    int ring;
};

struct csiDetect {
    int channel;
    double energy;
};

struct sdDetect {
    int channel;
    int energyRaw;
    float energy;
};

struct suDetect {
    int channel;
    int energyRaw;
};

bool sortYy(const yyDetect& lhs, const yyDetect& rhs) {return lhs.energy > rhs.energy;};
bool sortCsI(const csiDetect& lhs, const csiDetect& rhs) {return lhs.energy > rhs.energy;};
bool sortSd(const sdDetect& lhs, const sdDetect& rhs) {return lhs.energy > rhs.energy;};
bool sortSu(const suDetect& lhs, const suDetect& rhs) {return lhs.energyRaw > rhs.energyRaw;};

void decode (int run_num) {

    string f_name = Form("/home/reactions/RawRootFiles/12Be/tree%i.root", run_num); //open the raw root file you want to convert to the format used by my codes; put the correct path and file name
    TFile* f_in = TFile::Open(f_name.c_str(), "READ");
    TTree* tr_in = (TTree*) f_in->Get("Iris"); //reading the Iris tree;
    //only the ADC leafs of the tree are used in the code. The rest of the data has a vector-array configuration or something like that and presents reliability which introduces unnecessary uncertanties

    IDet* det = new IDet();

    Int_t* TYdMul = &(det->TYdMul);
    std::vector<Int_t>* TYdChannel = &(det->TYdChannel);
    std::vector<Double_t>* TYdEnergy = &(det->TYdEnergy);
    std::vector<Int_t>* TYdADC = &(det->TYdADC);
    std::vector<Int_t>* TYdNo = &(det->TYdNo);
    std::vector<Int_t>* TYdRing = &(det->TYdRing);
    std::vector<Double_t>* TYdTheta = &(det->TYdTheta);

    Int_t* TCsI1Mul = &(det->TCsI1Mul);
    std::vector<Int_t>* TCsI1Channel = &(det->TCsI1Channel);
    std::vector<Double_t>* TCsI1Energy = &(det->TCsI1Energy);
    std::vector<Double_t>* TCsI1Phi = &(det->TCsI1Phi);
    std::vector<Double_t>* TCsI1ADC = &(det->TCsI1ADC);

    Int_t* TCsI2Mul = &(det->TCsI2Mul);
    std::vector<Int_t>* TCsI2Channel = &(det->TCsI2Channel);
    std::vector<Double_t>* TCsI2Energy = &(det->TCsI2Energy);
    std::vector<Double_t>* TCsI2Phi = &(det->TCsI2Phi);
    std::vector<Double_t>* TCsI2ADC = &(det->TCsI2ADC);

    Double_t* TYdCsI1ETot = &(det->TYdCsI1ETot);
    Double_t* TYdCsI2ETot = &(det->TYdCsI2ETot);

    Int_t* TSSBADC = &(det->TSSBADC);
    Double_t* TSSBEnergy = &(det->TSSBEnergy);

    Int_t* TScADC = &(det->TScADC);
    Double_t* TScEnergy = &(det->TScEnergy);

    std::vector<Int_t>* TTrADC = &(det->TTrADC);
    std::vector<Double_t>* TTrEnergy = &(det->TTrEnergy);

    std::vector<Int_t>* TICChannel = &(det->TICChannel);
    std::vector<Double_t>* TICEnergy = &(det->TICEnergy);
    std::vector<Double_t>* TICADC = &(det->TICADC);

    Int_t* TSd1rMul = &(det->TSd1rMul);
    std::vector<Int_t>* TSd1rChannel = &(det->TSd1rChannel);
    std::vector<Double_t>* TSd1rEnergy = &(det->TSd1rEnergy);
    std::vector<Int_t>* TSd1rADC = &(det->TSd1rADC);

    Int_t* TSd1sMul = &(det->TSd1sMul);
    std::vector<Int_t>* TSd1sChannel = &(det->TSd1sChannel);
    std::vector<Double_t>* TSd1sEnergy = &(det->TSd1sEnergy);
    std::vector<Int_t>* TSd1sADC = &(det->TSd1sADC);

    Int_t* TSd2rMul = &(det->TSd2rMul);
    std::vector<Int_t>* TSd2rChannel = &(det->TSd2rChannel);
    std::vector<Double_t>* TSd2rEnergy = &(det->TSd2rEnergy);
    std::vector<Int_t>* TSd2rADC = &(det->TSd2rADC);
    Double_t* TSd2rEnergyCal = &(det->TSd2rEnergyCal);

    Int_t* TSd2sMul = &(det->TSd2sMul);
    std::vector<Int_t>* TSd2sChannel = &(det->TSd2sChannel);
    std::vector<Double_t>* TSd2sEnergy = &(det->TSd2sEnergy);
    std::vector<Int_t>* TSd2sADC = &(det->TSd2sADC);

    Double_t* TSdETot = &(det->TSdETot);
    std::vector<Double_t>* TSd1Theta = &(det->TSd1Theta);
    std::vector<Double_t>* TSd2Theta = &(det->TSd2Theta);
    Double_t* TSdThetaCM = &(det->TSdThetaCM);
    std::vector<Double_t>* TSd1Phi = &(det->TSd1Phi);
    std::vector<Double_t>* TSd2Phi = &(det->TSd2Phi);

    Int_t* TYuMul = &(det->TYuMul);
    std::vector<Int_t>* TYuChannel = &(det->TYuChannel);
    std::vector<Double_t>* TYuEnergy = &(det->TYuEnergy);
    std::vector<Int_t>* TYuADC = &(det->TYuADC);
    std::vector<Int_t>* TYuNo = &(det->TYuNo);
    std::vector<Int_t>* TYuRing = &(det->TYuRing);
    std::vector<Double_t>* TYuTheta = &(det->TYuTheta); // Yd theta angle

    Int_t* TSurMul = &(det->TSurMul);
    std::vector<Int_t>* TSurChannel = &(det->TSurChannel);
    std::vector<Double_t>* TSurEnergy = &(det->TSurEnergy);
    std::vector<Int_t>* TSurADC = &(det->TSurADC);
    Double_t* TSurEnergyCal = &(det->TSurEnergyCal);

    Int_t* TSusMul = &(det->TSusMul);
    std::vector<Int_t>* TSusChannel = &(det->TSusChannel);
    std::vector<Double_t>* TSusEnergy = &(det->TSusEnergy);
    std::vector<Int_t>* TSusADC = &(det->TSusADC);
    std::vector<Double_t>* TSuTheta = &(det->TSuTheta);
    std::vector<Double_t>* TSuPhi = &(det->TSuPhi);

    tr_in->SetBranchAddress("det", &det);

    // Output variables
    // YYD detector
    int YdMul;
    int YdChannel[128];
    int YdEnergyRaw[128];
    float YdEnergy[128];
    int YdRing[128];
    int YdSector[128];

    // YYU detector
    int YuMul;
    int YuChannel[128];
    int YuEnergyRaw[128];
    float YuEnergy[128];
    int YuRing[128];
    int YuSector[128];

    // Purpose of two readouts is to increase the range of operation (multiply it by 2 in total)
    // CsI 1 detector
    int CsI1Mul;
    int CsI1Channel[128];
    float CsI1EnergyRaw[128];

    // CsI 2 detector
    int CsI2Mul;
    int CsI2Channel[128]; // Channel corresponds to one christal - 2 means that the gain is set to a different value than 1
    float CsI2EnergyRaw[128];

    // IC chamber
    int ICChannel;
    float ICEnergyRaw;

    // Sd1r detector (Sd1 Rings)
    int Sd1rMul;
    int Sd1rChannel[128];
    int Sd1rEnergyRaw[128];
    float Sd1rEnergy[128];

    // Sd1s detector (Sd1 Sectors)
    int Sd1sMul;
    int Sd1sChannel[128];
    int Sd1sEnergyRaw[128];
    float Sd1sEnergy[128];

    // Sd2r detector (Sd2 Rings)
    int Sd2rMul;
    int Sd2rChannel[128];
    int Sd2rEnergyRaw[128];
    float Sd2rEnergy[128];

    // Sd2s detector (Sd2 Sectors)
    int Sd2sMul;
    int Sd2sChannel[128];
    int Sd2sEnergyRaw[128];
    float Sd2sEnergy[128];

    // Sur detector (Su Rings)
    int SurMul;
    int SurChannel[128];
    int SurEnergyRaw[128];

    // Sus detector (Su Sectors)
    int SusMul;
    int SusChannel[128];
    int SusEnergyRaw[128];

    ////////////////////////////
    // Yu Sector and Ring Map //
    ////////////////////////////
    std::map<int, std::pair<int, int> > YuSectorRingMap;
    for(int i = 0; i < 128; i++) {
        int sector = i/16;
        int ring = i%16;
        YuSectorRingMap[i] = std::make_pair(sector, ring);
    }

    //////////////////////////
    // Open the output file //
    //////////////////////////

    string fOut_name = Form("/home/reactions/Decoded_root/Be/NewDecode_%i.root", run_num); //Open the output file; put the correct path and file name
    TFile* f_out = new TFile(fOut_name.c_str(), "RECREATE");
    TTree* tr_out = new TTree("AutoTree", "AutoTree"); //create the new root tree

    //Branches created in the new file; the leafs are removed from the three as well as arrays of vectors or whatever that was
    tr_out->Branch("YdMul", &YdMul, "YdMul/I");
    tr_out->Branch("YdChannel", &YdChannel, "YdChannel[YdMul]/I");
    tr_out->Branch("YdEnergyRaw", &YdEnergyRaw, "YdEnergyRaw[YdMul]/I");
    tr_out->Branch("YdEnergy", &YdEnergy, "YdEnergy[YdMul]/F");
    tr_out->Branch("YdRing", &YdRing, "YdRing[YdMul]/I");
    tr_out->Branch("YdSector", &YdSector, "YdSection[YdMul]/I");

    tr_out->Branch("YuMul", &YuMul, "YuMul/I");
    tr_out->Branch("YuChannel", &YuChannel, "YuChannel[YuMul]/I");
    tr_out->Branch("YuEnergyRaw", &YuEnergyRaw, "YuEnergyRaw[YuMul]/I");
    tr_out->Branch("YuEnergy", &YuEnergy, "YuEnergy[YuMul]/F");
    tr_out->Branch("YuRing", &YuRing, "YuRing[YuMul]/I");
    tr_out->Branch("YuSector", &YuSector, "YuSector[YuMul]/I");

    tr_out->Branch("CsI1Mul", &CsI1Mul, "CsI1Mul/I");
    tr_out->Branch("CsI1Channel", &CsI1Channel, "CsI1Channel[CsI1Mul]/I");
    tr_out->Branch("CsI1EnergyRaw", &CsI1EnergyRaw, "CsI1EnergyRaw[CsI1Mul]/F");

    tr_out->Branch("CsI2Mul", &CsI2Mul, "CsI2Mul/I");
    tr_out->Branch("CsI2Channel", &CsI2Channel, "CsI2Channel[CsI2Mul]/I");
    tr_out->Branch("CsI2EnergyRaw", &CsI2EnergyRaw, "CsI2EnergyRaw[CsI2Mul]/F");

    tr_out->Branch("ICChannel", &ICChannel, "ICChannel/I");
    tr_out->Branch("ICEnergyRaw",&ICEnergyRaw, "ICEnergyRaw/F");

    tr_out->Branch("Sd1rMul", &Sd1rMul, "Sd1rMul/I");
    tr_out->Branch("Sd1rChannel", &Sd1rChannel, "Sd1rChannel[Sd1rMul]/I");
    tr_out->Branch("Sd1rEnergyRaw", &Sd1rEnergyRaw, "Sd1rEnergyRaw[Sd1rMul]/I");
    tr_out->Branch("Sd1rEnergy", &Sd1rEnergy, "Sd1rEnergy[Sd1rMul]/F");

    tr_out->Branch("Sd1sMul", &Sd1sMul, "Sd1sMul/I");
    tr_out->Branch("Sd1sChannel", &Sd1sChannel, "Sd1sChannel[Sd1sMul]/I");
    tr_out->Branch("Sd1sEnergyRaw", &Sd1sEnergyRaw, "Sd1sEnergyRaw[Sd1sMul]/I");
    tr_out->Branch("Sd1sEnergy", &Sd1sEnergy, "Sd1sEnergy[Sd1sMul]/F");

    tr_out->Branch("Sd2rMul", &Sd2rMul, "Sd2rMul/I");
    tr_out->Branch("Sd2rChannel", &Sd2rChannel, "Sd2rChannel[Sd2rMul]/I");
    tr_out->Branch("Sd2rEnergyRaw", &Sd2rEnergyRaw, "Sd2rEnergyRaw[Sd2rMul]/I");
    tr_out->Branch("Sd2rEnergy", &Sd2rEnergy, "Sd2rEnergy[Sd2rMul]/F");

    tr_out->Branch("Sd2sMul", &Sd2sMul, "Sd2sMul/I");
    tr_out->Branch("Sd2sChannel", &Sd2sChannel, "Sd2sChannel[Sd2sMul]/I");
    tr_out->Branch("Sd2sEnergyRaw", &Sd2sEnergyRaw, "Sd2sEnergyRaw[Sd2sMul]/I");
    tr_out->Branch("Sd2sEnergy", &Sd2sEnergy, "Sd2sEnergy[Sd2sMul]/F");

    tr_out->Branch("SurMul", &SurMul, "SurMul/I");
    tr_out->Branch("SurChannel", &SurChannel, "SurChannel[SurMul]/I");
    tr_out->Branch("SurEnergyRaw", &SurEnergyRaw, "SurEnergyRaw[SurMul]/I");

    tr_out->Branch("SusMul", &SusMul, "SusMul/I");
    tr_out->Branch("SusChannel", &SusChannel, "SusChannel[SusMul]/I");
    tr_out->Branch("SusEnergyRaw", &SusEnergyRaw, "SusEnergyRaw[SusMul]/I");


    ///////////////////////////////////
    // Open and set the calibrations //
    ///////////////////////////////////

    //Opening the Yd calibration file, reading it and assigning the relevant variables to the internaly defined variables
    double stripYd[128], Ydp0[128],Ydp1[128]; //definition of gain and offset used for the Yd calibration

    ifstream Alpha_peaksYd;
    Alpha_peaksYd.open("Yd_4815Ped.txt"); //Open the calibration file; put the correct path and file name
    if(Alpha_peaksYd.is_open()) {
        Alpha_peaksYd.ignore(256, '\n'); //Skips the first line in the .txt file to go and read the acctual parameters
        for(int i = 0; i < 128; i++) {
            Alpha_peaksYd >> stripYd[i] >> Ydp0[i] >> Ydp1[i] ;
            // cout << " Strip " <<  stripYd[i] << " P0 " << Ydp0[i] << " p1 " << Ydp1[i] << endl;
        }
    } else {
        cout << "No calib file for Yd ditector "<< endl;
    }
    Alpha_peaksYd.close();

    //Opening and reading the Yu calibration file, same as for Yd
    double stripYu[128], Yup0[128],Yup1[128], Yup2[128];

    ifstream Alpha_peaksYu;
    Alpha_peaksYu.open("Yu5225_4815Ped.txt"); //put the correct path and file name
    if(Alpha_peaksYu.is_open()) {
        Alpha_peaksYu.ignore(256, '\n');
        for(int i = 0; i < 128; i++) {
            Alpha_peaksYu >> stripYu[i] >> Yup0[i] >> Yup1[i] ;
            // cout << " Strip " <<  stripYu[i] << " P0 " << Yup0[i] << " p1 " << Yup1[i] << endl;
        }
    } else {
        cout << " No calib file for Yu detector "<< endl;
    }
    Alpha_peaksYu.close();

    //Sd1 detector
    double stripSd1r[24], stripSd1s[32], Sd1rp0[24], Sd1rp1[24], Sd1sp0[32], Sd1sp1[32];
    double stripSd2r[24], stripSd2s[32], Sd2rp0[24], Sd2rp1[24], Sd2sp0[32], Sd2sp1[32];

    ifstream Alpha_Sd1r, Alpha_Sd1s;
    ifstream Alpha_Sd2r, Alpha_Sd2s;

    // Opening and reading the calibration files for the Sd1 and Sd2 detectors
    Alpha_Sd1r.open("Sd1r_4815Ped.txt"); // put the correct path and file name
    if(Alpha_Sd1r.is_open()) {
        Alpha_Sd1r.ignore(256, '\n');
        for(int i = 0; i < 24; i++) {
            Alpha_Sd1r >> stripSd1r[i] >> Sd1rp0[i] >> Sd1rp1[i] ;
        }
    } else {cout << "No calib file for Sd1r "<< endl;}

    Alpha_Sd1s.open("Sd1s_4815Ped.txt"); // put the correct path and file name
    if(Alpha_Sd1s.is_open()) {
        Alpha_Sd1s.ignore(256, '\n');
        for(int i = 0; i < 32; i++) {
            Alpha_Sd1s >> stripSd1s[i] >> Sd1sp0[i] >> Sd1sp1[i] ;
        }
    } else {cout << "No calib file for Sd1s "<< endl;}

    Alpha_Sd2r.open("Sd2r_4815Ped.txt"); // put the correct path and file name
    if(Alpha_Sd2r.is_open()) {
        Alpha_Sd2r.ignore(256, '\n');
        for(int i = 0; i < 24; i++) {
            Alpha_Sd2r >> stripSd2r[i] >> Sd2rp0[i] >> Sd2rp1[i] ;
        }
    } else {cout << "No calib file for Sd2r "<< endl;}

    Alpha_Sd2s.open("Sd2s_4815Ped.txt"); // put the correct path and file name
    if(Alpha_Sd2s.is_open()) {
        Alpha_Sd2s.ignore(256, '\n');
        for(int i = 0; i < 32; i++) {
            Alpha_Sd2s >> stripSd2s[i] >> Sd2sp0[i] >> Sd2sp1[i] ;
        }
    } else {cout << "No calib file for Sd2s "<< endl;}
    // end of calibration files


    /////////////////////////
    // Loop through events //
    /////////////////////////

    long totalEvents = tr_in->GetEntries(); // get the total number of events in the tree
    cout << "Total number of events: " << totalEvents << endl;

    for(long ev = 0; ev < totalEvents; ev++) {
        if(ev % 10000 == 0) cout << "Current event: " << ev << endl;
        tr_in->GetEntry(ev); // pulling up the current event from the Iris tree

        YdMul = 0;
        YuMul = 0;
        CsI1Mul = 0;
        CsI2Mul = 0;
        Sd1rMul = 0;
        Sd1sMul = 0;
        Sd2rMul = 0;
        Sd2sMul = 0;
        SurMul = 0;
        SusMul = 0;
        for(int i = 0; i < 128; i++) {
            YdChannel[i] = 0;
            YdEnergyRaw[i] = 0;
            YdEnergy[i] = 0;
            YdRing[i] = 0;
            YdSector[i] = 0;
            YuChannel[i] = 0;
            YuEnergyRaw[i] = 0;
            YuEnergy[i] = 0;
            YuRing[i] = 0;
            YuSector[i] = 0;
            CsI1Channel[i] = 0;
            CsI1EnergyRaw[i] = 0;
            CsI2Channel[i] = 0;
            CsI2EnergyRaw[i] = 0;
            Sd1rChannel[i] = 0;
            Sd1rEnergyRaw[i] = 0;
            Sd1rEnergy[i] = 0;
            Sd1sChannel[i] = 0;
            Sd1sEnergyRaw[i] = 0;
            Sd1sEnergy[i] = 0;
            Sd2rChannel[i] = 0;
            Sd2rEnergyRaw[i] = 0;
            Sd2rEnergy[i] = 0;
            Sd2sChannel[i] = 0;
            Sd2sEnergyRaw[i] = 0;
            Sd2sEnergy[i] = 0;
            SurChannel[i] = 0;
            SurEnergyRaw[i] = 0;
            SusChannel[i] = 0;
            SusEnergyRaw[i] = 0;
        }

        ICChannel = -1;
        ICEnergyRaw = -1;


        /////////////////
        // Yd Detector //
        /////////////////
        std::vector<yyDetect> ydDetect_;
        for(size_t i = 0; i < det->TYdADC.size(); i++) { // Loop over the size of the YdADC vector, this is always 128
            int energyRaw = det->TYdADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Ydp0[i] + Ydp1[i];
                std::pair<int, int> region = YuSectorRingMap[i];
                yyDetect hit = {static_cast<int>(i), energyRaw, energy, region.first, region.second};
                ydDetect_.push_back(hit);
            }
        }


        /////////////////
        // Yu Detector //
        /////////////////
        std::vector<yyDetect> yuDetect_;
        for(size_t i = 0; i < det->TYuADC.size(); i++) { // Loop over the size of the YuADC vector, this is always 128
            int energyRaw = det->TYuADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Yup0[i] + Yup1[i];
                std::pair<int, int> region = YuSectorRingMap[i];
                yyDetect hit = {static_cast<int>(i), energyRaw, energy, region.first, region.second};
                yuDetect_.push_back(hit);
            }
        }


        ///////////////////
        // CsI Detectors //
        ///////////////////
        std::vector<csiDetect> csi1Detect_;
        for(size_t i = 0; i < det->TCsI1ADC.size(); i++) {
            float energyRaw = det->TCsI1ADC.at(i);
            if(energyRaw > 0) {
                csiDetect hit = {static_cast<int>(i), energyRaw};
                csi1Detect_.push_back(hit);
            }
        }
        std::vector<csiDetect> csi2Detect_;
        for(size_t i = 0; i < det->TCsI2ADC.size(); i++) {
            float energyRaw = det->TCsI2ADC.at(i);
            if(energyRaw > 0) {
                csiDetect hit = {static_cast<int>(i), energyRaw};
                csi2Detect_.push_back(hit);
            }
        }


        ////////
        // IC //
        ////////
        if(TICChannel->size() > 0) {
            ICChannel = 15;
            ICEnergyRaw = TICADC->at(15);
        }


        //////////////////
        // Sd Detectors //
        //////////////////

        // Treating the Sd1r part of the Sd1 detector. There are 2 Sd detectors. Sd1 is the one making the dE layer in the dE-E telechope. It has 24 rings and 32 sectors. The rings are oriented toward the beam.
        // The data is gathered in a way that the sectors and rings are registered separately which gives the rings as Sd1r. The Sd2 detector is treated the same as Sd1.

        //////////
        // Sd1r //
        //////////
        std::vector<sdDetect> sd1rDetect_;
        for(size_t i = 0; i < det->TSd1rADC.size(); i++) {
            int energyRaw = det->TSd1rADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Sd1rp0[i] + Sd1rp1[i];
                sdDetect hit = {static_cast<int>(i), energyRaw, energy};
                sd1rDetect_.push_back(hit);
            }
        }

        //////////
        // Sd1s //
        //////////
        std::vector<sdDetect> sd1sDetect_;
        for(size_t i = 0; i < det->TSd1sADC.size(); i++) {
            int energyRaw = det->TSd1sADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Sd1sp0[i] + Sd1sp1[i];
                sdDetect hit = {static_cast<int>(i), energyRaw, energy};
                sd1sDetect_.push_back(hit);
            }
        }

        //////////
        // Sd2r //
        //////////
        std::vector<sdDetect> sd2rDetect_;
        for(size_t i = 0; i < det->TSd2rADC.size(); i++) {
            int energyRaw = det->TSd2rADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Sd2rp0[i] + Sd2rp1[i];
                sdDetect hit = {static_cast<int>(i), energyRaw, energy};
                sd2rDetect_.push_back(hit);
            }
        }

        //////////
        // Sd2s //
        //////////
        std::vector<sdDetect> sd2sDetect_;
        for(size_t i = 0; i < det->TSd2sADC.size(); i++) {
            int energyRaw = det->TSd2sADC.at(i);
            if(energyRaw > 0) {
                float energy = static_cast<float>(energyRaw)*Sd2sp0[i] + Sd2sp1[i];
                sdDetect hit = {static_cast<int>(i), energyRaw, energy};
                sd2sDetect_.push_back(hit);
            }
        }


        //////////////////
        // Su Detectors //
        //////////////////

        // Another S3 detector positioned upstream in front of the Yu detector
        // The ring side of this detector is facing the target and its NOT facing the beam

        /////////
        // Sur //
        /////////
        std::vector<suDetect> surDetect_;
        for(size_t i = 0; i < det->TSurADC.size(); i++) {
            int energyRaw = det->TSurADC.at(i);
            if(energyRaw > 0) {
                suDetect hit = {static_cast<int>(i), energyRaw};
                surDetect_.push_back(hit);
            }
        }

        /////////
        // Sus //
        /////////
        std::vector<suDetect> susDetect_;
        for(size_t i = 0; i < det->TSusADC.size(); i++) {
            int energyRaw = det->TSusADC.at(i);
            if(energyRaw > 0) {
                suDetect hit = {static_cast<int>(i), energyRaw};
                susDetect_.push_back(hit);
            }
        }

        //////////////////////////////////////////////////////////////////
        // Sort vectors so first entry is highest energy (or energyRaw) //
        //////////////////////////////////////////////////////////////////

        std::sort(ydDetect_.begin(), ydDetect_.end(), sortYy);
        std::sort(yuDetect_.begin(), yuDetect_.end(), sortYy);
        std::sort(csi1Detect_.begin(), csi1Detect_.end(), sortCsI);
        std::sort(csi2Detect_.begin(), csi2Detect_.end(), sortCsI);
        std::sort(sd1rDetect_.begin(), sd1rDetect_.end(), sortSd);
        std::sort(sd1sDetect_.begin(), sd1sDetect_.end(), sortSd);
        std::sort(sd2rDetect_.begin(), sd2rDetect_.end(), sortSd);
        std::sort(sd2sDetect_.begin(), sd2sDetect_.end(), sortSd);
        std::sort(surDetect_.begin(), surDetect_.end(), sortSu);
        std::sort(susDetect_.begin(), susDetect_.end(), sortSu);


        //////////////////////////
        // Fill arrays for tree //
        //////////////////////////

        YdMul = 0;
        for(size_t i = 0; i < ydDetect_.size(); i++) {
            YdChannel[YdMul] = ydDetect_[i].channel;
            YdEnergyRaw[YdMul] = ydDetect_[i].energyRaw;
            YdEnergy[YdMul] = ydDetect_[i].energy/1000.;
            YdRing[YdMul] = ydDetect_[i].ring;
            YdSector[YdMul] = ydDetect_[i].sector;
            YdMul++;
        }

        YuMul = 0;
        for(size_t i = 0; i < yuDetect_.size(); i++) {
            YuChannel[YuMul] = yuDetect_[i].channel;
            YuEnergyRaw[YuMul] = yuDetect_[i].energyRaw;
            YuEnergy[YuMul] = yuDetect_[i].energy/1000.;
            YuRing[YuMul] = yuDetect_[i].ring;
            YuSector[YuMul] = yuDetect_[i].sector;
            YuMul++;
        }

        CsI1Mul = 0;
        for(auto csi: csi1Detect_) {
            CsI1Channel[CsI1Mul] = csi.channel;
            CsI1EnergyRaw[CsI1Mul] = csi.energy;
            CsI1Mul++;
        }

        CsI2Mul = 0;
        for(auto csi: csi1Detect_) {
            CsI2Channel[CsI2Mul] = csi.channel;
            CsI2EnergyRaw[CsI2Mul] = csi.energy;
            CsI2Mul++;
        }

        Sd1rMul = 0;
        for(auto sd: sd1rDetect_) {
            Sd1rChannel[Sd1rMul] = sd.channel;
            Sd1rEnergyRaw[Sd1rMul] = sd.energyRaw;
            Sd1rEnergy[Sd1rMul] = sd.energy/1000.;
            Sd1rMul++;
        }

        Sd1sMul = 0;
        for(auto sd: sd1sDetect_) {
            Sd1sChannel[Sd1sMul] = sd.channel;
            Sd1sEnergyRaw[Sd1sMul] = sd.energyRaw;
            Sd1sEnergy[Sd1sMul] = sd.energy/1000.;
            Sd1sMul++;
        }

        Sd2rMul = 0;
        for(auto sd: sd2rDetect_) {
            Sd2rChannel[Sd2rMul] = sd.channel;
            Sd2rEnergyRaw[Sd2rMul] = sd.energyRaw;
            Sd2rEnergy[Sd2rMul] = sd.energy/1000.;
            Sd2rMul++;
        }

        Sd2sMul = 0;
        for(auto sd: sd2sDetect_) {
            Sd2sChannel[Sd2sMul] = sd.channel;
            Sd2sEnergyRaw[Sd2sMul] = sd.energyRaw;
            Sd2sEnergy[Sd2sMul] = sd.energy/1000.;
            Sd2sMul++;
        }

        SurMul = 0;
        for(auto su: surDetect_) {
            SurChannel[SurMul] = su.channel;
            SurEnergyRaw[SurMul] = su.energyRaw;
            SurMul++;
        }

        SusMul = 0;
        for(auto su: susDetect_) {
            SusChannel[SusMul] = su.channel;
            SusEnergyRaw[SusMul] = su.energyRaw;
            SusMul++;
        }


        tr_out->Fill();
    } // End of the main while loop

    f_out->cd(); // Access the output file
    tr_out->Write(); // Write the output file
    f_out->Close(); // Close the output file

    return;
}
