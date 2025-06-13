#include <TH1D.h>
#include <cmath>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <set>
#include <iostream>

//global definitions
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Double_t omega_tel = 0.040 ;//sr 


Float_t histIntegralError(TH1D * histo);

Float_t totalCharge(vector<Int_t> runs, bool verbose = true, string infofile =  "/mnt/medley/LucasAnalysis/2023/runlist.csv");

std::vector<std::vector<TH1D*>> GetDDX(
    float Ea,
    float Eb,
    const std::vector<float>& angles,
    const std::vector<char>& particles,
    const std::vector<int>& runsForward,
    const std::vector<int>& runsBackward
)




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to calculate the error in the integral of a histogram
Float_t histIntegralError(TH1D * histo){
    //Function to calculate the error in the integral of a histogram
    Float_t sum = 0;
    Int_t nbins = histo->GetNbinsX();
    for (Int_t i = 1; i <= nbins; i++) {
        if (histo->GetBinContent(i) > 0) {
            Float_t error = histo->GetBinError(i);
            sum += error * error * histo->GetBinWidth(i) * histo->GetBinWidth(i);
        }
    }

    return sqrt(sum);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to calculate the total Charge of a selection of runs 
//returns the total charge in µC
//  -> 'runs' is the vector of runs to be summed
//  -> 'verbose' is a boolean to print the runs and their charges
//  -> 'infofile' is the path to the runlist file containing the charge information
Float_t totalCharge(vector<Int_t> runs, bool verbose = true, string infofile =  "/mnt/medley/LucasAnalysis/2023/runlist.csv") {
    
    Float_t TotalCharge = 0;
    TTree *InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(infofile.c_str(), "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F");

    if(verbose) {
        cout << "---------------------------------------------------" << endl;
        cout << "Total Charge Calculation" << endl;
        cout << "---------------------------------------------------" << endl;
        cout << "Info file: " << infofile << endl;
        cout << "Total Charge will be calculated in µC." << endl;
        cout << "Number of runs: " << runs.size() << endl;
        if (runs.empty()) {
            cout << "No runs provided. Returning 0." << endl;
            return 0;
        }
        cout << "Calculating total charge for runs: ";  
        for (const auto& run : runs) {
            cout << run << " ";
        }
        cout << endl;
    }

    Int_t RunN;
    Float_t ChargeFaraday;

    InfoTree->SetBranchAddress("RunN", &RunN);
    InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

 // Coloque os runs em um set para busca rápida
    std::set<Int_t> runSet(runs.begin(), runs.end());

    int nEntries = InfoTree->GetEntries();
    for(int i = 0; i < nEntries; i++) {
        InfoTree->GetEntry(i);
        if(runSet.find(RunN) != runSet.end()) {
            if(verbose) {
                cout << "RunN: " << RunN << ", ChargeFaraday: " << ChargeFaraday << " µC" << endl;
            }
            // Adiciona a carga Faraday ao total
            TotalCharge += ChargeFaraday;
        }
    }

    if(verbose) cout<<"---------------------------------------------------"<<endl;
    return TotalCharge;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to get the differential cross-section for a given energy range and angles
// Ea: lower energy limit in MeV for the neutron energy
// Eb: upper energy limit in MeV for the neutron energy
// angles: vector of angles in degrees for which the differential cross-section will be calculated
// particles: vector of particles for which the differential cross-section will be calculated
// runsForward: vector of run numbers for the forward angle
// runsBackward: vector of run numbers for the backward angle
std::vector<std::vector<TH1D*>> GetDDX(
    float Ea = 25,
    float Eb = 26,
    const std::vector<float>& angles,
    const std::vector<char>& particles,
    const std::vector<int>& runsForward,
    const std::vector<int>& runsBackward
){

    //First thing, calculate the total charge for the runs forward and backward: 
    //the vectors runsForward(Backward) contains intervals, so they whould be read two by two
    //reading charges:: 
    Int_t intervals = runsForward.size()/2; //Should always be pair
        if(runsForward.size() % 2 != 0 || runsBackward.size() % 2 != 0) {
        std::cerr << "Error: The number of runs in runsForward and runsBackward must be even." << std::endl;
        return {};
    } 
    
    Float_t TotalChargeForward = 0;
    Float_t TotalChargeBackward = 0;
    
    //Disclaimer: 
    // runsForward and runsBackward are vectors containing intervals of runs:
    //  std::vector<int> runsForward = {12, 12, 15, 15, 18, 21, 44, 44, 56, 59}; -->12, 15, 18 to 21, 44 to 59
    //to input this into totalCharge, we need to iterate over the vector two by two:
    // std::vector<int> runs = {12, 15, 18, 19, 20, 21, 44, 56, 57, 58, 59};
    
        std::vector<int> runs;

        for (size_t i = 0; i < runsForward.size(); i += 2) {
            int start = runsForward[i];
            int end = runsForward[i + 1];

            for (int run = start; run <= end; ++run) {
                runs.push_back(run);
            }
        }
    TotalChargeForward = totalCharge(runs, true);

    runs.clear();

    for (size_t i = 0; i < runsBackward.size(); i += 2) {
        int start = runsBackward[i];
        int end = runsBackward[i + 1];

        for (int run = start; run <= end; ++run) {
            runs.push_back(run);
        }
    }
    TotalChargeBackward = totalCharge(runs, true);
    //Now we have the total charge for the forward and backward runs, so we will print them for the user
    std::cout << "Total Charge Forward: " << TotalChargeForward << " µC" << std::endl;
    std::cout << "Total Charge Backward: " << TotalChargeBackward << " µC" << std::endl;

    return {};

}