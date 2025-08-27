// #include <iostream>
// #include <fstream>
// #include <string>
// #include "TChain.h"
// #include "TTree.h"
// #include "TFile.h"
// #include "TStopwatch.h"
#include "KVMaterial.h" // Ensure this header is included

// #include "/mnt/medley/LucasAnalysis/useful.h"

//bool considerangle = false,
TH1D *newE(
    TH1D *hh = nullptr, 
    string hname = "h", 
    string htitle = "htitle", 
    Double_t tSi1 = 53.4, 
    Int_t Z = 1,
    Int_t A = 1
){
    TStopwatch timer;
    timer.Start();
    
    if(!hh){
        cerr<<"Error: The input histogram is null."<<endl;
        return nullptr;
    }
    
    KVMaterial *det1 = new KVMaterial("Si");
    det1->SetThickness(tSi1 * KVUnits::um);
    cout << "Thickness of Si1: " << tSi1 << " µm" << endl;
    Float_t maxdE = det1->GetMaxDeltaE(Z,A); //condition for the max dE in Si1


    Int_t nbins = hh->GetNbinsX();
    Int_t ncounts = hh->GetEntries();
    TH1D *h = new TH1D(hname.c_str(), htitle.c_str(), nbins, hh->GetBinLowEdge(0),  hh->GetBinLowEdge(nbins+1));


    // TRandom2 *rn = new TRandom2();
    // std::time_t seed = std::time(nullptr);
    // rn->SetSeed(seed);
    
    bool considerangle = false;
    vector<Double_t> angles;
    if(considerangle) vector<Double_t> angles = giveMeTheAngle2(ncounts, 20);
    
    
    
    Double_t si1, newE,angle_factor;

    angle_factor = 1;
    Int_t not_calc = 0;
    for(Int_t i=0;i<ncounts;i++){
        si1 = hh->GetRandom();
        
        if(si1 > maxdE){
            not_calc++;
            continue; // Skip this iteration if si1 exceeds max dE
        }

        if(considerangle){
            angle_factor = 1 /cos(angles[i]);
            det1->SetThickness(tSi1 * angle_factor * KVUnits::um);    
        }
        

        newE = si1 + det1->GetEResFromDeltaE(Z, A, si1)/KVUnits::MeV;
        h->Fill(newE);
    }

    cout << "Number of entries not calculated due to max dE condition: " << not_calc << endl;
    cout << "Total entries processed: " << ncounts << endl;
    timer.Print();
    return h;

}
