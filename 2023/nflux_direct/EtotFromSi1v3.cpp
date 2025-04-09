//
// Script for reconstructing the proton spectra 
//
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <fstream>
#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002

#include "TStopwatch.h"
#include <time.h>

TTree *tx;

const int MAX_TEL = 8;


int parse_distances_data(int year, double *distances) {
    std::string filename = "/home/pi/ganil/kalscripts/PARAMETERS_TELESCOPES/distances" + std::to_string(year) + ".dat";
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    //double distances[MAX_TEL] = {0};
    int count = 0;
    std::string line;
    while (std::getline(file, line) && count < MAX_TEL) {
        if (!line.empty() && line.find("//") != 0) {
            size_t pos = line.find("//");
            if (pos != std::string::npos) {
                line = line.substr(0, pos); // Remove comments
            }
            std::istringstream iss(line);
            double value;
            if (iss >> value) {
                distances[count++] = value;
            }
        }
    }
    return 0;
}


void EtotFromSi1v3(){
    
TStopwatch timer;
double distances[MAX_TEL] = {0}; // Declare distances array
parse_distances_data(2023, distances);

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
Double_t tSi1 = 53.4; //µm
//Double_t tSi2 = 1015; //µm
//Double_t tCsI = 5*1e4; //µm

//definition of materials
KVMaterial *det1 = new KVMaterial("Si");//,tSi1*KVUnits::um);
// KVMaterial *det2 = new KVMaterial("Si",tSi2*KVUnits::um);
// KVMaterial *det3 = new KVMaterial("CsI",tCsI*KVUnits::um);
//KVMaterial *det2_dead = new KVMaterial("Si");//,dead_th*KVUnits::um);


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


//TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/reducedv61/370.root","READ");
TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/reducedv61/388.root","READ");
//ff = new TFile("/home/dearruda/ganil/medley_2023/reduced/370v5.root","READ");


tx = (TTree*) ff->Get("M");

Long64_t nEntries = tx->GetEntries();
cout<<"number of events: "<<floor(nEntries)<<"k."<<endl;


//vector<Double_t> angles = giveMeTheAngle(nEntries);

Int_t pid;
Float_t si1,si2,csi, ang, ENN, tofn, tof_measured, tof_correction;

tx->SetBranchAddress("si1",&si1);
tx->SetBranchAddress("si2",&si2);
tx->SetBranchAddress("csi",&csi);
tx->SetBranchAddress("ENN",&ENN);
tx->SetBranchAddress("tofn",&tofn);
tx->SetBranchAddress("tof_measured",&tof_measured);
tx->SetBranchAddress("tof_correction",&tof_correction);

tx->SetBranchAddress("ang",&ang);
tx->SetBranchAddress("PID",&pid);

double progress; 
Double_t Etot, newE, n_energy; 

TH1D *horig = new TH1D("horig","horig",400,0,40);
TH1D *hnew = new TH1D("hnew","hnew",400,0,40);

TH1D *horigD = new TH1D("horigD","horigD",400,0,40);
TH1D *hnewD = new TH1D("hnewD","hnewD",400,0,40);

TH2D *h2 = new TH2D("h2","h2",400,0,40,400,0,40);
TH2D *h2energies = new TH2D("h2energies","h2energies",400,0,40,400,0,40);
TH2D *h2D = new TH2D("h2D","h2D",400,0,40,400,0,40);

//nEntries = 1e5;
det1->SetThickness(tSi1 * KVUnits::um);
//for(Long64_t l = 0;l<nEntries/1.0;l++){
for(Long64_t l = 0;l<nEntries;l++){
    tx->GetEntry(l);
    Etot = si1+si2+csi;
    if(pid==1 && ang == 20.0 && si1>0 && si2>0&& csi>0){
       
        //cout<<"si1 = "<<si1<<", si2 = "<<si2<<", csi = "<<csi<<endl;
        
        horig->Fill(Etot);
        //det1->SetThickness(tSi1/cos(angles[l]) * KVUnits::um);
        newE = si1 + det1->GetEResFromDeltaE(1,1,si1);
        hnew->Fill(newE);

        n_energy = En(tof_measured - ToFparticle(newE,'p',distances[0]/1e3));
        
        //cout<<"\n#"<<l<<"\n(Enn_old, Enn):\n"<<ENN<<" , "<<n_energy<<endl;

        //cout<<"\n#"<<l<<"\n(tofn, tof_measured, tofpart, tofpart*, tofcorrection):\n"<<
        //tofn<<" , "<<tof_measured<<" , "<<ToFparticle(si1+si2+csi,'p',distances[0]*1E-3)<<" , "<<
        //ToFparticle(newE,'p',distances[0]*1E-3)<<" , "<<tof_correction<<endl;
        
        h2->Fill(n_energy,newE);//newE
        h2energies->Fill(newE,si1+si2+csi);
    }else if(pid==2 && ang == 20.0 && si1>0 && si2>0&& csi>0){
    /*
    
        horigD->Fill(Etot);
        //det1->SetThickness(tSi1/cos(angles[l]) * KVUnits::um);
        newE = si1 + det1->GetEResFromDeltaE(1,2,si1);
        hnewD->Fill(newE);
        n_energy = En(tof_measured - ToFparticle(newE,'d',distances[0]));
        //h2D->Fill(ENN,newE);
        h2D->Fill(n_energy,newE);
    */

    }

    progress = 100.0*l/(nEntries);

        // Print the progress
    std::cout << "Progress: " << progress << "%\r";
    std::cout.flush();
}



horig->Draw();
hnew->SetLineColor(kRed);
hnew->Draw("same");

TCanvas *cv2 = new TCanvas();
h2energies->Draw("colz");
// horigD->Draw();
// hnewD->SetLineColor(kRed);
// hnewD->Draw("same");

new TCanvas();
h2->Draw();
new TCanvas();
h2D->Draw();

cout<<"COMPLETED!!!"<<endl;
timer.Print();
return;


}