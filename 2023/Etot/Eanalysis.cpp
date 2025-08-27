// script to calculate the energies deposited in the different parts of te detector
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
#include "/mnt/medley/LucasAnalysis/useful.h" //version 2024.11.25.001

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


std::string extractname(string name, string what_to_inset = "_dE.root"){//verified 
    
    std::string newFilename;
    size_t lastSlash = name.find_last_of('/');
    size_t lastDot = name.find_last_of('.');
    
    if (lastSlash != std::string::npos && lastDot != std::string::npos && lastDot > lastSlash) {
        // Extrair o nome do arquivo sem a extensão
        std::string baseName = name.substr(lastSlash + 1, lastDot - lastSlash - 1);
        
        // Construir o novo nome de arquivo
        newFilename = baseName + what_to_inset;
        
        std::cout << "Novo nome do arquivo: " << newFilename << std::endl;
    } else {
        std::cerr << "Erro: Não foi possível extrair o nome base." << std::endl;
    }

    return newFilename;
}

void Eanalysis(string runfile = "/mnt/medley/LucasAnalysis/2023/reducedv61/370.root", bool savehistos = false, string savingname = "/mnt/medley/LucasAnalysis/2023/Etot/protons370_energies.root"){
    
TStopwatch timer;
double distances[MAX_TEL] = {0}; // Declare distances array
parse_distances_data(2023, distances);

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
Double_t tSi1 = 53.4; //µm
Double_t tSi2 = 1015; //µm
//Double_t tCsI = 5*1e4; //µm

//definition of materials
KVMaterial *det1 = new KVMaterial("Si",tSi1*KVUnits::um);
KVMaterial *det2 = new KVMaterial("Si",tSi2*KVUnits::um);
KVMaterial *det1n2 = new KVMaterial("Si",(tSi1+tSi2)*KVUnits::um);
KVMaterial *det3 = new KVMaterial("CsI");
//KVMaterial *det2_dead = new KVMaterial("Si");//,dead_th*KVUnits::um);


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


TFile *ff = new TFile(runfile.c_str(),"READ");

tx = (TTree*) ff->Get("M");
tx->SetBranchStatus("*",0); // Desabilita todos os branches
tx->SetBranchStatus("si1",1);
tx->SetBranchStatus("si2",1);
tx->SetBranchStatus("csi",1);
tx->SetBranchStatus("ang",1);
tx->SetBranchStatus("PID",1);
tx->SetBranchStatus("ENN",1);

Long64_t nEntries = tx->GetEntries();
cout<<"number of events: "<<floor(nEntries)<<"k."<<endl;

// Criar novo arquivo e árvore =======================================================
string newFileName = extractname(runfile, "_dE.root");
cout<<"New file name: "<<newFileName<<endl;
TFile *fn = new TFile(newFileName.c_str(), "RECREATE");
TTree *newTree = tx->CloneTree(0,"fast"); // Clona a estrutura, mas não copia os dados

//vector<Double_t> angles = giveMeTheAngle(nEntries);

Int_t pid;
Float_t si1,si2,csi, ang, ENN, tofn, tof_measured, tof_correction;

tx->SetBranchAddress("si1",&si1);
tx->SetBranchAddress("si2",&si2);
tx->SetBranchAddress("csi",&csi);
tx->SetBranchAddress("ENN",&ENN);
tx->SetBranchAddress("ang",&ang);
tx->SetBranchAddress("PID",&pid);

double progress; 
// Float_t newValue;
//     newTree->Branch("newValue", &newValue, "newValue/F");
Float_t    Etot, 
            Etot_dE1, 
            Etot_dE2,
            Etot_dE1_2,
            Etot_Eres; 

newTree->Branch("Etot", &Etot, "Etot/F");
newTree->Branch("Etot_dE1", &Etot_dE1, "Etot_dE1/F");
newTree->Branch("Etot_dE2", &Etot_dE2, "Etot_dE2/F");
newTree->Branch("Etot_dE1_2", &Etot_dE1_2, "Etot_dE1_2/F");
newTree->Branch("Etot_Eres", &Etot_Eres, "Etot_Eres/F");

Float_t si1l,csil, det12l;
Long64_t nanCounter = 0;

for(Long64_t l = 0;l<nEntries;l++){
    tx->GetEntry(l);
    Etot = si1+si2+csi;

    if(pid==1 && ang == 20.0 && si1>0. && si2>0&& csi>0){
       
        //cout<<"si1 = "<<si1<<", si2 = "<<si2<<", csi = "<<csi<<endl;
        //det1->SetThickness(tSi1/cos(angles[l]) * KVUnits::um);

        Etot_dE1 = si1 + det1->GetEResFromDeltaE(1,1,si1);

            csil = det2->GetEResFromDeltaE(1,1,si2);
            
            if(isnan(csil) || csil<0.05){
                csil = 0.0;
            }
            

            si1l = det1->GetDeltaEFromERes(1,1,si2 + csil);

        Etot_dE2 = si2 + si1l + csil;
        

        det12l = det1n2->GetEResFromDeltaE(1,1,si1 + si2 );
        if(isnan(det12l) || det12l<0.05){
            det12l = 0.0;
        }
        Etot_dE1_2 = si1 + si2 + det12l;
        
        
        if(csi<0.05){
            Etot_Eres = si1+si2;    
        }else{
            Etot_Eres = det1n2->GetDeltaEFromERes(1,1,csi) + csi;    
        }
        

    }else if(pid==2 && ang == 20.0 && si1>0.05 && si2>0.05&& csi>0.0){
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
    if(isnan(Etot_dE1)||isnan(Etot_dE2)||isnan(Etot_dE1_2)||isnan(Etot_Eres)){
        nanCounter++; 
        cout<<"\n.\n!\n.nanCounter = "<<nanCounter<<" - - - "<< 100.0*nanCounter/nEntries<<endl;
        cout<<"si1 = "<<si1<<", si2 = "<<si2<<", csi = "<<csi<<endl;
        cout<<"Etot_dE1 = "<<Etot_dE1<<", Etot_dE2 = "<<Etot_dE2<<", Etot_dE1_2 = "<<Etot_dE1_2<<", Etot_Eres = "<<Etot_Eres<<endl;
        cout<<"\n:::\n:::\n:::"<<endl;
        
        if(isnan(Etot_dE1)){
            Etot_dE1 = -1.0;
        }
        if(isnan(Etot_dE2)){
            Etot_dE2 = -1.0;
        }
        if(isnan(Etot_dE1_2)){
            Etot_dE1_2 = -1.0;
        }
        if(isnan(Etot_Eres)){
            Etot_Eres = -1.0;
        }

        if(Etot_dE1<0.01){
            Etot_dE1 =  -9;
        }
        if(Etot_dE2<0.01){
            Etot_dE2 = -9;
        }
        if(Etot_dE1_2<0.01){
            Etot_dE1_2 = -9;
        }
        if(Etot_Eres<0.01){
            Etot_Eres = -9;
        }
        //Etot_dE1 = Etot_dE2 = Etot_dE1_2 = Etot_Eres = 0.0;
        
        
        newTree->Fill(); // Preenche a nova árvore com os dados antigos e novos
    }else{
        newTree->Fill(); // Preenche a nova árvore com os dados antigos e novos
    }
    

    progress = 100.0*l/(nEntries);

        // Print the progress
    std::cout << "Progress: " << progress << "%\r";
    std::cout.flush();
}
newTree->Write(); // Escreve a nova árvore no arquivo
fn->Close();
ff->Close();


cout<<"COMPLETED!!!"<<endl;
timer.Print();



return;


}