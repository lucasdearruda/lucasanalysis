//
// Script for transforming the data for a proper format, so we can deal with it in a good way
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
#include "TStopwatch.h"
#include <time.h>

#include <regex>


void transform_data(string namefile = "n.xP_26.5MeV_Slypen2.csv", 
                    string datasetname = "Slypen+(2000)",
                    string particlename = "p", 
                    string formatfile = "EN:E:Ang:XS:XSerr", 
                    string reference = "22718.002"
                    ){

formatfile = std::regex_replace(formatfile, std::regex(":"), "/F:");
formatfile += "/F";
cout<<"File to be read: "<<namefile<<", format: "<<formatfile<<endl;

TTree *txs = new TTree("txs","txs");
txs->ReadFile(namefile.c_str(),formatfile.c_str());
//new TCanvas();

float EN_read, E_read, Ang_read, XS_read, XSerr_read;

txs->SetBranchAddress("EN", &EN_read);
txs->SetBranchAddress("E", &E_read);
txs->SetBranchAddress("Ang", &Ang_read);
txs->SetBranchAddress("XS", &XS_read);
txs->SetBranchAddress("XSerr", &XSerr_read);

/////////////////////////////////////////////////////////////////////////////////////////:
// Caminho do arquivo ROOT onde os dados ficarão armazenados

    string rootfile = "myformat/datasets.root";

    // Abrir (ou criar) o arquivo ROOT principal
    TFile *file = TFile::Open(rootfile.c_str(), "UPDATE"); // "UPDATE" mantém os dados existentes
    if (!file || file->IsZombie()) {
        cerr << "Erro ao abrir " << rootfile << endl;
        return;
    }

    // Verificar se a árvore final já existe
    TTree *final_tree = (TTree*) file->Get("DataTree");
    bool tree_exists = (final_tree != nullptr);

    Float_t EN, E, Ang, XS, XSerr;
    // Declare a char array to hold the string values
    char particle[5];  // Adjust size as needed for the string length
    char dataset[200];


    //dataset and particle    
    strncpy(particle, particlename.c_str(), sizeof(particle) - 1);
    particle[sizeof(particle) - 1] = '\0';  // Garantir a terminação nula

    strncpy(dataset, datasetname.c_str(), sizeof(dataset) - 1);
    dataset[sizeof(dataset) - 1] = '\0';  // Garantir a terminação nula

    // Agora as variáveis 'particle' e 'dataset' têm os conteúdos das strings
    std::cout << "Particle: " << particle << std::endl;
    std::cout << "Dataset: " << dataset << std::endl;

    char existing_dataset[200];  // Para armazenar o dataset existente no TTree
    bool dataset_exists = false;
    // Criar a árvore principal se ainda não existir
    if (!tree_exists) {
        final_tree = new TTree("DataTree", "Cross-section data");


        final_tree->Branch("EN", &EN, "EN/F");
        final_tree->Branch("E", &E, "E/F");
        final_tree->Branch("Ang", &Ang, "Ang/F");
        final_tree->Branch("XS", &XS, "XS/F");
        final_tree->Branch("XSerr", &XSerr, "XSerr/F");
        final_tree->Branch("Particle", &particle, "particle/C");
        final_tree->Branch("Dataset", &dataset,"dataset/C");
        

    }else{ // se existir settamos a coisa 
        final_tree->SetBranchAddress("EN", &EN);
        final_tree->SetBranchAddress("E", &E);
        final_tree->SetBranchAddress("Ang", &Ang);
        final_tree->SetBranchAddress("XS", &XS);
        final_tree->SetBranchAddress("XSerr", &XSerr);
        final_tree->SetBranchAddress("Particle", &particle);
        final_tree->SetBranchAddress("Dataset", &existing_dataset);
    }


    

    int nEntries_final_tree = final_tree->GetEntries();
    for (int i = 0; i < nEntries_final_tree; i++) {
        final_tree->GetEntry(i);
        

        // Depuração: Imprimir as strings antes da comparação
        //std::cout << "Comparando o dataset existente: '" << existing_dataset << "' com o dataset atual: '" << dataset << "'" << std::endl;
        
        // Comparação entre as strings de char com strcmp
        if (strcmp(existing_dataset, dataset) == 0) {
            dataset_exists = true;
            break;  // Se encontrar o dataset, não precisa continuar a busca
        }
    }

    if (dataset_exists) {
        std::cout << "Dataset encontrado!" << std::endl;
    } else {
        std::cout << "Dataset não encontrado!" << std::endl;
        final_tree->SetBranchAddress("Dataset", &dataset);
    }


int nEntries = txs->GetEntries();
for (int i = 0; i < nEntries; i++) {
    // Carregue os dados para a entrada i
    if(dataset_exists) break;
    txs->GetEntry(i);
    EN = EN_read;
    E = E_read;
    Ang = Ang_read;
    XS = XS_read;
    XSerr = XSerr_read;


    final_tree->Fill();
}

    if(!dataset_exists){
        file->Write();  // Salva as mudanças no arquivo ROOT
        cout << "Dados adicionados com sucesso para o dataset '" << dataset << "'." << endl;
    }else{
        cout << "Dados já existem para o dataset '" << dataset << "'." << endl;
    }
    


    return;
}