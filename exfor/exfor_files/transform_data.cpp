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
#include <set>
#include <map>


void showData(TString filename = "/mnt/medley/LucasAnalysis/exfor/exfor_files/myformat/datasets.root", bool discriminate_anlges= true) {
    //function mostly constructed by chatGpt
    // Abrindo o arquivo ROOT
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }
    
    // Pegando a TTree
    TTree *tree = (TTree*)file->Get("DataTree");
    if (!tree) {
        std::cerr << "Erro: TTree 'DataTree' não encontrada no arquivo." << std::endl;
        file->Close();
        return;
    }

    // Definição de variáveis para leitura
    float energy, angle;
    char particle[50], dataset[50];
    double reference;

    tree->SetBranchAddress("EN", &energy);
    tree->SetBranchAddress("Ang", &angle);
    tree->SetBranchAddress("Particle", &particle);
    tree->SetBranchAddress("Dataset", &dataset);
    tree->SetBranchAddress("Reference", &reference);
    // Mapa para armazenar valores únicos de dataset e ângulos correspondentes
    std::map<std::tuple<float, std::string, std::string, std::string>, std::set<float>> datasetAngles;
    
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        datasetAngles[{energy, std::string(dataset), std::string(particle), std::to_string(reference)}].insert(angle);
    }
    
    // Printando os valores únicos
    std::cout << "Datasets:" << std::endl;
    for (const auto &[key, angles] : datasetAngles) {
        auto [en, ds, pt, ref] = key;
        std::cout << en << "\t" << ds << "\t" << pt << "\t" << ref << std::endl;
        if(discriminate_anlges){
            
            std::cout << "angles: ";
        
            for (const auto &ang : angles) {
                std::cout <<"  -> "<< ang << endl;
            }
            std::cout << std::endl;
        }
        
    }
    
    // Fechando o arquivo
    file->Close();
}

void transform_data(std::string namefile = "n.xP_26.5MeV_Slypen2.csv", 
                    std::string datasetname = "Slypen+(2000)",
                    std::string particlename = "p", 
                    double ref = 22718.002,  
                    std::string rootfile = "datasets.root") {

    // Caminhos completos
    namefile = std::string(getenv("PWD")) + "/" + namefile;
    rootfile = std::string(getenv("PWD")) + "/" + rootfile;

    // TTree temporária para ler CSV
    TTree *txs = new TTree("txs","txs");
    txs->ReadFile(namefile.c_str(), "EN/F:E/F:Ang/F:XS/F:XSerr/F");

    // Variáveis para leitura
    float EN_read, E_read, Ang_read, XS_read, XSerr_read;
    txs->SetBranchAddress("EN", &EN_read);
    txs->SetBranchAddress("E", &E_read);
    txs->SetBranchAddress("Ang", &Ang_read);
    txs->SetBranchAddress("XS", &XS_read);
    txs->SetBranchAddress("XSerr", &XSerr_read);

    // Abre ou cria arquivo ROOT
    TFile *file = TFile::Open(rootfile.c_str(), "UPDATE");
    if (!file || file->IsZombie()) {
        std::cerr << "Erro ao abrir " << rootfile << std::endl;
        return;
    }

    TTree *final_tree = (TTree*) file->Get("DataTree");
    bool tree_exists = (final_tree != nullptr);

    // Variáveis para branch
    float EN, E, Ang, XS, XSerr;
    double reference = ref;

    // Variáveis para preencher nova TTree
    TString particle(particlename.c_str());
    TString dataset(datasetname.c_str());

    // Ponteiros para TTree existente
    TString *existing_particle = nullptr;
    TString *existing_dataset = nullptr;
    double existing_reference = 0;
    bool dataset_exists = false;

    if (!tree_exists) {
        std::cout << "TTree não existe: criando..." << std::endl;
        final_tree = new TTree("DataTree", "Cross-section data");

        final_tree->Branch("EN", &EN, "EN/F");
        final_tree->Branch("E", &E, "E/F");
        final_tree->Branch("Ang", &Ang, "Ang/F");
        final_tree->Branch("XS", &XS, "XS/F");
        final_tree->Branch("XSerr", &XSerr, "XSerr/F");
        final_tree->Branch("Particle", &particle);
        final_tree->Branch("Dataset", &dataset);
        final_tree->Branch("Reference", &reference, "Reference/D");
    } else {
        std::cout << "TTree existe: completando..." << std::endl;

        // Conecta ponteiros para TTree existente
        final_tree->SetBranchAddress("EN", &EN);
        final_tree->SetBranchAddress("E", &E);
        final_tree->SetBranchAddress("Ang", &Ang);
        final_tree->SetBranchAddress("XS", &XS);
        final_tree->SetBranchAddress("XSerr", &XSerr);
        final_tree->SetBranchAddress("Particle", &existing_particle);
        final_tree->SetBranchAddress("Dataset", &existing_dataset);
        final_tree->SetBranchAddress("Reference", &existing_reference);

        // Verifica se dataset já existe
        Long64_t nEntries_final_tree = final_tree->GetEntries();
        for (Long64_t i = 0; i < nEntries_final_tree; i++) {
            final_tree->GetEntry(i);
            if (*existing_dataset == datasetname.c_str() && existing_reference == reference) {
                dataset_exists = true;
                break;
            }
        }
    }

    if (dataset_exists) {
        std::cout << "Dataset já existe!" << std::endl;
    } else {
        std::cout << "Adicionando novo dataset..." << std::endl;
        Long64_t nEntries = txs->GetEntries();
        for (Long64_t i = 0; i < nEntries; i++) {
            txs->GetEntry(i);
            EN = EN_read;
            E = E_read;
            Ang = Ang_read;
            XS = XS_read;
            XSerr = XSerr_read;
            particle = particlename.c_str();
            dataset = datasetname.c_str();
            reference = ref;

            final_tree->Fill();
        }
        file->Write();
        std::cout << "Dados adicionados com sucesso para o dataset '" << datasetname << "'." << std::endl;
    }

    file->Close();
    std::cout << "Arquivo ROOT fechado." << std::endl;
}