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
    char particle[50], dataset[50], reference[50];
    
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
        datasetAngles[{energy, std::string(dataset), std::string(particle), std::string(reference)}].insert(angle);
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

void transform_data(string namefile = "n.xP_26.5MeV_Slypen2.csv", 
                    string datasetname = "Slypen+(2000)",
                    string particlename = "p", 
                    double ref = 22718.002,  // Alterado para 'double'
                    string formatfile = "EN:E:Ang:XS:XSerr",
                    string rootfile = "datasets.root"
                    ){

    formatfile = std::regex_replace(formatfile, std::regex(":"), "/F:");
    formatfile += "/F";
    cout<<"File to be read: "<<namefile<<", format: "<<formatfile<<endl;

    namefile = std::string(getenv("PWD"))+"/"+namefile;

    TTree *txs = new TTree("txs","txs");
    txs->ReadFile(namefile.c_str(),formatfile.c_str());

    float EN_read, E_read, Ang_read, XS_read, XSerr_read;

    txs->SetBranchAddress("EN", &EN_read);
    txs->SetBranchAddress("E", &E_read);
    txs->SetBranchAddress("Ang", &Ang_read);
    txs->SetBranchAddress("XS", &XS_read);
    txs->SetBranchAddress("XSerr", &XSerr_read);

    
    rootfile =std::string(getenv("PWD"))+"/"+rootfile;

    TFile *file = TFile::Open(rootfile.c_str(), "UPDATE"); // "UPDATE" mantém os dados existentes
    if (!file || file->IsZombie()) {
        cerr << "Erro ao abrir " << rootfile << endl;
        return;
    }

    TTree *final_tree = (TTree*) file->Get("DataTree");
    bool tree_exists = (final_tree != nullptr);

    double EN, E, Ang, XS, XSerr;  // Usando 'double'

    // Substituindo 'char' para 'double'
    double reference = ref;
    char particle[5];  // Ajuste o tamanho conforme necessário
    char dataset[200];

    //dataset and particle    
    strncpy(particle, particlename.c_str(), sizeof(particle) - 1);
    particle[sizeof(particle) - 1] = '\0';  // Garantir a terminação nula

    strncpy(dataset, datasetname.c_str(), sizeof(dataset) - 1);
    dataset[sizeof(dataset) - 1] = '\0';  // Garantir a terminação nula

    char existing_reference;
    bool dataset_exists = false;

    if (!tree_exists) {
        cout<<"TTree does not exist: creating it..."<<endl;
        final_tree = new TTree("DataTree", "Cross-section data");

        final_tree->Branch("EN", &EN, "EN/F");
        final_tree->Branch("E", &E, "E/F");
        final_tree->Branch("Ang", &Ang, "Ang/F");
        final_tree->Branch("XS", &XS, "XS/F");
        final_tree->Branch("XSerr", &XSerr, "XSerr/F");
        final_tree->Branch("Particle", &particle, "particle/C");
        final_tree->Branch("Dataset", &dataset,"dataset/C");
        final_tree->Branch("Reference", &reference,"reference/D");  // 'D' para 'double'

    } else {
        cout<<"TTree does exist: completing it..."<<endl;
        final_tree->SetBranchAddress("EN", &EN);
        final_tree->SetBranchAddress("E", &E);
        final_tree->SetBranchAddress("Ang", &Ang);
        final_tree->SetBranchAddress("XS", &XS);
        final_tree->SetBranchAddress("XSerr", &XSerr);
        final_tree->SetBranchAddress("Particle", &particle);
        final_tree->SetBranchAddress("Dataset", &dataset);
        final_tree->SetBranchAddress("Reference", &existing_reference);
    }

    int nEntries_final_tree = final_tree->GetEntries();
    for (int i = 0; i < nEntries_final_tree; i++) {
        final_tree->GetEntry(i);
        
        // Comparação entre double
        if (existing_reference == reference) {
            dataset_exists = true;
            break;
        }
    }

    if (dataset_exists) {
        std::cout << "Dataset encontrado!" << std::endl;
    } else {
        std::cout << "Dataset não encontrado!" << std::endl;
    }

    int nEntries = txs->GetEntries();
    for (int i = 0; i < nEntries; i++) {
        txs->GetEntry(i);
        EN = EN_read;
        E = E_read;
        Ang = Ang_read;
        XS = XS_read;
        XSerr = XSerr_read;

        final_tree->Fill();
    }

    if (!dataset_exists) {
        file->Write();  // Salva as mudanças no arquivo ROOT
        cout << "Dados adicionados com sucesso para o dataset '" << dataset << "', '" <<reference<<"'."<< endl;
    } else {
        cout << "Dados já existem para o dataset '" << dataset << "', '" <<reference<<"'."<< endl;
    }

    if (file) {
        file->Close();  // Fecha o arquivo ROOT
        std::cout << "Arquivo ROOT fechado com sucesso." << std::endl;
    } else {
        std::cerr << "Erro ao tentar fechar o arquivo ROOT." << std::endl;
    }

    return;
}
