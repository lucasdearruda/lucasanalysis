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


void createTree(){
    TRandom2 *rn = new TRandom2();

    TFile *f = new TFile("test.root", "RECREATE");
    TTree *tx = new TTree("tx", "tx");
    
    Float_t value, other;
    tx->Branch("value", &value, "value/F"); // Forma correta de criar o branch
    tx->Branch("other", &other, "other/F"); // Forma correta de criar o branch
    
    for (int i = 0; i < 25; ++i) {
        value = rn->Gaus();
        other = 10 + -0.5*rn->Gaus();
        std::cout << value << std::endl;
        tx->Fill();
    }

    tx->Write();
    f->Close();
    
    delete rn;
    delete f; // Boa prática para evitar memory leaks

    return;
}

void addbranch(){

    TRandom2 *rn = new TRandom2();

    TFile *f = new TFile("test.root", "UPDATE");
    TTree *tx = (TTree*) f->Get("tx");
    if (!tx) {
        std::cerr << "Error: Tree not found!" << std::endl;
        return;
    }
    Float_t Value;
    tx->SetBranchAddress("value", &Value); // Corrigido para usar o mesmo nome do branch

    Float_t newValue;
    tx->Branch("newValue", &newValue, "newValue/F"); // Corrigido para usar o mesmo nome do branch
    Long64_t nEntries = tx->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tx->GetEntry(i);
        newValue = std::abs(Value +1.0); // Exemplo de operação
        cout<<"Value: " << Value << ", newValue: " << newValue << endl;
        tx->Fill();
        /* code */
    }
    
    tx->Write();
    f->Close();
}
void addbranchNewFile() {
    TRandom2 *rn = new TRandom2();

    // Abrir o arquivo original
    TFile *f = new TFile("test.root", "READ");
    TTree *tx = (TTree*) f->Get("tx");
    if (!tx) {
        std::cerr << "Error: Tree not found!" << std::endl;
        return;
    }

    // Criar novo arquivo e árvore
    TFile *fn = new TFile("test_new.root", "RECREATE");
    TTree *newTree = tx->CloneTree(0); // Clona a estrutura

    // Criar variável para a branch original
    Float_t Value;
    tx->SetBranchAddress("value", &Value);

    // Criar nova branch no novo arquivo
    Float_t newValue;
    newTree->Branch("newValue", &newValue, "newValue/F");

    // Loop sobre as entradas e preencher a nova árvore
    Long64_t nEntries = tx->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tx->GetEntry(i);
        newValue = std::abs(Value + 1.0); // Operação exemplo
        std::cout << "Value: " << Value << ", newValue: " << newValue << std::endl;
        
        newTree->Fill(); // Preenche a nova árvore com os dados antigos e novos
    }

    // Escrever e fechar os arquivos
    newTree->Write();
    fn->Close();
    f->Close();

}


void addbranchNewFile2() {
    TRandom2 *rn = new TRandom2();

    // Abrir o arquivo original
    TFile *f = new TFile("test.root", "READ");
    TTree *tx = (TTree*) f->Get("tx");
    if (!tx) {
        std::cerr << "Error: Tree not found!" << std::endl;
        return;
    }



    tx->SetBranchStatus("*",0); // Desabilita todos os branches
    tx->SetBranchStatus("value",1);
    

    // Criar novo arquivo e árvore
    TFile *fn = new TFile("test_new.root", "RECREATE");
    TTree *newTree = tx->CloneTree(0, "fast"); // Clona a estrutura

    // Criar variável para a branch original
    Float_t Value;
    tx->SetBranchAddress("value", &Value);

    // Criar nova branch no novo arquivo
    Float_t newValue;
    newTree->Branch("newValue", &newValue, "newValue/F");

    // Loop sobre as entradas e preencher a nova árvore
    Long64_t nEntries = tx->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tx->GetEntry(i);
        newValue = std::abs(Value + 1.0); // Operação exemplo
        std::cout << "Value: " << Value << ", newValue: " << newValue << std::endl;
        
        newTree->Fill(); // Preenche a nova árvore com os dados antigos e novos
    }

    // Escrever e fechar os arquivos
    newTree->Write();
    fn->Close();
    f->Close();

}

void test() {
    
    createTree();
    return;
}