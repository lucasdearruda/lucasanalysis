#include <iostream>
#include <optional>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>

// Definição do struct para armazenar os dados
struct MedleyData {
    Int_t RunN;
    std::string MedleyConfig;
    std::string Target;
    Int_t RunTime;
    Float_t RunCharge;
};

// Função para obter os dados de uma run específica
std::optional<MedleyData> GetRunData(const std::string& filepath, Int_t requestedRunN) {
    // Criar o TTree e ler o arquivo
    TTree* InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(filepath.c_str(), "RunN/I:MedleyConfig/C:Target/C:RunTime/I:RunCharge/F");

    // Declarar variáveis para leitura dos dados
    Int_t RunN;
    char MedleyConfig[256];
    char Target[256];
    Int_t RunTime;
    Float_t RunCharge;

    // Associar as variáveis às branches do TTree
    InfoTree->SetBranchAddress("RunN", &RunN);
    InfoTree->SetBranchAddress("MedleyConfig", MedleyConfig);
    InfoTree->SetBranchAddress("Target", Target);
    InfoTree->SetBranchAddress("RunTime", &RunTime);
    InfoTree->SetBranchAddress("RunCharge", &RunCharge);

    // Loop sobre as entradas no TTree
    Long64_t nEntries = InfoTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        InfoTree->GetEntry(i);

        // Verificar se o RunN atual é o solicitado
        if (RunN == requestedRunN) {
            MedleyData entry;
            entry.RunN = RunN;
            entry.MedleyConfig = MedleyConfig;
            entry.Target = Target;
            entry.RunTime = RunTime;
            entry.RunCharge = RunCharge;

            delete InfoTree; // Liberar memória
            return entry; // Retorna os dados da run encontrada
        }
    }

    // Caso não encontre a run, liberar memória e retornar vazio
    delete InfoTree;
    return std::nullopt;
}

// Exemplo de uso
void test(Int_t runToSearch = 27) {
    std::string infofile = "/media/dearruda/Elements/LucasAnalysis/2023/MEDLEY2023.csv";
    //Int_t runToSearch = 27; // RunN que queremos buscar

    // Buscar dados da RunN específica
    auto result = GetRunData(infofile, runToSearch);

    // Exibir os resultados
    if (result) {
        std::cout << "RunN: " << result->RunN
                  << ", MedleyConfig: " << result->MedleyConfig
                  << ", Target: " << result->Target
                  << ", RunTime: " << result->RunTime
                  << ", RunCharge: " << result->RunCharge
                  << std::endl;
    } else {
        std::cout << "RunN " << runToSearch << " não encontrado no arquivo." << std::endl;
    }

    return;
}
