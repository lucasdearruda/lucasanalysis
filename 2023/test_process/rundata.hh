#include <iostream>
#include <optional>
#include <string>
#include <TTree.h>

using namespace std;
// Definição do struct para armazenar os dados
struct MedleyData {
    Int_t RunN;
    std::string MedleyConfig;
    std::string Target;
    Int_t RunTime;
    Float_t RunCharge;
};


std::string getCurrentTime() {
    char cur_time[128];
    time_t t = time(NULL);
    struct tm* ptm = localtime(&t);

    // Formata o tempo no formato desejado: "YYYY-MM-DD_HH:MM:SS"
    strftime(cur_time, sizeof(cur_time), "%Y-%m-%d_%H:%M:%S", ptm);

    return std::string(cur_time);
}


class RunDataHandler {
private:
    std::string filepath;
    TTree* infoTree;

public:
    // Construtor da classe
    RunDataHandler(const std::string& file) : filepath(file), infoTree(nullptr) {}

    // Destruidor para liberar recursos
    ~RunDataHandler() {
        if (infoTree) {
            delete infoTree;
        }
    }

    // Função para obter os dados de uma run específica
    std::optional<MedleyData> GetRunData(Int_t requestedRunN) {
        // Criar o TTree e ler o arquivo, se necessário
        if (!infoTree) {
            infoTree = new TTree("InfoTree", "InfoTree");
            infoTree->ReadFile(filepath.c_str(), "RunN/I:MedleyConfig/C:Target/C:RunTime/I:RunCharge/F");
        }

        // Declarar variáveis para leitura dos dados
        Int_t RunN;
        char MedleyConfig[256];
        char Target[256];
        Int_t RunTime;
        Float_t RunCharge;

        // Associar as variáveis às branches do TTree
        infoTree->SetBranchAddress("RunN", &RunN);
        infoTree->SetBranchAddress("MedleyConfig", MedleyConfig);
        infoTree->SetBranchAddress("Target", Target);
        infoTree->SetBranchAddress("RunTime", &RunTime);
        infoTree->SetBranchAddress("RunCharge", &RunCharge);

        // Loop sobre as entradas no TTree
        Long64_t nEntries = infoTree->GetEntries();
        for (Long64_t i = 0; i < nEntries; ++i) {
            infoTree->GetEntry(i);

            // Verificar se o RunN atual é o solicitado
            if (RunN == requestedRunN) {
                MedleyData entry;
                entry.RunN = RunN;
                entry.MedleyConfig = MedleyConfig;
                entry.Target = Target;
                entry.RunTime = RunTime;
                entry.RunCharge = RunCharge;
                return entry;  // Retorna os dados da run encontrada
            }
        }

        // Caso não encontre a run, retorna std::nullopt
        return std::nullopt;
    }
};