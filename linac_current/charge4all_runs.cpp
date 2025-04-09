#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <TFile.h>
#include <TGraph.h>
#include "calculateCharge.cpp"

void addChargeHeader(std::ofstream& outfile) {
    // Escreve a linha de cabeçalho com uma coluna adicional para Charge (µC)
    outfile << "#RUN\tDAY\tmonth\thour\tminute\tDAY\tmonth\thour\tminute\tCharge (µC)" << std::endl;
}

int charge4all_runs() {
    // Abre o arquivo ROOT e obtém o TGraph
    TFile *file = TFile::Open("total.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Erro ao abrir total.root" << std::endl;
        return 1;
    }
    TGraph *graph = (TGraph*)file->Get("ACCT11_temps");
    if (!graph) {
        std::cerr << "Erro ao obter TGraph de total.root" << std::endl;
        return 1;
    }

    // Abre o arquivo runs_summary.txt
    std::ifstream infile("runs_summary.txt");
    if (!infile.is_open()) {
        std::cerr << "Erro ao abrir runs_summary.txt" << std::endl;
        return 1;
    }

    // Abre o arquivo runs_results.txt para escrever os resultados
    std::ofstream outfile("runs_results.txt");
    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir runs_results.txt" << std::endl;
        return 1;
    }

    // Adiciona o cabeçalho com a coluna adicional para Charge (µC)
    addChargeHeader(outfile);

    std::string line;
    // Lê cada linha do arquivo de entrada e processa-a
    while (std::getline(infile, line)) {
        // Trim dos espaços em branco à esquerda e à direita
        line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        
        if (line.empty()) {
            continue;  // Ignora linhas vazias
        }

        // Extrai os valores da linha diretamente
        std::istringstream iss(line);
        int run_number, day1, month1, hour1, min1, day2, month2, hour2, min2;
        if (!(iss >> run_number >> day1 >> month1 >> hour1 >> min1 >> day2 >> month2 >> hour2 >> min2)) {
            std::cerr << "Erro ao ler linha: " << line << std::endl;
            continue;  // Avança para a próxima linha em caso de erro
        }

        // Calcula a carga usando a função calculateCharge
        Double_t charge = calculateCharge(graph, month1, day1, hour1, min1, month2, day2, hour2, min2);

        // Escreve os resultados no arquivo de saída
        outfile << run_number << "\t" << day1 << "\t" << month1 << "\t" << hour1 << "\t" << min1 << "\t"
                << day2 << "\t" << month2 << "\t" << hour2 << "\t" << min2 << "\t" << charge << std::endl;
    }

    // Fecha os arquivos
    infile.close();
    outfile.close();
    file->Close();

    return 0;
}
