// loadcals.cc
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>   // std::replace
#include <cctype>
#include <TString.h>   // ROOT: para Form()
#include <Rtypes.h>    // ROOT: para Float_t (se necessário)

using namespace std;
Float_t g1[8], g2[8], g3[8];

void loadCals(const char* filename = "src/calibrations.dat") { // created with chatgpt;; benchmarket by Lucas
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Erro ao abrir arquivo de calibracao: " << filename << endl;
        return;
    }

    string line;
    int idx = 0;

    while (getline(file, line) && idx < 8) {
        // Ignora comentários e linhas vazias
        if (line.empty() || line[0] == '#') continue;

        istringstream iss(line);
        if (!(iss >> g1[idx] >> g2[idx] >> g3[idx])) {
            cerr << "Formato inválido na linha (indice " << idx << "): " << line << endl;
            break;
        }
        idx++;
    }
}

void printCals(){
    for(int i=0;i<8;i++) cout << Form("%d | %10.9f, %10.9f, %10.9f",i+1,g1[i],g2[i],g3[i])<<endl;
}

// = = = = = = Benchmark = = = = = = Benchmark = = = = = = Benchmark = = = = = = Benchmark = = = = = =
//    ------------------------------------------------------------------
//   | Welcome to ROOT 6.32.08                        https://root.cern |
//   | (c) 1995-2024, The ROOT Team; conception: R. Brun, F. Rademakers |
//   | Built for linuxx8664gcc on May 01 2025, 14:12:33                 |
//   | From tags/6-32-08@6-32-08                                        |
//   | With c++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0                   |
//   | Try '.help'/'.?', '.demo', '.license', '.credits', '.quit'/'.q'  |
//    ------------------------------------------------------------------
//
// root [0] .L loadcals.cc 
// root [1] loadCals("calibrations.dat")
// root [2] printCals()
// 1 | 0.000402000, 0.002275000, 0.001694278
// 2 | 0.000357333, 0.002198070, 0.001178630
// 3 | 0.000411000, 0.002224000, 0.000822500
// 4 | 0.000441000, 0.002247000, 0.000911600
// 5 | 0.001735000, 0.002292000, 1.000000000
// 6 | 0.001888000, 0.002310000, 1.000000000
// 7 | 0.001826000, 0.002232000, 1.000000000
// 8 | 0.001872000, 0.002250000, 1.000000000
// root [3] 
// - - - b - - - e - - - n - - - c - - - h - - - m - - - a - - - r - - - k  - - - E - - - n - - - d - - -