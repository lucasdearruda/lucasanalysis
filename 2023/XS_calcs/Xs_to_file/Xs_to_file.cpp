#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "TH1D.h"
#include "TFile.h"
#include "TROOT.h"
#include "TApplication.h"

// Inclua sua função aqui
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_allAngles2.cpp"

void Xs_to_file(float Ea, float Eb, const std::string& filename) {
    std::vector<TH1D*> hists = Fe_allAngles2(Ea, Eb);
    std::ofstream outfile;

    // Abre para append, cria se não existir
    outfile.open(filename, std::ios::out | std::ios::app);

    if (!outfile.is_open()) {
        std::cerr << "Erro ao abrir arquivo: " << filename << std::endl;
        return;
    }

    float Eavg = (Ea + Eb) / 2.0;
    outfile << Eavg;

    for (int i = 0; i < hists.size(); ++i) {
        float integral = hists[i]->Integral("width");
        outfile << "," << integral;
    }

    outfile << std::endl;
    outfile.close();
}
