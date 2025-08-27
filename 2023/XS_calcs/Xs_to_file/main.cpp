#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TGraph2D.h"

using namespace std;

// Prototype
void Xs_to_file(float Ea, float Eb, const std::string& filename);

int main(int argc, char** argv) {
    const std::string outFile = "XS_output.csv";

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " Ea Eb" << std::endl;
        return 1;
    }

    float Ea = std::stof(argv[1]);
    float Eb = std::stof(argv[2]);

    Xs_to_file(Ea, Eb, outFile);

    std::cout << "Arquivo CSV gerado: " << outFile << std::endl;
    return 0;
}
