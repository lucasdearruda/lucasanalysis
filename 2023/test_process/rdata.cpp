
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TRandom2.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TGraph.h"
#include "TFile.h"
#include "TGraphErrors.h"

#include "rundata.hh"


#include <assert.h>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>  // Para std::cout e std::endl
#include <string>    // Para std::string


int main(int argc, char* argv[]) {

    

    string cur_time = getCurrentTime();
    clock_t tStart = clock();
    
    bool tof_correction_bool = kTRUE;

    cout<<"# OF PROVIDED ARGUMENTS = "<<argc<<endl;
    for(int i=1; i<argc;i++)cout<<"argument #"<<i<<": "<<argv[i]<<endl;
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

        std::string infofile = "/mnt/medley/LucasAnalysis/2023/MEDLEY2023.csv";
        Int_t nrun = 39;  // Exemplo de número da run

        // Criar uma instância da classe RunDataHandler
        RunDataHandler runDataHandler(infofile);

        // Obter os dados da run desejada
        auto result = runDataHandler.GetRunData(nrun);

        // Verificar se os dados foram encontrados
        if (result) {
            const MedleyData& data = result.value();
            std::cout << "Run Number: " << data.RunN << std::endl;
            std::cout << "Config: " << data.MedleyConfig << std::endl;
            std::cout << "Target: " << data.Target << std::endl;
            std::cout << "Runtime: " << data.RunTime << std::endl;
            std::cout << "Charge: " << data.RunCharge << std::endl;
        } else {
            std::cout << "Run not found." << std::endl;
        }


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
    //timer.Print();

    cout<<"Compiled from processRun_v6.1d.cpp, version 6.2024-12-06.1"<<endl;
    return 0;

}
