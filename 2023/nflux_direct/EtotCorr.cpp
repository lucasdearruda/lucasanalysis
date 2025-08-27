//
// Script for reconstructing the proton spectra 
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
#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002

#include "TStopwatch.h"
#include <time.h>

TH1D *EtotCorr( TH1D *h = nullptr, Double_t th_si1 = 53.4 ){//thickness of Si1 in µm
TStopwatch timer;


iff(h == nullptr){
    std::cerr << "Error: The input histogram is null." << std::endl;
    return nullptr;
}

KVMaterial *det1 = new KVMaterial("Si");



timer.Print();
}