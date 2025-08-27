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

TTree *tx;

void plot(){

    TFile *ff = new TFile("370_dE.root","READ");
    TTree *tx = (TTree*) ff->Get("M");
    tx->SetAlias("goodEvents", "Etot_dE1>0 && Etot_dE2>0 && Etot_dE1_2>0 && Etot_Eres>0");

    TCanvas *c = new TCanvas();
    c->Divide(2,2);

    c->cd(1);
    tx->Draw("Etot_dE1:Etot>>h(600,0,60,600,0,60)","goodEvents","colz");
    
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    c->cd(2);
    tx->Draw("Etot_dE2:Etot>>h2(600,0,60,4000,0,400)","goodEvents","colz");
    

    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    c->cd(3);
    tx->Draw("Etot_dE1_2:Etot>>h3(600,0,60,700,0,70)","goodEvents","colz");
    
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    c->cd(4);
    tx->Draw("Etot_Eres:Etot>>h4(600,0,60,400,0,40)","goodEvents","colz");

    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();
    
    return;
}