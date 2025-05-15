//
// Script for evaluation the ddx for Fe in 2023 runs
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
//#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002
#include "Fe_ddx.cpp" // for the energy calibration
#include "plotMeAxs.cpp" // for the energy calibration

#include "TStopwatch.h"
#include <time.h>

TH1D* hddx;
//Int_t which;
TCanvas* plotDDX(Float_t E = 20, char particle = 'p', Float_t angle = 20.0, Float_t dE = 1, Int_t which = 1){

    hddx = Fe_ddx(E - dE/2,
                     E + dE/2,
                            which, //with TTC but no dead correction
                        angle,
                     particle);

                     
    TCanvas *cv = new TCanvas(Form("cv_%.1fMeV_%c_%.1fdeg",E,particle,angle),Form("Canvas EN = %.1f MeV, %c, %.1fdeg",E,particle,angle));
    hddx->SetTitle(Form("DDX %c, %.1fdeg",particle,angle));
    hddx->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
    hddx->GetXaxis()->SetTitle(Form("E_{%c} (MeV)",particle));
    hddx->SetName("hddx");


    hddx->Draw();
    gPad->SetGridx();
    gPad->SetGridy();
    plotMeAxs();
    
    
    
    return cv;
}

void makeTL(TH1D*h = (TH1D*)gROOT->FindObject("hddx"), TGraph *gr = (TGraph*)gROOT->FindObject("gr"), string nameH = "reaction", string name1 = "Medley", string name2 = "Talys 2.0"){
    TLegend *tl = new TLegend(0.4,0.53,0.88,0.86);
    tl->SetHeader(nameH.c_str());
    tl->AddEntry(h,name1.c_str(),"lpf");
    tl->AddEntry(gr,name2.c_str(),"lpf");
    tl->Draw();
}