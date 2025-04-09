

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
#include "TLeaf.h" // Make sure to include the necessary headers
#include "TChain.h"

#include "TGraph.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"

#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

#include "/home/e802/Analysis/LucasAnalysis/useful.h"
//#include "functions.h"
///
//
//
//Function to calculate the subtraction factor::
double sub_fac(double Qa, double Qb, double ma, double mb){// A is CH2 and B is C

float mmch2 = 14.0; 
float mmc = 12.0; 

return (mmc/mmch2)*(Qa/Qb)*(ma/mb);
}


//TH1D *NeutronFromToF(TH1D *tofHisto, bool verbose = false){//verbose put some info regarding the function working 

//run 407: 
    //double Q_CH2 = 0.444e6 ;//µC
//run 406: 
    //double Q_C = 0.256e6 ;//µC

//run 370: 
    double Q_CH2 = 0.239e6 ;//µC
//run 388: 
    double Q_C = 0.306e6 ;//µC


int eval_flux(){

    TStopwatch timer;
    string namefile_CH2,namefile_C,namefile_cut;

    //namefile_CH2 = "run_407_specific.root";
    namefile_CH2 = "run_370_specific.root";
    //namefile_C = "run_406_specific.root";
    namefile_C = "run_388_specific.root";
    namefile_cut = "cut_p_tof.root";


    TTree *tx, *tx_C;
    TCutG *p_cut;

    TFile *file_CH2, *file_C, *file_cut;

    file_CH2 = new TFile(namefile_CH2.c_str(),"READ");
    file_C = new TFile(namefile_C.c_str(),"READ");
    file_cut = new TFile(namefile_cut.c_str(),"READ");

    file_cut->cd();
    p_cut = (TCutG*)file_cut->Get("cut_p_tof");

    //p_cut->Draw();

    file_CH2->cd();
    tx = (TTree*)file_CH2->Get("M");
    cout<<"Entries CH2: "<<tx->GetEntries()<<"."<<endl;

    file_C->cd();
    tx_C = (TTree*)file_C->Get("M");
    cout<<"Entries C: "<<tx_C->GetEntries()<<"."<<endl;
 //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//selection of events

    TCanvas *TCutsCanvas = new TCanvas("TCutsCanvas","Getting rid of ghosts...", 160,80,1600,550);
    TCutsCanvas->Divide(2,1);
    TCutsCanvas->cd(1);
    
    gPad->SetGridx();
    gPad->SetGridy();
    
    tx->Draw("si2:tofn>>h(800,-200,600,500,0,16)","PID==1","colz");
    p_cut->Draw("same");


    TCutsCanvas->cd(2);

    gPad->SetGridx();
    gPad->SetGridy();

    tx_C->Draw("si2:tofn>>hC(800,-200,600,500,0,16)","PID==1","colz");
    p_cut->SetVarY("si2");
    p_cut->Draw("same");
 //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
// tof plots 


    TCanvas *TTofsCanvas = new TCanvas("TTofsCanvas","Taking a look at the ToF plots...", 260,180,800,800);
    TH1D *htofn = new TH1D("htofn","htofn",210,40,250);
    TH1D *htofn_C = new TH1D("htofn_C","htofn_C",210,40,250);
    
    //elastic scattered events tof histo:
    TH1D *htofEE = new TH1D("htofEE","htofEE",210,40,250);
    double subfac = sub_fac(Q_CH2,Q_C,0.0237,0.0408);
    cout<<"subtraction factor = "<<subfac<<endl;

    tx->Draw("tofn>>htofn","PID==1 && cut_p_tof");
    tx_C->Draw("tofn>>htofn_C","PID==1 && cut_p_tof","same");

    for(Int_t b=1;b<=htofn->GetNbinsX();b++){
        htofEE->SetBinContent(b,htofn->GetBinContent(b) - subfac*htofn_C->GetBinContent(b));
    }

    htofEE->SetLineWidth(2);
    htofEE->SetLineColor(kRed);
    htofEE->Draw("same");

    gPad->SetGridx();
    gPad->SetGridy();

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    TH1D *hn = NeutronFromToF(htofEE,true);
    
    //creating a new histo multiplying the X axis by 1/cos²(20°)

    TH1D *hneutrons = ScaleXhisto(hn,pow(1.0/cos(TMath::Pi()*20/180),2),"hneutrons");


    TCanvas *NNCanvas = new TCanvas("NNCanvas","Neutrons...", 360,380,800,800);

    //load XS
    TGraph *xs = new TGraph("/home/e802/Analysis/LucasAnalysis/fluxMedley_script/XS/nn.org_np_2to40MeV_LABrf_20.40.60.80deg.csv","%lg %lg %*lg %*lg %*lg",",");


    //stuff for multiplicative factor: 
        //number of H atoms
    double NatH = 2*0.0237*6.02e23/14.0; //2 ⋅ [mass of CH2 ⋅ NA / MM(CH2)]
    double L = 464.72; //cm
    double cm2_per_barn = 1e-24 ;//cm²
    double omega_tel = 0.040 ;//sr 
    double omega_tgt = 0.00002272935814;// sr
    

    //multiplicative factor
    //double factor_mult =sin(TMath::Pi()*20/180)*pow(L,2)*barn_in_cm2/NatH/omega_tel/Q;
    double factor_mult = TMath::Pi()*2*sin(TMath::Pi()*20/180)*pow(L,2)/(pow(cos(TMath::Pi()*20/180),2)*NatH*omega_tel*Q_CH2*cm2_per_barn );

    double bincontent, bincounts, bincenter,bin_width, xs_value ;
    for(Int_t b=1;b<=hneutrons->GetNbinsX();b++){

        bincounts = hneutrons->GetBinContent(b);
        bincenter = hneutrons->GetBinCenter(b);

        xs_value = xs->Eval(bincenter)/1e3;// in b
        bin_width = hneutrons->GetBinWidth(b); // in MeV
        
        bincontent = factor_mult*bincounts/(xs_value*bin_width); //1e3 because xs is given in mb
        
        hneutrons->SetBinContent(b,bincontent); 
    }
    hneutrons->GetYaxis()->SetTitle("");
    hneutrons->Draw();



//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    timer.Print();
    return 0;
}