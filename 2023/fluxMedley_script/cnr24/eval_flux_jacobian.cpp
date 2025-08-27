/////////////////
/////// Script using for testing the tof transform with jacobian  -- version 2025

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

#include "/mnt/medley/LucasAnalysis/useful.h"


//IBM palette
const Int_t nColors = 5;
Int_t colors[nColors] = {
        TColor::GetColor("#648FFF"),  // Blue
        TColor::GetColor("#785EF0"),  // Purple
        TColor::GetColor("#DC267F"),  // Salmon
        TColor::GetColor("#FE6100"),  // Orange
        TColor::GetColor("#FFB000"),  // Gold
};

//#include "functions.h"
//
//Function to calculate the subtraction factor::
double sub_fac(double Qa, double Qb, double ma, double mb){// A is CH2 and B is C

float mmch2 = 14.0; 
float mmc = 12.0; 

return (mmc/mmch2)*(Qa/Qb)*(ma/mb);
}


//run 370: 
     //integrateur:
    //     double Q_CH2 = 239266 ;//µC
    double Q_CH2 = 222542 ;//µC
//run 388: 
//integrateur    
      //  double Q_C =305744 ;//µC
    double Q_C = 292294 ;//µC


int eval_flux_jacobian(){

    TStopwatch timer;
    string namefile_CH2,namefile_C,namefile_cut;

    //namefile_CH2 = "run_407_specific.root";
    namefile_CH2 = "/mnt/medley/LucasAnalysis/2023/reducedv61/370.root";
    //namefile_C = "run_406_specific.root";
    namefile_C = "/mnt/medley/LucasAnalysis/2023/reducedv61/388.root";
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


 tx->SetAlias("myprotons", "PID==1 && ang == 20");
 tx_C->SetAlias("myprotons", "PID==1 && ang == 20");
 //selection of events

    TCanvas *TCutsCanvas = new TCanvas("TCutsCanvas","Getting rid of ghosts...", 160,80,1600,550);
    TCutsCanvas->Divide(2,1);
    TCutsCanvas->cd(1);
    
    gPad->SetGridx();
    gPad->SetGridy();
    

    tx->Draw("si2:tofn>>h(800,-200,600,500,0,16)","myprotons","colz");
    p_cut->Draw("same");


    TCutsCanvas->cd(2);

    gPad->SetGridx();
    gPad->SetGridy();

    tx_C->Draw("si2:tofn>>hC(800,-200,600,500,0,16)","myprotons","colz");
    p_cut->SetVarY("si2");
    p_cut->Draw("same");
 //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
// tof plots 


    TCanvas *TTofsCanvas = new TCanvas("TTofsCanvas","Taking a look at the ToF plots...", 260,180,1600,800);
    
    
    TH1D *htofn = new TH1D("htofn","htofn",210,40,250);
    TH1D *htofn_C = new TH1D("htofn_C","htofn_C",210,40,250);
    
    TH1D *htofn_unc = new TH1D("htofn_unc","htofn_unc",210,40,250);
    TH1D *htofn_C_unc = new TH1D("htofn_C_unc","htofn_C_unc",210,40,250);
    

    //elastic scattered events tof histo:
    TH1D *htofEE = new TH1D("htofEE","htofEE",210,40,250);
    TH1D *htofEE_unc = new TH1D("htofEE_unc","htofEE_unc",210,40,250);
    double subfac = sub_fac(Q_CH2,Q_C,0.0237,0.0408);
    cout<<"subtraction factor = "<<subfac<<endl;

    TTofsCanvas->Divide(2,1);

    TTofsCanvas->cd(1);
    tx->Draw("tofn>>htofn","myprotons && cut_p_tof");
    tx_C->Draw("tofn>>htofn_C","myprotons && cut_p_tof","same");

    TTofsCanvas->cd(2);
    tx->Draw("tofn - tof_correction>>htofn_unc","myprotons && cut_p_tof");
    tx_C->Draw("tofn- tof_correction>>htofn_C_unc","myprotons && cut_p_tof","same");

    for(Int_t b=1;b<=htofn->GetNbinsX();b++){
        htofEE->SetBinContent(b,htofn->GetBinContent(b) - subfac*htofn_C->GetBinContent(b));
        htofEE_unc->SetBinContent(b,htofn_unc->GetBinContent(b) - subfac*htofn_C_unc->GetBinContent(b));
    }

    TTofsCanvas->cd(1);
    htofEE->SetLineWidth(2);
    htofEE->SetLineColor(kRed);
    htofEE->Draw("same");
    
    TTofsCanvas->cd(2);
    htofEE_unc->SetLineWidth(2);
    htofEE_unc->SetLineColor(kRed);
    htofEE_unc->Draw("same");

    gPad->SetGridx();
    gPad->SetGridy();

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    TH1D *hn = NeutronFromToFJacobian(htofEE,0);
    TH1D *hn_unc = NeutronFromToFJacobian(htofEE_unc,0);
    //I need to correct back to counts...


    //creating a new histo multiplying the X axis by 1/cos²(20°)

    TH1D *hneutrons = ScaleXhisto(hn,pow(1.0/cos(TMath::Pi()*20/180),2),"hneutrons");
    TH1D *hneutrons_unc = ScaleXhisto(hn_unc,pow(1.0/cos(TMath::Pi()*20/180),2),"hneutrons_unc");
    
    TCanvas *NNCanvas = new TCanvas("NNCanvas","Neutrons...", 360,380,790,652);
    //load XS
    TGraph *xs = new TGraph("/mnt/medley/LucasAnalysis/2023/fluxMedley_script/XS/nn.org_np_2to40MeV_LABrf_20.40.60.80deg.csv","%lg %lg %*lg %*lg %*lg",",");


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
    double bincontent_unc, bincounts_unc, bincenter_unc,bin_width_unc, xs_value_unc ;
    for(Int_t b=1;b<=hneutrons->GetNbinsX();b++){

        bincounts = hneutrons->GetBinContent(b);
        bincenter = hneutrons->GetBinCenter(b);

        xs_value = xs->Eval(bincenter)/1e3;// in b
        bin_width = hneutrons->GetBinWidth(b); // in MeV
        
        bincontent = factor_mult*bincounts/(xs_value*bin_width); //1e3 because xs is given in mb
        
        hneutrons->SetBinContent(b,bincontent); 

        //-------------------------------------------------------

        bincounts_unc = hneutrons_unc->GetBinContent(b);
        bincenter_unc = hneutrons_unc->GetBinCenter(b);

        xs_value_unc = xs->Eval(bincenter_unc)/1e3;// in b
        bin_width_unc = hneutrons_unc->GetBinWidth(b); // in MeV
        
        bincontent_unc = factor_mult*bincounts_unc/(xs_value_unc*bin_width_unc); //1e3 because xs is given in mb
        
        hneutrons_unc->SetBinContent(b,bincontent_unc); 

    }
    //hneutrons->GetYaxis()->SetTitle("");
    
    gStyle->SetOptStat(0);
    gStyle->SetTextSize(0.04);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYSize(0.04);
    
    NNCanvas->cd()->SetRightMargin(0.03);
    NNCanvas->cd()->SetTopMargin(0.08);
    NNCanvas->cd()->SetLeftMargin(0.14);
    NNCanvas->cd()->SetBottomMargin(0.12);
    
    hneutrons_unc->GetXaxis()->SetRangeUser(0.,45.);
    hneutrons_unc->GetYaxis()->SetTitle("n  sr^{-1} 1-MeV^{-1} #muC^{-1}");
    hneutrons_unc->GetXaxis()->SetTitle("E (MeV)");
    hneutrons_unc->SetTitle("");
    
    hneutrons_unc->GetXaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetYaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetXaxis()->SetTitleSize(0.06);
    hneutrons_unc->GetYaxis()->SetTitleSize(0.06);
    
    hneutrons->SetFillColor(colors[0]);
    hneutrons->SetFillStyle(3005);
    hneutrons->SetLineColor(colors[0]);
    hneutrons_unc->SetLineColor(colors[2]);

    hneutrons->SetLineWidth(2);
    hneutrons_unc->SetLineWidth(2);

    hneutrons_unc->Draw();
    hneutrons->Draw("same");

    gPad->SetGridx();
    gPad->SetGridy();

    TFile *f = new TFile("run111.root","READ");
    TGraphErrors *nspec = (TGraphErrors *)f->Get("Nspectrum");
    nspec->SetLineWidth(2);
    nspec->Draw("same");

    TLegend *tl = new TLegend(0.481663,0.55574,0.96,0.888519);
    //tl->AddEntry(hneutrons_unc,"uncorrected TOF", "lpf");
    //tl->AddEntry(hneutrons,"corrected TOF", "lpf");
    tl->SetTextSize(0.05);
    tl->AddEntry(hneutrons_unc,"TOF without correction", "lpf");
    tl->AddEntry(hneutrons,"TOF with correction", "lpf");
    tl->AddEntry(nspec,"PPACs", "lpf");
    tl->Draw();
    
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    timer.Print();
    return 0;
}