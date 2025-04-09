//
// Script for reconstructing the neutron flux from the direct method 
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
#include "/home/pi/ganil/kalscripts/eloss/GetPPlot.cpp"


#include "TStopwatch.h"
#include <time.h>

TH1D *pc, *pch2;

Double_t qC = 305744;//µC run 388
Double_t qCH2 = 222542;//µC run 370


Double_t mch2 = 0.0237;//g 
Double_t mc = 0.0408;//g 


Double_t m_n = 939.56542052;
Double_t m_p = 938.27208816;

//small function to calculate subtraction factor
double sub_fac(double Qch2, double Qc, double mch2, double mc){
// A is CH2 and B is C

float mmch2 = 14.0; 
float mmc = 12.0; 

return (Qch2/Qc)* ((mch2/mmch2) / (mc/mmc));
}


TH1D* flux_reconstruction(string ch2file = "protonsCH2_corrected.root", string cfile = "protonsC_corrected.root", string name = "hexp"){
//TH1D* flux_reconstruction(string ch2file = "protonsCH2.root", string cfile = "protonsC.root", string name = "hexp"){
TFile *ff = new TFile(ch2file.c_str(), "READ");

TCanvas *c1 = new TCanvas("CH2_cv","CH2_cv");

ff->Print();
pch2 =  (TH1D*)ff->Get(name.c_str());
pch2->SetNameTitle("hch2","hch2");
pch2->Draw();

TCanvas *c2 = new TCanvas("C_cv","C_cv");
ff = new TFile(cfile.c_str(), "READ");
ff->Print();
pc =  (TH1D*)ff->Get(name.c_str());
pc->SetNameTitle("hc","hc");
pc->Draw();



TH1D *hh = (TH1D*)pch2->Clone();
hh->SetNameTitle("hh","hh");
double subfac = sub_fac(qCH2,qC,mch2,mc);
cout<<"subtraction facto = "<< subfac<<endl;
hh->Add(pc,-subfac);

c1->cd();
hh->SetLineColor(kBlack);

double facX = pow(m_n + m_p,2)/(4*m_n*m_p*pow(cos(TMath::Pi()*20/180),2));

hh = ScaleXhisto(hh, facX,"hh_scaled");
hh->Draw("same");

//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.
// calculating N flux

TGraph *xs = new TGraph("/mnt/medley/LucasAnalysis/cross_secs/XS/nn.org_np_2to40MeV_LABrf_20.40.60.80deg.csv","%lg %lg %*lg %*lg %*lg",",");
double NA = 6.023e23;
double NatH = 2*mch2*NA/14.0; // [mass of CH2 ⋅ NA / MM(CH2)] //TAKE A LOOK HERE
double L = 464.72; //cm
double cm2_per_barn = 1e-24 ;//cm²
double omega_tel = 0.040 ;//sr 
double omega_tgt = 0.00002272935814;// sr

double factor_mult = TMath::Pi()*2*sin(TMath::Pi()*20/180)*pow(L,2)/(NatH*omega_tel*qCH2*cm2_per_barn );

TH1D *nn = (TH1D*)hh->Clone();
nn->SetNameTitle("nn","nn");

double bincontent, bincounts, bincenter,bin_width, xs_value ;

for(Int_t b=1;b<=nn->GetNbinsX();b++){
    bincounts = nn->GetBinContent(b);
    bincenter = nn->GetBinCenter(b);

    xs_value = xs->Eval(bincenter)/1e3;// in b
    bin_width = nn->GetBinWidth(b); // in MeV

    bincontent = factor_mult*bincounts/(xs_value*bin_width); //1e3 because xs is given in mb
    nn->SetBinContent(b,bincontent); 
}



TGraph *gF = GetPfactor();
TCanvas *cF = new TCanvas("Pfactor","Pfactor");
gF->Draw("ALP");

TCanvas *c3 = new TCanvas("Nflux_cv","Nflux_cv");


nn->Draw();

TH1D *nn2 = (TH1D*)nn->Clone();
nn2->SetNameTitle("nn2","nn2");

for(Int_t b=1;b<=nn2->GetNbinsX();b++){
    bincounts = nn2->GetBinContent(b);
    bincenter = nn2->GetBinCenter(b);

    bincontent = bincounts*gF->Eval(bincenter); 
    nn2->SetBinContent(b,bincontent); 

    
}
nn2->SetLineColor(kRed);
nn2->Draw("same");



TFile *fflux = new TFile("run111.root","READ");
TGraphErrors *nspec = (TGraphErrors *)fflux->Get("Nspectrum");
nspec->SetLineWidth(2);

nspec->Draw("same");




return nn;
}