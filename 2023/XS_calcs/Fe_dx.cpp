//
// Script for evaluation the cross section for protons in Carbon__ 2023 runs
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
#include "/mnt/medley/LucasAnalysis/2023/applyTTC/correctSpec2.cpp" // for correcting the spectrum
#include "/mnt/medley/LucasAnalysis/2023/Etot/newE.hh" // for the energy calibration


#include "TStopwatch.h"
#include <time.h>

Double_t mc = 0.0182;//g 
Double_t M = 55.845; //g/mol
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Double_t omega_tel = 0.040 ;//sr 


Int_t Ta;
//thin 389 395
void Fe_dx(Double_t Ea = 31, Double_t Eb = 32 ,Int_t runa = 396,Int_t runb = 405){
TStopwatch timer;


TChain *tx = new TChain("M");


string processed_runs = "/mnt/medley/LucasAnalysis/2023/reducedv61";

string name;
for(Int_t i=runa;i<=runb;i++)
{
    Bool_t fExist = true;
    name = Form("%s/%03d.root",processed_runs.c_str(),i);
    cout << name << endl;
    ifstream mfile;
    mfile.open(name);
    if(mfile)
    {
        mfile.close();
    }
    else
        fExist=false;
    if(fExist)
    {	  
        cout << "Adding " << name << endl;
        tx->Add(name.c_str());
        Ta = (Int_t) (tx->GetEntries());
        cout << "Entries " << Ta/1000 << "k "  << endl;
    }
    
}

cout<<"---------------------------------------------------"<<endl;
TTree *InfoTree = new TTree("InfoTree", "InfoTree");
InfoTree->ReadFile("/mnt/medley/LucasAnalysis/2023/runlist.csv", "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F");

Int_t RunN;
Int_t nbins = 100;
Float_t ChargeFaraday;
Float_t TotalChargeFaraday = 0;

InfoTree->SetBranchAddress("RunN", &RunN);
InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

for(int i=0;i<InfoTree->GetEntries();i++)
{
    InfoTree->GetEntry(i);
    //cout<< i<<", nrun = "<<RunN<< endl;
    if (RunN<= runb && RunN>=runa)
    {
        //InfoTree->GetEntry(i);
        cout << "RunN: " << RunN << " ChargeFaraday: " << ChargeFaraday << endl;
        TotalChargeFaraday += ChargeFaraday;
    }
}
cout<<"---------------------------------------------------"<<endl;
cout<< "Total Charge Faraday: " << TotalChargeFaraday <<" µC"<< endl;


TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", "READ");
TH1D *nflux = (TH1D*)ff->Get("nflux");
TCanvas * cflux = new TCanvas("cflux","cflux");
nflux->Draw();

Double_t nflux_En = nflux->GetBinContent(nflux->FindBin((Ea+Eb)/2.0));
cout << "nflux_En ("<<(Ea+Eb)/2.0<<" MeV): " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;


TCanvas * hp_cv = new TCanvas("hp","hp");
TH1D *hp = new TH1D("hp","hp",nbins,0,40);
TH1D *hp_si1 = new TH1D("hp_si1","hp_si1",10*nbins,0,40);


tx->Draw("E>>hp",Form("ENN>%f && ENN<%f && PID==2 && ang == 20",Ea, Eb));

tx->Draw("si1>>hp_si1",Form("ENN>%f && ENN<%f && PID==2 && ang == 20",Ea, Eb),"goff");
//TH1D *hpC = newE(hp_si1,"hpC",false);
TH1D *hpC = newE(hp_si1,"hpC",false,53.4,1,2);
hpC->SetLineColor(kRed);
hpC->SetLineWidth(2);
hpC->Rebin(10);
hpC->Draw("same");
//calculate corrected spectrum
TH1D *hTTC = correctSpec(hpC,true,"/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick/v7/eloss_d_20.0deg_025.0um.root","Fe_thick_Medley",25,'d',false);
hTTC->Draw("same");

TCanvas * xs_cv = new TCanvas("XS","XS");

Double_t L = 464.72; //cm
Double_t A = TMath::Pi()*pow(2.5/2,2); //cm^2
Double_t cm2_to_barn = 1e+24 ;//cm²
Double_t Factor = L*L/(Nc*omega_tel*TotalChargeFaraday*nflux_En);
Factor = Factor*cm2_to_barn;
cout<<"Factor: "<<Factor<<endl;

TH1D *xsp = new TH1D("xsp","xsp",nbins,0,40);
TH1D *xspC = new TH1D("xspC","xspC",nbins,0,40);



for(int i=1;i<=hpC->GetNbinsX();i++)
{
    xsp->SetBinContent(i,1e3*Factor*hpC->GetBinContent(i)/hpC->GetBinWidth(i));//1e3 changes to mb
    xspC->SetBinContent(i,1e3*Factor*hTTC->GetBinContent(i)/hTTC->GetBinWidth(i));//1e3 changes to mb
}
xsp->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
xsp->GetXaxis()->SetTitle("E_{deuterons} (MeV)");
xsp->Draw();

timer.Print();
}

