//
// Script for evaluation the cross section for protons in Carbon__ 2023 runs
//
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TCutG.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <fstream>
//#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002
#include "/mnt/medley/LucasAnalysis/2023/applyTTC/correctSpec2.cpp" // for correcting the spectrum



#include "TStopwatch.h"
#include <time.h>

Double_t mc = 0.0408;//g 
Double_t M = 12; //g/mol
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Double_t omega_tel = 0.040 ;//sr 


Int_t Ta;

void C_px2(Double_t Ea = 26.0, Double_t Eb = 27.0 ,Int_t runa = 35,Int_t runb = 39){
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
cout << "nflux_En: " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;


TCanvas * hp_cv = new TCanvas("hp","hp");

TCutG *tofcut;
cout << Form("\n- - -   LOADING NP contamination cut   - - -\n")<< "... " ;

gROOT->ProcessLine(".L  /mnt/medley/LucasAnalysis/2023/reducedv61/goodtof.C");
tofcut = (TCutG *)gROOT->GetListOfSpecials()->FindObject("goodtof");
if(tofcut!=NULL) cout << "  --> OK."<< endl;

tofcut->Print();

TH1D *hp = new TH1D("hp","hp",100,0,40);

tx->Draw("E>>hp",Form("ENN>%f && ENN<%f && PID==1 && ang == 20 &!goodtof",Ea, Eb));

//calculate corrected spectrum
TH1D *hTTC = correctSpec(hp,true,"/home/pi/ganil/kalscripts/eloss/results/UniformZ/MedleyCarbon/v7/eloss_p_20.0deg_075.0um.root","MedleyCarbon",75,'p',false);
hTTC->Draw("same");

TCanvas * xs_cv = new TCanvas("XS","XS");

Double_t L = 464.72; //cm
Double_t A = TMath::Pi()*pow(2.5/2,2); //cm^2
Double_t cm2_to_barn = 1e+24 ;//cm²
Double_t Factor = L*L/(Nc*omega_tel*TotalChargeFaraday*nflux_En);
Factor = Factor*cm2_to_barn;
cout<<"Factor: "<<Factor<<endl;

TH1D *xsp = new TH1D("xsp","xsp",100,0,40);
TH1D *xspC = new TH1D("xspC","xspC",100,0,40);



for(int i=1;i<=hp->GetNbinsX();i++)
{
    xsp->SetBinContent(i,1e3*Factor*hp->GetBinContent(i)/hp->GetBinWidth(i));//1e3 changes to mb
    xspC->SetBinContent(i,1e3*Factor*hTTC->GetBinContent(i)/hTTC->GetBinWidth(i));//1e3 changes to mb
}
xsp->GetYaxis()->SetRangeUser(0,1.1*xsp->GetMaximum());
xsp->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
xsp->GetXaxis()->SetTitle("E_{protons} (MeV)");
//xsp->Draw();



TChain *xsdata = new TChain("DataTree");
xsdata->Add("/mnt/medley/LucasAnalysis/exfor/exfor_files/myformat/datasets_C.root");




TGraphErrors *gSlypen = new TGraphErrors();
TGraphErrors *gSubramanian = new TGraphErrors();
xspC->SetLineColor(kBlue);
xsp->SetLineWidth(2);
xsp->SetFillColor(kBlue);
xsp->SetFillStyle(3004);
xsp->Draw("same");


xspC->SetLineColor(kRed);
xspC->SetLineWidth(3);
xspC->SetFillColor(kRed);
xspC->SetFillStyle(3005);
xspC->Draw("same");


gPad->SetGridx();
#include "TCutG.h"
gPad->SetGridy();

TLegend *tl = new TLegend(0.54,0.63,0.85,0.85);
tl->AddEntry(xsp,"Medley (2023)","lpf");
tl->AddEntry(xspC,"Medley TTC (2023)","lpf");
// tl->AddEntry(graph,"Slypen+(2000), 20deg","lpf");
// tl->AddEntry(g40,"Slypen+(2000), 40deg","lpf");
// tl->AddEntry(g60,"Slypen+(2000), 60deg","lpf");
// tl->AddEntry(g80,"Slypen+(2000), 80deg","lpf");

tl->Draw();

timer.Print();
}

