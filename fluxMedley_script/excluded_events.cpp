
//reconstruction of neutron flux
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TLatex.h"
#include "TStopwatch.h"
#include <string>
#include <iostream>
#include <fstream>
#include "/home/e802/Analysis/LucasAnalysis/useful.h"


const double L_Be_Medley = 464.72; //cm
const double omega_telescope = 0.0182538;// ± 0.0000498 sr,  calculated with SACALC3v14

const double massePol = 0.0237; //g
const double masseC = 0.0408; //g

const double masse_molairePol = 14.0266; //g/mol
const double masse_molaireC = 12.011; //g/mol

const double NA = 6.02e23;

const double number_C_atoms_Pol = (massePol/masse_molairePol)*NA; 
const double number_C_atoms_C =  (masseC/masse_molaireC)*NA; 

const double number_H_atoms_Pol = 2*number_C_atoms_Pol; 

double mn =  939.56542194; //MeV/c² - CODATA https://physics.nist.gov/cgi-bin/cuu/Value?mnc2mev
double mp = 938.27208943; //MeV/c² - CODATA https://physics.nist.gov/cgi-bin/cuu/Value?mpc2mev


///////////////////////////////////////////////////////////////////////////////////////////
//  on cherche la section efficace doublement differentielle pour Ein et cos(theta)cm
Double_t sigma_diff(Double_t E,Double_t cos)

{

TFile *f=new TFile("/home/e802/Analysis/LucasAnalysis/fluxMedley_script/XS/dsigma_np_online.root","readonly");

//printf("E=%f   cos=%f\n",E,cos);
TH2F *hh=(TH2F *)f->Get("hh");
Int_t nx=hh->GetXaxis()->GetNbins();
Double_t x[nx];
hh->GetXaxis()->GetLowEdge(x);
//for(Int_t i=0;i<nx;i++){printf("x=%f\n",x[i]);}

Int_t ny=hh->GetYaxis()->GetNbins();
Double_t Ener[ny];
hh->GetYaxis()->GetLowEdge(Ener);

Int_t j0=0;
for(Int_t j=1;j<ny;j++){
	if(E>=Ener[j-1] &&E <Ener[j]) j0=j;   // on cherche la ligne y qui correspond à l'energie E
	if(E>=Ener[ny-1]) j0=ny-1;
}

// on cree le TGraph    dsig/dOmega = f(cos(theta)) pour energie E
Double_t y[nx];
for(Int_t i=0;i<nx;i++){
	y[i]=hh->GetBinContent(i+1,j0+1);
}
TGraph *g=new TGraph(nx,x,y);

Double_t sig=g->Eval(cos);
//printf("E=%f   J0=%i     sig=%f\n",E,j0,sig);
f->Close();
return sig;
}


void excluded_events(double ang = 20.0){
Int_t NumberOfBins = 400;

TStopwatch timer;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//defining paths and file names
string namefilePol, namefileCarbon, namefilePPACs, path_runs;


path_runs = "/home/e802/Analysis/LucasAnalysis/reducedData";


namefilePPACs = "run111.root";

//pour 407/406:
// namefilePol  = "407.root";
// namefileCarbon  = "406.root";
// const double chargeC = 0.256*1e6; //µC
// const double chargePol = 0.444*1e6; //µC
// cout<<"Processing for runs 407 (CH2) and 406 (C)..."<<endl;


// //pour 370/388:
namefilePol  = "370.root";
namefileCarbon  = "388.root";
const double chargeC = 0.305*1e6; //µC
const double chargePol = 0.239*1e6; //µC
cout<<"Processing for runs 370 (CH2) and 388 (C)..."<<endl;


// //declaring files
TFile* filePPACs, *filePol, *fileCarbon;


//extracting and preparing PPACs flux - - - - - - - - - - - - - - - -
filePPACs = new TFile(namefilePPACs.c_str(),"READ");
//getting ppacs flux
TGraphErrors *PPACsFlux = (TGraphErrors*)filePPACs->Get("Nspectrum");
double D2 = 1/pow(L_Be_Medley,2);

//multiply PPACs spectrum for 1/D² so now it will be given in /cm²
for(int t=0;t<PPACsFlux->GetN();t++)
{
   PPACsFlux->SetPointY(t, PPACsFlux->GetPointY(t)*D2);
    PPACsFlux->SetPointError(t, 0.0,PPACsFlux->GetErrorY(t)*D2);
}
PPACsFlux->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
PPACsFlux->SetLineColor(kRed);
PPACsFlux->SetMarkerColor(kRed);
// finished extracting and preparing PPACs flux - - - - - - - - - - - - - - - -
PPACsFlux->GetYaxis()->SetMaxDigits(3);
PPACsFlux->Draw();
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Load the cut for ToF -- to remove the ghosts of bad HF

gROOT->ProcessLine(".L cut_ProtonsToF.C");
TCutG * cut_tof = (TCutG*)gROOT->GetListOfSpecials()->FindObject("cut_ProtonsToF");
cut_tof->SetLineWidth(2);
cut_tof->SetLineColor(kRed);
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


// // // opening and processing Pol file
TCanvas *cpol = new TCanvas("cpol","pol canvas protons",100,100,1400,600);
cpol->SetLeftMargin(0.14);
cpol->Divide(2,1);
cpol->cd(1);
gPad->SetGridx();
gPad->SetGridy();

filePol = new TFile(Form("%s/%s",path_runs.c_str(),namefilePol.c_str()),"READ");
TTree *trpol = (TTree *)filePol->Get("M");

//TH1D *protons_pol = new TH1D("protons_pol","protons_pol",NumberOfBins,0,40);
TH2D *protons_pol = new TH2D("protons_pol","protons_pol",NumberOfBins,0,40,NumberOfBins,0,40);
protons_pol->SetLineColor(kRed);
trpol->Draw("energy:ENN>>protons_pol", Form("PID==1 && ang ==%f ",ang),"colz");
protons_pol->GetXaxis()->SetTitle("E_{NN} (MeV)");
protons_pol->GetYaxis()->SetTitle("counts");
cut_tof->Draw("same");


cpol->cd(2);
gPad->SetGridx();
gPad->SetGridy();

cout<<"Counts histo Polyethylene: "<<protons_pol->GetEntries()<<", where "<<IntegralCut(cut_tof,protons_pol)<<"("<<IntegralCut(cut_tof,protons_pol)*1.0/protons_pol->GetEntries()<<") are valid."<<endl;

// // // opening and processing Carbon file
fileCarbon = new TFile(Form("%s/%s",path_runs.c_str(),namefileCarbon.c_str()),"READ");
TTree *trcarbon = (TTree *)fileCarbon->Get("M");


//TH1D *protons_carbon = new TH1D("protons_carbon","protons_carbon",NumberOfBins,0,40);
TH2D *protons_carbon = new TH2D("protons_carbon","protons_carbon",NumberOfBins,0,40,NumberOfBins,0,40);
protons_carbon->SetLineColor(kBlue);
trcarbon->Draw("energy:ENN>>protons_carbon", Form("PID==1 && ang ==%f  ",ang),"colz");
cut_tof->Draw("same");

cout<<"Counts histo Carbon: "<<protons_carbon->GetEntries()<<", where "<<IntegralCut(cut_tof,protons_carbon)<<"("<<IntegralCut(cut_tof,protons_carbon)*1.0/protons_carbon->GetEntries()<<") are valid."<<endl;

timer.Print("m");


}
