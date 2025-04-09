//calculation of scattered protons histo.
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
#include "/home/e802/Analysis/pre_analysis/libs/useful.h"


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

void scattered_protons(double ang = 20.0){
Int_t NumberOfBins = 400;

TStopwatch timer;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//defining paths and file names
string namefilePol, namefileCarbon, namefilePPACs, path_runs;


path_runs = "/home/e802/Analysis/LucasAnalysis/reducedData";


namefilePPACs = "run111.root";

// //pour 407/406:
// int runC=406;
// int runPol=407;
// namefilePol  = "407.root";
// namefileCarbon  = "406.root";
// // const double chargeC = 0.256*1e6; //µC
// // const double chargePol = 0.444*1e6; //µC
// cout<<"Processing for runs 407 (CH2) and 406 (C)..."<<endl;
// const double chargeC = 0.256*1e6*0.558404; //µC
// const double chargePol = 0.444*1e6*0.712411; //µC

//pour 370/388:
int runC=388;
int runPol=370;
namefilePol  = "370.root";
namefileCarbon  = "388.root";
const double chargeC = 0.305*1e6; //µC
const double chargePol = 0.239*1e6; //µC
cout<<"Processing for runs 370 (CH2) and 388 (C)..."<<endl;


// //declaring files
TFile *filePol, *fileCarbon;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Load the cut for ToF -- to remove the ghosts of bad HF

gROOT->ProcessLine(".L cut_ProtonsToF.C");
TCutG * cut_tof = (TCutG*)gROOT->GetListOfSpecials()->FindObject("cut_ProtonsToF");

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


// // // opening and processing Pol file
TCanvas *cpol = new TCanvas("cpol","pol canvas protons",100,100,800,600);
cpol->SetLeftMargin(0.14);
gPad->SetGridx();
gPad->SetGridy();

filePol = new TFile(Form("%s/%s",path_runs.c_str(),namefilePol.c_str()),"READ");
TTree *trpol = (TTree *)filePol->Get("M");

TH1D *protons_pol = new TH1D("protons_pol","protons_pol",NumberOfBins,0,40);
protons_pol->SetLineColor(kRed);
trpol->Draw("energy>>protons_pol", Form("PID==1 && ang ==%f  && cut_ProtonsToF",ang));
protons_pol->GetXaxis()->SetTitle("E_{NN} (MeV)");
protons_pol->GetYaxis()->SetTitle("counts");

// // // opening and processing Carbon file
fileCarbon = new TFile(Form("%s/%s",path_runs.c_str(),namefileCarbon.c_str()),"READ");
TTree *trcarbon = (TTree *)fileCarbon->Get("M");

TH1D *protons_carbon = new TH1D("protons_carbon","protons_carbon",NumberOfBins,0,40);
protons_carbon->SetLineColor(kBlue);
trcarbon->Draw("energy>>protons_carbon", Form("PID==1 && ang ==%f  && cut_ProtonsToF",ang),"same");

// //subtraction factor

double subfac = number_C_atoms_Pol*chargePol/(number_C_atoms_C*chargeC);

cout<<"\n::\n::\n:: The computed normatization factor (for C) is "<<subfac<<"\n::\n::"<<endl;
cout<<"\n -> Remember the operation is: CH2 - <Factor>*C, where <Factor> is show above.\n"<<endl;


TH1D *protons_spec = (TH1D*)protons_pol->Clone();
protons_spec->Add(protons_carbon,-subfac);
protons_spec->SetNameTitle("protons_spec","protons_spec");
protons_spec->SetLineColor(kBlack);
protons_spec->Draw("same");

TLegend *leg_protons = new TLegend(0.454887,0.644097,0.880952,0.885417);
leg_protons->AddEntry(protons_pol,"neutrons from CH2", "l");
leg_protons->AddEntry(protons_carbon,"neutrons from C", "l");
leg_protons->AddEntry(protons_spec,"resulting neutrons (excluding C)", "l");
leg_protons->Draw();

protons_spec->SaveAs(Form("scattered_protons_spectrum_runs%d-%d_%04.1fdeg.root",runPol,runC,ang));


timer.Print("m");


}
