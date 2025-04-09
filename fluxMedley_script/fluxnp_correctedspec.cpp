
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


void fluxnp_correctedspec(double ang = 20.0){
Int_t NumberOfBins = 400;

TStopwatch timer;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//defining paths and file names
string namefilePol, namefileCarbon, namefilePPACs;

namefilePPACs = "run111.root";

//pour 407/406:
// int runC=406;
// int runPol=407;
// namefilePol  = "protons_Tel1_run407_corrected.root";
// namefileCarbon  = "protons_Tel1_run406_corrected.root";
// const double chargeC = 0.256*1e6; //µC
// const double chargePol = 0.444*1e6; //µC
cout<<"Processing for ELOSS CORRECTED runs 407 (CH2) and 406 (C)..."<<endl;

//pour 370/388:
int runC=388;
int runPol=370;
// namefilePol  = "protons_Tel1_run370_corrected.root";
// namefileCarbon  = "protons_Tel1_run388_corrected.root";
namefilePol  = "altPID_pspecs/370_corrected.root";
namefileCarbon  = "altPID_pspecs/388_corrected.root";
const double chargeC = 0.305*1e6; //µC
const double chargePol = 0.239*1e6; //µC
cout<<"Processing for ELOSS CORRECTED runs 370 (CH2) and 388 (C)..."<<endl;


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


// // // opening and processing Pol file
TCanvas *cpol = new TCanvas("cpol","pol canvas protons",100,100,800,600);
cpol->SetLeftMargin(0.14);
gPad->SetGridx();
gPad->SetGridy();

filePol = new TFile(namefilePol.c_str(),"READ");

TH1D *protons_pol = (TH1D*)filePol->Get("hexp");
protons_pol->SetNameTitle("protons_pol","eloss corrected protons for Pol");
protons_pol->SetLineColor(kRed);
protons_pol->GetXaxis()->SetTitle("E_{protons} (MeV)");
protons_pol->GetYaxis()->SetTitle("counts");
protons_pol->SetFillStyle(0);
protons_pol->Draw();

// // // opening and processing Carbon file
fileCarbon = new TFile(namefileCarbon.c_str(),"READ");

TH1D *protons_carbon = (TH1D*)fileCarbon->Get("hexp");
protons_carbon->SetLineColor(kBlue);
protons_carbon->SetFillStyle(0);
protons_carbon->SetNameTitle("protons_carbon","eloss corrected protons for Carbon");
protons_carbon->Draw("same");
// //subtraction factor

double subfac = number_C_atoms_Pol*chargePol/(number_C_atoms_C*chargeC);

cout<<"\n::\n::\n:: The computed normatization factor (for C) is "<<subfac<<"\n::\n::"<<endl;
cout<<"\n -> Remember the operation is: CH2 - <Factor>*C, where <Factor> is show above.\n"<<endl;


TH1D *protons_spec = (TH1D*)protons_pol->Clone();
protons_spec->Add(protons_carbon,-subfac);
protons_spec->SetNameTitle("protons_spec","protons_spec");
protons_spec->SetLineColor(kBlack);
protons_spec->SetFillStyle(3004);
protons_spec->SetFillColor(kBlack);
protons_spec->Draw("same");

TLegend *leg_protons = new TLegend(0.454887,0.644097,0.880952,0.885417);
leg_protons->AddEntry(protons_pol,"protons from CH2", "l");
leg_protons->AddEntry(protons_carbon,"protons from C", "l");
leg_protons->AddEntry(protons_spec,"EE protons (excluding C)", "l");
leg_protons->Draw();

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

Double_t co2 = pow(cos(ang*TMath::Pi()/180),2);
//Double_t ee_factor = 4*mp*mn/pow(mp+mn,2);

//creating neutrons histogram
int nbins = protons_spec->GetNbinsX();
Double_t Xp[nbins+1];
Double_t Np[nbins]; //number of protons for each bin
Long64_t totalCounts = 0;
for(int i=1;i<=nbins+1;i++)
{
    Xp[i-1] = protons_spec->GetBinLowEdge(i)*1./co2; //vector with bin edges
    Np[i-1] = abs(protons_spec->GetBinContent(i)); //counts for each bin

		//cout<<Form("bin %04d [%04.1f,%04.1f] -->[%04.1f,%04.1f] : %08.0f", i,protons_spec->GetBinLowEdge(i),protons_spec->GetBinLowEdge(i+1),protons_spec->GetBinLowEdge(i)*1./co2,protons_spec->GetBinLowEdge(i+1)*1./co2,protons_spec->GetBinContent(i))<<endl;
		//totalCounts += (Long64_t)protons_spec->GetBinContent(i);
}


//cout<<"totalCounts = "<<totalCounts<<endl;

//histogram with proper binning
TH1D *neutron_spectrum = new TH1D("neutrons_spectrum","neutrons_spectrum",nbins,Xp);
//completing the histo
Double_t binWidth, En, Fn, dsig;
for(int i=1;i<=nbins;i++)
{
    En = neutron_spectrum->GetBinCenter(i);
    binWidth = neutron_spectrum->GetBinWidth(i);
    dsig = sigma_diff(En,ang*TMath::Pi()/180.)*1e-3*1e-24;// in cm²
    
    //Fn = Np/(Nat*dsig*Omega)/hFn->GetBinWidth(i); 

    Fn = Np[i-1]/(number_H_atoms_Pol*dsig*omega_telescope*chargePol)/binWidth;
    neutron_spectrum->SetBinContent(i,Fn);

}

TCanvas *cn = new TCanvas("cn","neutron spectrum",500,150,800,600);
cn->SetLeftMargin(0.14);
gPad->SetGridx();
gPad->SetGridy();


gStyle->SetOptStat(0);
neutron_spectrum->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
neutron_spectrum->GetXaxis()->SetTitle("E_{NN} (MeV)");
neutron_spectrum->GetYaxis()->SetMaxDigits(3);
neutron_spectrum->SetTitle(Form("runs %d (CH_{2}) minus %d (C)", runPol, runC));
neutron_spectrum->Draw();
gStyle->SetOptStat(0);
PPACsFlux->Draw("same");


TLegend *leg_flux = new TLegend(0.63,0.7,0.88,0.88);
leg_flux->AddEntry(neutron_spectrum,"Medley", "l");
leg_flux->AddEntry(PPACsFlux,"PPACs", "l");
leg_flux->SetTextSize(0.045);
leg_flux->Draw();

gStyle->SetOptStat(0);
cn->SaveAs(Form("nspec_r%d-%d.root", runC,runPol));
cn->SaveAs(Form("figures/nspec_r%d-%d.png", runC,runPol));
cpol->SaveAs(Form("figures/protonsSpec_r%d-%d.png", runC,runPol));

timer.Print("m");


}
