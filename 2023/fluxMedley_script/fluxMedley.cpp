
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

//TGraph * XS_fecther(int telN = 1,bool drawXS = false,bool verbose = false, string xsfile = "/home/e802/Analysis/pre_analysis/cross_secs/H1_MT2_angXS.csv"){

void fluxMedley(int telN = 1 ){

TStopwatch timer;


//loading XS 
TGraph *xs = XS_fecther(telN,1,0,"endf_B-VI.8_20deg.40deg.60deg.80deg.csv");


string namefilePol, namefileCarbon, namefilePPACs, path_runs;

namefilePPACs = "run111.root";
namefilePol  = "407.root";
namefileCarbon  = "406.root";

// namefilePol  = "370.root";
// namefileCarbon  = "388.root";


path_runs = "/home/e802/Analysis/LucasAnalysis/reducedData";

TFile* filePPACs, *filePol, *fileCarbon;



//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//CALCULATIONS FOR INTRODUCTION XS 

//number of protons detected for given energy:
// #(E) = Φ(E) * dσ/dΩ(E) * ΔΩ * Nat_H * Qrun

// The number of scattering nucleus Nat_H = (2mPol/MM)*NA, and dσ/dΩ(E) is tipically given in b/sr so we need a conversion factor [cm²/b] 
// note that NA = 0.602E+24 and b = 1E-24 cm², so this means we can switch NA*[cm²/b] = 0.602 

//using it we get: 

// #(E) = Φ * dσ/dΩ(E) * ΔΩ * (mpol/mm)* 1.204 * Qrun

//where mpol is the mass of CH2 and mm is its molar mass

//I will call Rm = (mpol/mm) and Cte = 1/(  ΔΩ * Rm * 1.204 * Qrun ), then 

//Φ(E) = #(E) * Cte /[dσ/dΩ(E)]

double mPol = 0.0237;//in grams
double QPol = 0.444*1E6;
double QC = 0.256*1E6;
// double QPol = 0.239*1E6;
//  double QC = 0.306*1E6;
double Rm = (mPol/14.0);//I consider mm[CH2] = 14g/mol
double omega_telescope = 0.0182538;// ± 0.0000498 sr,  calculated with SACALC3v14
double Cte = 1/(omega_telescope * Rm *1.204 * QPol);

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//  // //  // //  // //  // //  // //
//First: prepare PPACs measurement


filePPACs = new TFile(namefilePPACs.c_str(),"READ");
TCanvas *cppacs = new TCanvas("cppacs","cppacs",500,300,800,600);
cppacs->SetLeftMargin(0.16);

//retrieve TGraphErros from file
TGraphErrors *PPACsFlux = (TGraphErrors*)filePPACs->Get("Nspectrum");

//calculate 1/D² 
double D2 = 1/pow(L_Be_Medley,2);

//multiply PPACs spectrum for 1/D² so now it will be given in /cm² 
for(int t=0;t<PPACsFlux->GetN();t++)
{
   PPACsFlux->SetPointY(t, PPACsFlux->GetPointY(t)*D2);
    PPACsFlux->SetPointError(t, 0.0,PPACsFlux->GetErrorY(t)*D2);
}

//some plotting attributes 
PPACsFlux->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
PPACsFlux->SetLineColor(kRed);
PPACsFlux->SetMarkerColor(kRed);
PPACsFlux->Draw();
gPad->SetGridx();
gPad->SetGridy();
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Load the cut for ToF

gROOT->ProcessLine(".L cut_ProtonsToF.C");
TCutG * cut_tof = (TCutG*)gROOT->GetListOfSpecials()->FindObject("cut_ProtonsToF");

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
TCanvas *cpol = new TCanvas("cpol","pol canvas protons",100,100,800,600);

filePol = new TFile(Form("%s/%s",path_runs.c_str(),namefilePol.c_str()),"READ");
TTree *trpol = (TTree *)filePol->Get("M");

TH1D *protons_pol = new TH1D("protons_pol","protons_pol",400,0,40);
protons_pol->SetLineColor(kRed);
trpol->Draw("ENN>>protons_pol", Form("PID==1 && ang ==%f && ENN>2.5 && cut_ProtonsToF",telN*20.0),"colz");

//plotting this just to show without the tofcut 
trpol->Draw("ENN>>protons_pol0", Form("PID==1 && ang ==%f && ENN>2.5",telN*20.0),"colz same");

fileCarbon = new TFile(Form("%s/%s",path_runs.c_str(),namefileCarbon.c_str()),"READ");
TTree *trcarbon = (TTree *)fileCarbon->Get("M");

TH1D *protons_carbon = new TH1D("protons_carbon","protons_carbon",400,0,40);
protons_carbon->SetLineColor(kRed);
trcarbon->Draw("ENN>>protons_carbon", Form("PID==1 && ang ==%f && ENN>2.5 && cut_ProtonsToF",telN*20.0),"colz same");
trcarbon->Draw("ENN>>protons_carbon0", Form("PID==1 && ang ==%f && ENN>2.5",telN*20.0),"colz same");


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //======================================== Computation of the normalization constant =================================================
    //target information
    double mass[2], thickness[2], MM[2], NC[2], intensity[2], NF;
        //MM is molar mass and NC, numero des cibles (number of spallating centers)


    mass[0] = 0.0408;//grams, for carbon;
    mass[1] = 0.0237;//grams, fr ch2;

    // thickness[0] = 0.075; //mm, for carbon
    // thickness[1] = 0.0237; //mm, for ch2

    MM[0] = 12; // g.mol^-1, for carbon
    MM[1] = 12+2; // g.mol^-1, for ch2

    NC[0] = mass[0]/MM[0]; //NA atoms
    NC[1] = mass[1]/MM[1]; //NA atoms

    intensity[0] = QC; //charge provided by integrator, for run 406
    intensity[1] = QPol;//charge provided by integrator, for run 407

    NF = (NC[1]/NC[0])*(intensity[1]/intensity[0]); //to be multiplied by the 0-th histo (carbon)

    cout<<"\n::\n::\n:: The computed normatization factor (for C) is "<<NF<<"\n::\n::"<<endl;
    cout<<"\n -> Remember the operation is: CH2 - <Factor>*C, where <Factor> is show above.\n"<<endl;
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


TCanvas *cTofs = new TCanvas("cTofcs", "TOF plots (correct procedure)", 630,330,800,600);

TH1D *htof_pol = new TH1D("htof_pol","htof_pol",400, 0 ,400);
TH1D *htof_carbon = new TH1D("htof_carbon","htof_carbon",400, 0 ,400);
TH1D *htof_np;


trpol->Draw("tofn>>htof_pol", Form("PID==1 && ang ==%f && ENN>2.0 && cut_ProtonsToF",telN*20.0));
trcarbon->Draw("tofn>>htof_carbon", Form("PID==1 && ang ==%f && ENN>2.0 && cut_ProtonsToF",telN*20.0),"same");

htof_np = (TH1D*)htof_pol->Clone();
htof_np->SetNameTitle("htof_np","htof_np");

htof_np->Add(htof_carbon,-NF);

htof_np->SetLineColor(kRed);
htof_np->Draw("same");
//transforming
TH1D *nflux_np = NeutronFromToF(htof_np, true);//verbose put some info regarding the function working 


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


TH1D *fluxpol = (TH1D*)protons_pol->Clone();
fluxpol->SetNameTitle("fluxpol","fluxpol");
fluxpol->SetLineColor(kBlack);




    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //printing some info 
    cout<<"flux Pol stats: "<<endl;
    //I print some info to check it is not calculating stuff wrong:
    //number of bins:
    cout<<"bins: "<<fluxpol->GetNbinsX();
    //covering energies:
    cout<<", from "<<fluxpol->GetBinLowEdge(1)<<" MeV to "<<fluxpol->GetBinLowEdge(fluxpol->GetNbinsX()+1)<<"MeV"<<endl;

    cout<<"proton_carbon histo stats: "<<endl;
    //I print some info to check it is not calculating stuff wrong:
    //number of bins:
    cout<<"bins: "<<protons_carbon->GetNbinsX();
    //covering energies:
    cout<<", from "<<protons_carbon->GetBinLowEdge(1)<<" MeV to "<<protons_carbon->GetBinLowEdge(protons_carbon->GetNbinsX()+1)<<"MeV"<<endl;
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//subtraction C

fluxpol->Add(protons_carbon,-NF);

//uncomment line below for debugging
// cout<<Form("bin # |  center (MeV) |   value  |  xs(b/sr)")<<endl;

//some variables 
double xscalc, bw, value ;


float transformation_factor = 2*TMath::Pi()*sin(telN*20*TMath::Pi()/180);
// float calcTF[] = {3.885,3.567,3.176,2.896};
// float transformation_factor = calcTF[telN-1];
cout<<"transformation_factor = "<<transformation_factor<<endl;

//for each bin in the neutron energy measure histogram:
for(int b=1; b<= fluxpol->GetNbinsX(); b++){
    //Get the corresponding XS: 
    xscalc = transformation_factor*xs->Eval(fluxpol->GetBinCenter(b));

    //get the bin-width from the neutrons' energy histo
    bw = fluxpol->GetBinWidth(b);

    //calculate: Φ(E) = Cte * #(E) /[dσ/dΩ(E)]
    //here I divided to binwidth to be comparible with the other results
    value = Cte*fluxpol->GetBinContent(b)/(xscalc*bw);

    //fill the value in fluxpol 
    fluxpol->SetBinContent(b,value);
}


gStyle->SetOptStat(0);

TCanvas *cFlux = new TCanvas("cFlux", "Flux", 630,330,800,600);
//Axis labels and drawing 
 fluxpol->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
 fluxpol->GetXaxis()->SetTitle("E_{NN} (MeV)");
 fluxpol->GetYaxis()->SetMaxDigits(3);
 fluxpol->SetLineWidth(2);
 fluxpol->SetFillStyle(3005);
 fluxpol->SetFillColor(kBlack);
 fluxpol->SetLineColor(kBlack);
 fluxpol->Draw();



PPACsFlux->Draw("same");

TLegend *tl = new TLegend(0.533835, 0.730903,0.899749,0.901042);
tl->AddEntry(fluxpol,"flux from np", "fl");
tl->AddEntry(PPACsFlux,"flux PPACs", "fl");
tl->Draw();

gPad->SetGridx();
gPad->SetGridy();

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//getting flux for transforming tof-->Energy after subtraction: 
//this would be the correct procedure


timer.Print("m");


}