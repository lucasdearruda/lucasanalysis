
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TLatex.h"
#include <string>
#include <iostream>
#include <fstream>

void integrating_PPACs_flux(){

string namefilePPACs;

namefilePPACs = "/home/e802/Analysis/pre_analysis/ppac2023/run111.root";

TFile* filePPACs = new TFile(namefilePPACs.c_str(),"READ");
TCanvas *cppacs = new TCanvas("cppacs","cppacs",500,300,800,600);
cppacs->SetLeftMargin(0.16);
TGraphErrors *PPACsFlux = (TGraphErrors*)filePPACs->Get("Nspectrum");

double D2 = 1/pow(464.72,2);
for(int t=0;t<PPACsFlux->GetN();t++)
{
   PPACsFlux->SetPointY(t, PPACsFlux->GetPointY(t)*D2);
    PPACsFlux->SetPointError(t, 0.0,PPACsFlux->GetErrorY(t)*D2);
}
PPACsFlux->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
PPACsFlux->SetLineColor(kRed);
PPACsFlux->SetMarkerColor(kRed);
PPACsFlux->Draw();
gPad->SetGridx();
gPad->SetGridy();

//integrating: 
double x, y1,y2, x1,x2,xmed, dx, A, dA;


//integral value: 
A = 0; 

for(int t=0;t<PPACsFlux->GetN()-1;t++)
{
   y2 = PPACsFlux->GetPointY(t+1);
   y1 = PPACsFlux->GetPointY(t);
   
   x2 = PPACsFlux->GetPointX(t+1);
   x1 = PPACsFlux->GetPointX(t);
   
   dx = x2 - x1;


   dA = (y1+y2)*dx/2;

   A =A+dA;

}

cout<<"==================================================================================="<<endl;
cout<<"|                                     SUMMARY                                     |"<<endl;
cout<<"-----------------------------------------------------------------------------------"<<endl;
cout<<"..................................................................................."<<endl;
std::scientific;
cout<<"ϕ: "<<A<<" neutrons/µC/cm²."<<endl;
cout<<"for 2.5 cm diameter target: "<<A*TMath::Pi()*pow(2.5/2,2)<<" neutrons."<<endl;

}
