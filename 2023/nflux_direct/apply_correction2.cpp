#include "KVMaterial.h"
#include "KVUnits.h"
#include "TF1.h"
#include "TRandom2.h"
#include <TMath.h>
#include "TVector3.h"
#include "TMatrixT.h"
#include "/home/pi/ganil/libraries/LArruda.h"
#include "TStopwatch.h"
#include "TFile.h"
#include <iostream>
//#include "/mnt/medley/LucasAnalysis/2023/nflux_direct"

string ganil_folder= "/home/pi/ganil/"; 


// function to give the correction for dead particles 
TGraph *DeadParticlesCorrection(
  TTree* elossT = nullptr,  
  string pSpecFileName = "curve.root",
  bool savecurve = false
){
  TH1D *h = new TH1D("h","h",400,0,40);
  //F = #produced/ (#produced - #dead)
  //elossT->Draw("")



  return nullptr;
}


void apply_correction2(
  string pSpecFileName = "protonsCH2.root",
  string pElossFileName = Form("%s/kalscripts/eloss/results/UniformZ/CH2/v7/eloss_p_20deg_050.0um.root",ganil_folder.c_str()),  
  string target_mat = "CH2",
  Double_t th = 50
  ){


TStopwatch timer;

char particle = 'p';
float ang = 20;

TFile *fEloss = new TFile(pElossFileName.c_str(),"READ"); 

// Check if the file was successfully opened
if (!fEloss || fEloss->IsZombie()) {
    std::cerr << "Warning: The file " << pElossFileName << " does not exist or could not be opened." << std::endl;
    if (fEloss) fEloss->Close(); // Close the file if it is open but in a bad state
    delete fEloss; // Clean up
    return; // Exit the function if the file does not exist or could not be opened
}

//create the TTree to the loss file 
TTree *ElTree = (TTree*)fEloss->Get("SIM");


//defining the response function:
Int_t nbins = 200;
TH2D *h= new TH2D("h","h",nbins, 0,40,nbins,0,40);
TCanvas *c0 = new TCanvas("response_func_canvas","response_func_canvas",50,50,800,600);
c0->SetRightMargin(0.16);
ElTree->Draw("E:Erem>>h(200,0,40,200,0,40)","","goff");
gStyle->SetOptStat("e");
h->SetTitle(Form("protons leaving %s target, %4.1f#mum",target_mat.c_str(),th));
h->GetXaxis()->SetTitle("Measured proton energy (MeV)");
h->GetYaxis()->SetTitle("Initial proton energy (MeV)");
h->Draw("colz");
h->Draw("colz");
gPad->SetGridx();
gPad->SetGridy();


//TH1D *projY= h->ProjectionY("df",16,16);

TCanvas *cc = new TCanvas();

TH1D *projY[nbins+1];
for (int b = 1; b <=nbins; b++)
{
    projY[b] =  h->ProjectionY(Form("hproj%d",b),b,b);
}

///////////////////////////////////////////////////////////////////
//get new energy
//First: create the random generator 
TRandom2 *rand = new TRandom2();
std::time_t seed = std::time(nullptr);
rand->SetSeed(seed);


//open experimental file 
TFile *expFile = new TFile(pSpecFileName.c_str(),"READ");
TH1D *expSpec = (TH1D*)expFile->Get("hnew");

/*
TCanvas *c1 = new TCanvas("imaginary_p_Exp","imaginary_p_Exp",100,200,800,600);
//Lets construct an imaginary exp protons spec (monoenergetic!)
c1->SetLeftMargin(0.14);
gStyle->SetOptStat(0);
expSpec->SetTitle(pSpecFileName.c_str());
//expSpec->Draw();
*/
TH1D *corrected_exp = new TH1D("hexp","hexp",400,0,40);

Int_t binProj=0;
for (Int_t bn = 1; bn <= expSpec->GetNbinsX(); bn++)
{
    //we need to find the profile corresponding to that bin value 
binProj = h->GetXaxis()->FindBin(expSpec->GetBinCenter(bn));
  for (Int_t i = 0; i < expSpec->GetBinContent(bn); i++)
  {
    corrected_exp->Fill(projY[binProj]->GetRandom(rand));
  }
  
}
corrected_exp->SetFillColor(kRed);
corrected_exp->SetFillStyle(3005);
corrected_exp->SetLineColor(kRed);
corrected_exp->SetLineWidth(2);
corrected_exp->Draw();

//expSpec->SetFillColor(kBlue);
expSpec->SetLineColor(kBlue);
expSpec->SetLineWidth(2);
expSpec->Draw("same");

TLegend *tl = new TLegend(0.590226,0.736111,0.889724,0.887153);

tl->AddEntry(expSpec,"Initial distribution","lf");
tl->AddEntry(corrected_exp,"Corrected distribution","lf");
tl->Draw();
gPad->SetGridx();
gPad->SetGridy();

string SaveFileName = lastnameNoExt(pSpecFileName)  + "_corrected.root";
string SaveFileNameImage = "figures/"+lastnameNoExt(pSpecFileName)  + "_corrected.png";
string SaveFileNameImageRF = "figures/"+lastnameNoExt(pSpecFileName)  + "RF_corrected.png";
corrected_exp->SaveAs(SaveFileName.c_str());

c1->SaveAs(SaveFileNameImage.c_str());
c0->SaveAs(SaveFileNameImageRF.c_str());




timer.Print();

}