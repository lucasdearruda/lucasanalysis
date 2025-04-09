/*
////////////////////////////////////////////////////////
Function to correct spectrum, mode adapted to the new eloss files

*/
#include "KVMaterial.h"
#include "KVUnits.h"
#include "TF1.h"
#include "TRandom2.h"
#include <TMath.h>
#include "TVector3.h"
#include "TMatrixT.h"
//#include "/home/pi/ganil/libraries/LArruda.h"
#include "/mnt/medley/LucasAnalysis/useful.h"
#include "TStopwatch.h"
#include "TFile.h"
#include <iostream>
#include "/home/pi/ganil/kalscripts/eloss/GetPfactor.cpp"

string ganil_folder= "/home/pi/ganil/"; 



void splitPath(const std::string& fullPath, std::string& directory, std::string& filename) {
    size_t found = fullPath.find_last_of("/\\"); // Encontra última barra
    directory = fullPath.substr(0, found + 1);  // Inclui a barra final
    filename = fullPath.substr(found + 1);      // Parte após a barra
}

TGraph *Fcorr(TTree *S, Int_t nbins = 400, bool drawit = true){

TH1D *hprod = new TH1D("hprod","hprod",nbins,0,40);
TH1D *hlost = new TH1D("hlost","hlost",nbins,0,40);

S->Draw("E>>hprod","","goff");
S->Draw("E>>hlost","!transmitted","goff");

TGraph *Ffunc = new TGraph();
Double_t F; 
for(int i = 1; i <= nbins; i++){
    F = hprod->GetBinContent(i)/(hprod->GetBinContent(i) - hlost->GetBinContent(i));
    //F = F/hprod->GetBinWidth(i);
    Ffunc->AddPoint(hprod->GetBinCenter(i),F);
}
Ffunc->GetXaxis()->SetTitle("Energy (MeV)");
Ffunc->GetYaxis()->SetTitle("F factor ");

if(drawit){
    TCanvas *corrFcv = new TCanvas("corrFactorCv","correction Factor -- lost particles",150,150,800,678);
    corrFcv->SetLeftMargin(0.14);

    Ffunc->Draw("ALP");
    gPad->SetGridx();
    gPad->SetGridy();
}

return Ffunc;

}


TH1D* correctSpec(TH1D *expSpec = nullptr,
                        bool correctDead = false,
                        string pElossFileName = "/home/pi/ganil/kalscripts/eloss/results/UniformZ/CH2/v7/eloss_p_20.0deg_050.0um.root", 
                        string target_mat = "CH2",
                        Double_t th = 50,
                        char particle = 'p',
                        bool saveCanvas = false,
                        string pSpecFileName = "myspec.root"
                        ){

TStopwatch timer;


//defining the particle
Int_t Afis, Zfis;
defineParticle(particle, &Zfis, &Afis);


//defining the angle
float ang = 20;

// //open experimental file 
// TFile *expFile = new TFile(pSpecFileName.c_str(),"READ");
// TH1D *expSpec = (TH1D*)expFile->Get(hist_name.c_str());


//reading file for energy loss:
TFile *fEloss = new TFile(pElossFileName.c_str(),"READ"); 

// Check if the file was successfully opened
if (!fEloss || fEloss->IsZombie()) {
    std::cerr << "Warning: The file " << pElossFileName << " does not exist or could not be opened." << std::endl;
    if (fEloss) fEloss->Close(); // Close the file if it is open but in a bad state
    delete fEloss; // Clean up
}

//create the TTree to the loss file 
TTree *ElTree = (TTree*)fEloss->Get("SIM");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//defining the response function:

//First of all, the number of bins will depend on the spectrum I want to correct:
Int_t nbins = expSpec->GetNbinsX();
cout<<"Definig the number of bins for the response function equal to NbinsX for input spec = "<<nbins<<endl;


//Now, I will create the response function
TH2D *h= new TH2D("h","h",nbins, expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1),nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));


//Here I will draw the response function
ElTree->Draw("E:Erem>>h","","goff");
gStyle->SetOptStat("e");
h->SetTitle(Form("protons leaving %s target, %4.1f#mum",target_mat.c_str(),th));
h->GetXaxis()->SetTitle("Measured proton energy (MeV)");
h->GetYaxis()->SetTitle("Initial proton energy (MeV)");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//creating the projections over Y axis
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

//Create a canvas to temporarly store the corrected spectrum
TCanvas *corrCv = new TCanvas("corrCv","corrCv",50,50,1598,678);
corrCv->SetLeftMargin(0.14);
corrCv->Divide(2,1);



int cvID = 1;
corrCv->cd(cvID);

gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
gStyle->SetTextSize(0.04);
gStyle->SetTitleYSize(0.04);

gPad->SetGridx();
gPad->SetGridy();
// Configurações do canvas
corrCv->cd(cvID)->SetRightMargin(0.14);
corrCv->cd(cvID)->SetTopMargin(0.08);
corrCv->cd(cvID)->SetLeftMargin(0.14);
corrCv->cd(cvID)->SetBottomMargin(0.14);

h->Draw("colz");
h->GetYaxis()->SetTitleSize(0.06);
h->GetYaxis()->SetLabelSize(0.06);
h->GetZaxis()->SetLabelSize(0.06);
h->GetXaxis()->SetTitleSize(0.06);
h->GetXaxis()->SetLabelSize(0.06);

//create the correct spectrum
TH1D *corrected_exp = new TH1D("corrected_exp","corrected_exp",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));

Int_t counts; 
for (Int_t bin = 1; bin <= nbins; bin++)
{
  counts = expSpec->GetBinContent(bin);
  for (Int_t i = 0; i < counts; i++)
  {
    corrected_exp->Fill(projY[bin]->GetRandom(rand));
    //corrected_exp->Fill(expSpec->GetBinCenter(bin));
  }
  
}
corrected_exp->SetFillColor(kRed);
corrected_exp->SetFillStyle(3005);
corrected_exp->SetLineColor(kRed);
corrected_exp->SetLineWidth(2);
expSpec->SetLineColor(kBlack);
expSpec->SetLineWidth(2);

cvID =2;
corrCv->cd(cvID);


// Configurações do canvas
corrCv->cd(cvID)->SetRightMargin(0.03);
corrCv->cd(cvID)->SetTopMargin(0.08);
corrCv->cd(cvID)->SetLeftMargin(0.14);
corrCv->cd(cvID)->SetBottomMargin(0.14);

expSpec->SetTitle("");
expSpec->Draw();

gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
gStyle->SetTextSize(0.04);
gStyle->SetTitleYSize(0.04);
expSpec->GetXaxis()->SetTitle("Energy (MeV)");
expSpec->GetXaxis()->SetTitleSize(0.06);
expSpec->GetXaxis()->SetLabelSize(0.06);

expSpec->GetYaxis()->SetTitle("counts");
expSpec->GetYaxis()->SetMaxDigits(2);
expSpec->GetYaxis()->SetTitleSize(0.06);
expSpec->GetYaxis()->SetLabelSize(0.06);


gPad->SetGridx();
gPad->SetGridy();
corrected_exp->Draw("same");

// Find the maximum counts between corrected_exp and expSpec
Double_t maxCounts = TMath::Max(corrected_exp->GetMaximum(), expSpec->GetMaximum());
std::cout << "The maximum counts between corrected_exp and expSpec is: " << maxCounts << std::endl;

expSpec->GetYaxis()->SetRangeUser(0,maxCounts*1.15);


//regarding the directory    
std::string directory = "figures";

// Verificar se o diretório existe usando gSystem->AccessPathName
if (gSystem->AccessPathName(directory.c_str())) { // Diretório NÃO existe
    if (gSystem->mkdir(directory.c_str(), true) != 0) { // Tentar criar o diretório
        std::cerr << "Error creating directory: " << directory << std::endl;
        return nullptr;
    }
    std::cout << "Directory created successfully: " << directory << std::endl;
} else {
    std::cout << "Directory already exists: " << directory << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//generate correction factor for lost particles
//Create a canvas to temporarly store the corrected spectrum

TH1D *c2_exp;
if(correctDead){
    TCanvas *corrFcv = new TCanvas("corrFactorCv","correction Factor -- lost particles",150,150,800,678);
    corrFcv->SetLeftMargin(0.14);
    cout<<"\n.\n.\n.nbins = "<<nbins<<endl;

    TH1D *hprod = new TH1D("hprod","hprod",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));
    TH1D *hlost = new TH1D("hlost","hlost",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));

    std::string directory, filename;

    splitPath(pElossFileName, directory, filename);

    cout<<"Debugging splitPath:"<<endl;
    cout<<"directory = "<<directory<<endl;
    cout<<"filename = "<<filename<<endl;

    TGraph *Ffunc = GetPfactor(filename,directory,800, false);//= Fcorr(ElTree,nbins,false);
    Ffunc->GetXaxis()->SetTitle("Energy (MeV)");
    Ffunc->GetYaxis()->SetTitle("F factor ");
    Ffunc->Draw("ALP");

    gPad->SetGridx();
    gPad->SetGridy();


    c2_exp = (TH1D*)corrected_exp->Clone();
    c2_exp->SetNameTitle("c2_exp","c2_exp");

    for(int i = 1; i <= c2_exp->GetNbinsX(); i++){
        c2_exp->SetBinContent(i,corrected_exp->GetBinContent(i)*Ffunc->Eval(corrected_exp->GetBinCenter(i)));
    }
    corrCv->cd(2);
    c2_exp->SetLineColor(kBlue);
    c2_exp->SetFillColor(kCyan);
    c2_exp->SetFillStyle(3004);
    c2_exp->SetLineWidth(2);
    c2_exp->Draw("same");
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(expSpec,"Experimental Spectrum","lpf");
    leg->AddEntry(corrected_exp,"Corrected Spectrum","lpf");
    leg->AddEntry(c2_exp,"Corrected Spectrum with F factor","lpf");
    leg->Draw();

}else{ // if correctDead is false...
    corrCv->cd(2);
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(expSpec,"Experimental Spectrum","lpf");
    leg->AddEntry(corrected_exp,"Corrected Spectrum","lpf");
    leg->Draw();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(saveCanvas){
    corrCv->SaveAs(Form("%s/%s_TTC.root",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.png",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.pdf",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    cout << "Saving canvas as: " << Form("%s/%s_TTC.root (.pdf and .png)",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()) << endl;
    
}


timer.Print();

if(correctDead)return c2_exp;
else return corrected_exp;

}