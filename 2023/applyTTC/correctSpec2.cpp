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

string ganil_folder= "/home/pi/ganil/"; 


TH1D* correctSpec(TH1D *expSpec = nullptr,
                        string pSpecFileName = "myspec.root",
                        string pElossFileName = "/home/pi/ganil/kalscripts/eloss/results/UniformZ/C/v6/eloss_p_0.3deg_075.0um.root", 
                        string target_mat = "C",
                        Double_t th = 75,
                        char particle = 'p',
                        bool saveCanvas = false
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

//I need to create a canvas to draw the response function
//TCanvas *c0 = new TCanvas("response_func_canvas","response_func_canvas",50,50,800,600);
//c0->SetRightMargin(0.16);

//Here I will draw the response function
ElTree->Draw("E:Erem>>h","","goff");
gStyle->SetOptStat("e");
h->SetTitle(Form("protons leaving %s target, %4.1f#mum",target_mat.c_str(),th));
h->GetXaxis()->SetTitle("Measured proton energy (MeV)");
h->GetYaxis()->SetTitle("Initial proton energy (MeV)");
//h->Draw("colz");
//gPad->SetGridx();
//gPad->SetGridy();

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
expSpec->SetLineColor(kBlue);
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


if(saveCanvas){
    corrCv->SaveAs(Form("%s/%s_TTC.root",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.png",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.pdf",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    cout << "Saving canvas as: " << Form("%s/%s_TTC.root (.pdf and .png)",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()) << endl;
    
}
return corrected_exp;
timer.Print();

}