#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TCutG.h"
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
//#include "useful.h"
#include <algorithm> //to use function swap


using namespace std; // so doesnt need to declare "std::"

void checkcuts_v5(int nruni=407,int tel = 1){
//plot simulations over the calibrated runs, for the chosen telescope

    //opening files
    char name[1000];
    TChain *tx = NULL;

    if (!tx)
        tx = new TChain("AD");

    vector<Int_t> entries;
    
    Bool_t fExist = true;
    for (Int_t j = 0; j <= 999 && fExist; j++)
    {
        ifstream mfile;
        sprintf(name, "/media/dearruda/Elements/RootA_2023/r%04d_%03da.root", nruni, j);
        mfile.open(name);
        if (mfile)
        {
            mfile.close();
        }
        else
        {
            fExist = false;
        }
        if (fExist)
        {
            cout << "Adding " << name << endl;
            tx->Add(name);
            entries.push_back((Int_t)(tx->GetEntries()));
            cout << "Entries: " << entries.back() / 1000 << "k " << endl;
        }
    }
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//Defining canvases
TCanvas *C0  = new TCanvas("c0","dE1 vs dE2", 150,10,800,800);
TCanvas *C1  = new TCanvas("c1","dE2 vs dE3", 800+150,10,800,800);


C0->cd()->SetRightMargin(0.16);
C0->cd()->SetLeftMargin(0.12);
C0->cd()->SetBottomMargin(0.12);


C1->cd()->SetRightMargin(0.16);
C1->cd()->SetLeftMargin(0.12);
C1->cd()->SetBottomMargin(0.12);

 gStyle->SetOptStat(0);

int nbinsx = 1024;
int nbinsy =1024;

//we will plot three things: dE1:dE2, dE2:dE3 and dE1:dTOT (dE1:dE1+dE2+dE3)

//declaring the maximum range
float Xmax[3] ,Ymax[3];

Xmax[0]=50;//dE1:dE2 max range in Mev
Xmax[1]=45;//dE2:dE3 max range in Mev
Xmax[2]=45;//dE1:dETOT max range in Mev

Ymax[0]=10;//dE1:dE2 max range in Mev
Ymax[1]=30;//dE2:dE3 max range in Mev
Ymax[2]=10;//dE1:dETOT max range in Mev

//declaring experimental histos
TH2D *hist[3];

hist[0] = new TH2D("dE_dE","dE1:dE2", nbinsx,0,Xmax[0],nbinsy,0,Ymax[0]);
hist[1] = new TH2D("dE_Eres","dE2:dE3", nbinsx,0,Xmax[1],nbinsy,0,Ymax[1]);


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//declaring cuts    
    string cuts_path =   "/media/dearruda/Elements/LucasAnalysis/2023/PIDv6";

    TCutG *cut_protons;
    TCutG *cut_deuterons;
    TCutG *cut_tritons;

    TCutG *cut_protonsPT;
    TCutG *cut_deuteronsPT;
    TCutG *cut_tritonsPT;

    TCutG *cut_protonsCsI;
    TCutG *cut_deuteronsCsI;
    TCutG *cut_tritonsCsI;

    TCutG *cut_protonsNPTCsI;
    TCutG *cut_deuteronsNPTCsI;
    
    TCutG *cut_he;
    TCutG *cut_alphas;

    //Opening cuts: 
    //he-3
    cout << "\n -> loading HELIUM3 cut:" << Form(".L %s/he3_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/he3_%d.C", cuts_path.c_str(),tel));
    cut_he = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("he3_%d",tel));
    if(cut_he!=NULL) cout << "  --> OK."<< endl;
    cut_he->SetLineWidth(2);
    cut_he->SetLineColor(kCyan);
    //alphas
    cout << "\n -> loading ALPHA cut:" << Form(".L %s/alphas%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/alphas%d.C", cuts_path.c_str(),tel));
    cut_alphas = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("alphas%d",tel));
    if(cut_alphas!=NULL) cout << "  --> OK."<< endl;
    cut_alphas->SetLineWidth(2);
    cut_alphas->SetLineColor(kMagenta);


    //protons NPT
    cout << "\n -> loading NPT PROTONS cut:" << Form(".L %s/p_npt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/p_npt_%d.C", cuts_path.c_str(),tel));
    cut_protons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_npt_%d",tel));
    if(cut_protons!=NULL) cout << "  --> OK."<< endl;
    cut_protons->SetLineWidth(2);
    cut_protons->SetLineColor(kRed);
    //deuterons NPT
    cout << "\n -> loading NPT DEUTERONS cut:" << Form(".L %s/d_npt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/d_npt_%d.C", cuts_path.c_str(),tel));
    cut_deuterons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_npt_%d",tel));
    if(cut_deuterons!=NULL) cout << "  --> OK."<< endl;
    cut_deuterons->SetLineWidth(2);
    cut_deuterons->SetLineColor(kGreen);
    //tritons NPT
    cout << "\n -> loading NPT PROTONS cut:" << Form(".L %s/t_npt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/t_npt_%d.C", cuts_path.c_str(),tel));
    cut_tritons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_npt_%d",tel));
    if(cut_tritons!=NULL) cout << "  --> OK."<< endl;
    cut_tritons->SetLineWidth(2);
    cut_tritons->SetLineColor(kBlack);



    //protons PT
    cout << "\n -> loading PT PROTONS cut:" << Form(".L %s/p_pt_dE_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/p_pt_dE_%d.C", cuts_path.c_str(),tel));
    cut_protonsPT = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_dE_%d",tel));
    if(cut_protonsPT!=NULL) cout << "  --> OK."<< endl;
    cut_protonsPT->SetLineWidth(2);
    cut_protonsPT->SetLineColor(kRed);
    //deuterons PT
    cout << "\n -> loading PT DEUTERONS cut:" << Form(".L %s/d_pt_dE_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/d_pt_dE_%d.C", cuts_path.c_str(),tel));
    cut_deuteronsPT = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_dE_%d",tel));
    if(cut_deuteronsPT!=NULL) cout << "  --> OK."<< endl;
    cut_deuteronsPT->SetLineWidth(2);
    cut_deuteronsPT->SetLineColor(kGreen);
    //tritons PT
    cout << "\n -> loading PT PROTONS cut:" << Form(".L %s/t_pt_dE_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/t_pt_dE_%d.C", cuts_path.c_str(),tel));
    cut_tritonsPT = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_dE_%d",tel));
    if(cut_tritonsPT!=NULL) cout << "  --> OK."<< endl;
    cut_tritonsPT->SetLineWidth(2);
    cut_tritonsPT->SetLineColor(kBlack);

    //protons CSI
    cout << "\n -> loading CsI PROTONS cut:" << Form(".L %s/p_pt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/p_pt_%d.C", cuts_path.c_str(),tel));
    cut_protonsCsI = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_%d",tel));
    if(cut_protonsCsI!=NULL) cout << "  --> OK."<< endl;
    cut_protonsCsI->SetLineWidth(2);
    cut_protonsCsI->SetLineColor(kRed);
    //deuterons CSI
    cout << "\n -> loading CsI DEUTERONS cut:" << Form(".L %s/d_pt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/d_pt_%d.C", cuts_path.c_str(),tel));
    cut_deuteronsCsI = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_%d",tel));
    if(cut_deuteronsCsI!=NULL) cout << "  --> OK."<< endl;
    cut_deuteronsCsI->SetLineWidth(2);
    cut_deuteronsCsI->SetLineColor(kGreen);
    //tritons CSI
    cout << "\n -> loading CsI PROTONS cut:" << Form(".L %s/t_pt_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/t_pt_%d.C", cuts_path.c_str(),tel));
    cut_tritonsCsI = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_%d",tel));
    if(cut_tritonsCsI!=NULL) cout << "  --> OK."<< endl;
    cut_tritonsCsI->SetLineWidth(2);
    cut_tritonsCsI->SetLineColor(kBlack);


    //protons NPT CSI
    cout << "\n -> loading CsI PROTONS cut:" << Form(".L %s/pnptcsi_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/pnptcsi_%d.C", cuts_path.c_str(),tel));
    cut_protonsNPTCsI = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("pnptcsi_%d",tel));
    if(cut_protonsNPTCsI!=NULL) cout << "  --> OK."<< endl;
    cut_protonsNPTCsI->SetLineWidth(2);
    cut_protonsNPTCsI->SetLineColor(kRed);
    //deuterons NPT CSI
    cout << "\n -> loading CsI DEUTERONS cut:" << Form(".L %s/dnptcsi_%d.C", cuts_path.c_str(),tel)<< "... " << endl;
    gROOT->ProcessLine(Form(".L %s/dnptcsi_%d.C", cuts_path.c_str(),tel));
    cut_deuteronsNPTCsI = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("dnptcsi_%d",tel));
    if(cut_deuteronsNPTCsI!=NULL) cout << "  --> OK."<< endl;
    cut_deuteronsNPTCsI->SetLineWidth(2);
    cut_deuteronsNPTCsI->SetLineColor(kGreen);


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


//float finetunedE1 = 1.00;
//float finetunedE2 = 1.00;
C0->cd();
tx->Draw(Form("Medley_%d_dE1:Medley_%d_dE2>>dE_dE",tel,tel),"","colz");
gPad->SetLogz();
cut_he->Draw("same");
cut_alphas->Draw("same");
cut_protons->Draw("same");
cut_deuterons->Draw("same");
cut_tritons->Draw("same");

cut_protonsPT->Draw("same");
cut_deuteronsPT->Draw("same");
cut_tritonsPT->Draw("same");


C1->cd();
tx->Draw(Form("Medley_%d_dE2:Medley_%d_Eres>>dE_Eres",tel,tel),"","colz");
gPad->SetLogz();
cut_protonsCsI->Draw("same");
cut_deuteronsCsI->Draw("same");
cut_tritonsCsI->Draw("same");


cut_protonsNPTCsI->Draw("same");
cut_deuteronsNPTCsI->Draw("same");




}
