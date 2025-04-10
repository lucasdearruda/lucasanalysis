#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
//#include "useful.h"
#include <algorithm> //to use function swap
#include <time.h> //for the time

using namespace std; // so doesnt need to declare "std::"

string pathRuns = "/mnt/medley/runswithMM_2022";

void CalOverData( int nruni, int nrunf,int tel = 1, float gap = 0.0){
//plot simulations over the calibrated runs, for the chosen telescope

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //variables to insert data-time to output file---------------
    char cur_time[128];

    time_t      t;
    struct tm*  ptm;

    t = time(NULL);
    ptm = localtime(&t);

    strftime(cur_time, 128, "%Y-%m-%d_%H:%M:%S", ptm);
    //start clock
    clock_t tStart = clock();
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //opening files
    char name[1000];
    TChain *tx = NULL;

    if (!tx)
        tx = new TChain("RD");

    vector<Int_t> entries;

    for (int i = nruni; i <= nrunf; i++)
    {
        Bool_t fExist = true;
        for (Int_t j = 0; j <= 999 && fExist; j++)
        {
            ifstream mfile;
            sprintf(name, "%s/r%04d_%03dr.root",pathRuns.c_str(), i, j);
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
    }
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//opening kaliveda SIMs

string filename[6];
string pathCal = "/mnt/medley/kaliveda_results/Telescopes_dE_E_results2023";
cout<<"running for telescope "<<tel<<"..."<<endl;
filename[1] = Form("%s/kal_p_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //p production plot
filename[2] = Form("%s/kal_d_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //d production plot
filename[3] = Form("%s/kal_t_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //t production plot
filename[4] = Form("%s/kal_h_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //h production plot
filename[5] = Form("%s/kal_a_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //a production plot

//Defining canvases
TCanvas *C0  = new TCanvas("c0","dE1 vs dE2", 150,10,800,800);
TCanvas *C1  = new TCanvas("c1","dE2 vs dE3", 350,10,800,800);
TCanvas *C2  = new TCanvas("c2","dE1 vs dTOT", 550,10,800,800);


C0->cd()->SetRightMargin(0.16);
C0->cd()->SetLeftMargin(0.12);
C0->cd()->SetBottomMargin(0.12);


C1->cd()->SetRightMargin(0.16);
C1->cd()->SetLeftMargin(0.12);
C1->cd()->SetBottomMargin(0.12);


C2->cd()->SetRightMargin(0.16);
C2->cd()->SetLeftMargin(0.12);
C2->cd()->SetBottomMargin(0.12);


 gStyle->SetOptStat(0);

 TFile *g[6];
 TTree *TR[6];
 Long64_t nEntries[6];

//Tfile and TTree attribution
for(int i=1; i<6; i++){
    g[i] = new TFile(filename[i].c_str());
        TR[i] = (TTree*)g[i]->Get("SIM");
        nEntries[i] = TR[i]->GetEntries();
        cout<<Form("Entries in the tree SIM (%s): ",filename[i].c_str())<<nEntries[i]<< endl;

}

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

hist[0] = new TH2D("exp0","dE1:dE2", nbinsx,0,Xmax[0],nbinsy,0,Ymax[0]);
hist[1] = new TH2D("exp1","dE2:dE3", nbinsx,0,Xmax[1],nbinsy,0,Ymax[1]);
hist[2] = new TH2D("exp2","dE1:dE1dE2dE3", nbinsx,0,Xmax[2],nbinsy,0,Ymax[2]);

//-------------------------------
//Drawing stuff

TGraph *gp1,*gp2,*gp3; //protons
TGraph *gd1,*gd2,*gd3; //deuterons
TGraph *gt1,*gt2,*gt3; //tritons
TGraph *gh1,*gh2,*gh3; //helium
TGraph *ga1,*ga2,*ga3; //alphas

gp1 = new TGraph();
gp2 = new TGraph();
gp3 = new TGraph();


gd1 = new TGraph();
gd2 = new TGraph();
gd3 = new TGraph();


gt1 = new TGraph();
gt2 = new TGraph();
gt3 = new TGraph();



gh1 = new TGraph();
gh2 = new TGraph();
gh3 = new TGraph();



ga1 = new TGraph();
ga2 = new TGraph();
ga3 = new TGraph();

float finetunedE1 = 1.00;
float finetunedE2 = 1.00;

tx->Draw(Form("Medley_%d_dE1*%f:Medley_%d_dE2*%f>>exp0",tel,finetunedE1,tel,finetunedE2),"","colz");
tx->Draw(Form("Medley_%d_dE2*%f:Medley_%d_Eres>>exp1",tel,finetunedE2,tel),"","colz");
tx->Draw(Form("Medley_%d_dE1:Medley_%d_dE1+Medley_%d_dE2+Medley_%d_Eres>>exp2",tel,tel,tel,tel),"","colz");



//--------------------------------
C0->cd();
hist[0]->GetXaxis()->SetTitle("#Delta E_{2} (MeV)");
hist[0]->GetYaxis()->SetTitle("#Delta E_{1} (MeV)");
hist[0]->Draw("colz");
TR[1]->Draw("si1:si2>>gp1","","same");
TR[2]->Draw("si1:si2>>gd1","","same");
TR[3]->Draw("si1:si2>>gt1","","same");
TR[4]->Draw("si1:si2>>gh1","","same");
TR[5]->Draw("si1:si2>>ga1","","same");
gPad->SetGridx();
gPad->SetGridy();

gPad->SetLogz();

C1->cd();
hist[1]->GetXaxis()->SetTitle("#Delta E_{3} (MeV)");
hist[1]->GetYaxis()->SetTitle("#Delta E_{2} (MeV)");
hist[1]->Draw("colz");

TR[1]->Draw(Form("si2:csi+%f>>gp2",gap),"","same");
TR[2]->Draw(Form("si2:csi+%f>>gd2",gap),"","same");
TR[3]->Draw(Form("si2:csi+%f>>gt2",gap),"","same");
TR[4]->Draw(Form("si2:csi+%f>>gh2",gap),"","same");
TR[5]->Draw(Form("si2:csi+%f>>ga2",gap),"","same");

gPad->SetGridx();
gPad->SetGridy();

gPad->SetLogz();

C2->cd();
hist[2]->GetXaxis()->SetTitle("#Delta E_{1}+#Delta E_{2}+#Delta E_{3} (MeV)");
hist[2]->GetYaxis()->SetTitle("#Delta E_{1} (MeV)");
hist[2]->Draw("colz");

TR[1]->Draw(Form("si1:si1+si2+csi+%f>>gp3",gap),"","same");
TR[2]->Draw(Form("si1:si1+si2+csi+%f>>gd3",gap),"","same");
TR[3]->Draw(Form("si1:si1+si2+csi+%f>>gt3",gap),"","same");
TR[4]->Draw(Form("si1:si1+si2+csi+%f>>gh3",gap),"","same");
TR[5]->Draw(Form("si1:si1+si2+csi+%f>>ga3",gap),"","same");

gPad->SetGridx();
gPad->SetGridy();

gPad->SetLogz();

gp1->SetLineColor(kBlack);
gp2->SetLineColor(kBlack);
gp3->SetLineColor(kBlack);
gp1->SetLineWidth(kBlack);
gp2->SetLineWidth(kBlack);
gp3->SetLineWidth(kBlack);


gd1->SetLineColor(kRed);
gd2->SetLineColor(kRed);
gd3->SetLineColor(kRed);
gd1->SetLineWidth(kRed);
gd2->SetLineWidth(kRed);
gd3->SetLineWidth(kRed);

gt1->SetLineColor(kGreen);
gt2->SetLineColor(kGreen);
gt3->SetLineColor(kGreen);
gt1->SetLineWidth(kGreen);
gt2->SetLineWidth(kGreen);
gt3->SetLineWidth(kGreen);

gh1->SetLineColor(kBlue);
gh2->SetLineColor(kBlue);
gh3->SetLineColor(kBlue);
gh1->SetLineWidth(kBlue);
gh2->SetLineWidth(kBlue);
gh3->SetLineWidth(kBlue);

ga1->SetLineColor(kPink);
ga2->SetLineColor(kPink);
ga3->SetLineColor(kPink);
ga1->SetLineWidth(kPink);
ga2->SetLineWidth(kPink);
ga3->SetLineWidth(kPink);


}
