//
// Script for reconstructing the proton spectra 
//
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <fstream>
#include "/media/dearruda/Elements/LucasAnalysis/useful.h" //version6.10.2024.0002

#include "TStopwatch.h"
#include <time.h>

TTree *tx;


void checkTOF(){



    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    Double_t tSi1 = 53.4; //µm
    Double_t tSi2 = 1015; //µm
    //Double_t tCsI = 5*1e4; //µm

    //definition of materials
    KVMaterial *det1 = new KVMaterial("Si");//,tSi1*KVUnits::um);
    KVMaterial *det2 = new KVMaterial("Si",tSi2*KVUnits::um);
    // KVMaterial *det3 = new KVMaterial("CsI",tCsI*KVUnits::um);
    //KVMaterial *det2_dead = new KVMaterial("Si");//,dead_th*KVUnits::um);

    det1->SetThickness(tSi1 * KVUnits::um);
    det1->SetThickness(tSi2 * KVUnits::um);
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


    string points_file = "points.csv";
    TTree* pointsTree = new TTree("pointsTree", "pointsTree");
    pointsTree->ReadFile(points_file.c_str(), "X/F:Y/F:use_csi/I");

    // Declarar variáveis para leitura dos dados
    Int_t use_csi;
    Float_t X;
    Float_t Y;

    // Associar as variáveis às branches do TTree
    pointsTree->SetBranchAddress("use_csi", &use_csi);
    pointsTree->SetBranchAddress("X", &X);
    pointsTree->SetBranchAddress("Y", &Y);

    vector<Double_t> Xp; 
    vector<Double_t> Yp; 
    vector<Int_t> ucsi; 

    Int_t Np = pointsTree->GetEntries();

    cout<<"#: (si2,si1)"<<endl;
    for(int i=0; i<Np;i++){
        pointsTree->GetEntry(i);
        Xp.push_back(X);
        Yp.push_back(Y);
        ucsi.push_back(use_csi);
        cout<<i+1<<": "<<X<<" "<<Y<<" use csi? "<<use_csi<<endl;
    }

    TFile *tf = new TFile("/media/dearruda/Elements/LucasAnalysis/2023/reducedv61/028.root","READ");
    tx = (TTree*)tf->Get("M");

    TCanvas *Cpoints = new TCanvas("cpoints","cpoints",1000,500);

    tx->Draw("si1:si2>>h(400,0,16,250,0,2.5)","PID==1 && ang ==20", "colz");
    
    gPad->SetGridx();
    gPad->SetGridy();

    Double_t rad = 0.1;//Mev
    TEllipse *circles[Np];
    TLatex *tlx[Np];

    TCanvas *cv[Np];

    TH1D *htof_all[Np];
    TH1D *hsi[Np];
    TH1D *hsicsi[Np];
    TH1D *htof_point[Np];


    //get gamma flash:
    TNamed* gammaflash = (TNamed*)tf->Get("gFlash_T1(ns)");
    const char* value = gammaflash->GetTitle();
    Double_t gflash = std::atof(value);
    Double_t tgamma = provideTgama(1);

    //
    //rawtof 40 a 250 
    tx->SetAlias("rt", Form("%f - rawtof + %f", gflash,tgamma ));
    cout <<"  --> alias defined: 'rt' == '"<<  tx->GetAlias("rt") <<"'"<<endl;
    

    //drawing points
    for(int i=0; i<Np;i++){
        Cpoints->cd();
        circles[i] = new TEllipse(Xp[i],Yp[i],rad,rad);
        circles[i]->SetLineColor(kBlack);
        circles[i]->SetFillColor(kBlack);
        circles[i]->SetFillStyle(3005);
        circles[i]->SetLineWidth(2);
        circles[i]->Draw("same"); 


        if(i<11){
            tlx[i] = new TLatex(Xp[i] +0.15, Yp[i] +0.11,Form("%d",i+1));    
        }else{
            tlx[i] = new TLatex(Xp[i] +0.16, Yp[i] -0.15,Form("%d",i+1));
        }
        tlx[i]->SetTextColor(kBlack);
        tlx[i]->Draw("same");

        htof_all[i] = new TH1D(Form("htof_all%d",i+1),Form("htof_all%d",i+1),350,50,400);
        
        htof_all[i]->SetLineColor(kBlue);
        htof_all[i]->SetLineWidth(2);

        htof_point[i] = new TH1D(Form("htof_p%d",i+1),Form("htof_p%d",i+1),350,50,400);

        htof_point[i]->SetLineColor(kRed);
        htof_point[i]->SetFillColor(kRed);
        htof_point[i]->SetFillStyle(3005);
        htof_point[i]->SetLineWidth(2);

        hsi[i] = new TH1D(Form("hsi_%d",i+1),Form("hsi_%d",i+1),1600,0,40);
        hsicsi[i] = new TH1D(Form("hsicsi_%d",i+1),Form("hsicsi_%d",i+1),1600,0,40);

        tx->SetAlias(Form("circle%d",i+1),Form("PID ==1 && ang == 20 && pow(si2 - %f,2)+pow(si1 - %f,2)<=%f",Xp[i],Yp[i],rad*rad));


        cv[i] = new TCanvas(Form("Point_%d", i+1),Form("Point_%d", i+1),50,50,1000,500);
        
        
        //division
        if(!ucsi[i]){
            cv[i]->Divide(2,1);
        }else{
            cv[i]->Divide(3,1);
        }

        //common part:
         cv[i]->cd(1);
            
            hsi[i]->SetLineWidth(2);
            hsi[i]->SetLineColor(kRed);
            hsi[i]->SetFillColor(kRed);
            hsi[i]->SetFillStyle(3005);

            hsicsi[i]->SetLineWidth(2);
            hsicsi[i]->SetLineColor(kBlack);
            hsicsi[i]->SetFillColor(kBlack);
            hsicsi[i]->SetFillStyle(3004);
            

            tx->Draw(Form("si1+si2>>hsi_%d",i+1), Form("circle%d",i+1));
            tx->Draw(Form("si1+si2+csi>>hsicsi_%d",i+1), Form("circle%d",i+1),"same");

            int lastBinNonZero = hsicsi[i]->GetNbinsX();
            while(hsicsi[i]->GetBinContent(lastBinNonZero)<=1){
                lastBinNonZero--;
                if(lastBinNonZero == 0) break;
            }

            cout<<"point #"<<i+1<<": lastBinNonzero = "<<lastBinNonZero<<", centered at "<<hsicsi[i]->GetBinCenter(lastBinNonZero)<<" MeV."<<endl;
            hsi[i]->GetXaxis()->SetRangeUser(Xp[i]+Yp[i] -1.0,  floor(hsicsi[i]->GetBinCenter(lastBinNonZero)+1.0));

            gPad->SetGridx();
            gPad->SetGridy();
            
            cv[i]->cd(2);

            tx->Draw(Form("rt>>htof_all%d",i+1),"PID==1 && ang == 20");
            
            if(ucsi[i]) tx->Draw(Form("rt>>htof_p%d",i+1),Form("circle%d && si1+si2+csi>=11.5",i+1),"same");
            else    tx->Draw(Form("rt>>htof_p%d",i+1),Form("circle%d",i+1),"same");
            gPad->SetLogy();
            gPad->SetGridx();
            gPad->SetGridy();
        if(ucsi[i]){
            cout<<"point #"<<i+1<<"part 3..."<<endl;
            //cv[i]->cd(3);
           // tx->Draw(Form("%f:%f>>h2e%d",,,i+1),Form("circle%d",i+1),"same");
            

            

        }
    }


}