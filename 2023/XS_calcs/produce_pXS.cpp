#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_allAngles2.cpp" 

#include <vector>
#include <iostream>
#include "TStopwatch.h"

void produce_pXS(Float_t Ea = 25, Float_t Eb = 26, Float_t binMeV = 1, string filename = "Fe_pXS.root"){ 
    
    TStopwatch timer;
    TCanvas *cc = new TCanvas("cc","cc");

    //string of histograms; 
    std::vector<TH1D*> h;


    TGraph2D * XSgr = new TGraph2D();
    Float_t intg;

    Int_t nEbins = 0;
    for (Float_t i = Ea; i < Eb; i+=binMeV) {
        nEbins++;
    }//calculate the number of histograms

    TFile *ff = new TFile(filename.c_str(), "RECREATE");

    //histogram[angle][Energy_bin]
    TH1D *hist[8][nEbins];

    TGraph *gAngle[nEbins];
    TGraph *gCosAngle[nEbins];

    //for each bin, create a graph
    for(int i = 0; i < nEbins; i++) {
        gAngle[i] = new TGraph();
        gCosAngle[i] = new TGraph();
    }



    Float_t cos_angle, angle, En;

    TCanvas *Cv[nEbins];
    for(int i = 0; i < nEbins; i++) {
        //Mean energy of that bin:
        En = Ea + i * binMeV/2;
        Cv[i] = new TCanvas(Form("Cv_%d",i),Form("Cv_%02.1f_MeV",En),100+i*5,100+i*5,1200,700);
        Cv[i]->Divide(3,3);
    }

    
    // for each energy bin, create histograms and fill the graphs
    for (int i = 0; i < nEbins; i++) {
 
        //receive string of histograms;
        h = Fe_allAngles2(Ea+i*binMeV, Ea+(i+1)*binMeV);
        hist[0][i] = h[0]; //20deg
        hist[0][i]->SetNameTitle(Form("h_20_%d",i),Form("h_20_%d",i));
        hist[0][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[0][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[1][i] = h[1]; //40deg
        hist[1][i]->SetNameTitle(Form("h_40_%d",i),Form("h_40_%d",i));
        hist[1][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[1][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[2][i] = h[2]; //60deg
        hist[2][i]->SetNameTitle(Form("h_60_%d",i),Form("h_60_%d",i));
        hist[2][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[2][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[3][i] = h[3]; //80deg
        hist[3][i]->SetNameTitle(Form("h_80_%d",i),Form("h_80_%d",i));
        hist[3][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[3][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[4][i] = h[4];  //100deg
        hist[4][i]->SetNameTitle(Form("h_100_%d",i),Form("h_100_%d",i));
        hist[4][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[4][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[5][i] = h[5];  //120deg
        hist[5][i]->SetNameTitle(Form("h_120_%d",i),Form("h_120_%d",i));
        hist[5][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[5][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[6][i] = h[6];  //140deg
        hist[6][i]->SetNameTitle(Form("h_140_%d",i),Form("h_140_%d",i));
        hist[6][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[6][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");
        hist[7][i] = h[7];  //160deg
        hist[7][i]->SetNameTitle(Form("h_160_%d",i),Form("h_160_%d",i));
        hist[7][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[7][i]->GetYaxis()->SetTitle("d#sigma/d#OmegadE (mb/sr/MeV)");

        for (int l = 0; l < 8; l++)
        {
            ff->cd();
            Cv[i]->cd(l+1);
            hist[l][i]->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            hist[l][i]->Write();
        }

        //calculating the energy:
        En = Ea + i * binMeV/2;
        
        for(int j=0;j<8;j++){
            angle = 20 + 20 * j;
            cos_angle=  cos(  angle*TMath::Pi()/180.0 ) ;

            intg = h[j]->Integral("width");
            
            XSgr->AddPoint(angle,En,intg);
            
            gAngle[i]->SetPoint(j, angle,intg);
            gCosAngle[i]->SetPoint(j, cos_angle,intg);
            
            std::cout << "Angle: " << 20 + 20 * j << " MeV: " << i + 2 << " Value: " << intg << std::endl;
        }

        gAngle[i]->SetTitle(Form("d#sigma/d#Omega (mb/sr) for %02.1f MeV",En));
        gAngle[i]->GetXaxis()->SetTitle("Angle (deg)");
        gAngle[i]->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
        gAngle[i]->Write();
        gCosAngle[i]->SetTitle(Form("d#sigma/d#Omega (mb/sr) for %02.1f MeV",En));
        gCosAngle[i]->GetXaxis()->SetTitle("cos(#theta)");
        gCosAngle[i]->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
        gCosAngle[i]->Write();
        Cv[i]->Write();
    }
    ff->cd();
    XSgr->Write();
    TCanvas *cXS = new TCanvas("cXS","cXS");
    XSgr->Draw();

    ff->Close();

    timer.Print();
}