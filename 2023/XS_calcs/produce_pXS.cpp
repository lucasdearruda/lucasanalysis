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
    Float_t Ea_zero, Eb_zero;
    for (int i = 0; i < nEbins; i++) {
        //calculating the energy:
        En = Ea + (i+0.5)*binMeV;
        Ea_zero = Ea + i * binMeV;
        Eb_zero = Ea + (i+1) * binMeV;
        std::cout << "Energy bin: " << i << " En: " << En << " Ea: " << Ea_zero << " Eb: " << Eb_zero << std::endl;
        //receive string of histograms;
        h = Fe_allAngles2(Ea+i*binMeV, Ea+(i+1)*binMeV);
        hist[0][i] = h[0]; //20deg
        hist[0][i]->SetNameTitle(Form("h_20deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 20 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[0][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[0][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[1][i] = h[1]; //40deg
        hist[1][i]->SetNameTitle(Form("h_40deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 40 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[1][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[1][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[2][i] = h[2]; //60deg
        hist[2][i]->SetNameTitle(Form("h_60deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 60 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[2][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[2][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[3][i] = h[3]; //80deg
        hist[3][i]->SetNameTitle(Form("h_80deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 80 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[3][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[3][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[4][i] = h[4];  //100deg
        hist[4][i]->SetNameTitle(Form("h_100deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 100 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[4][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[4][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[5][i] = h[5];  //120deg
        hist[5][i]->SetNameTitle(Form("h_120deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 120 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[5][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[5][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[6][i] = h[6];  //140deg
        hist[6][i]->SetNameTitle(Form("h_140deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 140 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[6][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[6][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[7][i] = h[7];  //160deg
        hist[7][i]->SetNameTitle(Form("h_160deg_En_%.1f_%.1f_MeV_bin%d",Ea_zero,Eb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 160 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[7][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[7][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");

        for (int l = 0; l < 8; l++)
        {
            ff->cd();
            Cv[i]->cd(l+1);
            hist[l][i]->Draw();
            gPad->SetGridx();
            gPad->SetGridy();
            //hist[l][i]->Write();
        }

        
        
        for(int j=0;j<8;j++){
            angle = 20 + 20 * j;
            cos_angle=  cos(  angle*TMath::Pi()/180.0 ) ;

            intg = h[j]->Integral("width");
            
            XSgr->AddPoint(angle,En,intg);
            
            gAngle[i]->SetPoint(j, angle,intg);
            gCosAngle[i]->SetPoint(j, cos_angle,intg);
            
            std::cout << "Angle: " << 20 + 20 * j << " MeV: " << i + 2 << " Value: " << intg << std::endl;
        }

        gAngle[i]->SetTitle(Form("Fe(n,Xp) - d#sigma/d#Omega (mb/sr) for %02.1f MeV",En));
        gAngle[i]->SetName(Form("ddOmega_%02.1f_MeV",En));
        gAngle[i]->GetXaxis()->SetTitle("Angle (deg)");
        gAngle[i]->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
        //gAngle[i]->Write();
        gCosAngle[i]->SetTitle(Form("Fe(n,Xp) - d#sigma/d#Omega (mb/sr) for %02.1f MeV",En));
        gCosAngle[i]->SetName(Form("ddOmega_%02.1f_MeV_cos",En));
        gCosAngle[i]->GetXaxis()->SetTitle("cos(#theta)");
        gCosAngle[i]->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
        //gCosAngle[i]->Write();
    }
    for (int i = 0; i < nEbins; i++) {
        Cv[i]->Write();
    }
    for (int i = 0; i < nEbins; i++) {
        gAngle[i]->Write();
        gCosAngle[i]->Write();
    }
    for (int i = 0; i < nEbins; i++) {
        for (int l = 0; l < 8; l++)
        {
            hist[l][i]->Write();
        }
    }

    ff->cd();
    XSgr->Write();
    TCanvas *cXS = new TCanvas("cXS","cXS");
    XSgr->Draw();

    ff->Close();

    timer.Print();
}