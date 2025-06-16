#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_allAngles2.cpp" 
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/include/bin_width.hh" 
#include <vector>
#include <iostream>
#include "TStopwatch.h"
#include "include/xs_functions.hh"
//xs_functions contains: 
/// -> total Charge calculation
/// -> global def of Nc, omega_tel, and NA


//target dependent defintions:
Double_t mc = 0.0851;//g 
Double_t M = 55.845; //g/mol

//define binMeV in a HH file 

void produce_pXSv3(Float_t Ea = 25, Float_t Eb = 26, string filename = "Fe_p_XSv3.root"){ 

    TStopwatch timer;
    
    TCanvas *cc = new TCanvas("cc","cc");

    //string of histograms; 
    std::vector<TH1D*> h;

    //Create the 2D graph for cross-section data
    TGraph2DErrors * XSgr = new TGraph2DErrors();
    Float_t dsigma_dOmega, dsigma_dOmegaEr;
    

    Int_t nEbins = 0;
    //Add some warning for safety 

    if( Ea + binMeV > Eb) {
            std::cerr << "Warning: The last energy bin will not be filled completely. Adjust the binning or energy range." << std::endl;
            return;
    }

    for (Float_t i = Ea; i < Eb; i+=binMeV) {
        nEbins++;
    }//calculate the number of histograms

    TFile *ff = new TFile(filename.c_str(), "RECREATE");

    //histogram[angle][Energy_bin]
    TH1D *hist[8][nEbins];

    TGraphErrors *gAngle[nEbins];
    TGraphErrors *gCosAngle[nEbins];

    //for each bin, create a graph
    for(int i = 0; i < nEbins; i++) {
        gAngle[i] = new TGraphErrors();
        gCosAngle[i] = new TGraphErrors();
    }



    Float_t cos_angle, angle, En, angleEr, cosAngleEr;
    Float_t Ea_zero, Eb_zero;
    Int_t intEa_zero, intEb_zero,decEa_zero, decEb_zero,intEn, decEn;

    TCanvas *Cv[nEbins];
    for(int i = 0; i < nEbins; i++) {
        //Mean energy of that bin:
        En = Ea + i * binMeV/2;
        Cv[i] = new TCanvas(Form("Cv_%d",i),Form("Cv_%02.1f_MeV",En),100+i*5,100+i*5,1200,700);
        Cv[i]->Divide(3,3);
    }

    string target_details = "Fe_thick";
    // for each energy bin, create histograms and fill the graphs
    Int_t npoint=0;
    for (int i = 0; i < nEbins; i++) {
        //calculating the energy:
        En = Ea + (i+0.5)*binMeV;
        Ea_zero = Ea + i * binMeV;
        Eb_zero = Ea + (i+1) * binMeV;

        //Convert parts to mode X.Y --> XpY, so it gets better to load up the names of the objects
        intEn = static_cast<int>(En); // Convert to integer for binning
        decEn = static_cast<int>(En * 10) % 10;

        intEa_zero = static_cast<int>(Ea_zero); // Convert to integer for binning
        decEa_zero = static_cast<int>(Ea_zero * 10) % 10;
        
        intEb_zero = static_cast<int>(Eb_zero); // Convert to integer for binning
        decEb_zero = static_cast<int>(Eb_zero * 10) % 10;
        
        
        std::cout<<"########################## debugging info ##########################"<<std::endl;
        std::cout<< "Processing energy bin: " << i << " with En: " << En << ", Ea_zero: " << Ea_zero << ", Eb_zero: " << Eb_zero << std::endl;
        std::cout << "intEn: " << intEn << " decEn: " << decEn << " intEa_zero: " << intEa_zero << " decEa_zero: " << decEa_zero << " intEb_zero: " << intEb_zero << " decEb_zero: " << decEb_zero << std::endl;
        std::cout<<" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "<<std::endl;
        std::cout << "Energy bin: " << i << " En: " << En << " Ea: " << Ea_zero << " Eb: " << Eb_zero << std::endl;
        //receive string of histograms;
        h = Fe_allAngles2(Ea+i*binMeV, Ea+(i+1)*binMeV);
        hist[0][i] = h[0]; //20deg
        hist[0][i]->SetNameTitle(Form("%s_20deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 20 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[0][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[0][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[1][i] = h[1]; //40deg
        hist[1][i]->SetNameTitle(Form("%s_40deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 40 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[1][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[1][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[2][i] = h[2]; //60deg
        hist[2][i]->SetNameTitle(Form("%s_60deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 60 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[2][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[2][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[3][i] = h[3]; //80deg
        hist[3][i]->SetNameTitle(Form("%s_80deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 80 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[3][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[3][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[4][i] = h[4];  //100deg
        hist[4][i]->SetNameTitle(Form("%s_100deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 100 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[4][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[4][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[5][i] = h[5];  //120deg
        hist[5][i]->SetNameTitle(Form("%s_120deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 120 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[5][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[5][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[6][i] = h[6];  //140deg
        hist[6][i]->SetNameTitle(Form("%s_140deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 140 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[6][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[6][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");
        hist[7][i] = h[7];  //160deg
        hist[7][i]->SetNameTitle(Form("%s_160deg_En_%dp%d_%dp%d_MeV_bin%d",target_details.c_str(),intEa_zero,decEa_zero,intEb_zero,decEb_zero,i),Form("Fe(n,Xp) #theta_{LAB} = 160 deg, E_{NN} = %.1f to %.1f MeV",Ea_zero, Eb_zero));
        hist[7][i]->GetXaxis()->SetTitle("E_p (MeV)");
        hist[7][i]->GetYaxis()->SetTitle("d^{2}#sigma/d#OmegadE (mb/sr/MeV)");

        for (int l = 0; l < 8; l++)
        {
            ff->cd();
            Cv[i]->cd(l+1);
            hist[l][i]->Draw("E1 HIST");
            gPad->SetGridx();
            gPad->SetGridy();
            //hist[l][i]->Write();
        }

        
        
        for(int j=0;j<8;j++){
            angle = 20 + 20 * j;
            angleEr = 1;
            cos_angle=  cos(  angle*TMath::Pi()/180.0 ) ;
            //cosAngleEr =  cos(  angleEr*TMath::Pi()/180.0 ) ;
            //we have to propagate the error in the angle to the cos(angle)
            cosAngleEr =  sin(  angle*TMath::Pi()/180.0 )*angleEr*TMath::Pi()/180.0 ;

            intg = h[j]->Integral("width");
            intgEr = histIntegralError(h[j]); // Get the error in the integral
            XSgr->SetPoint(npoint,En,angle,intg);

            XSgr->SetPointError(npoint++,0,angleEr,intgEr); // Add error in the integral
            
            gAngle[i]->SetPoint(j, angle,intg);
            gCosAngle[i]->SetPoint(j, cos_angle,intg);
            
            gAngle[i]->SetPointError(j, 1,intgEr );
            gCosAngle[i]->SetPointError(j, cosAngleEr,intgEr);
            

            std::cout << "Angle: (" << 20 + 20 * j << "#pm"<< angleEr<<")^{o}, " << En << " MeV. Value: " << intg <<" #pm "<< intgEr<<" mb/sr"<< std::endl;
        }

        gAngle[i]->SetTitle(Form("Fe(n,Xp) - d#sigma/d#Omega (mb/sr) for %02.1f #pm %.1f  MeV",En,binMeV));
        gAngle[i]->SetName(Form("ddOmega_%s_%dp%d_MeV",target_details.c_str(),intEn,decEn));
        gAngle[i]->GetXaxis()->SetTitle("Angle (deg)");
        gAngle[i]->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
        //gAngle[i]->Write();
        gCosAngle[i]->SetTitle(Form("Fe(n,Xp) - d#sigma/d#Omega (mb/sr) for %02.1f #pm %.1f  MeV",En,binMeV));
        gCosAngle[i]->SetName(Form("ddOmega_%s_%dp%d_MeV_cos",target_details.c_str(),intEn,decEn));
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
    XSgr->SetName(Form("%s_DDX", target_details.c_str()));
    XSgr->SetTitle("Fe(n,Xp) - d^{2}#sigma/d#OmegadE (mb/sr/MeV) vs Angle and Energy");
    XSgr->GetYaxis()->SetTitle("Angle (deg)");
    XSgr->GetXaxis()->SetTitle("Energy (MeV)");
    XSgr->GetZaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
    XSgr->Write();
    TCanvas *cXS = new TCanvas("cXS","cXS");
    XSgr->Draw();

    ff->Close();

    timer.Print();
}