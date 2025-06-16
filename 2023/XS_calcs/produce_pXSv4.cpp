#include <vector>
#include <iostream>
#include "TStopwatch.h"
#include "include/xs_function.cxx"

//define binMeV in a HH file 
//std::vector<std::vector<std::vector<TH1D*>>> allHistos;

void produce_pXSv4(Float_t Ea = 28, Float_t Eb = 29, string filename = "Fe_p_DDX.root"){ 


    TStopwatch timer;
    std::vector<int> runsForward = {396,405};
    std::vector<int> runsBackward = {410,414};



    // std::vector<std::vector<TH1D*>> GetDDX(
    //     float Ea = 25,
    //     float Eb = 26,
    //     const std::vector<float>& angles =  {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0},
    //     const std::vector<char>& particles = {'p', 'd', 't', 'h', 'a'},
    //     const std::vector<int>& runsForward = {},
    //     const std::vector<int>& runsBackward = {},
    //     bool Match_CORR = true,
    //     bool TTC_CORR = true,
    //     string target = "Fe",
    //     bool DrawFlux = false
    // );

   //std::vector<std::vector<TH1D*>> hp = GetDDX(Ea,Eb,{20,40,60,80,100,120,140,160},{'p'},runsForward,runsBackward,true);
   
   std::vector<float> angles = {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0};
   std::vector<std::vector<TH1D*>> hp = GetDDX(Ea,Eb,angles,{'p'},runsForward,runsBackward,true);

   TFile *ff = new TFile(filename.c_str(), "UPDATE");

   TGraphErrors *gAngle[hp.size()];
   TGraphErrors *gCosAngle[hp.size()];
   Double_t integral;// error;

    //Since hp is given by hp[particle][angle], we need to loop over he elements to Writhe the histograms to the file
    for (size_t i = 0; i < hp.size(); ++i) {//particle
        //for each particle we will have a graph as a function of the angles 
        gAngle[i] = new TGraphErrors();
        gCosAngle[i] = new TGraphErrors();
        for (size_t j = 0; j < hp[i].size(); ++j) {//angles
            TH1D* h = hp[i][j];
            
            //integral = hp[i][j]->Integral("width");
            Double_t maxint, minint, errorint; //calculating manually the max and min integral and the error   
            maxint = minint = errorint = 0;
            for(int bin = 1; bin <= hp[i][j]->GetNbinsX(); bin++){
                maxint+= (hp[i][j]->GetBinContent(bin)+hp[i][j]->GetBinError(bin) )* hp[i][j]->GetBinWidth(bin);
                minint+= (hp[i][j]->GetBinContent(bin)-hp[i][j]->GetBinError(bin) )* hp[i][j]->GetBinWidth(bin);
            }
            errorint = fabs(maxint - minint)/2; //error is half the difference between max and min integral

            //integral = hp[i][j]->IntegralAndError(1,hp[i][j]->GetNbinsX(),error,"width");
            integral = hp[i][j]->Integral("width");
            
            // cout<<"Comparison between errors :"<<endl;
            // cout<<"Error from IntegralAndError: " << error << endl;
            // cout<<"Error from manual calculation: " << errorint << endl;
            // cout<<"Integral: " << integral << endl;
            
            // std::cout << "Pressione ENTER para continuar...";
            // std::cin.get(); 

            gAngle[i]->SetPoint(j, angles[j], integral);
            gAngle[i]->SetPointError(j, 5, errorint);// Assuming a 5deg error on the angle

            gCosAngle[i]->SetPoint(j, cos(TMath::Pi()*angles[j]/180), integral);
            gCosAngle[i]->SetPointError(j, sin(TMath::Pi()*angles[j]/180)*TMath::Pi()*5/180, errorint);// Assuming a 5deg error on the angle
            //couting the plotted values for both plots:
            std::cout << "Particle: " << i << ", Angle: " << angles[j] << " -> Integral: " << integral << ", Error: " << errorint << std::endl;


            if (h) {
                h->Write();
                std::cout << "Wrote histogram: " << h->GetName() << std::endl;
                std::cout << " - - - - - - - -> [" << i << "][" << j << "]"<< std::endl;
            } else {
                std::cout << "Histogram at index [" << i << "][" << j << "] is null." << std::endl;
            }
        }
        

        TCanvas *cc = new TCanvas("cAngle","Differential Cross Section vs Angle",150,150,1600,600);
        if(i == 0 ){
            //cout<<"\n::\n::\n::\n::\n::\n::\n::\n::\n::\n:: Im HERE"<<endl;
            
            cc->Divide(2,1);
            cc->cd(1);
            gAngle[i]->Draw();
            cc->cd(2);
            gCosAngle[i]->Draw();
        }else{
            cc->cd(1);
            gAngle[i]->Draw("same");
            cc->cd(2);
            gCosAngle[i]->Draw("same");
        }
        
        gAngle[i]->GetXaxis()->SetTitle("Angle (degrees)");
        gAngle[i]->GetYaxis()->SetTitle("Differential Cross Section (mb/sr)");
        gCosAngle[i]->GetXaxis()->SetTitle("cos(Angle)");
        gCosAngle[i]->GetYaxis()->SetTitle("Differential Cross Section (mb/sr)");
        gAngle[i]->SetName(Form("gAngle_particle%zu_ENN_%s_%s_MeV", i,DoubleToXpY(Ea).c_str(),DoubleToXpY(Eb).c_str()));
        gCosAngle[i]->SetName(Form("gCosAngle_particle%zu_ENN_%s_%s_MeV", i,DoubleToXpY(Ea).c_str(),DoubleToXpY(Eb).c_str()));
        gAngle[i]->Write();
        gCosAngle[i]->Write();

    }

    
    ff->Close();

    timer.Print();
}