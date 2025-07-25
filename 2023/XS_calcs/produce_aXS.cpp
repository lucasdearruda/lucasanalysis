#include <vector>
#include <iostream>
#include "TStopwatch.h"
#include "include/xs_function.cxx"
void produce_aXS(Float_t Ea = 28, Float_t Eb = 29, string filename = "Fe_a_DDX.root");

//define binMeV in a HH file 
//std::vector<std::vector<std::vector<TH1D*>>> allHistos;
void CloseAllRootFiles() {
    TIter next(gROOT->GetListOfFiles());
    TObject* obj;
    while ((obj = next())) {
        TFile* file = dynamic_cast<TFile*>(obj);
        if (file && file->IsOpen()) {
            file->Close();
        }
    }
    gROOT->GetListOfFiles()->Clear();  // Optional
}

void RunSeveral_vBin(){
    cout<<"Running RunSeveral with variable binning..."<<endl;


    // Lista de valores
    std::vector<Float_t> vals = {
        3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5,
        14.5, 15.5, 16.5, 17.6, 18.8, 20.1, 21.6, 23.2, 25, 27,
        29.3, 32, 35, 38.4, 40
    };

    for (size_t i = 0; i < vals.size() - 1; i++) {
        Float_t Ea = vals[i];
        Float_t Eb = vals[i+1];

        cout << "Calling produce_aXS with Ea = " << Ea << ", Eb = " << Eb << endl;

        // Se quiser, pode modificar o filename com os valores, ou deixar fixo
        produce_aXS(Ea, Eb, "nightowl/Fe_a_DDX.root");
    }

}

void produce_aXS(Float_t Ea = 28, Float_t Eb = 29, string filename = "Fe_a_DDX.root"){ 


    TStopwatch timer;
    //std::vector<int> runsForward = {396,405};
    std::vector<int> runsBackward = {410,412,383,384};
    std::vector<int> runsForward = {396,405};
   // std::vector<int> runsBackward = {410,412};


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
   std::vector<std::vector<TH1D*>> hp = GetDDX(Ea,Eb,angles,{'a'},runsForward,runsBackward,true,false);

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
            integral = maxint = minint = errorint = 0;
            for(int bin = 1; bin <= hp[i][j]->GetNbinsX(); bin++){
                maxint+= (hp[i][j]->GetBinContent(bin)+hp[i][j]->GetBinError(bin) )* hp[i][j]->GetBinWidth(bin);
                minint+= (hp[i][j]->GetBinContent(bin)-hp[i][j]->GetBinError(bin) )* hp[i][j]->GetBinWidth(bin);
                integral += hp[i][j]->GetBinContent(bin) * hp[i][j]->GetBinWidth(bin);
            }
            errorint = fabs(maxint - minint)/2; //error is half the difference between max and min integral

            //integral = hp[i][j]->IntegralAndError(1,hp[i][j]->GetNbinsX(),error,"width");
            
            
    //integral = hp[i][j]->Integral("width");
            
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

    string cur_time = getCurrentTime();
    TNamed *processing_info= new TNamed("processed by /mnt/medley/LucasAnalysis/2023/XS_calcs/produce_aXSv4.cpp:version4.2025-06-17.0",cur_time);
    processing_info->Write();
    ff->Close();
    CloseAllRootFiles();
    timer.Print();
}