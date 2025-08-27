// script to calculate the energies deposited in the different parts of te detector
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
#include "/mnt/medley/LucasAnalysis/useful.h" //version 2024.11.25.001

#include "TStopwatch.h"
#include <time.h>



Double_t Gaussian(Double_t *x, Double_t *par) {
    return par[0] *  TMath::Exp( -0.5 * TMath::Power( (x[0] - par[1])/par[2] ,2)  + par[3] );
}

TTree *tx;



void Calibration_3alpha_test(){

    string namefile = "cloned_branch.root";

    TFile *f =new TFile(namefile.c_str(),"READ");
    TTree *tx = (TTree*) f->Get("AD");


    TH1D *h = new TH1D("h","h",420,4.8,6);
    tx->Draw("Medley_1_dE2>>h");



    float threshold = 0.25;
    

    int N;
    
    Double_t *x;
    Double_t *y;



    TSpectrum *s = new TSpectrum();
    s->Search(h,2,"",threshold);
    N = s->GetNPeaks();
    x = s->GetPositionX();
    y = s->GetPositionY();

    cout<<Form("N = %d", N)<<endl;
    cout <<"center: "<<  x[0] << " , " << x[1] << " , " << x[2] << endl; 
    cout <<"height: "<<  y[0] << " , " << y[1] << " , " << y[2] << endl<<endl; 
    float bin_threshold = 50;
    cout<< "bin threshold: "<< bin_threshold <<endl;


    Int_t bin_center;
    Int_t bn = 0;
    Double_t upper_limit = 0;
    Double_t lower_limit = 0;

    TF1 *func[N];

    for(int i=0;i<N;i++){

        cout << Form("for peak #%d: - - - - - - - - - - -  \n\n",i+1)<<endl;

        cout<<Form("x[%d] = %f",i,x[i])<<endl;
        cout<<Form("y[%d] = %f",i,y[i])<<endl;
        bin_center = h->FindBin(x[i]);
        
        
        //finding upper limit
        bn = bin_center;
        cout << Form("bin center: %d", bin_center)<<endl;
        while(h->GetBinContent(bn)>bin_threshold){
            bn++;
            if(bn == h->GetNbinsX()){
                cout << "histo never below threshold!"<<endl;
                break;
            }
        }

        cout<<"upper limit for this peak: "<< h->GetBinCenter(bn) <<endl;
        upper_limit = h->GetBinCenter(bn);
        
        //finding lower limit 85% 
        bn = bin_center;
        while(h->GetBinContent(bn)>y[i]*0.75){
            bn--;
            if(bn == 1){
                cout << "histo never below threshold!"<<endl;
                break;
            }
        }
        cout<<"lower limit for this peak: "<< h->GetBinCenter(bn) <<endl;
        lower_limit = h->GetBinCenter(bn);
        
        func[i] = new TF1(Form("func%d",i),Gaussian,lower_limit,upper_limit ,4);
        func[i]->SetParameters(y[i],x[i],0.05,0);

        h->Fit(func[i],"MR+");

        cout << " - - - - - - - - - - - - - "<<endl;

    }



    return;


}