//function to get Pfactor from the TTC files 

#include <sys/stat.h>


TGraph * GetPfactor(
    string filename = "eloss_p_20.0deg_050.0um.root", 
    string path = "/home/pi/ganil/kalscripts/eloss/results/UniformZ/CH2/v7/",
    Float_t nbins = 800,
    bool drawit = false
    ){
    

    TFile *f = new TFile(Form("%s/%s",path.c_str(),filename.c_str()),"READ");
    TTree *t = (TTree*)f->Get("SIM");
    Float_t Pfactor;

    TH1D *hE = new TH1D("hE","hEnergy",nbins,0,40);
    TH1D *hEl = new TH1D("hEl","hEnergyLostParticles",nbins,0,40);
    
    
    if(drawit){
        TCanvas *cc = new TCanvas("cc","cc",1600,600);
        cc->Divide(2,1);
        cc->cd(1);

        t->Draw("E>>hE");
        hE->GetYaxis()->SetRangeUser(0,1.1*hE->GetMaximum());
        t->Draw("E>>hEl","detected<1", "same");
        cc->cd(2);

    }else{
        t->Draw("E>>hE","","goff");
        t->Draw("E>>hEl","detected<1", "goff");

    }
    

    t->Draw("E>>hE","","goff");
    hE->GetYaxis()->SetRangeUser(0,1.1*hE->GetMaximum());
    t->Draw("E>>hEl","detected<1", "goff");


    TGraph *gr = new TGraph();

    for (int i = 0; i < nbins; i++) {
        if(hE->GetBinContent(i) > hEl->GetBinContent(i)){
            Pfactor = hE->GetBinContent(i)/(hE->GetBinContent(i) - hEl->GetBinContent(i));
            gr->AddPoint(hE->GetBinCenter(i),Pfactor);
            //if(drawit) cout<<"Ei = "<<hE->GetBinCenter(i)<<"MeV, Pfactor = "<<Pfactor<<endl;
        }
    }

    gr->SetMarkerStyle(20);
    if(drawit){
        gr->Draw("ALP");
    }
    

return gr;
}

