#include "src/somefunctions.cc"
//#include "/mnt/medley/LucasAnalysis/useful.h"
#include "src/newE.cc"
#include "include/somefunctions.hh"
#include "apply_correction2.cc"
#include <any>

//////
// This script is in its vertions version1.2025-07-14.0

void Medleyflux(string filename = "runs_035-038_CH2.root", 
                string pspecFilename = "pspec_runs_035-038_CH2.root",
                bool savePspec = true,
                string pElossFileName ="/home/pi/ganil/kalscripts/eloss/results/UniformZ/CH2/v7/eloss_p_20.0deg_050.0um.root",  
                string target_mat = "CH2",
                Double_t th = 50
            ) 
{
  

    TFile *runsff = new TFile(filename.c_str(),"READ");

    TChain *tx = (TChain*)runsff->Get("M");
    tx->Print();

    Long64_t nentries = tx->GetEntries();
  


    //Int_t particle; 
            
        //     Float_t energy; 
            
        //     Float_t si1; 
        //     Float_t si2; 
        //     Float_t csi; 
        //     Float_t rawtof; // raw ToF
        //     Float_t tof_measured;
        //     Float_t tof_correction; 
        //     Float_t tofnn; 
        //     Float_t tofpart; 
        //     Float_t tofnn_woc;// without correction 
            
        //     Float_t ENN;
        //     Float_t ENNwoc; 

        //    // TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
        //     TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
            
        //     TBranch * rawtof_Branch = newtree->Branch("rawtof", &rawtof, "rawtof/F");
        //     TBranch * tof_measured_Branch = newtree->Branch("tof_measured", &tof_measured, "tof_measured/F");
        //     TBranch * tof_correction_Branch  = newtree->Branch("tof_correction", &tof_correction, "tof_correction/F");
            
            
        //     TBranch * tofnn_Branch  = newtree->Branch("tofnn", &tofnn, "tofnn/F");
        //     TBranch * tofpart_Branch  = newtree->Branch("tofpart", &tofpart, "tofpart/F");
        //     TBranch * tofnn_woc_Branch  = newtree->Branch("tofnn_woc", &tofnn_woc, "tofnn_woc/F");
            
        //     TBranch * EnBranch  = newtree->Branch("ENN", &ENN, "ENN/F");
        //     TBranch * EnWocBranch  = newtree->Branch("ENNwoc", &ENNwoc, "ENNwoc/F");
            
        //     //debbugging part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        //     TBranch * branch_si1  = newtree->Branch("si1", &si1, "si1/F");
        //     TBranch * branch_si2  = newtree->Branch("si2", &si2, "si2/F");
        //     TBranch * branch_csi  = newtree->Branch("csi", &csi, "csi/F");
            

    TGraph *xs = new TGraph("/mnt/medley/LucasAnalysis/2023/fluxMedley_script/XS/nn.org_np_2to40MeV_LABrf_20.40.60.80deg.csv","%lg %lg %*lg %*lg %*lg",",");

    //we already have Qtot, 
    Double_t Omegatel = 0.0094; //SACALCV3 sr
    //Double_t cos_alpha = cos(TMath::Pi()*45/180);
    double NatH = 2*0.0237*6.02e23/14.0; //2 ⋅ [mass of CH2 ⋅ NA / MM(CH2)]
    double L = 504.4; //cm
    double mb_to_cm2 = 1e-27 ;//cm²


    Int_t nbins = 160;
    TH1D *flux = new TH1D("flux","flux",nbins,0,40);
    TH1D *flux_woc = new TH1D("flux_woc","flux_woc",nbins,0,40);
    

    tx->Draw("ENNwoc>>flux_woc");
    tx->Draw("ENN>>flux","","same");

    flux_woc->SetLineColor(kRed);
    flux->SetLineColor(kBlue);

    TNamed* charge = (TNamed*) runsff->Get("charge_in_uC");
    TNamed* duration = (TNamed*) runsff->Get("duration_of_the_run_in_sec");
    
    double charge_val = atof(charge->GetTitle());
    double duration_val = atof(duration->GetTitle());

    // Agora você tem os valores nas variáveis charge_val e duration_val
    cout << "Charge (uC): " << charge_val << endl;
    cout << "Duration (s): " << duration_val << endl;

    Double_t factor = pow(L,2)/(NatH*mb_to_cm2*Omegatel*charge_val);
    cout<<"Factor = "<<factor<<endl;


    for(Int_t i=1;i<=nbins;i++){
        


        flux->SetBinContent(i,factor*flux->GetBinContent(i)/(flux->GetBinWidth(i)*xs->Eval(flux->GetBinCenter(i))));
        flux_woc->SetBinContent(i,factor*flux_woc->GetBinContent(i)/(flux_woc->GetBinWidth(i)*xs->Eval(flux_woc->GetBinCenter(i))));

    }

    TCanvas *newCv = new TCanvas("ncv","ncv");
    TH1D *fluxEn = new TH1D("fluxEn","fluxEn",80*nbins,0,40);
    tx->Draw("si1>>fluxEn");

    TCanvas *newCvE = new TCanvas("ncvE","ncvE");
    TH1D *totE = newE(fluxEn, "totE", false, 53.4,1,1);
    totE->Rebin(80);
    totE->SetLineColor(kRed);
    totE->SetTitle("Energy distribution after Si1");
    totE->GetXaxis()->SetTitle("Energy (MeV)");
    totE->GetYaxis()->SetTitle("Counts");
    totE->Draw();

    double mp = 938.27208816;
    double mn = 939.56542052;
    double factX = pow(mp+mn,2)/(4*mn*mp*pow(cos(20.0*TMath::Pi()/180),2));
    cout<<"factX = "<<factX<<endl;
    TH1D *totEx = ScaleXhisto(totE, factX, "totEx", true);
    totEx->SetLineColor(kBlue);
    totEx->Draw("same");

    TH1D *totExTTC = apply_correction2(totEx,
                                        pElossFileName,
                                        target_mat,
                                        th);
    totExTTC->SetName("totExTTC");
    totExTTC->SetLineColor(kGreen);
    newCvE->cd();
    totExTTC->Draw("same");

    if(savePspec){
        TFile *pspecFile = new TFile(pspecFilename.c_str(),"RECREATE");
        totExTTC->Write();
        pspecFile->Close();
        cout<<"Saved the energy spectrum to " << pspecFilename << endl;
    } 

    



    TH1D * totExFlux = (TH1D *)totExTTC->Clone("totExFlux");
    totExFlux->SetTitle("Flux after Si1");
    totExFlux->GetXaxis()->SetTitle("Energy (MeV)");
    
    nbins = totEx->GetNbinsX();
    
    for(Int_t i=1;i<=nbins;i++){
    
        totExFlux->SetBinContent(i,factor*totExFlux->GetBinContent(i)/(totExFlux->GetBinWidth(i)*xs->Eval(totExFlux->GetBinCenter(i))));

    }

    TCanvas *cFlux = new TCanvas("cFlux","cFlux");
    totExFlux->Draw();
    //tx->Draw()



}

TH1D *subtractFlux(TH1D *hch2, TH1D *hC, Float_t Qch2 = 445900,Float_t Qc = 764700) {
    if (!hch2 || !hC) {
        cerr << "Error: One or both histograms are null." << endl;
        return nullptr;
    }

    if (hch2->GetNbinsX() != hC->GetNbinsX()) {
        cerr << "Error: Histograms have different number of bins." << endl;
        return nullptr;
    }

    double NatH = 2*0.0237*6.02e23/14.0; //2 ⋅ [mass of CH2 ⋅ NA / MM(CH2)]
    double NaC = 0.0408*6.02e23/12.0; //  mass of C * NA / MM(C)

    double factor = (NatH*Qch2)/(Qc*NaC);
    cout<< "Factor for subtraction: " << factor << endl;

    TH1D *result = (TH1D *)hch2->Clone("result");
    result->Add(hC, -factor); // Subtract h2 from h1

    return result;
}