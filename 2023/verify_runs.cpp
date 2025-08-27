//verify_runs.cpp, version 2025-08-25.1
// Script for veryfing the runs for 2023, so I can use in my calculations 

#include <iostream>     // for std::cout, std::cerr, std::endl
#include <fstream>      // for std::ifstream
#include <string>       // for std::string
#include <filesystem> // C++17
namespace fs = std::filesystem;

//for the time:
#include <time.h>
//for the vectors and pairs:
#include <vector>
#include <utility> 
#include "TVector3.h"

using namespace std;

#include "/mnt/medley/LucasAnalysis/2023/src/loadcals.cc" //for loading calibrations:


int maxTel = 8; //max telescope number
//cuts for npt particles
TCutG *cut_protons[8], *cut_deuterons[8], *cut_tritons[8]; 
TCutG *cut_he[8], *cut_alphas[8], *cut_he34[8], *cut_dt[8], *cutH[8];

//cuts for pt particles
TCutG *cut_protonsPT[8], *cut_deuteronsPT[8], *cut_tritonsPT[8];
TCutG *cut_protonsCsI[8], *cut_deuteronsCsI[8], *cut_tritonsCsI[8]; 

//cut npt csi 
TCutG *cut_protonsNPTCsI[8], *cut_deuteronsNPTCsI[8];

Float_t GetGflash(TTree *tx, Int_t telN = 1, Int_t runN=0,Float_t sigma= 4,Float_t threshold= 0.002,float guess = 400,float tofmin = 300,float tofmax = 500, bool closecanvas = false, bool saveIt=true){

	TCanvas *Ctof = new TCanvas("ctof",Form("Time of flight"),50,50,600,600);
    string branchname=Form("Medley_%d_dE2_ToF",telN);

    TH1D *htof;
    htof = new TH1D("htof","htof",200,tofmin,tofmax);

    tx->Draw(Form("%s>>htof",branchname.c_str()), Form("%s>%f && %s<%f ",branchname.c_str(),tofmin,branchname.c_str(),tofmax));

    double tof_peak[2];
    TSpectrum *spec = new TSpectrum();
    //Int_t npeaks = spec->Search(htof, 6,"",0.1);//originally 0.3
	Int_t npeaks = spec->Search(htof, sigma,"",threshold);//originally 0.3
    Double_t *xpeaks = spec->GetPositionX();
    cout<<"ToF plot: "<<npeaks<<" peak found.\n";
    for(int t =0;t<npeaks;t++){
        cout<<"position peak"<<t<<": "<<xpeaks[t]<<" ns.\n";
    }

    if (npeaks >= 2) {
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = xpeaks[npeaks -1] - 10;
        double fitRangeMax = xpeaks[npeaks -1] + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[npeaks -1], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }else{ //if did not find the peak
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = guess - 10;
        double fitRangeMax = guess  + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[0], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }

    if(!closecanvas){
        TLatex *latex = new TLatex();
        latex->SetNDC();
        latex->SetTextColor(kRed);
        latex->SetTextSize(0.04);
        latex->DrawLatex(0.5, 0.8, Form("ToF peak: %.2f ns", tof_peak[0]));


        // Sigma em preto, logo abaixo
        TLatex *latexSigma = new TLatex();
        latexSigma->SetNDC();
        latexSigma->SetTextColor(kBlack);
        latexSigma->SetTextSize(0.035);
        latexSigma->DrawLatex(0.5, 0.75, Form("sigma = %.2f", sigma));

        // Threshold em preto, abaixo do sigma
        TLatex *latexThreshold = new TLatex();
        latexThreshold->SetNDC();
        latexThreshold->SetTextColor(kBlack);
        latexThreshold->SetTextSize(0.035);
        latexThreshold->DrawLatex(0.5, 0.70, Form("threshold = %.4f", threshold));


		gStyle->SetOptStat("e");
		gPad->SetGridx();
		gPad->SetGridy();
		htof->GetYaxis()->SetMaxDigits(3);
		htof->GetYaxis()->SetTitle("counts");
        htof->SetTitle(Form("Run%d - tel %d", runN, telN));
		htof->GetXaxis()->SetTitle("raw (inv)ToF (ns)");
    }
    if(saveIt){
        Ctof->SaveAs(Form("verifications/r%d_t%d.png",runN,telN));
        Ctof->SaveAs(Form("verifications/PDFs/r%d_t%d.pdf",runN,telN));
        Ctof->SaveAs(Form("verifications/ROOTs/r%d_t%d.root",runN,telN));
    }
	if(closecanvas) Ctof->Close();
	
	return tof_peak[0];
}


std::string getCurrentTime() {
    char cur_time[128];
    time_t t = time(NULL);
    struct tm* ptm = localtime(&t);

    // Formata o tempo no formato desejado: "YYYY-MM-DD_HH:MM:SS"
    strftime(cur_time, sizeof(cur_time), "%Y-%m-%d_%H:%M:%S", ptm);

    return std::string(cur_time);
}

void plotPIDall(TTree *Tree, int runN = 0, bool saveCanvas = true){
    loadCals("/mnt/medley/LucasAnalysis/2023/src/calibrations.dat");
    TCanvas *cPID = new TCanvas("cPID", "cPID", 50, 50, 1400, 700);
    cPID->Divide(2, 1);
    
    TH2D *hde12 = new TH2D("hde12", "dE1 vs dE2", 400, 0, 30, 400, 0, 10);  
    hde12->GetYaxis()->SetTitle("dE1 (MeV)");
    hde12->GetXaxis()->SetTitle("dE2 (MeV)");

    TH2D *hde23 = new TH2D("hde23", "dE2 vs Eres", 400, 0, 40, 400, 0, 18);  
    hde23->GetYaxis()->SetTitle("dE2 (MeV)");
    hde23->GetXaxis()->SetTitle("Eres (MeV)");

    for (int telN = 1; telN <= 4; telN++) {
        cPID->cd(1);
        Tree->Draw(Form("Medley_%d_dE1*%f:Medley_%d_dE2*%f>>hde12", telN,g1[telN-1], telN,g2[telN-1]), Form("Medley_%d_dE1>0 && Medley_%d_dE2>0", telN, telN),"colz");
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogz();

        if(telN<=4){
            cut_protons[telN]->Draw("same");
            cut_deuterons[telN]->Draw("same");
            cut_tritons[telN]->Draw("same");
            cut_he[telN]->Draw("same");
            cut_alphas[telN]->Draw("same");
            cut_protonsPT[telN]->Draw("same");
            cut_deuteronsPT[telN]->Draw("same");
            cut_tritonsPT[telN]->Draw("same");
        
        }    // }else{
        //     cut_he34[telN]->Draw("same");
        //     if(telN != 7){
        //         cut_protons[telN]->Draw("same");
        //         cut_dt[telN]->Draw("same");
        //     }else{
        //         cutH[telN]->Draw("same");    
        //     }
            
        // }
        

        cPID->cd(2);
        Tree->Draw(Form("Medley_%d_dE2*%f:Medley_%d_Eres*%f>>hde23", telN,g2[telN-1], telN,g3[telN-1]), Form("Medley_%d_dE2>0 && Medley_%d_Eres>0", telN, telN),"colz");
        gPad->SetGridx();
        gPad->SetGridy();
        gPad->SetLogz();
        
        if(telN<=4){
            cut_protonsCsI[telN]->Draw("same");
            cut_deuteronsCsI[telN]->Draw("same");
            cut_tritonsCsI[telN]->Draw("same");
            cut_protonsNPTCsI[telN]->Draw("same");
            cut_deuteronsNPTCsI[telN]->Draw("same");
        }


        if (saveCanvas) {
            cPID->SaveAs(Form("verifications/r%d_PID_dE_t%d.png",runN, telN));
            cPID->SaveAs(Form("verifications/PDFs/r%d_PID_dE_t%d.pdf",runN, telN));
            cPID->SaveAs(Form("verifications/ROOTs/r%d_PID_dE_t%d.root",runN, telN));
        }
        
    }


}


void plotPID(TTree *Tree, int telN = 1, int runN = 0, bool saveCanvas = true){
    loadCals("/mnt/medley/LucasAnalysis/2023/src/calibrations.dat");
    TCanvas *cPID = new TCanvas("cPID", "cPID", 50, 50, 1400, 700);
    cPID->Divide(2, 1);
    
    TH2D *hde12 = new TH2D("hde12", "dE1 vs dE2", 400, 0, 30, 400, 0, 10);  
    hde12->GetYaxis()->SetTitle("dE1 (MeV)");
    hde12->GetXaxis()->SetTitle("dE2 (MeV)");

    TH2D *hde23 = new TH2D("hde23", "dE2 vs Eres", 400, 0, 40, 400, 0, 18);  
    hde23->GetYaxis()->SetTitle("dE2 (MeV)");
    hde23->GetXaxis()->SetTitle("Eres (MeV)");


    cPID->cd(1);
    Tree->Draw(Form("Medley_%d_dE1*%f:Medley_%d_dE2*%f>>hde12", telN,g1[telN-1], telN,g2[telN-1]), Form("Medley_%d_dE1>0 && Medley_%d_dE2>0", telN, telN),"colz");
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    if(telN<=4){
        cut_protons[telN]->Draw("same");
        cut_deuterons[telN]->Draw("same");
        cut_tritons[telN]->Draw("same");
        cut_he[telN]->Draw("same");
        cut_alphas[telN]->Draw("same");
        cut_protonsPT[telN]->Draw("same");
        cut_deuteronsPT[telN]->Draw("same");
        cut_tritonsPT[telN]->Draw("same");
    
    }    // }else{
    //     cut_he34[telN]->Draw("same");
    //     if(telN != 7){
    //         cut_protons[telN]->Draw("same");
    //         cut_dt[telN]->Draw("same");
    //     }else{
    //         cutH[telN]->Draw("same");    
    //     }
        
    // }
    

    cPID->cd(2);
    Tree->Draw(Form("Medley_%d_dE2*%f:Medley_%d_Eres*%f>>hde23", telN,g2[telN-1], telN,g3[telN-1]), Form("Medley_%d_dE2>0 && Medley_%d_Eres>0", telN, telN),"colz");
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();
    
    if(telN<=4){
        cut_protonsCsI[telN]->Draw("same");
        cut_deuteronsCsI[telN]->Draw("same");
        cut_tritonsCsI[telN]->Draw("same");
        cut_protonsNPTCsI[telN]->Draw("same");
        cut_deuteronsNPTCsI[telN]->Draw("same");
    }


    if (saveCanvas) {
        cPID->SaveAs(Form("verifications/r%d_PID_dE_t%d.png",runN, telN));
        cPID->SaveAs(Form("verifications/PDFs/r%d_PID_dE_t%d.pdf",runN, telN));
        cPID->SaveAs(Form("verifications/ROOTs/r%d_PID_dE_t%d.root",runN, telN));
    }
    
}

void loadCuts(){
    cout<< "Loading cuts for PID verification..." << endl;

string cuts_path =   "/mnt/medley/LucasAnalysis/2023/PIDv6";
bool warn_tels58 = false;
for(int telN = 1; telN<=8; telN++){//just up tel 4 for now 2025-08-13
    if(telN<=4){
        
        //NPT CUTS dEdE (upper bands) =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
            //protons npt
            cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
            cout << Form("\n->[tel %d] loading NPT PROTONS CUT: %s/p_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/p_npt_%d.C", cuts_path.c_str(),telN));
            cut_protons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_npt_%d",telN));
            if(cut_protons[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_protons[telN]->SetName(Form("p_npt_%d",telN));
            cut_protons[telN]->SetLineWidth(2);
            cut_protons[telN]->SetLineColor(kRed);
        
            //deuterons npt
            cout << Form("\n->[tel %d] loading NPT DEUTERONS CUT: %s/d_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/d_npt_%d.C", cuts_path.c_str(),telN));
            cut_deuterons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_npt_%d",telN));
            if(cut_deuterons[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_deuterons[telN]->SetName(Form("d_npt_%d",telN));
            cut_deuterons[telN]->SetLineWidth(2);
            cut_deuterons[telN]->SetLineColor(kGreen);

            //tritons npt
            cout << Form("\n->[tel %d] loading NPT TRITONS CUT: %s/t_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/t_npt_%d.C", cuts_path.c_str(),telN));
            cut_tritons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_npt_%d",telN));
            if(cut_tritons[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_tritons[telN]->SetName(Form("t_npt_%d",telN));
            cut_tritons[telN]->SetLineWidth(2);
            cut_tritons[telN]->SetLineColor(kBlack);
        
        //NPT CUTS dEdE (upper bands) =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
            //helium 
            cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH (he,alpha) CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
            cout << Form("\n->[tel %d] loading NPT HE3 CUT: %s/he3_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/he3_%d.C", cuts_path.c_str(),telN));
            cut_he[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("he3_%d",telN));
            if(cut_he[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_he[telN]->SetName(Form("he3_%d",telN));
            cut_he[telN]->SetLineWidth(2);
            cut_he[telN]->SetLineColor(kRed);
        
            //alphas 
            cout << Form("\n->[tel %d] loading NPT ALPHAS CUT: %s/alphas%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/alphas%d.C", cuts_path.c_str(),telN));
            cut_alphas[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("alphas%d",telN));
            if(cut_alphas[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_alphas[telN]->SetName(Form("alphas%d",telN));
            cut_alphas[telN]->SetLineWidth(2);
            cut_alphas[telN]->SetLineColor(kGreen);

        //PT CUTS in dEdE (bottom band).=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
            //protons
            cout << Form("\n- - -   LOADING PUNCH-THROUGH CUTS FOR dE1:dE2 (BOTTOM BRANCH)   - - -\n")<< "... " ;
            cout << Form("\n->[tel %d] loading PT PROTONS CUT: %s/p_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/p_pt_dE_%d.C", cuts_path.c_str(),telN));
            cut_protonsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_dE_%d",telN));
            if(cut_protonsPT[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_protonsPT[telN]->SetName(Form("p_pt_dE_%d",telN));
            cut_protonsPT[telN]->SetLineWidth(2);
            cut_protonsPT[telN]->SetLineColor(kRed);
            //deuterons:
            cout << Form("\n->[tel %d] loading PT DEUTERONS CUT: %s/d_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/d_pt_dE_%d.C", cuts_path.c_str(),telN));
            cut_deuteronsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_dE_%d",telN));
            if(cut_deuteronsPT[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_deuteronsPT[telN]->SetName(Form("d_pt_dE_%d",telN));
            cut_deuteronsPT[telN]->SetLineWidth(2);
            cut_deuteronsPT[telN]->SetLineColor(kGreen);
            //tritons:
            cout << Form("\n->[tel %d] loading PT TRITONSCUT: %s/t_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/t_pt_dE_%d.C", cuts_path.c_str(),telN));
            cut_tritonsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_dE_%d",telN));
            if(cut_tritonsPT[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_tritonsPT[telN]->SetName(Form("t_pt_dE_%d",telN));
            cut_tritonsPT[telN]->SetLineWidth(2);
            cut_tritonsPT[telN]->SetLineColor(kBlack);

        // //PT csi CUTS =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
            //protons
            cout << Form("\n- - -   LOADING PUNCH-THROUGH CUTS FOR dE2:Eres (CsI)   - - -\n")<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/p_pt_%d.C", cuts_path.c_str(),telN));
            cut_protonsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_%d",telN));
            if(cut_protonsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_protonsCsI[telN]->SetName(Form("p_pt_%d",telN));
            cut_protonsCsI[telN]->SetLineWidth(2);
            cut_protonsCsI[telN]->SetLineColor(kRed);
            //deuterons:
            cout << Form("\n->[tel %d] loading PT DEUTERONS CUT: %s/d_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/d_pt_%d.C", cuts_path.c_str(),telN));
            cut_deuteronsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_%d",telN));
            if(cut_deuteronsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_deuteronsCsI[telN]->SetName(Form("d_pt_%d",telN));
            cut_deuteronsCsI[telN]->SetLineWidth(2);
            cut_deuteronsCsI[telN]->SetLineColor(kGreen);
            //tritons:
            cout << Form("\n->[tel %d] loading PT TRITONS CUT: %s/d_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/t_pt_%d.C", cuts_path.c_str(),telN));
            cut_tritonsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_%d",telN));
            if(cut_tritonsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_tritonsCsI[telN]->SetName(Form("t_pt_%d",telN));
            cut_tritonsCsI[telN]->SetLineWidth(2);
            cut_tritonsCsI[telN]->SetLineColor(kGreen);

        // //NPT csi CUTS =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
            //protons
            cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE2:Eres (CsI)   - - -\n")<< "... " ;
            cout << Form("\n->[tel %d] loading NPT PROTONS (CSI) CUT: %s/pnptcsi_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/pnptcsi_%d.C", cuts_path.c_str(),telN));
            cut_protonsNPTCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("pnptcsi_%d",telN));
            if(cut_protonsNPTCsI[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_protonsNPTCsI[telN]->SetName(Form("pnptcsi_%d",telN));
            cut_protonsNPTCsI[telN]->SetLineWidth(2);
            cut_protonsNPTCsI[telN]->SetLineColor(kRed);
            //deuterons:
            cout << Form("\n->[tel %d] loading NPT DEUTERONS CUT: %s/dnptcsi_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/dnptcsi_%d.C", cuts_path.c_str(),telN));
            cut_deuteronsNPTCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("dnptcsi_%d",telN));
            if(cut_deuteronsNPTCsI[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_deuteronsNPTCsI[telN]->SetName(Form("dnptcsi_%d",telN));
            cut_deuteronsNPTCsI[telN]->SetLineWidth(2);
            cut_deuteronsNPTCsI[telN]->SetLineColor(kGreen);
    }else{
        //     if(!warn_tels58){
        //         cout << "\n . . .  - - - - - - - - - - - - - -   . . . \n . . . - - - TELESCOPES 5 - 8  - - -  . . . "<<"\n . . .  - - - - - - - - - - - - - -   . . . "<<endl;
        //         warn_tels58 = true;
        //     }
        //     if(telN != 7 ){
        //         //cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
        //         cout << Form("\n->[tel %d] loading NPT PROTONS CUT: %s/p_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        //         gROOT->ProcessLine(Form(".L  %s/p_npt_%d.C", cuts_path.c_str(),telN));
        //         cut_protons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_npt_%d",telN));
        //         if(cut_protons[telN]!=NULL) cout << "  --> OK."<< endl;
        //         cut_protons[telN]->SetName(Form("p_npt_%d",telN));
        //         cut_protons[telN]->SetLineWidth(2);
        //         cut_protons[telN]->SetLineColor(kRed);

        //         cout << Form("\n->[tel %d] loading NPT DEUTERON-TRITONS CUT: %s/dt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        //         gROOT->ProcessLine(Form(".L  %s/dt_%d.C", cuts_path.c_str(),telN));
        //         cut_dt[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("dt_%d",telN));
        //         if(cut_dt[telN]!=NULL) cout << "  --> OK."<< endl;
        //         cut_dt[telN]->SetName(Form("dt_%d",telN));
        //         cut_dt[telN]->SetLineWidth(2);
        //         cut_dt[telN]->SetLineColor(kBlue);

        //     }else{
        //         cout << Form("\n->[tel %d] loading H(1, 2 or 3) CUT: %s/h_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        //         gROOT->ProcessLine(Form(".L  %s/h_%d.C", cuts_path.c_str(),telN));
        //         cutH[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("h_%d",telN));
        //         if(cutH[telN]!=NULL) cout << "  --> OK."<< endl;
        //         cutH[telN]->SetName(Form("h_%d",telN));
        //         cutH[telN]->SetLineWidth(2);
        //         cutH[telN]->SetLineColor(kRed);
        //     }
        //     cout << Form("\n->[tel %d] loading He(3 or 4) CUTS: %s/he_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        //     gROOT->ProcessLine(Form(".L  %s/he_%d.C", cuts_path.c_str(),telN));
        //     cut_he34[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("he_%d",telN));
        //     if(cut_he34[telN]!=NULL) cout << "  --> OK."<< endl;
        //     cut_he34[telN]->SetName(Form("he_%d",telN));
        //     cut_he34[telN]->SetLineWidth(2);
        //     cut_he34[telN]->SetLineColor(8);
         cout<<". 5 6 7 8 ."<<endl;
         }
}
}

void verify_runs(int runa = -1,int runb = -1, float sigma = 4, float threshold = 0.002, float guess = 400, float tofmin = 300, float tofmax = 560){
    
    string cur_time = getCurrentTime();
    clock_t tStart = clock();
    
    loadCals();
    loadCuts();
    //This function is used to verify the runs and their charges
    //It will print the runs and their charges to the console
    string path_to_runs = "/mnt/medley/RootA_2023";

    
    cout << "Verifying runs..." << endl;
    std::regex run_regex("r([0-9]+)_.*\\.root"); // pega o número da run



for (const auto &entry : fs::directory_iterator(path_to_runs)) {
    if (entry.is_regular_file() && entry.path().extension() == ".root") {
        string filename = entry.path().string();
        string basename = entry.path().filename().string();
        std::smatch match;

        if (std::regex_match(basename, match, run_regex)) {
            int runNumber = stoi(match[1].str());

            // Verifica se está no intervalo
            bool processFile = false;
            if (runa == -1 && runb == -1) {
                processFile = true; // processa todos
            } else if (runa > 0 && runb >= runa && runNumber >= runa && runNumber <= runb) {
                processFile = true; // processa somente dentro do intervalo
            }

            if (!processFile) continue;

            cout << "Processing run: " << runNumber << " (" << filename << ")" << endl;

            // Aqui você abre o TFile, pega o TTree e chama GetGflash pros 8 telescopios
            TFile *file = TFile::Open(filename.c_str(), "READ");
            if (!file || file->IsZombie()) {
                cerr << "Error opening file: " << filename << endl;
                continue;
            }

            TTree *tx = (TTree*)file->Get("AD");
            if (!tx) {
                cerr << "No tree found in file: " << filename << endl;
                file->Close();
                continue;
            }

            for (int tel = 1; tel <= 8; ++tel) {
                GetGflash(tx, tel, runNumber, sigma, threshold, guess, tofmin,tofmax);
                plotPID(tx, tel, runNumber, true);
            }
            //plotPID(tx, runNumber, true);

            file->Close();
        }
    }
}
cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;



}