//process the run for the reduced run. It is based on the version 2;
//in this version, we use the classic PID (ΔE1 . ΔE2 + ΔE2.Eres ) and also add info for the measured ToF as well as the correction needed.
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSpectrum.h"
#include "TRandom2.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"
#include "TGraph.h"
#include "TFile.h"
#include "TGraphErrors.h"
//#include "TStopwatch.h"

#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

// Inclua o functions.hh, mas o root-config cuida das libs do ROOT, então nada redundante
#include "/home/e802/Analysis/LucasAnalysis/include/functions.hh"

int main(int argc, char* argv[]) {

    

    string cur_time = getCurrentTime();
    clock_t tStart = clock();

    cout<<"Starting execution in "<<cur_time<<endl;

    // Check if the correct number of arguments are provided
    if (argc < 2 || argc > 3) {
    	std::cerr << "Usage: " << argv[0] << " nrun [percentage]" << std::endl;
    	return 1; // Return error code
	}

	// Assuming the first argument is the run number
	int nrun = std::stoi(argv[1]);
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    // Check if the optional percentage argument is provided
    double percentage = 100.0; // Default value
    if (argc == 3) {
        try {
            percentage = std::stod(argv[2]);
            if (percentage <= 0 || percentage > 100) {
                std::cerr << "Warning: Percentage must be between 0 and 100. Defaulting to 100." << std::endl;
                percentage = 100.0; // Set to default value
            }
        } catch (const std::exception& e) {
            std::cerr << "Invalid percentage argument: " << e.what() << ". Defaulting to 100." << std::endl;
            percentage = 100.0; // Set to default value
        }
    }

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    //Create the TGraphs for correction:
    TGraphErrors *correctionP, *correctionD;

    string correctionsPath = "/home/e802/Analysis/pre_analysis/runs/Corrections";
    
    TFile *file = new TFile(Form("%s/dcorrection.root", correctionsPath.c_str()), "READ");
    TFile *fileP = new TFile(Form("%s/pcorrection.root", correctionsPath.c_str()), "READ");

    correctionD = (TGraphErrors *)file->Get("dcorrection");
    correctionP = (TGraphErrors *)fileP->Get("pcorrection");

    double Ppoint = correctionP->GetPointX(correctionP->GetN()-1);
    double Dpoint = correctionD->GetPointX(correctionD->GetN()-1);

    cout<<"Ppoint = "<<Ppoint<<" MeV. "<<endl;
    cout<<"Dpoint = "<<Dpoint<<" MeV. "<<endl;
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //CREATE THE REDUCED TTREE

    const char * outputFileName  = Form("/home/e802/Analysis/LucasAnalysis/reducedData/%03ddeb_v3.root",nrun);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
    TTree *newtree = new TTree("M","Medley's reduced data TTree");

    Int_t particle; 
    
    Float_t energy; 
    
    Float_t si1; 
    Float_t si2; 
    Float_t csi; 
    Float_t tof_measured;
    Float_t tof_correction; 
    Float_t tofn; 
    Float_t ang; 
    Float_t ENN; 

    TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");

    TBranch * tof_measured_Branch = newtree->Branch("tof_measured", &tof_measured, "tof_measured/F");
    TBranch * tof_correction_Branch  = newtree->Branch("tof_correction", &tof_correction, "tof_correction/F");

    TBranch * tofBranch  = newtree->Branch("tofn", &tofn, "tofn/F");
    TBranch * angleBranch  = newtree->Branch("ang", &ang, "ang/F");
    TBranch * EnBranch  = newtree->Branch("ENN", &ENN, "ENN/F");
    
    //debbugging part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    TBranch * branch_si1  = newtree->Branch("si1", &si1, "si1/F");
    TBranch * branch_si2  = newtree->Branch("si2", &si2, "si2/F");
    TBranch * branch_csi  = newtree->Branch("csi", &csi, "csi/F");
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //open run
    char name[1000];
    TChain *tx = NULL;

    if (!tx)
        tx = new TChain("AD");

    vector<Int_t> entries;
    Bool_t fExist = true;
    for (Int_t j = 0; j <= 999 && fExist; j++)
    {
        ifstream mfile;
        sprintf(name, "/data/e802X/e802/acquisition/RootA/r%04d_%03da.root", nrun, j);
        mfile.open(name);
        if (mfile)
        {
            mfile.close();
        }
        else
        {
            fExist = false;
        }
        if (fExist)
        {
            cout << "Adding " << name << endl;
            tx->Add(name);
            entries.push_back((Int_t)(tx->GetEntries()));
            cout << "Entries: " << entries.back() / 1000 << "k " << endl;
        }
    }

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //OPENING CUTS  FOR PID:

    //cuts for npt particles
    TCutG *cut_protons[5], *cut_deuterons[5], *cut_tritons[5]; 
    TCutG *cut_he[5], *cut_alphas[5];
    
    //cuts for pt particles
    TCutG *cut_protonsPT[5], *cut_deuteronsPT[5], *cut_tritonsPT[5];
    TCutG *cut_protonsCsI[5], *cut_deuteronsCsI[5], *cut_tritonsCsI[5]; 

    //cut npt csi 
    TCutG *cut_protonsNPTCsI[5], *cut_deuteronsNPTCsI[5];

    string cuts_path =   "/home/e802/Analysis/LucasAnalysis/PIDcuts_debugging2";
    
    for(int telN = 1; telN<=4; telN++){
        
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
        cut_tritonsPT[telN]->SetLineColor(kGreen);

    // //PT csi CUTS =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //protons
        cout << Form("\n- - -   LOADING PUNCH-THROUGH CUTS FOR dE2:Eres (CsI)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading PT PROTONS (CSI) CUT: %s/p_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
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
    
    }


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    ///events processing

    Float_t Gflash[5];
    Float_t tgamma[5];

    for(int n=1;n<=4;n++){
        //Fitting Gamma flash 

        Gflash[n] = GetGflash(tx,Form("Medley_%d_dE2_ToF",n),500,true); // changing it to false, it will let the fitting plot opened to check 
        tgamma[n] = provideTgama(n);
    }


    Float_t L[9];
    //distances Medley target --> Telescope
    L[1] = 157.6*1e-3;
    L[2] = 157.8*1e-3;
    L[3] = 159.0*1e-3;
    L[4] = 160.5*1e-3;
    L[5] = 158.5*1e-3;
    L[6] = 155.8*1e-3;
    L[7] = 156.5*1e-3;
    L[8] = 157.4*1e-3;

    // 1: proton
    // 2: deuteron
    // 3: triton
    // 4: helium-3
    // 5: alpha
    // 6: neutron -- not used of course (but the map follows it in the function)

    Float_t invTof;
    Float_t de2,de1,Eres;
    de2 = 0;
    de1 = 0;
    Eres = 0;

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    Long64_t nEntries = tx->GetEntries();

    cout<<" Total entries : "<<nEntries<<endl; 

    for (Long64_t entry = 0; entry < nEntries*percentage/100.0; entry++) {
        //get entry 
        tx->GetEntry(entry);
        particle = 0; //particle starts unidentified
        tofn = 0; //tof starts zero too 

        for(int telN = 1; telN <= 4; telN++){ // FOR EACH TELESCOPE:
            
            //critical verifications: 
            invTof = tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue();
            
            if(invTof > Gflash[telN]) break; //if the event is faster than gamma... go to next telecope

            de1 =  tx->GetLeaf(Form("Medley_%d_dE1",telN))->GetValue();
            de2 =  tx->GetLeaf(Form("Medley_%d_dE2",telN))->GetValue();
            Eres =  tx->GetLeaf(Form("Medley_%d_Eres",telN))->GetValue();
            energy = de1+de2+Eres;

            si1 = de1;
            si2 = de2;
            csi = Eres;
            if( de1 >= 0.02 && de2 >= 0.02 ){ // if it is not trash -- energy threshold

                if( 
                    (cut_protons[telN]->IsInside(de2,de1) && cut_protonsNPTCsI[telN]->IsInside(Eres,de2))
                    ||
                    (cut_protonsPT[telN]->IsInside(de2,de1) && cut_protonsCsI[telN]->IsInside(Eres,de2)) 
                ){
                    //it is a proton!
                    particle = 1;

                }else if(
                    (cut_deuterons[telN]->IsInside(de2,de1) && cut_deuteronsNPTCsI[telN]->IsInside(Eres,de2))
                    ||
                    (cut_deuteronsPT[telN]->IsInside(de2,de1)&&cut_deuteronsCsI[telN]->IsInside(Eres,de2)) 
                ){
                    //it is a deuteron!
                    particle = 2;

                }else if(
                    cut_tritons[telN]->IsInside(de2,de1)
                    ||
                    ( cut_tritonsPT[telN]->IsInside(de2,de1)   &&   cut_tritonsCsI[telN]->IsInside(Eres,de2) )
                ){
                    //it is a triton!
                    particle = 3;

                }else if(cut_he[telN]->IsInside(de2,de1)){
                    //it is a helium!
                    particle = 4;

                }else if(cut_alphas[telN]->IsInside(de2,de1)){
                    //it is an alpha!
                    particle = 5;

                }
            }

            //if a particle was identified

            if(particle){
                
                ang = 20.*telN;
                //evaluate tof normally
                tof_measured = Gflash[telN] - invTof + tgamma[telN];
                tofn = tof_measured - ToFparticleNumber(energy,particle,L[telN]);
                tof_correction = 0;
                //if needed you update for the corrected value 
                if((particle == 1) && (energy <=Ppoint)){// if it is proton and is inside the correction region
                    //cout<<"\nproton. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionP->Eval(energy)<<" ns."<<endl;
                    tof_correction = correctionP->Eval(energy);
                    tofn = tofn + tof_correction;
                    //cout<<"proton. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }else if((particle == 2)  && (energy <=Dpoint)){// if it is deuteron and is inside the correction region
                    //cout<<"\ndeuteron. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionD->Eval(energy)<<" ns."<<endl;
                    tof_correction =correctionD->Eval(energy);
                    tofn = tofn + tof_correction;
                    //cout<<"deuteron. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }

                ENN = En(tofn);
                break; //continue to next entry
            }


        }//finished looping over telescopes

        if(particle){ // if it was identified we write it down, otherwise we do not do anything 
            
            newtree->Fill();
        }


        // Print the progress
        std::cout << "Progress ["<< nrun <<"]: " << 100.0*entry/(nEntries*percentage/100.0) << "%\r";
        std::cout.flush();
        
    }



    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


    //some information to add
        TNamed *percentage_proc= new TNamed("processed (percentage): ",to_string(percentage));
        TNamed *processing_info= new TNamed("processed by processRun_debugging3",cur_time);
        TNamed *charge_in_uC = NULL;
        TNamed *MedleyTarget= NULL;
        TNamed *duration_of_the_run_in_sec= NULL;
        TNamed *configuration_telescopes= NULL;
        TNamed *gamma_flash_in_ns_T1 = new TNamed("gFlash_T1(ns)",to_string(Gflash[1]).c_str());
        TNamed *gamma_flash_in_ns_T2 = new TNamed("gFlash_T2(ns)",to_string(Gflash[2]).c_str());
        TNamed *gamma_flash_in_ns_T3 = new TNamed("gFlash_T3(ns)",to_string(Gflash[3]).c_str());
        TNamed *gamma_flash_in_ns_T4 = new TNamed("gFlash_T4(ns)",to_string(Gflash[4]).c_str());


        TNamed *dist_Medleytarget_T1_in_mm = new TNamed("dist_Medleytarget_T1(mm)",to_string(1e3*L[1]).c_str());
        TNamed *dist_Medleytarget_T2_in_mm= new TNamed("dist_Medleytarget_T2(mm)",to_string(1e3*L[2]).c_str());
        TNamed *dist_Medleytarget_T3_in_mm= new TNamed("dist_Medleytarget_T3(mm)",to_string(1e3*L[3]).c_str());
        TNamed *dist_Medleytarget_T4_in_mm= new TNamed("dist_Medleytarget_T4(mm)",to_string(1e3*L[4]).c_str());

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


    string infofile = "/home/e802/Analysis/LucasAnalysis/reducedData/MEDLEY2023.csv";

    TTree *InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(infofile.c_str(), "RunN/F:MedleyConfig/C:Target/C:RunTime/I:RunCharge/F");

        // Declare variables to hold branch data
        Float_t RunN;
        char MedleyConfig[256]; // Assuming maximum length of MedleyConfig is 255 characters
        char Target[256]; // Assuming maximum length of Target is 255 characters
        Int_t RunTime;
        Float_t RunCharge;


        // Set branch addresses
        InfoTree->SetBranchAddress("RunN", &RunN);
        InfoTree->SetBranchAddress("MedleyConfig", &MedleyConfig);
        InfoTree->SetBranchAddress("Target", &Target);
        InfoTree->SetBranchAddress("RunTime", &RunTime);
        InfoTree->SetBranchAddress("RunCharge", &RunCharge);


        // Loop over entries in the tree
        Long64_t nEntriesInfo = InfoTree->GetEntries();
        for (Long64_t i = 0; i < nEntriesInfo; i++) {
            InfoTree->GetEntry(i);

            if(RunN == nrun){
                if(RunCharge>0){
                    charge_in_uC = new TNamed("RunCharge(C)",Form("%s",to_string(RunCharge*percentage/100.).c_str()));
                }else{
                    charge_in_uC = new TNamed("RunCharge(C)","not_available");
                }
                configuration_telescopes = new TNamed("MedleyConfig",Form("%s",MedleyConfig));
                MedleyTarget = new TNamed("MedleyTarget",Form("%s",Target));
                if(RunTime>0){
                    duration_of_the_run_in_sec = new TNamed("RunTime(s)",to_string(RunTime*percentage/100.).c_str());
                }else{
                    duration_of_the_run_in_sec = new TNamed("RunTime(s)","not_available");
                }
                
            }

            // Print or use the information from each entry
            // std::cout << "Entry " << i << ":\n";
            // std::cout << "RunN: " << RunN << "\n";
            // std::cout << "MedleyConfig: " << MedleyConfig << "\n";
            // std::cout << "Target: " << Target << "\n";
            // std::cout << "RunTime: " << RunTime << "\n";
            // std::cout << "RunCharge: " << RunCharge << "\n";
        }
    
    percentage_proc->Write();
    processing_info->Write();
    charge_in_uC->Write();
    MedleyTarget->Write();
    duration_of_the_run_in_sec->Write();
    configuration_telescopes->Write();
    gamma_flash_in_ns_T1->Write();
    gamma_flash_in_ns_T2->Write();
    gamma_flash_in_ns_T3->Write();
    gamma_flash_in_ns_T4->Write();


    dist_Medleytarget_T1_in_mm->Write();
    dist_Medleytarget_T2_in_mm->Write();
    dist_Medleytarget_T3_in_mm->Write();
    dist_Medleytarget_T4_in_mm->Write();


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    TH1D * hp[5],*hd[5],*ht[5],*hh[5],*ha[5];
    TH2D * hpn[5],*hdn[5],*htn[5],*hhn[5],*han[5];


    for(int telN = 1; telN <= 4; telN++){
        hp[telN] = new TH1D(Form("protonsE%d",telN),Form("total energy for protons Tel %d", telN), 1000,0,50);
        hd[telN] = new TH1D(Form("deuteronsE%d",telN),Form("total energy for deuterons Tel %d", telN), 1000,0,50);
        ht[telN] = new TH1D(Form("tritonsE%d",telN),Form("total energy for tritons Tel %d", telN), 1000,0,50);
        hh[telN] = new TH1D(Form("heliumE%d",telN),Form("total energy for helium3 Tel %d", telN), 1000,0,50);
        ha[telN] = new TH1D(Form("alphasE%d",telN),Form("total energy for alphas Tel %d", telN), 1000,0,50);


        hpn[telN] = new TH2D(Form("protonsE_En%d",telN),Form("Eprotons:ENN_Tel %d", telN), 1000,0,50, 1000,0,50);
        hdn[telN] = new TH2D(Form("deuteronsE_En%d",telN),Form("Edeuterons:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
        htn[telN] = new TH2D(Form("tritonsE_En%d",telN),Form("Etritons:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
        hhn[telN] = new TH2D(Form("heliumE_En%d",telN),Form("Ehelium:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
        han[telN] = new TH2D(Form("alphasE_En%d",telN),Form("Ealphas:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
    
        newtree->Draw(Form("energy>>protonsE%d",telN),Form("PID==1&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy>>deuteronsE%d",telN),Form("PID==2&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy>>tritonsE%d",telN),Form("PID==3&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy>>heliumE%d",telN),Form("PID==4&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy>>alphasE%d",telN),Form("PID==5&&ang==%.1f",telN*20.0));

        newtree->Draw(Form("energy:ENN>>protonsE_En%d",telN),Form("PID==1&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy:ENN>>deuteronsE_En%d",telN),Form("PID==2&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy:ENN>>tritonsE_En%d",telN),Form("PID==3&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy:ENN>>heliumE_En%d",telN),Form("PID==4&&ang==%.1f",telN*20.0));
        newtree->Draw(Form("energy:ENN>>alphasE_En%d",telN),Form("PID==5&&ang==%.1f",telN*20.0));

    }

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    outputFile->Write();
    //
    // Close the output file
    outputFile->Close();


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
    //timer.Print();

    cout<<"Compiled from processRun_debugging3.cpp"<<endl;
    return 0;

}
