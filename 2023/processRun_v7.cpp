//version 7.2025-08-29.0
// We updated the ifs for the angles; adding an error message also if the configuration is not defined in the run database

//version 7.2025-08-27.0
//I've changed the GetGflash function to include the runN and the param[] to it


//version 7.2025-08-21.2
//I updated the ToF correction for protons and deuterons
// removed the Dpoint and Ppoint variables, they were not used

//version 7.2025-08-20.2
//Diego spotted an error in the angles attribution, it was fixed today
//We updated some printing information regarding the progress

//version 7.2025-08-05.2
// Fixed the information from the .csv
//version 7.2025-08-05.1
// Fixed the TNamed for the processing info to include the correct file name and time.
//version 7.2025-08-05.0
// I added the correction for the TNamed 'info_tofNN' to include the tof_correction if it is enabled.

//version 7.2025-07-23.0
//I added the processing for Tels 5-8...

//version 7.2025-07-21.0
//The difference is that you have a 'continue' instead of 'break' in the loop that processes the telescopes; when you verify the invToF...

//version 6.2024-12-06.1
// The difference with 6.0 is that you write the event as soon as you find a particle (you skip looping over the other telesocpes)
//The difference with the 5.0 is in the use of GetRunData function, the adaptation of the angles according with the run (i.e., if its a 'R' run the T1-4 goes automatically to 100-160° range)
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
#include "TStyle.h"
#include "TROOT.h"

#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

// Inclua o functions.hh, mas o root-config cuida das libs do ROOT, então nada redundante
#include "/mnt/medley/LucasAnalysis/2023/include/functions2.hh"

int main(int argc, char* argv[]) {

    

    string cur_time = getCurrentTime();
    clock_t tStart = clock();
    
    bool tof_correction_bool = kTRUE;

    cout<<"# OF PROVIDED ARGUMENTS = "<<argc<<endl;
    for(int i=1; i<argc;i++)cout<<"argument #"<<i<<": "<<argv[i]<<endl;
    // Check if the correct number of arguments are provided
    if (argc < 2 || argc > 4) {
    	std::cerr << "Usage: " << argv[0] << " nrun [percentage] [tofcorrection: 0 for disable or 1 for enabled]" << std::endl;
    	return 1; // Return error code
	}

	// Assuming the first argument is the run number
	int nrun = std::stoi(argv[1]);
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    
    string infofile = "/mnt/medley/LucasAnalysis/2023/MEDLEY2023v7.csv";
    cout<<"Obtaining run data from "<<infofile<<"..."<<endl;

    //obtaining runData: 
    
    std::optional<MedleyData> result = GetRunData(infofile, nrun); //it was 'auto' before

    if(result){
        cout<<". . . runData info properly obtained. "<<endl;
        result->Print();
    }else{
        cout<<". . . runData info NOT properly obtained. "<<endl;
    } 
    // Check if the optional percentage argument is provided
    double percentage = 100.0; // Default value
    if (argc >= 3) {// --> First condition for tritons: t_npt_1 AND csi < 0.5 MeV)
//process the run for the reduced run. It is based on the version 2;
//in this version, we use the classic PID (ΔE1 . ΔE2 + ΔE2.Eres ) and also add info for the measured ToF as well as the correction needed.

        percentage = std::stod(argv[2]);
        
        if (percentage <= 0.0 || percentage > 100.0) {
            std::cerr << "Warning: Percentage must be between 0 and 100. Defaulting to 100." << std::endl;
            percentage = 100.0; // Set to default value
        }
        
    }
    cout<<"Processing "<<percentage<<"% of the run."<<endl;
    if (argc == 4) {
     
     tof_correction_bool = std::stoi(argv[3]);
     

    }
    if(tof_correction_bool) cout<<"TOF CORRECTION enabled!("<<tof_correction_bool<<")"<<endl;
    else cout<<"TOF CORRECTION disabled!("<<tof_correction_bool<<")"<<endl;

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    TGraphErrors *correctionP, *correctionD;
    double Ppoint, Dpoint;
    if(tof_correction_bool){

    //Create the TGraphs for correction:
        //TGraphErrors *correctionP, *correctionD;

        string correctionsPath = "/mnt/medley/LucasAnalysis/2023/Corrections";
        
        TFile *file = new TFile(Form("%s/newDcorr2023.root", correctionsPath.c_str()), "READ");
        TFile *fileP = new TFile(Form("%s/newPcorr2023.root", correctionsPath.c_str()), "READ");

        correctionD = (TGraphErrors *)file->Get("dcorrection");
        correctionP = (TGraphErrors *)fileP->Get("pcorrection");

        Ppoint = correctionP->GetPointX(correctionP->GetN()-1);
        Dpoint = correctionD->GetPointX(correctionD->GetN()-1);

        cout<<"Ppoint = "<<Ppoint<<" MeV. "<<endl;
        cout<<"Dpoint = "<<Dpoint<<" MeV. "<<endl;
        

    }
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //CREATE THE REDUCED TTREE

    const char * outputFileName  = Form("/mnt/medley/LucasAnalysis/2023/reducedv7/%03d.root",nrun);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
    TTree *newtree = new TTree("M","Medley's reduced data TTree");

    Int_t particle; 
    
    Float_t energy; 
    
    Float_t si1; 
    Float_t si2; 
    Float_t csi; 
    Float_t tof_measured;
    Float_t rtof;
    Float_t tof_correction; 
    Float_t tofn; 
    Float_t ang; 
    Float_t ENN; 

    TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
    
    TBranch * rawtof = newtree->Branch("rawtof", &rtof, "rawtof/F");
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
        sprintf(name, "/mnt/medley/RootA_2023/r%04d_%03da.root", nrun, j);
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

    int maxTel = 8; //max telescope number
    //cuts for npt particles
    TCutG *cut_protons[8], *cut_deuterons[8], *cut_tritons[8]; 
    TCutG *cut_he[8], *cut_alphas[8], *cut_he34[8], *cut_dt[8], *cutH[8];
    
    //cuts for pt particles
    TCutG *cut_protonsPT[8], *cut_deuteronsPT[8], *cut_tritonsPT[8];
    TCutG *cut_protonsCsI[8], *cut_deuteronsCsI[8], *cut_tritonsCsI[8]; 

    //cut npt csi 
    TCutG *cut_protonsNPTCsI[8], *cut_deuteronsNPTCsI[8];

    string cuts_path =   "/mnt/medley/LucasAnalysis/2023/PIDv6";
    bool warn_tels58 = false;
    for(int telN = 1; telN<=8; telN++){
        
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
                cut_tritonsPT[telN]->SetLineColor(kGreen);

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
            if(!warn_tels58){
                cout << "\n . . .  - - - - - - - - - - - - - -   . . . \n . . . - - - TELESCOPES 5 - 8  - - -  . . . "<<"\n . . .  - - - - - - - - - - - - - -   . . . "<<endl;
                warn_tels58 = true;
            }
            if(telN != 7 ){
                //cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
                cout << Form("\n->[tel %d] loading NPT PROTONS CUT: %s/p_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
                gROOT->ProcessLine(Form(".L  %s/p_npt_%d.C", cuts_path.c_str(),telN));
                cut_protons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_npt_%d",telN));
                if(cut_protons[telN]!=NULL) cout << "  --> OK."<< endl;
                cut_protons[telN]->SetName(Form("p_npt_%d",telN));
                cut_protons[telN]->SetLineWidth(2);
                cut_protons[telN]->SetLineColor(kRed);

                cout << Form("\n->[tel %d] loading NPT DEUTERON-TRITONS CUT: %s/dt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
                gROOT->ProcessLine(Form(".L  %s/dt_%d.C", cuts_path.c_str(),telN));
                cut_dt[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("dt_%d",telN));
                if(cut_dt[telN]!=NULL) cout << "  --> OK."<< endl;
                cut_dt[telN]->SetName(Form("dt_%d",telN));
                cut_dt[telN]->SetLineWidth(2);
                cut_dt[telN]->SetLineColor(kBlue);

            }else{
                cout << Form("\n->[tel %d] loading H(1, 2 or 3) CUT: %s/h_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
                gROOT->ProcessLine(Form(".L  %s/h_%d.C", cuts_path.c_str(),telN));
                cutH[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("h_%d",telN));
                if(cutH[telN]!=NULL) cout << "  --> OK."<< endl;
                cutH[telN]->SetName(Form("h_%d",telN));
                cutH[telN]->SetLineWidth(2);
                cutH[telN]->SetLineColor(kRed);
            }
            cout << Form("\n->[tel %d] loading He(3 or 4) CUTS: %s/he_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
            gROOT->ProcessLine(Form(".L  %s/he_%d.C", cuts_path.c_str(),telN));
            cut_he34[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("he_%d",telN));
            if(cut_he34[telN]!=NULL) cout << "  --> OK."<< endl;
            cut_he34[telN]->SetName(Form("he_%d",telN));
            cut_he34[telN]->SetLineWidth(2);
            cut_he34[telN]->SetLineColor(8);

        }

    }


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    ///events processing

    Float_t Gflash[maxTel];
    Float_t tgamma[maxTel];
    float param[maxTel][4];
    for(int i=0;i<5;i++) cout<<"/ : /"<<i<<" /"<<endl; // just to spot it on the terminal 
    for(int n=1;n<=8;n++){
        //Fitting Gamma flash 



// Float_t GetGflash(TTree *tx, 
//                     Int_t telN = 1, 
//                     Int_t runN=0,
//                     Float_t sigma= 4,
//                     Float_t threshold= 0.002,
//                     float guess = 400,
//                     float param[] = nullptr,
//                     float tofmin = 300,
//                     float tofmax = 500, 
//                     bool closecanvas = false, 
//                     bool saveIt=true){

        Gflash[n] = GetGflash(tx,n,nrun,4,0.002,500,param[n-1],300,550); // changing it to false, it will let the fitting plot opened to check 
        tgamma[n] = provideTgama(n);

        cout<<".\n.\n.\nFitting Gamma flash for telescope "<<n<<": Gflash = "<<Gflash[n]<<" ns. tgamma = "<<tgamma[n]<<" ns."<<endl;

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

        for(int telN = 1; telN <= 8; telN++){ // FOR EACH TELESCOPE:
            
            //critical verifications: 
            invTof = tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue();
            rtof = invTof;

            if(invTof > Gflash[telN]) continue; //if the event is faster than gamma... go to next telecope:: THIS IS RIGHT NOW (before it was wrong, with break) 2025-07-21

            de1 =  tx->GetLeaf(Form("Medley_%d_dE1",telN))->GetValue();
            de2 =  tx->GetLeaf(Form("Medley_%d_dE2",telN))->GetValue();
            Eres =  tx->GetLeaf(Form("Medley_%d_Eres",telN))->GetValue();
            energy = de1+de2+Eres;

            si1 = de1;
            si2 = de2;
            csi = Eres;
            if( de1 >= 0.02 && de2 >= 0.02 ){ // if it is not trash -- energy threshold

                if(telN<=4){                    
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
                        cut_tritons[telN]->IsInside(de2,de1) && csi<0.5
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
                }else{//particle classification for telescopes 5-8::
                    
                    if(
                        cut_he34[telN]->IsInside(de2,de1)
                    ){
                        particle = 7; // He-3 or alpha
                    // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ its an He-(3 or 4); tels 5 6 7 8 
                    }else{
                        if(telN != 7){ // telescopes 5, 6  or 8 
                            if(cut_protons[telN]->IsInside(de2,de1)){
                                //it is a proton!
                                particle = 1;
                            }else if(cut_dt[telN]->IsInside(de2,de1)){
                                //it is a deuteron or triton!
                                particle = 6;
                            }
                            // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ can be a proton or a D/T  for tels 5 6 8
                        }else{ //telescope 7:
                            if(cutH[telN]->IsInside(de2,de1)){
                                //it is a proton, deuteron or triton!
                                particle = 6;
                            }
                            // _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ can be a H(1, 2 or 3)  for tels 7
                        }

                    }

                }
            }

            //if a particle was identified

            if(particle){
                
                if( (result->MedleyConfig) == 'R'){
                    ang = 180-20.*telN; // [1, 2, 3, 4] --> [160, 140, 120, 100] //IMPORTANT CORRECTION 2025-08-20 (spotted by Diego)    
                                        // [5, 6, 7, 8] --> [80, 60, 40, 20]
                }else if( (result->MedleyConfig) == 'N'){
                    ang = 20.*telN; // [1, 2, 3, 4] --> [20, 40, 60, 80] and [5, 6, 7, 8] --> [100, 120, 140, 160]
                }else{
                    ang = -1; //undefined configuration
                    cerr << "\n\nWARNING: Medley configuration is not defined in the run database!! Check it out!!\n\n"<<endl;                
                }
                //evaluate tof normally
                tof_measured = Gflash[telN] - invTof + tgamma[telN];
                tofn = tof_measured - ToFparticleNumber(energy,particle,L[telN]);
                tof_correction = 0;
                //if needed you update for the corrected value 
                if((particle == 1) && (energy <=Ppoint) && (telN <=4)){// if it is proton and is inside the correction region
                    //cout<<"\nproton. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionP->Eval(energy)<<" ns."<<endl;
                    if(tof_correction_bool ) tof_correction = correctionP->Eval(energy);
                    tofn = tofn + tof_correction;
                    //cout<<"proton. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }else if((particle == 2)  && (energy <=Dpoint)&& (telN <=4)){// if it is deuteron and is inside the correction region
                    //cout<<"\ndeuteron. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionD->Eval(energy)<<" ns."<<endl;
                    if(tof_correction_bool )  tof_correction =correctionD->Eval(energy);
                    tofn = tofn + tof_correction;
                    //cout<<"deuteron. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }

                ENN = En(tofn);

                newtree->Fill();
                break; // Go to the next entry

            }


        }//finished looping over telescopes



        // Print the progress
        if(entry % 10000 == 0){
            std::cout << "Progress ["<< nrun <<"]: " << 100.0*entry/(nEntries*percentage/100.0) << "%\r";
            std::cout.flush();
        }
        
        
    }



    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


    //some information to add
        TNamed *percentage_proc= new TNamed("processed (percentage): ",to_string(percentage));
        TNamed *processing_info= new TNamed("processed by processRun_v7.cpp:",Form("version 7.2025-08-29.0 in %s",cur_time.c_str()));
        TNamed *charge_in_uC = NULL;
        TNamed *MedleyTarget= NULL;
        TNamed *duration_of_the_run_in_sec= NULL;
        TNamed *configuration_telescopes= NULL;
        TNamed *gamma_flash_in_ns_T1 = new TNamed("gFlash_T1(ns)",to_string(Gflash[1]).c_str());
        TNamed *gamma_flash_in_ns_T2 = new TNamed("gFlash_T2(ns)",to_string(Gflash[2]).c_str());
        TNamed *gamma_flash_in_ns_T3 = new TNamed("gFlash_T3(ns)",to_string(Gflash[3]).c_str());
        TNamed *gamma_flash_in_ns_T4 = new TNamed("gFlash_T4(ns)",to_string(Gflash[4]).c_str());
        TNamed *gamma_flash_in_ns_T5 = new TNamed("gFlash_T5(ns)",to_string(Gflash[5]).c_str());
        TNamed *gamma_flash_in_ns_T6 = new TNamed("gFlash_T6(ns)",to_string(Gflash[6]).c_str());
        TNamed *gamma_flash_in_ns_T7 = new TNamed("gFlash_T7(ns)",to_string(Gflash[7]).c_str());
        TNamed *gamma_flash_in_ns_T8 = new TNamed("gFlash_T8(ns)",to_string(Gflash[8]).c_str());

        TNamed *dist_Medleytarget_T1_in_mm = new TNamed("dist_Medleytarget_T1(mm)",to_string(1e3*L[1]).c_str());
        TNamed *dist_Medleytarget_T2_in_mm= new TNamed("dist_Medleytarget_T2(mm)",to_string(1e3*L[2]).c_str());
        TNamed *dist_Medleytarget_T3_in_mm= new TNamed("dist_Medleytarget_T3(mm)",to_string(1e3*L[3]).c_str());
        TNamed *dist_Medleytarget_T4_in_mm= new TNamed("dist_Medleytarget_T4(mm)",to_string(1e3*L[4]).c_str());
        TNamed *dist_Medleytarget_T5_in_mm= new TNamed("dist_Medleytarget_T5(mm)",to_string(1e3*L[5]).c_str());
        TNamed *dist_Medleytarget_T6_in_mm= new TNamed("dist_Medleytarget_T6(mm)",to_string(1e3*L[6]).c_str());
        TNamed *dist_Medleytarget_T7_in_mm= new TNamed("dist_Medleytarget_T7(mm)",to_string(1e3*L[7]).c_str());
        TNamed *dist_Medleytarget_T8_in_mm= new TNamed("dist_Medleytarget_T8(mm)",to_string(1e3*L[8]).c_str());


    
        TNamed *info_rawtof= new TNamed("raw_tof","raw tof from raw data (inv tof without gamma flash)");
        TNamed *info_tof_measured= new TNamed("tof_measured","Gflash - invtof + tgamma == tofNN+ tof_part");
        TNamed *info_tof_correction= new TNamed("tof_correction","correction for proton and deuteron particles, if needed");
        TNamed *info_tofNN= new TNamed("tofNN","tofNN = tof_measured - ToFparticleNumber(energy,particle,L[telN])");
        if(tof_correction_bool)info_tofNN = new TNamed("tofNN","tofNN = tof_measured - ToFparticleNumber(energy,particle,L[telN]) + tof_correction");
        
        TNamed *info_ENN= new TNamed("ENN","ENN = En(tofNN) -- energy in MeV from ToFNN");
        TNamed *info_si1= new TNamed("si1","energy deposited in si1 (dE1)");
        TNamed *info_si2= new TNamed("si2","energy deposited in si2 (dE2)");
        TNamed *info_csi= new TNamed("csi","energy deposited in CsI (Eres)");

        if (result) {
            charge_in_uC = new TNamed("RunCharge(C)",Form("%s",to_string((result->RunCharge)*percentage/100.).c_str()));
            configuration_telescopes = new TNamed("MedleyConfig",Form("%s",(result->MedleyConfig).c_str()));
            MedleyTarget = new TNamed("MedleyTarget",Form("%s",(result->Target).c_str()));
            duration_of_the_run_in_sec = new TNamed("RunTime(s)",to_string((result->RunTime)*percentage/100.).c_str());
        } else {
            charge_in_uC = new TNamed("RunCharge(C)","not_available");
            configuration_telescopes = new TNamed("MedleyConfig","not_available");
            MedleyTarget = new TNamed("MedleyTarget","not_available");
            duration_of_the_run_in_sec = new TNamed("RunTime(s)","not_available");
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
    gamma_flash_in_ns_T5->Write();
    gamma_flash_in_ns_T6->Write();
    gamma_flash_in_ns_T7->Write();
    gamma_flash_in_ns_T8->Write();


    dist_Medleytarget_T1_in_mm->Write();
    dist_Medleytarget_T2_in_mm->Write();
    dist_Medleytarget_T3_in_mm->Write();
    dist_Medleytarget_T4_in_mm->Write();
    dist_Medleytarget_T5_in_mm->Write();
    dist_Medleytarget_T6_in_mm->Write();
    dist_Medleytarget_T7_in_mm->Write();
    dist_Medleytarget_T8_in_mm->Write();

    info_rawtof->Write();
    info_tof_measured->Write();
    info_tof_correction->Write();   
    info_tofNN->Write();
    info_ENN->Write();
    info_si1->Write();
    info_si2->Write();
    info_csi->Write();

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

    cout<<"Compiled from processRun_v7.cpp, version 7.2025-08-29.0"<<endl;
    return 0;

}
