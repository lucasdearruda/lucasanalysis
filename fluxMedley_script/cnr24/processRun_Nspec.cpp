
//
//
//
// ROUTINE TO PROCESS RUNS SPECIFICALLY TO OBTAIN THE NEUTRON FLUX 
//   use single banana as selection 
//   runs only for telescope 1 
//   keep wrong ToF and corrected one
//
//
//

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
#include "TLeaf.h" // Make sure to include the necessary headers
#include "TChain.h"

#include "TGraph.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStopwatch.h"

#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

//take a look at the functions:
#include "functions.h"


int main(int argc, char* argv[]) {

        // Check if the correct number of arguments are provided
    if (argc < 2 || argc > 3) {
        std::cerr << "Script to process runs specially for neutron flux reconstruction for 2023." << std::endl;
        std::cerr << "Important notes:" << std::endl;
    	std::cerr << " --> Use single bananas." << std::endl;
        std::cerr << " --> Keep good and bad ToF (and Enn) stored." << std::endl;
        std::cerr << " --> Process only telescope 1!" << std::endl<< std::endl;
        std::cerr << "Usage: " << argv[0] << " nrun [percentage]" << std::endl;
    	return 1; // Return error code
	}

	// Assuming the first argument is the run number
	int nrun = std::stoi(argv[1]);
    
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

    //ok up to here
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


 

//int processRun_Nspec(int nrun = 407){
//Open TOF corrections

    TStopwatch timer;
    //Create the TGraphs for correction:
    TGraphErrors *correctionP, *correctionD;
    string correctionsPath = "/home/e802/Analysis/pre_analysis/runs/Corrections";


    TFile *file = new TFile(Form("%s/dcorrection.root", correctionsPath.c_str()), "READ");

    cout<< "reading correction from "<<Form("dcorrection.root")<<"..."<<endl;


    correctionD = (TGraphErrors *)file->Get("dcorrection");

    cout<< "creating canvases..."<<endl;
    TCanvas *Crr, *CC;

    Crr = new TCanvas("Crr","Loaded corrections",100,100,1600,800);
    Crr->Divide(2,1);

    Crr->cd(1);
    Crr->cd(1)->SetRightMargin(0.12);
    Crr->cd(1)->SetLeftMargin(0.16);
    cout<< "plotting correction..."<<endl;

    correctionD->SetTitle("tof correction");
    correctionD->SetMarkerColor(kBlue);
    correctionD->SetLineColor(kBlue);
    correctionD->SetLineWidth(3);
    correctionD->SetMarkerStyle(21);
    correctionD->Draw("APL");
    gPad->SetGridx();
    gPad->SetGridy();

    Crr->cd(2);
    Crr->cd(2)->SetRightMargin(0.12);
    Crr->cd(2)->SetLeftMargin(0.16);


    TFile *fileP = new TFile(Form("%s/pcorrection.root", correctionsPath.c_str()), "READ");

    cout<< "reading correction from "<<Form("pcorrection.root")<<"..."<<endl;

    correctionP = (TGraphErrors *)fileP->Get("pcorrection");
    cout<< "plotting correction..."<<endl;

    correctionP->SetTitle("tof correction P");
    correctionP->SetMarkerColor(kRed);
    correctionP->SetLineColor(kRed);
    correctionP->SetLineWidth(3);
    correctionP->SetMarkerStyle(26);
    correctionP->Draw("APL");
    gPad->SetGridx();
    gPad->SetGridy();

    double Ppoint = correctionP->GetPointX(correctionP->GetN()-1);
    double Dpoint = correctionD->GetPointX(correctionD->GetN()-1);

    cout<<"maximum energy Ppoint = "<<Ppoint<<" MeV. "<<endl;
    cout<<"maximum energy Dpoint = "<<Dpoint<<" MeV. "<<endl;
    
    //ok up to here

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //CREATE THE REDUCED TTREE

    const char * outputFileName  = Form("run_%03d_specific.root",nrun);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    TTree *newtree = new TTree("M","Medley's reduced data TTree");

    Int_t particle; //particle type (p, d, t, h, a) = (1,2,3,4,5)
    Float_t energy; //total energy in MeV
    Float_t si1; //energy of  si1 in MeV
    Float_t si2;  //energy of  si2 in MeV
    Float_t csi;  //energy of csi in MeV
    Float_t tofn; //tof for neutron
    Float_t ENN; //energy for neutron
    Float_t tofn_unc; //uncorrected tof for neutron
    Float_t ENN_unc; //uncorrecter energy for neutron

    TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
    TBranch * branch_si1  = newtree->Branch("si1", &si1, "si1/F");
    TBranch * branch_si2  = newtree->Branch("si2", &si2, "si2/F");
    TBranch * branch_csi  = newtree->Branch("csi", &csi, "csi/F");
    TBranch * tofBranch  = newtree->Branch("tofn", &tofn, "tofn/F");
    TBranch * EnBranch  = newtree->Branch("ENN", &ENN, "ENN/F");

    TBranch * tofBranch_unc  = newtree->Branch("tofn_unc", &tofn_unc, "tofn_unc/F");
    TBranch * EnBranch_unc  = newtree->Branch("ENN_unc", &ENN_unc, "ENN_unc/F");
    
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
    //OPENING CUTS FOR PID:
    TCutG *cut_protons;
    TCutG *cut_deuterons;
    TCutG *cut_tritons;

    TCutG *cut_he; 
    TCutG *cut_alphas;
    
    string cuts_path = "/home/e802/Analysis/LucasAnalysis/PID_singleBananaCuts";
    
    Int_t telN = 1;
    //loading cuts

    //     //protons
    cout << Form("\n- - - - - - - - - LOADING CUTS  - - - - - - - - -\n")<< "... " ;
    cout << Form("\n->[tel %d] loading PROTONS CUT: %s/protonsT%d.C", telN,cuts_path.c_str(),telN)<< "... " ;
    gROOT->ProcessLine(Form(".L %s/protonsT%d.C", cuts_path.c_str(),telN));
    cut_protons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("protonsT%d",telN));
    if(cut_protons!=NULL) cout << "  --> OK."<< endl;
    cut_protons->SetLineWidth(2);
    cut_protons->SetLineColor(kRed);

    //deuterons
    cout << Form("\n->[tel %d] loading DEUTERONS CUT: %s/deuteronsT%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
    gROOT->ProcessLine(Form(".L  %s/deuteronsT%d.C", cuts_path.c_str(),telN));
    cut_deuterons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("deuteronsT%d",telN));
    if(cut_deuterons!=NULL) cout << "  --> OK."<< endl;
    cut_deuterons->SetLineWidth(2);
    cut_deuterons->SetLineColor(kRed);

    //tritons
    cout << Form("\n->[tel %d] loading TRITONS CUT: %s/tritonsT%d.C", telN, cuts_path.c_str(),telN)<< "... " ;
    gROOT->ProcessLine(Form(".L  %s/tritonsT%d.C", cuts_path.c_str(),telN));
    cut_tritons = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("tritonsT%d",telN));
    if(cut_tritons!=NULL) cout << "  --> OK."<< endl;
    cut_tritons->SetLineWidth(2);
    cut_tritons->SetLineColor(kRed);

    //helium
    cout << Form("\n->[tel %d] loading HELIUM CUT: %s/heliumT%d.C", telN, cuts_path.c_str(),telN)<< "... " ;
    gROOT->ProcessLine(Form(".L  %s/heliumT%d.C", cuts_path.c_str(),telN));
    cut_he = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("heliumT%d",telN));
    if(cut_he!=NULL) cout << "  --> OK."<< endl;
    cut_he->SetLineWidth(2);
    cut_he->SetLineColor(kRed);

    //alphas
    cout << Form("\n->[tel %d] loading ALPHAS CUT: %s/alphasT%d.C",telN,cuts_path.c_str(),telN)<< "... " ;
    gROOT->ProcessLine(Form(".L  %s/alphasT%d.C", cuts_path.c_str(),telN));
    cut_alphas = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("alphasT%d",telN));
    if(cut_alphas!=NULL) cout << "  --> OK."<< endl;
    cut_alphas->SetLineWidth(2);
    cut_alphas->SetLineColor(kRed);

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//testing cuts 
    TCanvas *testcanvas = new TCanvas();
    tx->Draw("Medley_1_dE1:Medley_1_dE2+Medley_1_Eres>>h(2048,0,40,2048,0,10)","","colz");
    cut_alphas->Draw("same");
    cut_he->Draw("same");
    cut_tritons->Draw("same");
    cut_deuterons->Draw("same");
    cut_protons->Draw("same");
    testcanvas->Write();
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    Float_t Gflash;
    Float_t tgamma;
    Float_t L = 0.1437;//meters

    Gflash = GetGflash(tx,"Medley_1_dE2_ToF",500,true); // changing it to false, it will let the fitting plot opened to check 
    tgamma = provideTgama(1);
    cout<<"Gflash = "<<Gflash<<" ns."<<endl;
    cout<<"tgamma = "<<tgamma<<" ns."<<endl;
    


    // 1: proton
    // 2: deuteron
    // 3: triton
    // 4: helium-3
    // 5: alpha
    // 6: neutron -- not used of course (but the map follows it in the function)
    // 7: unknown 


    Long64_t nEntries = tx->GetEntries();
    cout<<" Total entries : "<<nEntries<<endl; 

    Float_t invTof;
    //Float_t percentage = 1;
    for (Long64_t entry = 0; entry < nEntries*percentage/100.0; entry++) {

        tx->GetEntry(entry);
        
        //initializing variables
        particle = 0; //particle starts unidentified
        energy =  0;
        si1 = 0;
        si2 = 0;
        csi = 0;
        tofn = 0; //tof starts zero too 
        ENN = 0;
        tofn_unc = 0;
        ENN_unc = 0;

        invTof = tx->GetLeaf("Medley_1_dE2_ToF")->GetValue();

        if(invTof > Gflash) continue; //if the event is faster than gamma... forget it => next event

        si1 =  tx->GetLeaf("Medley_1_dE1")->GetValue();
        si2 =  tx->GetLeaf("Medley_1_dE2")->GetValue();
        csi =  tx->GetLeaf("Medley_1_Eres")->GetValue();
        energy = si1+si2+csi;
    
        tofn_unc = Gflash - invTof + tgamma; // some default value for unknown particles 
        
        if( si1 >= 0.02 && si2 >= 0.02 ){ // if it is not trash -- energy threshold
            
            if( cut_protons->IsInside(si2+csi,si1) ){
                //it is a proton!
                particle = 1;
                

            }else if(cut_deuterons->IsInside(si2+csi,si1)){
                //it is a deuteron!
                particle = 2;
                    

            }else if(cut_tritons->IsInside(si2+csi,si1)){
                //it is a triton!
                particle = 3;                if(tofn_unc != tofn) 
                    ENN_unc = En(tofn_unc);
                else    
                    ENN_unc = ENN;

            }else if(cut_he->IsInside(si2+csi,si1)){
                //it is a helium!
                particle = 4;
        
            }else if(cut_alphas->IsInside(si2+csi,si1)){
                //it is an alpha!
                particle = 5;
            }

            //if particle was identified
            if(particle){
                
                tofn = Gflash - invTof + tgamma - ToFparticleNumber(energy,particle,L);
                
                tofn_unc= tofn; //if there is no correction to be applied: tofn == tofn_unc

                //if needed you update for the corrected value 
                if((particle == 1) && (energy <=Ppoint)){// if it is proton and is inside the correction region
                    //cout<<"\nproton. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionP->Eval(energy)<<" ns."<<endl;
                    tofn = tofn + correctionP->Eval(energy);
                    //cout<<"proton. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }else if((particle == 2)  && (energy <=Dpoint)){// if it is deuteron and is inside the correction region
                    //cout<<"\ndeuteron. energy = "<<energy<<" MeV. ToF = "<<tof<<" ns. Correction"<< correctionD->Eval(energy)<<" ns."<<endl;
                    tofn = tofn + correctionD->Eval(energy);
                    //cout<<"deuteron. energy = "<<energy<<" MeV. NEW ToF = "<<tof<<" ns."<<endl;
                }

                ENN = En(tofn);
                
                if(tofn_unc != tofn) 
                    ENN_unc = En(tofn_unc);
                else    
                    ENN_unc = ENN;

            }
            
        
        }
        newtree->Fill();
        //if its out of energy ranges on si1 and Si2 and also out of the cuts, everything but energies and tofn_unc will be zero

        std::cout << "Progress ["<< nrun <<"]: " << 100.0*entry/(nEntries*percentage/100.0) << "%\r";
        std::cout.flush();
    }//loop over entries 


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    string infofile = "/home/e802/Analysis/LucasAnalysis/reducedData/MEDLEY2023.csv";

    //some information to add
    TNamed *charge_in_uC = NULL;
    TNamed *MedleyTarget= NULL;
    TNamed *duration_of_the_run_in_sec= NULL;
    TNamed *configuration_telescopes= NULL;
    TNamed *Processed_with = new TNamed("processed_with","processRun_Nspec.cpp");
    TNamed *gamma_flash_in_ns_T1 = new TNamed("gFlash_T1(ns)",to_string(Gflash).c_str());

    TNamed *dist_Medleytarget_T1_in_mm = new TNamed("dist_Medleytarget_T1(mm)",to_string(1e3*L).c_str());

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
                    charge_in_uC = new TNamed("RunCharge(C)",Form("%s",to_string(RunCharge).c_str()));
                }else{
                    charge_in_uC = new TNamed("RunCharge(C)","not_available");
                }
                configuration_telescopes = new TNamed("MedleyConfig",Form("%s",MedleyConfig));
                MedleyTarget = new TNamed("MedleyTarget",Form("%s",Target));
                if(RunTime>0){
                    duration_of_the_run_in_sec = new TNamed("RunTime(s)",to_string(RunTime).c_str());
                }else{
                    duration_of_the_run_in_sec = new TNamed("RunTime(s)","not_available");
                }
                
            }
        }


    charge_in_uC->Write();
    MedleyTarget->Write();
    duration_of_the_run_in_sec->Write();
    configuration_telescopes->Write();
    Processed_with->Write();
    gamma_flash_in_ns_T1->Write();


    dist_Medleytarget_T1_in_mm->Write();
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    outputFile->Write();
    // Close the output file
    outputFile->Close();
    timer.Print();
    return 0;
}


/*
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    outputFile->Write();
    // Close the output file
    outputFile->Close();
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    timer.Print();
    //cout<<"Compiled from processRun_Nspec.cpp"<<endl;

    return 0;
}

*/