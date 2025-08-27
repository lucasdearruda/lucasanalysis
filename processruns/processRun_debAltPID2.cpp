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

    const char * outputFileName  = Form("/home/e802/Analysis/LucasAnalysis/reducedData/altPID2/%03ddeb.root",nrun);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
    TTree *newtree = new TTree("M","Medley's reduced data TTree");

    Int_t particle; 
    
    Float_t energy; 
    
    Float_t si1; 
    Float_t si2; 
    Float_t csi; 

    Float_t tofn; 
    
    Float_t ang; 
    
    Float_t ENN; 

    TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
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
    TCutG *cut_protons[5];
    TCutG *cut_deuterons[5];
    TCutG *cut_tritons[5];

    TCutG *cut_he[5]; 
    TCutG *cut_alphas[5];
    
    string cuts_path = "/home/e802/Analysis/LucasAnalysis/PID_singleBananaCuts";

    for(int telN = 1; telN<=4; telN++){
    
    //     //protons
        cout << Form("\n- - - - - - - - - LOADING CUTS  - - - - - - - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading PROTONS CUT: %s/protonsT%d.C", telN,cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L %s/protonsT%d.C", cuts_path.c_str(),telN));
        cut_protons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("protonsT%d",telN));
        if(cut_protons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_protons[telN]->SetLineWidth(2);
        cut_protons[telN]->SetLineColor(kRed);
    
        //deuterons
        cout << Form("\n->[tel %d] loading DEUTERONS CUT: %s/deuteronsT%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/deuteronsT%d.C", cuts_path.c_str(),telN));
        cut_deuterons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("deuteronsT%d",telN));
        if(cut_deuterons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_deuterons[telN]->SetLineWidth(2);
        cut_deuterons[telN]->SetLineColor(kRed);
    
        //tritons
        cout << Form("\n->[tel %d] loading TRITONS CUT: %s/tritonsT%d.C", telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/tritonsT%d.C", cuts_path.c_str(),telN));
        cut_tritons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("tritonsT%d",telN));
        if(cut_tritons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_tritons[telN]->SetLineWidth(2);
        cut_tritons[telN]->SetLineColor(kRed);
    
        //helium
        cout << Form("\n->[tel %d] loading HELIUM CUT: %s/heliumT%d.C", telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/heliumT%d.C", cuts_path.c_str(),telN));
        cut_he[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("heliumT%d",telN));
        if(cut_he[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_he[telN]->SetLineWidth(2);
        cut_he[telN]->SetLineColor(kRed);
    
        //alphas
        cout << Form("\n->[tel %d] loading ALPHAS CUT: %s/alphasT%d.C",telN,cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/alphasT%d.C", cuts_path.c_str(),telN));
        cut_alphas[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("alphasT%d",telN));
        if(cut_alphas[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_alphas[telN]->SetLineWidth(2);
        cut_alphas[telN]->SetLineColor(kRed);
    
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

                if( cut_protons[telN]->IsInside(de2+csi,de1) ){
                    //it is a proton!
                    particle = 1;

                }else if(cut_deuterons[telN]->IsInside(de2+csi,de1)){
                    //it is a deuteron!
                    particle = 2;
                     

                }else if(cut_tritons[telN]->IsInside(de2+csi,de1)){
                    //it is a triton!
                    particle = 3;


                }else if(cut_he[telN]->IsInside(de2+csi,de1)){
                    //it is a helium!
                    particle = 4;
                  
                }else if(cut_alphas[telN]->IsInside(de2+csi,de1)){
                    //it is an alpha!
                    particle = 5;
      
                }
            }

            //if a particle was identified

            if(particle){
                
                ang = 20.*telN;
                //evaluate tof normally
                tofn = Gflash[telN] - invTof + tgamma[telN] - ToFparticleNumber(energy,particle,L[telN]);
                
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
        TNamed *processing_info= new TNamed("processed by processRun_debAltPID2",cur_time);
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

    

    cout<<"Compiled from processRun_debAltPID2.cpp"<<endl;
    return 0;

}
