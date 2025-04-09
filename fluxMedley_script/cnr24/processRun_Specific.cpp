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

//
//
//
//
// ROUTINE TO PROCESS RUNS SPECIFICALLY TO OBTAIN THE NEUTRON FLUX 
//   use single banana as selection 
//   runs only for telescope 1 
//   keep wrong ToF and corrected one 

using namespace std;

//PT energies from kaliveda files. ptinfo[1][1] = 12.6047, is the PT energy for protons first tel 
float ptinfo[6][9] = {
    // Linha em branco inicial
    {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    // p = 1
    {0.0, 12.6047, 12.6605, 12.6446, 12.6845, 8.66254, 8.80618, 8.70244, 8.81416},
    // d = 2
    {0.0, 16.9538, 17.0256, 17.0016, 17.0575, 11.5673, 11.7668, 11.6311, 11.7747},
    // t = 3
    {0.0, 20.1059, 20.1857, 20.1617, 20.2256, 13.65, 13.8815, 13.7219, 13.8974},
    // h = 4
    {0.0, 0.0, 0.0, 0.0, 0.0, 30.6794, 31.1821, 30.839, 31.214},
    // a = 5
    {0.0, 0.0, 0.0, 0.0, 0.0, 34.5975, 35.1641, 34.7731, 35.204}
};

double En(double tof, double L = 4.6472){//L in meters, 2023 reference
	//return neutron energy in MeV, given its time of fligt.
    //Special cases: returns ZERO if ToF<0 or gamma == NaN.
    //
    //

    //if tof is zero return 0
	if(tof<=0){ return 0;} 
	double v = L/tof;
	double c = 0.299792458; //m/ns
	double gamma = 1/sqrt(1-(v/c)*(v/c));
	if(gamma!=gamma) return 0; //case where gamma == nan;
	double mc2 = 939.56542052; //MeV
	return (gamma -1)*mc2; //MeV
}



double ToFparticleNumber(double Eparticle, int particle, double L = 4.6472)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'
// 1: proton
// 2: deuteron
// 3: triton
// 4: helium-3
// 5: alpha
// 6: neutron -- not used of course (but the map follows it in the function)


	double c= 0.299792458;
	double m;

	switch(particle){//data from codata
        	case 6://neutrons
        		m = 939.56542052;
        		break;
		case 1://protons
			m = 938.27208816;
			break;
		case 2://deuterons
			m = 1875.61294257;
			break;
		case 3://tritons 
			m = 2808.92113298;
			break;
		case 4: //helium3
			m =  2808.39160743;
			break;
		case 5://alphas
			m = 3727.3794066;
			break;
		default:
			m = 938.27208816;//proton
	}

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


//function for getting gamma flash position in ns
Float_t GetGflash(TTree *tx, string branchname="Medley_1_dE2_ToF", float guess = 500, bool closecanvas = false){


	TCanvas *Ctof = new TCanvas("ctof",Form("Time of flight"),50,50,600,600);

    TH1D *htof;
    htof = new TH1D("htof","htof",500,100,600);

    tx->Draw(Form("%s>>htof",branchname.c_str()));

    double tof_peak[2];
    TSpectrum *spec = new TSpectrum();
    Int_t npeaks = spec->Search(htof, 6,"",0.1);//originally 0.3
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


	if(closecanvas) Ctof->Close();
	
	return tof_peak[0];
}



double provideTgama(int telN, int year = 2023){//just to save some lines in my script 

	double tgamma;

	if(year == 2023){
		switch(telN){
			case 1:
				tgamma = 15.742;
				break;
			case 2:
				tgamma = 15.700;
				break;
			case 3:
				tgamma = 15.635;
				break;
			case 4:
				tgamma = 15.554;
				break;
			case 5:
				tgamma = 15.465;
				break;
			case 6:
				tgamma = 15.381;
				break;
			case 7:
				tgamma = 15.311;
				break;
			case 8:
				tgamma = 15.264;
				break;
			default:
				tgamma = 15.742;
		}
	}else{
		tgamma = 15;
		cout<<"****************************************************************"<<endl;
		cout<<"	year "<<year<<" nor available! Setting tgamma = 15 ns."<<endl;
		cout<<"****************************************************************"<<endl;
	}
	return tgamma;
}


int main(int argc, char* argv[]) {

    TStopwatch timer;

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
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //Open kaliveda simulations

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

/*
    TGraph *gt[5][4];
    Float_t Emin[5][4];

    for(int tel = 1; tel <= 4; tel++) {
        for(int particle = 1; particle <= 3; particle++) {
            TString particleName;
            if (particle == 1)
                particleName = "p";
            else if (particle == 2)
                particleName = "d";
            else if (particle == 3)
                particleName = "t";
            
            TString fileName = Form("/home/e802/Analysis/LucasAnalysis/kaliveda_Edep_results/E0si_plots/kal_%s_Tel%d_plot.root", particleName.Data(), tel);

            TFile *fileplot = new TFile(fileName, "READ");
            gt[tel][particle] = (TGraph*)fileplot->Get("Graph");

            // Clone the TGraph object to prevent it from being deleted when closing the file
            gt[tel][particle] = (TGraph*)gt[tel][particle]->Clone();
            // new TCanvas();
            // gt[tel][particle]->Draw();
            Emin[tel][particle] = gt[tel][particle]->Eval(TMath::MaxElement(gt[tel][particle]->GetN(),gt[tel][particle]->GetX()));
            //cout<<"Tel = "<<tel<<" particle = "<<particle<<" emin = "<<Emin[tel][particle]<<endl;
            fileplot->Close();
        }
    }
*/

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //Open TOF corrections

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
    ///events processing

    Float_t Gflash;
    Float_t tgamma;

    
    Gflash = GetGflash(tx,Form("Medley_%d_dE2_ToF",telN),500,true); // changing it to false, it will let the fitting plot opened to check 
    tgamma = provideTgama(telN);



    Float_t L;
    //distances Medley target --> Telescope
    L = 157.6*1e-3;
    
    // 1: proton
    // 2: deuteron
    // 3: triton
    // 4: helium-3
    // 5: alpha
    // 6: neutron -- not used of course (but the map follows it in the function)
    // 7: unknown 

    Float_t invTof;
    Float_t de2,de1,Eres;
    de2 = 0;
    de1 = 0;
    Eres = 0;

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //Loop over TTree

    //uncertainty stuff
    /*
        TRandom2 *rand = new TRandom2();
        double sigmaSi1 = 0.08;//0.5
        double sigmaSi2 = 0.1;//.07
        double sigmacsi = 0.2;//.07
        double si1_plus_si2;//.07


        double combined_sigma = sqrt(pow(sigmaSi1,2)+pow(sigmaSi2,2));
        cout<<"combined sigma= "<<combined_sigma<<endl;
    */

    Long64_t nEntries = tx->GetEntries();

    cout<<" Total entries : "<<nEntries<<endl; 

    // Loop through each entry again to calculate and set the new branch values
    for (Long64_t entry = 0; entry < nEntries*percentage/100.0; entry++) {
        //get entry 
        tx->GetEntry(entry);
        particle = 0; //particle starts unidentified
        tofn = 0; //tof starts zero too 


        //critical verifications: 
        invTof = tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue();
            
        if(invTof > Gflash) break; //if the event is faster than gamma... forget it => next event

        de1 =  tx->GetLeaf(Form("Medley_%d_dE1",telN))->GetValue();
        de2 =  tx->GetLeaf(Form("Medley_%d_dE2",telN))->GetValue();
        Eres =  tx->GetLeaf(Form("Medley_%d_Eres",telN))->GetValue();
        energy = de1+de2+Eres;


        si1 = de1;
        si2 = de2;
        csi = Eres;
    
        if( de1 >= 0.02 && de2 >= 0.02 ){ // if it is not trash -- energy threshold

            
            if( cut_protons->IsInside(de2+csi,de1) ){
                //it is a proton!
                particle = 1;
                

            }else if(cut_deuterons->IsInside(de2+csi,de1)){
                //it is a deuteron!
                particle = 2;
                    

            }else if(cut_tritons->IsInside(de2+csi,de1)){
                //it is a triton!
                particle = 3;

            }else if(cut_he->IsInside(de2+csi,de1)){
                //it is a helium!
                particle = 4;
        
            }else if(cut_alphas->IsInside(de2+csi,de1)){
                //it is an alpha!
                particle = 5;
            }

            //if a particle was identified

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

                //si1_plus_si2 = rand->Gaus(de1,sigmaSi1) + rand->Gaus(de2,sigmaSi2);
                //si1_plus_si2 = de1 + de2;
                newtree->Fill();

            }
        }

        // Print the progress
        std::cout << "Progress ["<< nrun <<"]: " << 100.0*entry/(nEntries*percentage/100.0) << "%\r";
        std::cout.flush();
        
    }

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


    //some information to add
        TNamed *charge_in_uC = NULL;
        TNamed *MedleyTarget= NULL;
        TNamed *duration_of_the_run_in_sec= NULL;
        TNamed *configuration_telescopes= NULL;
        TNamed *Processed_with = new TNamed("processed_with","processRun_Specific.cpp");
        TNamed *gamma_flash_in_ns_T1 = new TNamed("gFlash_T1(ns)",to_string(Gflash).c_str());
    
        TNamed *dist_Medleytarget_T1_in_mm = new TNamed("dist_Medleytarget_T1(mm)",to_string(1e3*L).c_str());
        

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
    TH1D * hp,*hd,*ht,*hh,*ha;
    TH2D * hpn,*hdn,*htn,*hhn,*han;


    hp = new TH1D(Form("protonsE%d",telN),Form("total energy for protons Tel %d", telN), 1000,0,50);
    hd = new TH1D(Form("deuteronsE%d",telN),Form("total energy for deuterons Tel %d", telN), 1000,0,50);
    ht = new TH1D(Form("tritonsE%d",telN),Form("total energy for tritons Tel %d", telN), 1000,0,50);
    hh = new TH1D(Form("heliumE%d",telN),Form("total energy for helium3 Tel %d", telN), 1000,0,50);
    ha = new TH1D(Form("alphasE%d",telN),Form("total energy for alphas Tel %d", telN), 1000,0,50);


    hpn = new TH2D(Form("protonsE_En%d",telN),Form("Eprotons:ENN_Tel %d", telN), 1000,0,50, 1000,0,50);
    hdn = new TH2D(Form("deuteronsE_En%d",telN),Form("Edeuterons:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
    htn = new TH2D(Form("tritonsE_En%d",telN),Form("Etritons:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
    hhn = new TH2D(Form("heliumE_En%d",telN),Form("Ehelium:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
    han = new TH2D(Form("alphasE_En%d",telN),Form("Ealphas:ENN_Tel%d", telN), 1000,0,50, 1000,0,50);
    
    newtree->Draw(Form("energy>>protonsE%d",telN),Form("PID==1"));
    newtree->Draw(Form("energy>>deuteronsE%d",telN),Form("PID==2"));
    newtree->Draw(Form("energy>>tritonsE%d",telN),Form("PID==3"));
    newtree->Draw(Form("energy>>heliumE%d",telN),Form("PID==4"));
    newtree->Draw(Form("energy>>alphasE%d",telN),Form("PID==5"));

    newtree->Draw(Form("energy:ENN>>protonsE_En%d",telN),Form("PID==1"));
    newtree->Draw(Form("energy:ENN>>deuteronsE_En%d",telN),Form("PID==2"));
    newtree->Draw(Form("energy:ENN>>tritonsE_En%d",telN),Form("PID==3"));
    newtree->Draw(Form("energy:ENN>>heliumE_En%d",telN),Form("PID==4"));
    newtree->Draw(Form("energy:ENN>>alphasE_En%d",telN),Form("PID==5"));

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    outputFile->Write();
    //
    // Close the output file
    outputFile->Close();


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
    timer.Print();
    cout<<"Compiled from processRun_Specific.cpp"<<endl;
    return 0;

}
