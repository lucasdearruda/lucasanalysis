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
//#include "useful.h"
//in this verion we add the kaliveda plots 

using namespace std;


double En(double tof, double L = 4.6472){//L in meters, 2023 reference
	//return neutron energy in MeV, given its time of fligt.

	//if(tof<=0){ cout<<"ERROR. Tof<=0 is not allowed. Check the other steps. I will return -1 for that."<<endl; return -1;}
	if(tof<=0){ return 0;}
	double v = L/tof;
	double c = 0.299792458; //m/ns
	double gamma = 1/sqrt(1-(v/c)*(v/c));
	if(gamma!=gamma) return 0; //case where gamma = nan;
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


//function for getting gamma flash
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
		cout<<"	year "<<year<<" nor available!"<<endl;
		cout<<"****************************************************************"<<endl;
	}
	return tgamma;
}


int main(int argc, char* argv[]) {

    TStopwatch timer;

    // Check if the correct number of arguments are provided
    if (argc < 2 || argc > 3) {
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

    cout<<"Ppoint = "<<Ppoint<<" MeV. "<<endl;
    cout<<"Dpoint = "<<Dpoint<<" MeV. "<<endl;
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //CREATE THE REDUCED TTREE

    const char * outputFileName  = Form("/home/e802/Analysis/LucasAnalysis/reducedData/%03ddeb.root",nrun);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    TTree *newtree = new TTree("M","Medley's reduced data TTree");

    Int_t particle; 
    Int_t didPT;//it punched through?  
    Int_t ptflag; //punch through flag -- true if the event in in the PT region (i.e. X_pt_dE_nTel && X_npt_nTel) 
    Float_t energy; 
    Float_t si1; 
    Float_t si2; 
    Float_t csi; 
    Float_t energyKal; 
    Float_t tofn; 
    Float_t ang; 
    Float_t ENN; 

    TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
    TBranch * tofBranch  = newtree->Branch("tofn", &tofn, "tofn/F");
    TBranch * angleBranch  = newtree->Branch("ang", &ang, "ang/F");
    TBranch * EnBranch  = newtree->Branch("ENN", &ENN, "ENN/F");
    TBranch * EnKalBranch  = newtree->Branch("energykal", &energyKal, "EnergyKal/F");
    
    //debbugging part 
    TBranch * branchPT  = newtree->Branch("didPT", &didPT, "didPT/I");
    TBranch * branch_si1  = newtree->Branch("si1", &si1, "si1/F");
    TBranch * branch_si2  = newtree->Branch("si2", &si2, "si2/F");
    TBranch * branch_csi  = newtree->Branch("csi", &csi, "csi/F");
    TBranch * branch_ptflag  = newtree->Branch("ptflag", &ptflag, "ptflag/I");
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
    TCutG *cut_protons[5];
    TCutG *cut_deuterons[5];
    TCutG *cut_tritons[5];

    TCutG *cut_he[5]; 
    TCutG *cut_alphas[5];
    
    //cuts for pt particles
    TCutG *cut_protonsPT[5];
    TCutG *cut_deuteronsPT[5];
    TCutG *cut_tritonsPT[5];

    TCutG *cut_protonsCsI[5];
    TCutG *cut_deuteronsCsI[5];
    TCutG *cut_tritonsCsI[5]; 

    //cut npt csi 
    TCutG *cut_protonsNPTCsI[5];
    TCutG *cut_deuteronsNPTCsI[5];


    string cuts_path =   "/home/e802/Analysis/LucasAnalysis/PIDcuts3";
    
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
    //Loop over TTree
    TRandom2 *rand = new TRandom2();
    double sigmaSi1 = 0.08;//0.5
    double sigmaSi2 = 0.1;//.07
    double sigmacsi = 0.2;//.07
    double si1_plus_si2;//.07



    Long64_t nEntries = tx->GetEntries();

    cout<<" Total entries : "<<nEntries<<endl; 

    // Loop through each entry again to calculate and set the new branch values

    double combined_sigma = sqrt(pow(sigmaSi1,2)+pow(sigmaSi2,2));
    cout<<"combined sigma= "<<combined_sigma<<endl;
    for (Long64_t entry = 0; entry < nEntries*percentage/100.0; entry++) {
        //get entry 
        tx->GetEntry(entry);
        particle = 0; //particle starts unidentified
        tofn = 0; //tof starts zero too 

        for(int telN = 1; telN <= 4; telN++){
            //critical verifications: 
            invTof = tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue();
            
            if(invTof > Gflash[telN]) break; //if the event is faster than gamma... 

            de1 =  tx->GetLeaf(Form("Medley_%d_dE1",telN))->GetValue();
            de2 =  tx->GetLeaf(Form("Medley_%d_dE2",telN))->GetValue();
            Eres =  tx->GetLeaf(Form("Medley_%d_Eres",telN))->GetValue();
            energy = de1+de2+Eres;
            energyKal = -1;

            si1 = de1;
            si2 = de2;
            csi = Eres;
            ptflag = 0;
            if( de1 >= 0.02 && de2 >= 0.02 ){ // if it is not trash -- energy threshold

                //Randomized the PT decision 
                if( (cut_protons[telN]->IsInside(de2,de1) && cut_protonsNPTCsI[telN]->IsInside(Eres,de2))||(cut_protonsPT[telN]->IsInside(de2,de1) && cut_protonsCsI[telN]->IsInside(Eres,de2)) ){
                    //it is a proton!
                    particle = 1;
                    if((cut_protonsPT[telN]->IsInside(de2,de1) && cut_protonsCsI[telN]->IsInside(Eres,de2))){
                    //notORIGINAL: if((cut_protonsPT[telN]->IsInside(rand->Gaus(de2,sigmaSi2),rand->Gaus(de1,sigmaSi1)) && cut_protonsCsI[telN]->IsInside(rand->Gaus(Eres,sigmacsi),rand->Gaus(de2,sigmaSi2)))){                        
                        didPT = 1;
                    }else{
                        didPT = 0;
                    }

                    //condition for ptflag
                    if((cut_protons[telN]->IsInside(de2,de1) && cut_protonsPT[telN]->IsInside(de2,de1))){                    
                        ptflag = 1;
                    }

                }else if((cut_deuterons[telN]->IsInside(de2,de1) && cut_deuteronsNPTCsI[telN]->IsInside(Eres,de2))||(cut_deuteronsPT[telN]->IsInside(de2,de1)&&cut_deuteronsCsI[telN]->IsInside(Eres,de2)) ){
                    //it is a deuteron!
                    particle = 2;
                     if((cut_deuteronsPT[telN]->IsInside(de2,de1)&&cut_deuteronsCsI[telN]->IsInside(Eres,de2))){
                    //notORIGINAL: if((cut_deuteronsPT[telN]->IsInside(rand->Gaus(de2,sigmaSi2),rand->Gaus(de1,sigmaSi1))&&cut_deuteronsCsI[telN]->IsInside(rand->Gaus(Eres,sigmacsi),rand->Gaus(de2,sigmaSi2)))){
                        didPT = 1;
                    }else{
                        didPT = 0;
                    }

                    //condition for ptflag
                    if((cut_deuterons[telN]->IsInside(de2,de1) && cut_deuteronsPT[telN]->IsInside(de2,de1))){                    
                        ptflag = 1;
                    }

                }else if(cut_tritons[telN]->IsInside(de2,de1)||(cut_tritonsPT[telN]->IsInside(de2,de1)&&cut_tritonsCsI[telN]->IsInside(Eres,de2))){
                    //it is a triton!
                    particle = 3;
                    if((cut_tritonsPT[telN]->IsInside(de2,de1)&&cut_tritonsCsI[telN]->IsInside(Eres,de2))){
                    //notORIGINAL: if((cut_tritonsPT[telN]->IsInside(rand->Gaus(de2,sigmaSi2),rand->Gaus(de1,sigmaSi1))&&cut_tritonsCsI[telN]->IsInside(rand->Gaus(Eres,sigmacsi),rand->Gaus(de2,sigmaSi2)))){
                         didPT = 1;
                    }else{
                        didPT = 0;
                    }

                    //condition for ptflag
                    if((cut_deuterons[telN]->IsInside(de2,de1) && cut_deuteronsPT[telN]->IsInside(de2,de1))){                    
                        ptflag = 1;
                    }

                }else if(cut_he[telN]->IsInside(de2,de1)){
                    //it is a helium!
                    particle = 4;
                    didPT = -1;
                }else if(cut_alphas[telN]->IsInside(de2,de1)){
                    //it is an alpha!
                    particle = 5;
                    didPT = -1;
                }
            }

            //if a particle was identified

            if(particle){
                
                ang = 80+20.*telN;
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
                
                si1_plus_si2 = rand->Gaus(de1,sigmaSi1) + rand->Gaus(de2,sigmaSi2);
                //si1_plus_si2 = de1 + de2;
                
                if(particle<=3 && didPT>0){
                    //energyKal = gt[telN][particle]->Eval(de1+de2);
                    energyKal = gt[telN][particle]->Eval(si1_plus_si2);
                    energyKal = rand->Gaus(energyKal,combined_sigma);
                }else{
                    energyKal = si1_plus_si2;
                    //energyKal = de1 + de2; 
                }




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
        TNamed *charge_in_uC = NULL;
        TNamed *MedleyTarget= NULL;
        TNamed *duration_of_the_run_in_sec= NULL;
        TNamed *configuration_telescopes= NULL;
        TNamed *processed_with= NULL;
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

            // Print or use the information from each entry
            // std::cout << "Entry " << i << ":\n";
            // std::cout << "RunN: " << RunN << "\n";
            // std::cout << "MedleyConfig: " << MedleyConfig << "\n";
            // std::cout << "Target: " << Target << "\n";
            // std::cout << "RunTime: " << RunTime << "\n";
            // std::cout << "RunCharge: " << RunCharge << "\n";
        }

    processed_with = new TNamed("processed_with","processRun_debugging_Reversed.cpp");
    processed_with->Write();
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
    
        newtree->Draw(Form("energy>>protonsE%d",telN),Form("PID==1&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy>>deuteronsE%d",telN),Form("PID==2&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy>>tritonsE%d",telN),Form("PID==3&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy>>heliumE%d",telN),Form("PID==4&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy>>alphasE%d",telN),Form("PID==5&&ang==%.1f",80+20.*telN));

        newtree->Draw(Form("energy:ENN>>protonsE_En%d",telN),Form("PID==1&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy:ENN>>deuteronsE_En%d",telN),Form("PID==2&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy:ENN>>tritonsE_En%d",telN),Form("PID==3&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy:ENN>>heliumE_En%d",telN),Form("PID==4&&ang==%.1f",80+20.*telN));
        newtree->Draw(Form("energy:ENN>>alphasE_En%d",telN),Form("PID==5&&ang==%.1f",80+20.*telN));

    }

    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    outputFile->Write();
    //
    // Close the output file
    outputFile->Close();


    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
    timer.Print();
    cout<<"Compiled from processRun_debugging_Reversed.cpp"<<endl;
    return 0;

}
