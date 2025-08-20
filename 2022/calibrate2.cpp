#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include <assert.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
//#include "useful.h"
#include <algorithm> //to use function swap
#include <time.h> //for the time
#include <vector>
#include "TMinuit.h"

using namespace std; // so doesnt need to declare "std::"

string pathRuns = "/mnt/medley/runswithMM_2022";

//TGraph

// Variáveis globais para a função de custo
std::vector<TGraph *> g_myplots;
TH2D* g_raw1;
double g_g3;
cout << ".\n.\n.\n.Functions to use:" << endl;

cout << "  loadCuts(telN)" << endl
     << "     loads the cuts for the given telescope number (1-5)." << endl << endl;

cout << "  plotCuts(telN, type, c)" << endl
     << "     plots the cuts for the given telescope number (1-5)." << endl
     << "     Type 1 for dE1:dE2 cuts and 2 for dE2:Eres cuts." << endl
     << "     c is an optional TCanvas pointer." << endl << endl;

cout << "  getInterpolatedYs(TGraph* g, double x_target)" << endl
     << "     returns the interpolated y values for the given x_target in the TGraph g." << endl << endl;

cout << "  Q(TH2D* myraw, std::vector<TGraph*> myplot," << endl
     << "    Double_t X, Double_t threshold = 0.05)" << endl
     << "     calculates the Q value for the given TH2D histogram and TGraph vector." << endl << endl;

cout << "  definePlots(opt)" << endl
     << "     returns a vector of TGraph* with the plots defined." << endl
     << "     opt = 1 for si1si2 and 2 for si2csi." << endl << endl;

cout << "  setGuesses(guess1, guess2, guess3)" << endl
     << "     sets the guesses for the function of cost." << endl
     << "     Default is 1, 1, 1." << endl << endl;

cout << "  loadRuns(nruni, nrunf)" << endl
     << "     loads the runs from nruni to nrunf." << endl
     << "     Default is 63 to 63." << endl << endl;

cout << "  plotRawDataG(tel)" << endl
     << "     plots the raw data for the given telescope number (1-5)." << endl << endl;

cout << "  setBins(nbins)" << endl
     << "     sets the number of bins for the histograms." << endl
     << "     Default is 1500." << endl << endl;

cout << "  plotPlots(plot_nr, c)" << endl
     << "     plots the defined plots." << endl
     << "     plot_nr = 1 for si1si2 and 2 for si2csi." << endl
     << "     c is an optional TCanvas pointer." << endl << endl;

cout << "  printCampaignCals()" << endl
     << "     print the calibrations for 2022 campaign." << endl;
     

cout << "\nExample:" << endl;
cout << "  int telN = 1;" << endl;
cout << "  loadRuns(69, 75);" << endl;
cout << "  setGuesses(0.000402, 0.002275, 0.0010260);" << endl;
cout << "  loadCuts(telN);" << endl;
cout << "  loadPlots(telN, 0);" << endl;
cout << "  plotRawDataG(telN);" << endl;
cout << "  plotCuts(telN, 1, (TCanvas*)rawC_1);" << endl;
cout << "  plotCuts(telN, 2, (TCanvas*)rawC_2);" << endl;


cout << "\n----------------------------------------------------\n";
cout << "Calibration settings, campaign 2022:\n";
cout << "1 | 0.000402, 0.002275, 0.0016942782*\n";
cout << "2 | 0.000357333, 0.00219807, 0.00117863**\n";
cout << "3 | 0.000411, 0.002224, 0.0008225\n";
cout << "4 | 0.000441, 0.002247, 0.0009116\n";
cout << "5 | 0.001735, 0.002292, 1\n";
cout << "6 | 0.001888, 0.002310, 1\n";
cout << "7 | 0.001826, 0.002232, 1\n";
cout << "8 | 0.001872, 0.002250, 1\n";
cout << "* | Updated values 2025-07-10\n";
cout << "**| Updated values 2025-08-13\n";



//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//CUTS VARIABLES FOR PID:

//cuts for npt particles
TCutG *cut_protons[5], *cut_deuterons[5], *cut_tritons[5]; 
TCutG *cut_he[5], *cut_alphas[5];

//cuts for pt particles
TCutG *cut_protonsPT[5], *cut_deuteronsPT[5], *cut_tritonsPT[5];
TCutG *cut_protonsCsI[5], *cut_deuteronsCsI[5], *cut_tritonsCsI[5]; 

//cut npt csi 
TCutG *cut_protonsNPTCsI[5], *cut_deuteronsNPTCsI[5];

string cuts_path =   "/mnt/medley/LucasAnalysis/2022/PID";

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//
//  FILES AND PLOTS DEFINITIONS
//

TFile *g[6];

TTree *TR[6];

//std::vector<TGraph*> myplots;
TChain *tx = NULL;

Long64_t nEntries[6];
string filename[6];
string pathCal = "/mnt/medley/kaliveda_results/Telescopes_dE_E_results2023";

TGraph *si1si2[6];
TGraph *si2csi[6];

TH2D *raw1, *raw2;
TH2D *h1, *h2;

Double_t g1,g2,g3;
g1 = g2 = g3 = -1;
Int_t Nbins=1500;
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//
//IBM palette definition
//

const Int_t nColors = 5;
Int_t colors[nColors] = {
        TColor::GetColor("#648FFF"),  // Blue
        TColor::GetColor("#785EF0"),  // Purple
        TColor::GetColor("#DC267F"),  // Salmon
        TColor::GetColor("#FE6100"),  // Orange
        TColor::GetColor("#FFB000"),  // Gold
};
//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


void printCampaignCals(){
    cout << "\n----------------------------------------------------\n";
    cout << "Calibration settings, campaign 2022:\n";
    //cout << "1 | 0.000402, 0.002275, 0.0010260\n";0.0016942782
    cout << "1 | 0.000402, 0.002275, 0.0016942782*\n";
    //cout << "2 | 0.000416, 0.002238, 0.0011229\n";
    cout << "2 | 0.000357333, 0.00219807, 0.00117863**\n";
    cout << "3 | 0.000411, 0.002224, 0.0008225\n";
    cout << "4 | 0.000441, 0.002247, 0.0009116\n";
    cout << "5 | 0.001735, 0.002292, 1\n";
    cout << "6 | 0.001888, 0.002310, 1\n";
    cout << "7 | 0.001826, 0.002232, 1\n";
    cout << "8 | 0.001872, 0.002250, 1\n";
    cout << "* | Updated values 2025-07-10\n";
    cout << "**| Updated values 2025-08-13\n";
    cout << "----------------------------------------------------\n";
    
    return;
}

void loadCuts(int telN=1){
    cout<<"\n\n --- --- --- LOADING FOR TELESCOPE "<<telN<<" --- --- ---\n"<<endl;
         
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
        cut_tritonsCsI[telN]->SetLineColor(kBlack);

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

void setEnergyAliases(int telN=1 ){
    cout <<"g's values: g1 = "<<g1<<", g2 = "<<g2<<", g3 = "<<g3<<endl;
    
    tx->SetAlias("dE1", Form("Medley_%d_dE1*%f",telN,g1));
    tx->SetAlias("dE2", Form("Medley_%d_dE2*%f",telN,g2));
    tx->SetAlias("Eres", Form("Medley_%d_Eres*%f",telN,g3));
    tx->SetAlias("E", "dE1 + dE2 + Eres");        
    cout<<"Aliases set for telescope "<<telN<<":\n";
    cout<<"dE1 = Medley_"<<telN<<"_dE1 * "<<g1<<endl;
    cout<<"dE2 = Medley_"<<telN<<"_dE2 * "<<g2<<endl;
    cout<<"Eres = Medley_"<<telN<<"_Eres * "<<g3<<endl;
    cout<<"E = dE1 + dE2 + Eres"<<endl;
}

void setParticlesAliases(int telN=1 ){
   
    loadCuts(telN);
    tx->SetAlias("protons", Form("(p_npt_%d && pnptcsi_%d)||(p_pt_dE_%d && p_pt_%d)",telN,telN,telN,telN));
    tx->SetAlias("deuterons", Form("!protons && ((d_npt_%d && dnptcsi_%d)||(d_pt_dE_%d && d_pt_%d))",telN,telN,telN,telN));
    tx->SetAlias("tritons", Form("!protons &&!deuterons ((t_npt_%d && Eres<0.5)||(t_pt_dE_%d && t_pt_%d))",telN,telN,telN));
    tx->SetAlias("he3", Form("he3_%d",telN));
    tx->SetAlias("alphas", Form("alphas%d",telN));
    
    cout<<"Aliases set for telescope "<<telN<<":\n";
    cout<<"protons = (p_npt_"<<telN<<" && pnptcsi_"<<telN<<")||(p_pt_dE_"<<telN<<" && p_pt_"<<telN<<")"<<endl;
    cout<<"deuterons = !protons && ((d_npt_"<<telN<<" && dnptcsi_"<<telN<<")||(d_pt_dE_"<<telN<<" && d_pt_"<<telN<<"))"<<endl;
    cout<<"tritons = !protons &&!deuterons ((t_npt_"<<telN<<" && Eres<0.5)||(t_pt_dE_"<<telN<<" && t_pt_"<<telN<<"))"<<endl;
    cout<<"he3 = he3_"<<telN<<endl;
    cout<<"alphas = alphas"<<telN<<endl;
    
    if(cut_protons[telN] != nullptr){
        cut_protons[telN]->SetVarY("dE1");
        cut_protons[telN]->SetVarX("dE2");
        cout<<"Cut for protons set to dE1:dE2"<<endl;
    }
    
    if(cut_protonsPT[telN] != nullptr){    
        cut_protonsPT[telN]->SetVarY("dE1");
        cut_protonsPT[telN]->SetVarX("dE2");
        cout<<"Cut for protons PT set to dE1:dE2"<<endl;
    }

    if(cut_protonsCsI[telN]!= nullptr){
        cut_protonsCsI[telN]->SetVarY("dE2");
        cut_protonsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for protons CsI set to dE2:Eres"<<endl;
    }

    if(cut_protonsNPTCsI[telN]!=nullptr){
        cut_protonsNPTCsI[telN]->SetVarY("dE2");
        cut_protonsNPTCsI[telN]->SetVarX("Eres");    
        cout<<"Cut for protons NPT CsI set to dE2:Eres"<<endl;
    }

    if(cut_deuterons[telN] != nullptr){
        cut_deuterons[telN]->SetVarY("dE1");
        cut_deuterons[telN]->SetVarX("dE2");
        cout<<"Cut for deuterons set to dE1:dE2"<<endl;
    }
    
    if(cut_deuteronsPT[telN] != nullptr){
        cut_deuteronsPT[telN]->SetVarY("dE1");
        cut_deuteronsPT[telN]->SetVarX("dE2");
        cout<<"Cut for deuterons PT set to dE1:dE2"<<endl;
    }

    if(cut_deuteronsCsI[telN] != nullptr){
        cut_deuteronsCsI[telN]->SetVarY("dE2");
        cut_deuteronsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for deuterons CsI set to dE2:Eres"<<endl;
    }
    if(cut_deuteronsNPTCsI[telN] != nullptr){
        cut_deuteronsNPTCsI[telN]->SetVarY("dE2");
        cut_deuteronsNPTCsI[telN]->SetVarX("Eres");
        cout<<"Cut for deuterons NPT CsI set to dE2:Eres"<<endl;
    }

    if(cut_tritons[telN] != nullptr){
        cut_tritons[telN]->SetVarY("dE1");
        cut_tritons[telN]->SetVarX("dE2");
        cout<<"Cut for tritons set to dE1:dE2"<<endl;
    }
    if(cut_tritonsPT[telN] != nullptr){
        cut_tritonsPT[telN]->SetVarY("dE1");
        cut_tritonsPT[telN]->SetVarX("dE2");
        cout<<"Cut for tritons PT set to dE1:dE2"<<endl;
    }
    if(cut_tritonsCsI[telN] != nullptr){
        cut_tritonsCsI[telN]->SetVarY("dE2");
        cut_tritonsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for tritons CsI set to dE2:Eres"<<endl;  
    }
    if(cut_he[telN] != nullptr){
        cut_he[telN]->SetVarY("dE1");
        cut_he[telN]->SetVarX("dE2");
        cout<<"Cut for he3 set to dE1:dE2"<<endl;
    }
    if(cut_alphas[telN] != nullptr){
        cut_alphas[telN]->SetVarY("dE1");
        cut_alphas[telN]->SetVarX("dE2");
        cout<<"Cut for alphas set to dE1:dE2"<<endl;
    }
    
    return;
}


void plotCuts(int telN=1,int type = 1, TCanvas* c = nullptr){

    if(c == nullptr){
        c = new TCanvas(Form("c_cuts_%d",telN),Form("Cuts for telescope %d",telN),800,600);
    }
    c->cd();

    if(type == 1){

        if(cut_protons[telN] != nullptr){
         cut_protons[telN]->Draw();
        }else{
            cout<<Form("cut_protons[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_deuterons[telN] != nullptr){
            cut_deuterons[telN]->Draw("same");
        }else{
            cout<<Form("cut_deuterons[%d] is null. Skipping draw.", telN)<<endl;
        }
         
        if(cut_tritons[telN] != nullptr){
            cut_tritons[telN]->Draw("same");
        }else{
            cout<<Form("cut_tritons[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_he[telN] != nullptr){
            cut_he[telN]->Draw("same");
        }else{
            cout<<Form("cut_he[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_alphas[telN] != nullptr){
            cut_alphas[telN]->Draw("same");
        }else{
            cout<<Form("cut_alphas[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_protonsPT[telN] != nullptr){
            cut_protonsPT[telN]->Draw("same");
        }else{
            cout<<Form("cut_protonsPT[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_deuteronsPT[telN] != nullptr){
            cut_deuteronsPT[telN]->Draw("same");
        }else{
            cout<<Form("cut_deuteronsPT[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_tritonsPT[telN] != nullptr){
            cut_tritonsPT[telN]->Draw("same");
        }else{
            cout<<Form("cut_tritonsPT[%d] is null. Skipping draw.", telN)<<endl;
        } 

        // cut_protons[telN]->Draw("same");
        // cut_protonsPT[telN]->Draw("same");
        // cut_deuterons[telN]->Draw("same");
        // cut_deuteronsPT[telN]->Draw("same");
        // cut_tritons[telN]->Draw("same");
        // cut_tritonsPT[telN]->Draw("same");
        // cut_he[telN]->Draw("same");
        // cut_alphas[telN]->Draw("same");

    }else if(type == 2){
        //verify it each cut is not null before drawing
        if(cut_protonsCsI[telN] != nullptr){
            cut_protonsCsI[telN]->Draw("same");
        }else{
            cout<<Form("cut_protonsCsI[%d] is null. Skipping draw.", telN)<<endl;
        } 

        if(cut_deuteronsCsI[telN] != nullptr){
            cut_deuteronsCsI[telN]->Draw("same");
        }else{
            cout<<Form("cut_deuteronsCsI[%d] is null. Skipping draw.", telN)<<endl;
        }
        if(cut_tritonsCsI[telN] != nullptr){
            cut_tritonsCsI[telN]->Draw("same");
        }else{
            cout<<Form("cut_tritonsCsI[%d] is null. Skipping draw.", telN)<<endl;
        }
        if(cut_protonsNPTCsI[telN] != nullptr){
            cut_protonsNPTCsI[telN]->Draw("same");
        }else{
            cout<<Form("cut_protonsNPTCsI[%d] is null. Skipping draw.", telN)<<endl;
        } 
        if(cut_deuteronsNPTCsI[telN] != nullptr){
            cut_deuteronsNPTCsI[telN]->Draw("same");
        }else{
            cout<<Form("cut_deuteronsNPTCsI[%d] is null. Skipping draw.", telN)<<endl;
        } 

    }else{

        cout<<"Invalid type. Use 1 for dE1:dE2 cuts and 2 for dE2:Eres cuts."<<endl;
        return;

    }

    return;
}




std::vector<double> getInterpolatedYs(TGraph* g, double x_target) {
    std::vector<double> results;

    int n = g->GetN();
    double x1, y1, x2, y2;

    for (int i = 0; i < n - 1; ++i) {
        g->GetPoint(i, x1, y1);
        g->GetPoint(i + 1, x2, y2);

        // Verifica se x_target está entre x1 e x2 (em qualquer ordem)
        if ((x1 <= x_target && x_target <= x2) || (x2 <= x_target && x_target <= x1)) {
            // Interpolação linear
            double t = (x_target - x1) / (x2 - x1);
            double y_interp = y1 + t * (y2 - y1);
            results.push_back(y_interp);
        }
    }

    return results;
}


double Q(TH2D* myraw, std::vector<TGraph*> myplot, Double_t X, Double_t threshold = 0.05) {
    int binX = myraw->GetXaxis()->FindBin(X);
    int nbinsY = myraw->GetYaxis()->GetNbins();
    int minBinY = myraw->GetYaxis()->FindBin(threshold);
    double Q = 0.0;

    for (int binY = minBinY; binY <= nbinsY; ++binY) {
        double binCenterY = myraw->GetYaxis()->GetBinCenter(binY);
        double binContent = myraw->GetBinContent(binX, binY);
        if (binContent == 0) continue; // pula bins vazios

        double minVal = std::numeric_limits<double>::max();
        bool hasValid = false;

        for (auto graph : myplot) {
            auto ys = getInterpolatedYs(graph,X);

            if (!ys.empty()) {
                double localMin = *std::min_element(ys.begin(), ys.end());
                if (localMin < minVal) {
                    minVal = localMin;
                    hasValid = true;
                }
            }
        }

        if (hasValid) {
            Q += minVal * binContent;  // pondera pelo conteúdo do bin
        }
    }

    return Q;

}
std::vector<TGraph*> definePlots(int opt = 1){
    std::vector<TGraph*> myplots;
    if(opt == 1){
        cout<< "Setting plots to si1si2"<<endl;
        for (int i = 0; i < 5; i++) {
            myplots.push_back(si1si2[i]);
        }
    }else{
        cout<< "Setting plots to si2csi"<<endl;
        for (int i = 0; i < 5; i++) {
            myplots.push_back(si2csi[i]);
        }
    }
    
    return myplots;
}

void setGuesses(Double_t guess1=1, Double_t guess2 =1, Double_t guess3 = 1){
    g1 = guess1;
    g2 = guess2;
    g3 = guess3;
    cout<<"Guesses successfully set to:"<<endl;
    cout<<g1<<endl;
    cout<<g2<<endl;
    cout<<g3<<endl;
    cout<<" - - -"<<endl;
    return;
}


void setNbins(Double_t nbins = 1500){
    Nbins = nbins;
    cout<<"Setting nbins to:"<<nbins<<endl;
    cout<<" - - -"<<endl;
    return;
}

char particles[] = {'p','d','t','h','a'};

int lw = 2; //line width


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

void loadPlots(int tel = 1, bool verbose = true){


    cout<<"running for telescope "<<tel<<"..."<<endl;
    filename[0] = Form("%s/kal_p_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //p production plot
    filename[1] = Form("%s/kal_d_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //d production plot
    filename[2] = Form("%s/kal_t_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //t production plot
    filename[3] = Form("%s/kal_h_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //h production plot
    filename[4] = Form("%s/kal_a_Tel%d_wToF_v20240422.root",pathCal.c_str(),tel); //a production plot

    //initializing the TGraph
    for(int i=0; i<5; i++){
        si1si2[i] = new TGraph();
        si1si2[i]->SetNameTitle(Form("si1si2_%c",particles[i]),Form("si1si2_%c",particles[i]));
    
        si1si2[i]->GetXaxis()->SetTitle("si2 (MeV)");
        si1si2[i]->GetYaxis()->SetTitle("si1 (MeV)");
        si1si2[i]->SetLineColor(colors[i]);
        si1si2[i]->SetLineWidth(lw);
        si1si2[i]->SetMarkerColor(colors[i]);

        si2csi[i] = new TGraph();
        si2csi[i]->SetNameTitle(Form("si2csi_%c",particles[i]),Form("si2csi_%c",particles[i]));
        si2csi[i]->GetXaxis()->SetTitle("csi (MeV)");
        si2csi[i]->GetYaxis()->SetTitle("si2 (MeV)");
        si2csi[i]->SetLineColor(colors[i]);
        si2csi[i]->SetLineWidth(lw);
        si2csi[i]->SetMarkerColor(colors[i]);
        
    }

    Double_t si1,si2,csi;
    //Tfile and TTree attribution
    for(int i=0; i<5; i++){ //for each file 
        cout<<"Reading .."<<i<<", from "<<filename[i]<<" :"<<endl;
        g[i] = new TFile(filename[i].c_str());
            TR[i] = (TTree*)g[i]->Get("SIM");
            nEntries[i] = TR[i]->GetEntries();
            if(verbose)cout<<Form("Entries in the tree SIM (%s): ",filename[i].c_str())<<nEntries[i]<< endl;

            TR[i]->SetBranchAddress("si1",&si1);
            TR[i]->SetBranchAddress("si2",&si2);
            TR[i]->SetBranchAddress("csi",&csi);

            for(Long64_t j=0; j<nEntries[i]; j++){
                TR[i]->GetEntry(j);
                si1si2[i]->SetPoint(j,si2,si1);
                si2csi[i]->SetPoint(j,csi,si2);
                //if(verbose) cout<<"Point "<<j<<" : "<<si1<<" , "<<si2<<" , "<<csi<<endl;
            }
        g[i]->Close(); // Close the file after reading
    }

}
void plotPlots(int plot_nr = 1, TCanvas* c = nullptr) {
    if (!c) {
        c = new TCanvas(Form("c_plot%d", plot_nr), "Auto-created Canvas", 800, 600);
    }
    c->cd();

    if (plot_nr == 1) {
        for (int i = 4; i >= 0; --i) {
            si1si2[i]->Draw(i == 4 ? "" : "same");
        }
    } else if (plot_nr == 2) {
        // Cria histograma frame vazio para definir eixos fixos de 0 a 40 em X e Y
        TH2F* frame = new TH2F("frame", "Frame for plot 2", 40, 0, 40, 40, 0, 40);
        frame->SetStats(0); // sem estatísticas no canto
        if(c != nullptr)frame->Draw("same");      // desenha o frame primeiro
        else frame->Draw();      // desenha o frame primeiro

        // Agora desenha todos os histogramas com "same"
        for (int i = 4; i >= 0; --i) {
            si2csi[i]->Draw("same");
        }
    } else {
        cout << "Invalid plot number. Please choose 1 or 2." << endl;
    }
}

TChain* loadRuns(int nruni =63 ,int nrunf =63){

    if(!nruni){
        cout<<"No run number given, using default run 111"<<endl;
        nruni = 111;
        nrunf = 111;
    }

    char name[1000];
    if (!tx)
        tx = new TChain("RD");

    vector<Int_t> entries;

    for (int i = nruni; i <= nrunf; i++)
    {
        Bool_t fExist = true;
        for (Int_t j = 0; j <= 999 && fExist; j++)
        {
            ifstream mfile;
            sprintf(name, "%s/r%04d_%03dr.root",pathRuns.c_str(), i, j);
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
    }

    return tx;

}

TCanvas *rawC = NULL;

void plotRawDataG(int tel = 1){

    
    if(!tx){
        cout<<"There are no runs loaded. Use loadRuns(runA,runB) to load some."<<endl;
        return;
    }
    if(g1*g2*g3 < 0){
        cout<<"Guesses not set, using default values: g1 = 1, g2 = 1, g3 = 1"<<endl;
        g1 = g2 = g3 = 1;
    }   


    if(!rawC) {
        rawC = new TCanvas("rawC","raw data", 150,10,1600,800);  
        rawC->Divide(2,1);
    }


    rawC->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();
    
    if (raw1) {
        delete raw1;
        raw1 = nullptr;
        gDirectory->Delete("raw1;*");  // Make sure to clear ROOT's internal copy too
    }
    if (raw2) {
        delete raw2;
        raw2 = nullptr;
        gDirectory->Delete("raw2;*");
    }

    
    raw1 = new TH2D("raw1","raw1", Nbins,0,40,Nbins,0,10);
    raw2 = new TH2D("raw2","raw2", Nbins,0,40,Nbins,0,25);

    tx->Draw(Form("Medley_%d_SI_DE1*%f:Medley_%d_SI_DE2*%f>>raw1", tel,g1,tel,g2),"","colz");


    rawC->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    tx->Draw(Form("Medley_%d_SI_DE2*%f:Medley_%d_Eres*%f>>raw2", tel,g2, tel,g3),"","colz");

}


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double g1 = par[0];
    double g2 = par[1];
    double g3 = g_g3;

    setGuesses(g1, g2, g3);
    loadPlots(1,0);
    g_myplots = definePlots();

    f = Q(g_raw1, g_myplots, 5)
      + Q(g_raw1, g_myplots, 10)
      + Q(g_raw1, g_myplots, 15)
      + Q(g_raw1, g_myplots, 20);
}


void runMinuitFit() {
    g_raw1 = raw1;  // seu TH2D*
    g_g3 = 0.0010260; // g3 fixo

    TMinuit minuit(2); // dois parâmetros: g1 e g2
    minuit.SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;
    double maxpct = 0.1; // tolerância de convergência
    double minpct = 0.1; // tolerância mínima de convergência

    // Definindo parâmetros iniciais conforme seus valores
    minuit.mnparm(0, "g1", 0.000402, 0.00001, 0.000402*(1-minpct), 0.000402*(1+maxpct), ierflg);
    minuit.mnparm(1, "g2", 0.002275, 0.00001, 0.000402*(1-minpct), 0.000402*(1+maxpct), ierflg);

    // Rodar a minimização (MIGRAD)
    arglist[0] = 1000; // máximo de iterações
    minuit.mnexcm("MIGRAD", arglist, 1, ierflg);

    double g1_fit, g2_fit, err1, err2;
    minuit.GetParameter(0, g1_fit, err1);
    minuit.GetParameter(1, g2_fit, err2);

    std::cout << "Fit results:" << std::endl;
    std::cout << "g1 = " << g1_fit << " ± " << err1 << std::endl;
    std::cout << "g2 = " << g2_fit << " ± " << err2 << std::endl;
    std::cout << "g3 (fixo) = " << g_g3 << std::endl;
}