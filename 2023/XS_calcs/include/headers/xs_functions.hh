#ifndef XS_FUNCTIONS_HH
#define XS_FUNCTIONS_HH

#include <TH1D.h>
#include <cmath>
#include <TGraph2DErrors.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TGraph.h>
#include <set>
#include <iostream>

// Declarar variáveis globais como extern para evitar múltiplas definições
extern Double_t mc;
extern Double_t M;

extern const Double_t NA;     // constante, pode ser const mesmo
extern Double_t Nc;           // essa depende de mc e M, precisa definir no .cxx
extern const Double_t omega_tel;

// Protótipos de funções...
// ex:
// Float_t histIntegralError(TH1D * histo);


//Declaring of the functions : _ _ : - - : - - : _ _ :  _ _ : - - : - - : _ _ :  _ _ : - - : - - : _ _ :  _ _ : - - : - - : _ _ :  _ _ : - - : - - : _ _ : 
Double_t GetSi1Thickness(Double_t angle = 20);
std::string DoubleToXpY(double value) ;
std::pair<int, int> get_ZA(char particle);
Int_t pCode(char particle = 'p');
Float_t histIntegralError(TH1D * histo);
Float_t totalCharge(vector<Int_t> runs, bool verbose = true, string infofile =  "/mnt/medley/LucasAnalysis/2023/runlist.csv");
std::vector<std::vector<TH1D*>> GetDDX(
    float Ea = 25,
    float Eb = 26,
    const std::vector<float>& angles =  {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0},
    const std::vector<char>& particles = {'p', 'd', 't', 'h', 'a'},
    const std::vector<int>& runsForward = {},
    const std::vector<int>& runsBackward = {},
    bool Match_CORR = true,
    bool TTC_CORR = true,
    string target = "Fe",
    bool DrawFlux = false
);
//Long64_t LoadRunsToChain(TChain* chain, const std::vector<int>& runIntervals, const std::string& dir, bool verbose = true);
TGraph * histoToGraph(TH1D *hh);
Float_t histIntegralError(TH1D * histo);
void Attribute_Target(string target_name = "Fe");
#endif