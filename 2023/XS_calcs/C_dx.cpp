//
// Script for evaluation the cross section for protons in Carbon__ 2023 runs
//
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <fstream>
//#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002
#include "/mnt/medley/LucasAnalysis/2023/applyTTC/correctSpec2.cpp" // for correcting the spectrum
#include "include/run_loader.hh"
#include "TStopwatch.h"
#include <time.h>
//#include "/home/pi/ganil/kalscripts/eloss/GetPPlot.cpp"


Double_t mc = 0.0408;//g 
Double_t M = 12; //g/mol
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Double_t omega_tel = 0.040 ;//sr 




void C_dx(Double_t Ea = 29.0, Double_t Eb = 30.0,Int_t runa = 35,Int_t runb = 39, Float_t Ereq = 32.5){
TStopwatch timer;

Float_t TotalChargeFaraday;
TTree *tx = loadRuns(runa,runb, &TotalChargeFaraday);


TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", "READ");
TH1D *nflux = (TH1D*)ff->Get("nflux");
TCanvas * cflux = new TCanvas("cflux","cflux");
nflux->Draw();

Double_t nflux_En = nflux->GetBinContent(nflux->FindBin((Ea+Eb)/2.0));
cout << "nflux_En: " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;


TCanvas * hp_cv = new TCanvas("hp","hp");
TH1D *hp = new TH1D("hp","hp",100,0,40);
tx->Draw("E>>hp",Form("ENN>%f && ENN<%f && PID==2 && ang == 20",Ea, Eb));

//calculate corrected spectrum
TH1D *hTTC = correctSpec(hp,true,"/home/pi/ganil/kalscripts/eloss/results/UniformZ/MedleyCarbon/v6/eloss_d_20.0deg_075.0um.root","Carbon",75,'d',false);
hTTC->Draw("same");

TCanvas * xs_cv = new TCanvas("XS","XS");

Double_t L = 464.72; //cm
Double_t A = TMath::Pi()*pow(2.5/2,2); //cm^2
Double_t cm2_to_barn = 1e+24 ;//cm²
Double_t Factor = L*L/(Nc*omega_tel*TotalChargeFaraday*nflux_En);
Factor = Factor*cm2_to_barn;
cout<<"Factor: "<<Factor<<endl;

TH1D *xsp = new TH1D("xsp","xsp",100,0,40);
TH1D *xspC = new TH1D("xspC","xspC",100,0,40);



for(int i=1;i<=hp->GetNbinsX();i++)
{
    xsp->SetBinContent(i,1e3*Factor*hp->GetBinContent(i)/hp->GetBinWidth(i));//1e3 changes to mb
    xspC->SetBinContent(i,1e3*Factor*hTTC->GetBinContent(i)/hTTC->GetBinWidth(i));//1e3 changes to mb
}
xsp->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
xsp->GetXaxis()->SetTitle("E_{protons} (MeV)");
//xsp->Draw();


TFile *fdata = new TFile("/mnt/medley/LucasAnalysis/exfor/exfor_files/myformat/datasets_C.root","READ");


TTree *txs = (TTree*)fdata->Get("DataTree");


float EN,E, Ang, XS, XSerr;
Char_t Particle[50], Dataset[50], Reference[50];

// 3. Use SetBranchAddress para associar as variáveis aos ramos da árvore
txs->SetBranchAddress("EN", &EN);
txs->SetBranchAddress("E", &E);
txs->SetBranchAddress("Ang", &Ang);
txs->SetBranchAddress("XS", &XS);
txs->SetBranchAddress("XSerr", &XSerr);
txs->SetBranchAddress("Particle", Particle);
txs->SetBranchAddress("Dataset", Dataset);
txs->SetBranchAddress("Reference", Reference);


// 4. Crie o TGraphErrors para armazenar os pontos e os erros
TGraphErrors *graph = new TGraphErrors();
TGraphErrors *g40 = new TGraphErrors();
TGraphErrors *g60 = new TGraphErrors();
TGraphErrors *g80 = new TGraphErrors();


string nameref;
// 5. Loop na árvore com a condição 'Ang == 20'
int nEntries = txs->GetEntries();
for (int i = 0; i < nEntries; i++) {
    // Carregue os dados para a entrada i
    txs->GetEntry(i);

    if(EN == Ereq){
        nameref = Reference;
        // Verifique a condição (Ang == 20)
        if (Ang == 20) {
            // Adicione o ponto ao gráfico, incluindo o erro
            graph->SetPoint(graph->GetN(), E, XS);             // Adiciona o ponto
            graph->SetPointError(graph->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
        }else if(Ang == 40){
            g40->SetPoint(g40->GetN(), E, XS);             // Adiciona o ponto
            g40->SetPointError(g40->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
        }else if(Ang == 60){
            g60->SetPoint(g60->GetN(), E, XS);             // Adiciona o ponto
            g60->SetPointError(g60->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
        }else if(Ang == 80){
            g80->SetPoint(g80->GetN(), E, XS);             // Adiciona o ponto
            g80->SetPointError(g80->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
        }
    }
    
}
graph->SetMarkerStyle(20);  // 20 é o código para círculos
graph->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
graph->SetMarkerColor(kGreen); 
graph->SetLineColor(kGreen); 
// 6. Desenhe o gráfico
graph->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
graph->GetXaxis()->SetTitle("E_{deuterons} (MeV)");
graph->Draw("ALP");  // "A" para eixos, "P" para pontos

g40->SetMarkerStyle(21);  // 20 é o código para círculos
g40->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
g40->SetMarkerColor(kBlack); 
g40->SetLineColor(kBlack); 

g60->SetMarkerStyle(22);  // 20 é o código para círculos
g60->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
g60->SetMarkerColor(kMagenta); 
g60->SetLineColor(kMagenta);

g80->SetMarkerStyle(47);  // 20 é o código para círculos
g80->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
g80->SetMarkerColor(kPink); 
g80->SetLineColor(kPink);

// g40->Draw("LP same ");  
// g60->Draw("LP same ");  
// g80->Draw("LP same ");  

graph->SetTitle("C(n,Xd), 20^{o}, 2023");

xspC->SetLineColor(kBlue);
xsp->SetLineWidth(2);
xsp->SetFillColor(kBlue);
xsp->SetFillStyle(3004);
xsp->Draw("same");


xspC->SetLineColor(kRed);
xspC->SetLineWidth(3);
xspC->SetFillColor(kRed);
xspC->SetFillStyle(3005);
xspC->Draw("same");


gPad->SetGridx();
gPad->SetGridy();

TLegend *tl = new TLegend(0.54,0.63,0.85,0.85);
tl->AddEntry(xsp,"Medley (2023)","lpf");
tl->AddEntry(xspC,"Medley TTC (2023)","lpf");
tl->AddEntry(graph,Form("%s, 20deg",nameref.c_str()),"lpf");
// tl->AddEntry(g40,Form("%s, 40deg",nameref.c_str()),"lpf");
// tl->AddEntry(g60,Form("%s, 60deg",nameref.c_str()),"lpf");
// tl->AddEntry(g80,Form("%s, 80deg",nameref.c_str()),"lpf");

tl->Draw();

timer.Print();
}

