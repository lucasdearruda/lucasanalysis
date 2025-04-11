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

#include "TStopwatch.h"
#include <time.h>

Double_t mc = 0.0182;//g 
Double_t M = 55.845; //g/mol
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Double_t omega_tel = 0.040 ;//sr 


Int_t Ta;
//thin 389 395
void Fe_px(Double_t Ea = 31, Double_t Eb = 32 ,Int_t runa = 396,Int_t runb = 405){
TStopwatch timer;


TChain *tx = new TChain("M");


string processed_runs = "/mnt/medley/LucasAnalysis/2023/reducedv61";

string name;
for(Int_t i=runa;i<=runb;i++)
{
    Bool_t fExist = true;
    name = Form("%s/%03d.root",processed_runs.c_str(),i);
    cout << name << endl;
    ifstream mfile;
    mfile.open(name);
    if(mfile)
    {
        mfile.close();
    }
    else
        fExist=false;
    if(fExist)
    {	  
        cout << "Adding " << name << endl;
        tx->Add(name.c_str());
        Ta = (Int_t) (tx->GetEntries());
        cout << "Entries " << Ta/1000 << "k "  << endl;
    }
    
}

cout<<"---------------------------------------------------"<<endl;
TTree *InfoTree = new TTree("InfoTree", "InfoTree");
InfoTree->ReadFile("/mnt/medley/LucasAnalysis/2023/runlist.csv", "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F");

Int_t RunN;
Float_t ChargeFaraday;
Float_t TotalChargeFaraday = 0;

InfoTree->SetBranchAddress("RunN", &RunN);
InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

for(int i=0;i<InfoTree->GetEntries();i++)
{
    InfoTree->GetEntry(i);
    //cout<< i<<", nrun = "<<RunN<< endl;
    if (RunN<= runb && RunN>=runa)
    {
        //InfoTree->GetEntry(i);
        cout << "RunN: " << RunN << " ChargeFaraday: " << ChargeFaraday << endl;
        TotalChargeFaraday += ChargeFaraday;
    }
}
cout<<"---------------------------------------------------"<<endl;
cout<< "Total Charge Faraday: " << TotalChargeFaraday <<" µC"<< endl;


TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", "READ");
TH1D *nflux = (TH1D*)ff->Get("nflux");
TCanvas * cflux = new TCanvas("cflux","cflux");
nflux->Draw();

Double_t nflux_En = nflux->GetBinContent(nflux->FindBin((Ea+Eb)/2.0));
cout << "nflux_En ("<<(Ea+Eb)/2.0<<" MeV): " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;


TCanvas * hp_cv = new TCanvas("hp","hp");
TH1D *hp = new TH1D("hp","hp",100,0,40);
tx->Draw("E>>hp",Form("ENN>%f && ENN<%f && PID==1 && ang == 20",Ea, Eb));

//calculate corrected spectrum
TH1D *hTTC = correctSpec(hp,true,"/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick/v7/eloss_p_20.0deg_025.0um.root","Fe_thick_Medley",25,'p',false);
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
xsp->Draw();
/*
TTree *txs = new TTree("txs","txs");

//E          ANG        DATA      1ERR-T     1DATA      2ERR-T     2
//MEV        ADEG       MB/SR/MEV  MB/SR/MEV  MB/SR/MEV  MB/SR/MEV
txs->ReadFile("/mnt/medley/LucasAnalysis/exfor/exfor_files/n.xP_26.5MeV_Slypen.csv","E/F:Ang/F:XS/F:XSerr/F");
//new TCanvas();

float E, Ang, XS, XSerr;

// 3. Use SetBranchAddress para associar as variáveis aos ramos da árvore
txs->SetBranchAddress("E", &E);
txs->SetBranchAddress("Ang", &Ang);
txs->SetBranchAddress("XS", &XS);
txs->SetBranchAddress("XSerr", &XSerr);

// // 4. Crie o TGraphErrors para armazenar os pontos e os erros
// TGraphErrors *graph = new TGraphErrors();
// TGraphErrors *g40 = new TGraphErrors();
// TGraphErrors *g60 = new TGraphErrors();
// TGraphErrors *g80 = new TGraphErrors();

// // 5. Loop na árvore com a condição 'Ang == 20'
// int nEntries = txs->GetEntries();
// for (int i = 0; i < nEntries; i++) {
//     // Carregue os dados para a entrada i
//     txs->GetEntry(i);

//     // Verifique a condição (Ang == 20)
//     if (Ang == 20) {
//         // Adicione o ponto ao gráfico, incluindo o erro
//         graph->SetPoint(graph->GetN(), E, XS);             // Adiciona o ponto
//         graph->SetPointError(graph->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 40){
//         g40->SetPoint(g40->GetN(), E, XS);             // Adiciona o ponto
//         g40->SetPointError(g40->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 60){
//         g60->SetPoint(g60->GetN(), E, XS);             // Adiciona o ponto
//         g60->SetPointError(g60->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 80){
//         g80->SetPoint(g80->GetN(), E, XS);             // Adiciona o ponto
//         g80->SetPointError(g80->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }
// }
// graph->SetMarkerStyle(20);  // 20 é o código para círculos
// graph->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// graph->SetMarkerColor(kGreen); 
// graph->SetLineColor(kGreen); 
// // 6. Desenhe o gráfico
// graph->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
// graph->GetXaxis()->SetTitle("E_{protons} (MeV)");
// graph->Draw("ALP");  // "A" para eixos, "P" para pontos

// g40->SetMarkerStyle(21);  // 20 é o código para círculos
// g40->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g40->SetMarkerColor(kBlack); 
// g40->SetLineColor(kBlack); 

// g60->SetMarkerStyle(22);  // 20 é o código para círculos
// g60->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g60->SetMarkerColor(kMagenta); 
// g60->SetLineColor(kMagenta);

// g80->SetMarkerStyle(47);  // 20 é o código para círculos
// g80->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g80->SetMarkerColor(kPink); 
// g80->SetLineColor(kPink);

// g40->Draw("LP same ");  
// g60->Draw("LP same ");  
// g80->Draw("LP same ");  

// graph->SetTitle("C(n,Xp), 20^{o}, 2023");

// xspC->SetLineColor(kBlue);
// xsp->SetLineWidth(2);
// xsp->SetFillColor(kBlue);
// xsp->SetFillStyle(3004);
// xsp->Draw("same");


// xspC->SetLineColor(kRed);
// xspC->SetLineWidth(3);
// xspC->SetFillColor(kRed);
// xspC->SetFillStyle(3005);
// xspC->Draw("same");


// gPad->SetGridx();
// gPad->SetGridy();

// TLegend *tl = new TLegend(0.54,0.63,0.85,0.85);
// tl->AddEntry(xsp,"Medley (2023)","lpf");
// tl->AddEntry(xspC,"Medley TTC (2023)","lpf");
// tl->AddEntry(gra
// // 4. Crie o TGraphErrors para armazenar os pontos e os erros
// TGraphErrors *graph = new TGraphErrors();
// TGraphErrors *g40 = new TGraphErrors();
// TGraphErrors *g60 = new TGraphErrors();
// TGraphErrors *g80 = new TGraphErrors();

// // 5. Loop na árvore com a condição 'Ang == 20'
// int nEntries = txs->GetEntries();
// for (int i = 0; i < nEntries; i++) {
//     // Carregue os dados para a entrada i
//     txs->GetEntry(i);

//     // Verifique a condição (Ang == 20)
//     if (Ang == 20) {
//         // Adicione o ponto ao gráfico, incluindo o erro
//         graph->SetPoint(graph->GetN(), E, XS);             // Adiciona o ponto
//         graph->SetPointError(graph->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 40){
//         g40->SetPoint(g40->GetN(), E, XS);             // Adiciona o ponto
//         g40->SetPointError(g40->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 60){
//         g60->SetPoint(g60->GetN(), E, XS);             // Adiciona o ponto
//         g60->SetPointError(g60->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }else if(Ang == 80){
//         g80->SetPoint(g80->GetN(), E, XS);             // Adiciona o ponto
//         g80->SetPointError(g80->GetN() - 1, 0, XSerr); // Adiciona o erro no X (XSerr) e no Y (0, pois não temos erro em Y)
//     }
// }
// graph->SetMarkerStyle(20);  // 20 é o código para círculos
// graph->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// graph->SetMarkerColor(kGreen); 
// graph->SetLineColor(kGreen); 
// // 6. Desenhe o gráfico
// graph->GetYaxis()->SetTitle("mb/sr#dot1-MeV");
// graph->GetXaxis()->SetTitle("E_{protons} (MeV)");
// graph->Draw("ALP");  // "A" para eixos, "P" para pontos

// g40->SetMarkerStyle(21);  // 20 é o código para círculos
// g40->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g40->SetMarkerColor(kBlack); 
// g40->SetLineColor(kBlack); 

// g60->SetMarkerStyle(22);  // 20 é o código para círculos
// g60->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g60->SetMarkerColor(kMagenta); 
// g60->SetLineColor(kMagenta);

// g80->SetMarkerStyle(47);  // 20 é o código para círculos
// g80->SetMarkerSize(1.5);  // Tamanho dos pontos (ajuste conforme necessário)
// g80->SetMarkerColor(kPink); 
// g80->SetLineColor(kPink);

// g40->Draw("LP same ");  
// g60->Draw("LP same ");  
// g80->Draw("LP same ");  

// graph->SetTitle("C(n,Xp), 20^{o}, 2023");

// xspC->SetLineColor(kBlue);
// xsp->SetLineWidth(2);
// xsp->SetFillColor(kBlue);
// xsp->SetFillStyle(3004);
// xsp->Draw("same");


// xspC->SetLineColor(kRed);
// xspC->SetLineWidth(3);
// xspC->SetFillColor(kRed);
// xspC->SetFillStyle(3005);
// xspC->Draw("same");


// gPad->SetGridx();
// gPad->SetGridy();

// TLegend *tl = new TLegend(0.54,0.63,0.85,0.85);
// tl->AddEntry(xsp,"Medley (2023)","lpf");
// tl->AddEntry(xspC,"Medley TTC (2023)","lpf");
// tl->AddEntry(graph,"Slypen+(2000), 20deg","lpf");
// tl->AddEntry(g40,"Slypen+(2000), 40deg","lpf");
// tl->AddEntry(g60,"Slypen+(2000), 60deg","lpf");
// tl->AddEntry(g80,"Slypen+(2000), 80deg","lpf");

// tl->Draw();ph,"Slypen+(2000), 20deg","lpf");
// tl->AddEntry(g40,"Slypen+(2000), 40deg","lpf");
// tl->AddEntry(g60,"Slypen+(2000), 60deg","lpf");
// tl->AddEntry(g80,"Slypen+(2000), 80deg","lpf");

// tl->Draw();
*/
timer.Print();
}

