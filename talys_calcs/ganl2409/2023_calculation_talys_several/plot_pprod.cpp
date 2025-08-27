#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include <string>
#include <iostream>
#include <fstream>
#include "/home/dearruda/ganil/medley_2023/pre_analysis/libs/useful.h"

#include <time.h>

void plot_pprod(){
    //Script to plot the proton production for a given energy
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    //variables to insert data-time to output file---------------
    char cur_time[128];

    time_t      t;
    struct tm*  ptm;
    
    t = time(NULL);
    ptm = localtime(&t);
    
    strftime(cur_time, 128, "%Y-%m-%d_%H:%M:%S", ptm); 
    //start clock
    clock_t tStart = clock(); 
    
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

    float energy = 16;
    //TFile *g = new TFile("experimental/p28thick.root");
    
    TFile *g = new TFile("experimental/p15Thick.root");
    
    //TH1D *histo = (TH1D*)g->Get("hps0_13");
    TH1D *histo = (TH1D*)g->Get("hps0_6");
    gStyle->SetOptStat(0);
    TCanvas *Exp_canvas = new TCanvas("Exp","Exp",2250-1920,60,1000,800);
    histo->Scale(1/histo->GetMaximum());
    

    histo->SetLineWidth(3);
    histo->GetYaxis()->SetTitle("arbitrary untis [a.u.]");
    histo->SetTitle("");
    histo->Draw("hist");
    //reading TALYS' results 

    double factorX = 0.8;
   TH1D *scaledhisto = new TH1D("scaledhisto","scaledhisto",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),factorX*histo->GetXaxis()->GetXmax());
    for(int i=1;i<histo->GetNbinsX();i++ ){
        scaledhisto->SetBinContent(i,histo->GetBinContent(i));
        cout<<"i= "<<i<<" bincontent = "<<histo->GetBinContent(i)<<endl;
    }
    scaledhisto->SetLineColor(kPink);
    scaledhisto->SetLineWidth(3);
    //scaledhisto->Draw("same");

    string talys_file = Form("/home/dearruda/Documents/2023_calculation_talys_several/resimulation/Fe/pspec%08.3f.tot",energy);
    //string talys_file = Form("/home/dearruda/Documents/2023_calculation_talys_several/resimulation/C0/pspec%08.3f.tot",energy);
    cout<<"--> Opening "<<lastname(talys_file).c_str()<<endl; 

// # 
// # # energies =   185
// #  E-out    Total       Direct    Pre-equil.  Mult. preeq  Compound
    TTree *Cs = new TTree("CS", "CSTTree");
    Cs->ReadFile(talys_file.c_str(), "Eout/D:total/D:direct/D:preeq/D:multP/D:compound/D");
    float Emin, Emax;
    Emin = Cs->GetMinimum("Eout");
    Emax = Cs->GetMaximum("Eout");
    Cs->Print();
    cout<<"Emin = "<<Emin<<" MeV\tEmax = "<<Emax<<" MeV"<<endl;
    double maxtotal = Cs->GetMaximum("total");
    cout<<"Maximum (total)= "<<maxtotal<<" au."<<endl;


    const char * branch2;

    double er, s;    
    //TCanvas *P_canvas = new TCanvas("P","P",2000-1920,40,1000,800);
    TGraph *g_total = new TGraph(Cs->GetEntries());
    TGraph *g_direct = new TGraph(Cs->GetEntries());
    TGraph *g_preeq = new TGraph(Cs->GetEntries());
    TGraph *g_multP = new TGraph(Cs->GetEntries());
    TGraph *g_compound = new TGraph(Cs->GetEntries());

    Double_t eout, total, direct, preeq, multP, compound;

    Cs->SetBranchAddress("Eout", &eout);
    Cs->SetBranchAddress("total", &total);
    Cs->SetBranchAddress("direct", &direct);
    Cs->SetBranchAddress("preeq", &preeq);
    Cs->SetBranchAddress("multP", &multP);
    Cs->SetBranchAddress("compound", &compound);

        
    for (Long64_t i = 0; i < Cs->GetEntries(); i++) {
        Cs->GetEntry(i);
        g_total->SetPoint(i, eout, total*1/maxtotal);
        g_direct->SetPoint(i, eout, direct*1/maxtotal);
        g_preeq->SetPoint(i, eout, preeq*1/maxtotal);
        g_multP->SetPoint(i, eout, multP*1/maxtotal);
        g_compound->SetPoint(i, eout, compound*1/maxtotal);
    }
    
    //smooth plots:
    // TSpline3 *gs_total = new TSpline3("gs_total",g_total->GetX(),g_total->GetY(),g_total->GetN());
    // TSpline3 *gs_direct = new TSpline3("gs_direct",g_direct->GetX(),g_direct->GetY(),g_direct->GetN());
    // TSpline3 *gs_preeq = new TSpline3("gs_preeq",g_preeq->GetX(),g_preeq->GetY(),g_preeq->GetN());
    // TSpline3 *gs_multP = new TSpline3("gs_multP",g_multP->GetX(),g_multP->GetY(),g_multP->GetN());
    // TSpline3 *gs_compound = new TSpline3("gs_compound",g_compound->GetX(),g_compound->GetY(),g_compound->GetN());

    int line_thickness = 5;
    g_total->GetXaxis()->SetTitle("E_{p} (MeV)");
    g_total->GetYaxis()->SetTitle("#sigma [a.u.]");
    g_total->SetLineColor(kBlack);
    g_total->SetLineWidth(line_thickness);
    g_total->SetTitle("");
    g_total->Draw("same");

    g_direct->SetLineStyle(2);
    g_direct->SetLineColor(kRed);
    g_direct->SetLineWidth(line_thickness);
    g_direct->Draw("same");

    g_preeq->SetLineColor(kBlue);
    g_preeq->SetLineStyle(3);
    g_preeq->SetLineWidth(line_thickness);
    g_preeq->Draw("same");

    g_multP->SetLineColor(kGreen+3);
    g_multP->SetLineStyle(8);
    g_multP->SetLineWidth(line_thickness);
    g_multP->Draw("same");
    
    g_compound->SetLineColor(kMagenta);
    g_compound->SetLineStyle(5);
    g_compound->SetLineWidth(line_thickness);
    g_compound->Draw("same");
    
    
    gPad->SetGridx();
    gPad->SetGridy();    

    TLegend *lg = new TLegend(0.56,0.55,0.9,0.9);
    lg->SetHeader(Form("^{nat}Fe(n,pX), E_{N} = %.1f MeV", energy),"C");
    TLegendEntry *header = (TLegendEntry*)lg->GetListOfPrimitives()->First();
    //header->SetTextAlign(22);
    //header->SetTextColor(2);
    header->SetTextSize(0.035);
    //lg->AddEntry(Form()"#^{nat}Fe(n,pX), E_{N} = %.1f MeV", energy));
    lg->AddEntry(histo,Form("Exp. E_{NN} = %d to %d MeV",(int)(energy -1),(int)(energy +1)),"l");
    lg->AddEntry(g_total,"Total","l");
    lg->AddEntry(g_direct,"Direct","l");
    lg->AddEntry(g_preeq,"Pre-eq.","l");
    lg->AddEntry(g_multP,"Mult. pre-eq.","l");
    lg->AddEntry(g_compound,"Compound","l");
    lg->Draw();

cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;


}
