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

#include "TStopwatch.h"
#include <time.h>

//function to get the first non zero entry in given branch 
Long64_t GetFirstNonZero(TTree *G, const char* branchName = "Medley_1_SI_DE2TS"){
    ULong64_t value, nE;
    value = 0;
    G->SetBranchAddress(branchName,&value);
    nE=G->GetEntries();
    for(Int_t i=0; i<nE;i++){
    	G->GetEntry(i);
	if(value){
	//cout<<"evt "<<i<<endl;
	break;

	}

	}
    return value;
}


Long64_t GetLasttNonZero(TTree *G, const char* branchName = "Medley_1_SI_DE2TS"){
    ULong64_t value, nE;
    value = 0;
    G->SetBranchAddress(branchName,&value);
    nE=G->GetEntries();
    for(Int_t i=nE-1; i>=0;i--){
    	G->GetEntry(i);
	if(value){
	//cout<<"evt "<<i<<endl;
	break;

	}

	}
    return value;
}

//function to get the beam live time 
TH1D *beamlivetime(TTree * tx, string branchname="Medley_1_SI_DE2TS",float *time_sec=NULL, float *avc = NULL, float *totaltime = NULL,float threshold = 20,string histoname="hc", bool drawit =false){

//Float_t firstvalue = GetFirstNonZero(tx)/1.e8;
Long64_t firstvalue = GetFirstNonZero(tx,branchname.c_str())/1.e8;
Long64_t lastvalue = GetLasttNonZero(tx,branchname.c_str())/1.e8;

Int_t nbinstime = round(lastvalue - firstvalue) +1;
if(totaltime) *totaltime = lastvalue - firstvalue;

cout << ">>>> 'beamlivetime' function report: (considering TS branches)-------"<<endl;
cout<<".\n. first value found: "<<firstvalue<<" seconds."<<endl;
cout<<". last value found: "<<lastvalue<<" seconds."<<endl;
cout<<". bins TS histo: "<<nbinstime<<" seconds.\n."<<endl;
cout << ">>>>-----------------------------------------------------------------"<<endl;


  TH1D *hc=new TH1D(histoname.c_str(),Form("counts over time: %s",histoname.c_str()),nbinstime,0,nbinstime);

    //cout<<"Maximum TT = "<<tx->GetMaximum
if(drawit){
    hc->GetXaxis()->SetTitle("time (s)");
    tx->Draw(Form("%s/1.e8 - %lld>>%s",branchname.c_str(),firstvalue,histoname.c_str()),"");	// Divide by 1.e8 to get it in s, because time stamp happens every 10 ns;
    TF1 *f1 = new TF1("f1","pol0",0,nbinstime);
    f1->SetParameter(0,threshold);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->Draw("same");
}else{
    tx->Draw(Form("%s/1.e8 - %lld>>%s",branchname.c_str(),firstvalue,histoname.c_str()),"", "goff");	// Divide by 1.e8 to get it in s, because time stamp happens every 10 ns;
}
  
   Int_t nbins=hc->GetNbinsX();
  //cout << "nbins= " << nbins << " limits: " << h1->GetXaxis()->GetXmin() << " " << h1->GetXaxis()->GetXmax() << endl;
  Int_t realbeamtime=0;

  Float_t beamthreshold=threshold; // For C experiment: beamthreshold=40.; for Cr experiment: beamthreshold=10.;
  Float_t cc=0.;
  Float_t CR = 0.;
  for(Int_t i=0;i<nbins;i++){
   cc=hc->GetBinContent(i+1);

   if(cc>=beamthreshold)
   {
 //  cout << "bin " << i << " cc= " << cc << endl;
     realbeamtime++;
     CR += cc;
   }

  }

  CR = 1.*CR/realbeamtime; //in secs
  cout << "\n\n>> Livetime summary -------------------------------------------------"<<endl;

  cout << ". Real beam time: " << realbeamtime <<  "(s)  = " << realbeamtime/3600. << " (h). Avg CR = " << CR<< "/s" << endl;
  if(time_sec)*time_sec = realbeamtime;
  if(avc)*avc = CR;
    cout << ".\n---------------------------------------------------------------------"<<endl;

return hc;



}



Int_t Ta;


TChain *t1=NULL;

TChain *openRuns(Int_t nb1 = 41,Int_t nb2 = 41){
TStopwatch timer;




if(!nb1 ||!nb2 ){
	cerr<<"You must provide nb1 and nb2... Please, try again."<<endl;
 	return nullptr;  
}
  cout << "Runs from " << nb1 << " to " << nb2 <<endl;
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1001111);
  gStyle->SetFillStyle(1000); 
  
  char name[200];
   
  //cout << t1 << endl;
  if(!t1)
    t1 =new TChain("AD");
  for(Int_t i=nb1;i<=nb2;i++)
    {
      Bool_t fExist = true;
        for(Int_t j=0;j<999&&fExist;j++)
        {

            sprintf(name,"/mnt/medley/RootA_2024/r%04d_%03da.root",i,j);
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
                t1->Add(name);
                Ta = (Int_t) (t1->GetEntries());
                cout << "Entries " << Ta/1000 << "k "  << endl;
            }
        }
    }
  

  beamlivetime(t1);
  timer.Print();

 return t1;
}
