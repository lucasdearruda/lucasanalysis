#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include <string>
#include <iostream>
#include <fstream>
#include "useful.h"

#include <time.h>

void myAna(int nruni, int nrunf){

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
    //opening files 
    char name[1000];
    TChain *tx = NULL;

    if (!tx)
        tx = new TChain("M");

    vector<Int_t> entries;

    for (int i = nruni; i <= nrunf; i++)
    {
      Bool_t fExist = true;
      
      ifstream mfile;
      sprintf(name, "/home/e802/Analysis/LucasAnalysis/reducedData/%03d.root", i);
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
          cout <<"i = "<< i<< " Entries: " << entries.back() / 1000 << "k " << endl;
      }
  
    }
    //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
    
    
TCanvas *c1=NULL;
TTreeViewer *v1=NULL;
//   //  new TBrowser();
  if(!v1)
    v1 =new TTreeViewer(tx);
  else
    {
      delete v1;
      v1 =new TTreeViewer(tx);
    }

 if(!c1)
   c1 = new TCanvas("c1","myAnalysis",400,10,900,700);
 c1->ToggleEventStatus();
 
 gPad->SetGridx();
 gPad->SetGridy();
 //gPad->ToggleToolBar();
 
 gPad->SetLogy(0);
 gPad->SetLogz(1);
    


cout <<"\nexecution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;


}
