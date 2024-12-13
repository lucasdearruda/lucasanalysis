#include <stdlib.h>
#include <stdio.h>
//#include "stream.h"
#include "../useful.h"

//TChain *t1=NULL;
TCanvas *c1=NULL;
TTreeViewer *v1=NULL;
Int_t Ta;

TChain *myAna(int nb1= 0, int nb2=0)
{

   TChain *t1=NULL;

  if(!nb1 ||!nb2 ){
	cerr<<"You must provide nb1 and nb2... Please, try again."<<endl;
 	return NULL;  
}
  cout << "Runs from " << nb1 << " to " << nb2 <<endl;
  gStyle->SetPalette(55);
  gStyle->SetOptStat(1001111);
  gStyle->SetFillStyle(1000); 
  
  char name[200];
   
  cout << t1 << endl;
  if(!t1)
    t1 =new TChain("M");
  for(Int_t i=nb1;i<=nb2;i++)
    {
        Bool_t fExist = true;
            
        sprintf(name,"/media/dearruda/Elements/LucasAnalysis/2023/reduced/%03dv5.root",i);
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
  

  //  new TBrowser();
  if(!v1)
    v1 =new TTreeViewer(t1);
  else
    {
      delete v1;
      v1 =new TTreeViewer(t1);
    }

 if(!c1)
   c1 = new TCanvas("c1","myAnalysis",400,10,522,400);
 c1->ToggleEventStatus();
 
 gPad->SetGridx();
 gPad->SetGridy();
 //gPad->ToggleToolBar();
 
 gPad->SetLogy(0);
 gPad->SetLogz(1);
 
 return t1;
}   
