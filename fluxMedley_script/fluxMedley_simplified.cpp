
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TLatex.h"
#include <string>
#include <iostream>
#include <fstream>
#include "/home/e802/Analysis/pre_analysis/libs/useful.h"



void fluxMedley_simplified(int telN = 1 ){// telescope is an input

  //---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
  // TIMING STUFF ...

      char cur_time[128];
      time_t      t;
      struct tm*  ptm;
      t = time(NULL);
      ptm = localtime(&t);
      strftime(cur_time, 128, "%Y-%m-%d_%H:%M:%S", ptm);
      //start clock
      clock_t tStart = clock();


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
 //Filanames section

string namefilePol, namefileCarbon, namefilePPACs;

//Pol_run (run407) 
namefilePol = "/data/e802X/e802/acquisition/RootA/PIDruns/r0407_000a_T14_1-4.root";

//PPACS_run to compare
namefilePPACs = "run111.root";

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//PPACs section

//open PPACS and create Canvas... 
TFile* filePPACs = new TFile(namefilePPACs.c_str(),"READ");
TCanvas *cppacs = new TCanvas("cppacs","cppacs",500,300,800,600);
cppacs->SetLeftMargin(0.16);

//retrieve TGraphErros from file
TGraphErrors *PPACsFlux = (TGraphErrors*)filePPACs->Get("Nspectrum");

//calculate 1/D² 
double D2 = 1/pow(464.72,2);

//multiply PPACs spectrum for 1/D² so now it will be given in /cm² 
for(int t=0;t<PPACsFlux->GetN();t++)
{
   PPACsFlux->SetPointY(t, PPACsFlux->GetPointY(t)*D2);
    PPACsFlux->SetPointError(t, 0.0,PPACsFlux->GetErrorY(t)*D2);
}

//some plotting attributes 
PPACsFlux->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
PPACsFlux->SetLineColor(kRed);
PPACsFlux->SetMarkerColor(kRed);
PPACsFlux->Draw();
gPad->SetGridx();
gPad->SetGridy();

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Loading experimental run 

Long64_t nEntries[2];

TFile* filePol = new TFile(namefilePol.c_str(),"READ");

cout << "Processing " << namefilePol.c_str() << endl;
TTree* txPol = (TTree*) filePol->Get("AD");

nEntries[0] = txPol->GetEntries();
cout << "Entries [CH2]: " << nEntries[0]/ 1000 << "k " << endl;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Gamma flash and Tgamma


//Fitting Gamma flash 
double Gflash_Pol = GetGflash(txPol,Form("Medley_%d_dE2_ToF",telN),500,true); // changing it to false, it will let the fitting plot opened to check 

//print result: 
cout<<"\ngflash CH2 = "<<Gflash_Pol<<" ns."<<endl;



//gemme tgamma for given telescope, it is about 15 ns
double tgamma = provideTgama(telN);

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//Aliases definitions:

//define ToF: tgamma + Gflash  - ToF_corrected
txPol->SetAlias(Form("TOF"),Form("%f+%f - Medley_%d_dE2_ToFc_protons",tgamma,Gflash_Pol,telN));


//Total energy is siply the sum of values for the respective branches: 
txPol->SetAlias("Etot",Form("Medley_%d_dE1+Medley_%d_dE2+Medley_%d_Eres",telN,telN,telN));

//define positive events, to avoid getting problematic stuff
txPol->SetAlias("allpositive",Form("Medley_%d_dE1>0.05 && Medley_%d_dE2>=0.0 && Medley_%d_Eres>=0.0",telN,telN,telN));

//good events are both POSITIVE + has ToF bigger than tgamma
txPol->SetAlias("good_events",Form("TOF>%f && allpositive",tgamma));//goodtof

//define tof for the neutron with generated given proton
txPol->SetAlias("tof_p",Form("TOF - ToFproton(Etot)") );//it is the TOF minus the tof of a proton with Etot energy

//Define the energy of the neutron which has tof_p as tof value 
txPol->SetAlias("Enn_p",Form("En(tof_p))") );

//select protons 
txPol->SetAlias("protons",Form("PID_T%d == 1", telN));

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//TOF spectrum: 

//create canvas
TCanvas *cv2 = new TCanvas("cv2", "ToF Pol", 900,100,800,600);

//create histo with desired binning: 300 bins from about 25 ns to 325 ns, just has to be large enough to fit all the selected events 
TH1D *htof_Pol = new TH1D("htof_Pol", "time of flight - Pol", 300, floor(tgamma)+10, floor(tgamma)+310);
htof_Pol->GetXaxis()->SetTitle("ToF (ns)");


//Draw the toF spectrum for the neutrons that generate the protons we collected (we dont subtract C in this script )
txPol->Draw("tof_p>>htof_Pol","protons&&tof_p<310&&tof_p>20");


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.

//calculations

//number of protons detected for given energy:
// #(E) = Φ(E) * dσ/dΩ(E) * ΔΩ * Nat_H * Qrun

// The number of scattering nucleus Nat_H = (2mPol/MM)*NA, and dσ/dΩ(E) is tipically given in b/sr so we need a conversion factor [cm²/b] 
// note that NA = 0.602E+24 and b = 1E-24 cm², so this means we can switch NA*[cm²/b] = 0.602 

//using it we get: 

// #(E) = Φ * dσ/dΩ(E) * ΔΩ * (mpol/mm)* 1.204 * Qrun

//where mpol is the mass of CH2 and mm is its molar mass

//I will call Rm = (mpol/mm) and Cte = 1/(  ΔΩ * Rm * 1.204 * Qrun ), then 

//Φ(E) = #(E) * Cte /[dσ/dΩ(E)]

double mPol = 0.0237;//in grams
double QPol = 0.444*1E6;
double Rm = (mPol/14.0);//I consider mm[CH2] = 14g/mol
double omega_telescope = 0.0182538;// ± 0.0000498 sr,  calculated with SACALC3v14


double Cte = 1/(omega_telescope * Rm *1.204 * QPol);


//Now I create a canvas to plot the neutron energy measured spectrum (still no cross section involved)
TCanvas *Energy_canvas = new TCanvas("Energy_canvas", "neutron energy - Form Subtracted ToF", 130,130,800,600);

//The neutronFromToF will conver an histo in ToF to neutron energy, conserving the entries 
TH1D *hEn = NeutronFromToF(htof_Pol);
hEn->SetNameTitle("hEn","hEn");
hEn->GetXaxis()->SetRangeUser(0,50);
//Draw it so we can see
hEn->Draw();


//Now I call the XS fetcher, that will load the cross section into a TGraph so its easier to work with 
TGraph *xs = XS_fecther(telN,1,0,"endf_B-VI.8_20deg.40deg.60deg.80deg.csv");


TCanvas *cFlux = new TCanvas("cFlux", "Flux", 630,330,800,600);
TH1D *hflux = (TH1D*)hEn->Clone();
hflux->SetNameTitle("hflux","hflux");
cout<<"flux stats: "<<endl;

//I print some info to check it is not calculating stuff wrong:

//number of bins:
cout<<"bins: "<<hEn->GetNbinsX();
//covering energies:
cout<<", from "<<hEn->GetBinLowEdge(1)<<" MeV to "<<hEn->GetBinLowEdge(hEn->GetNbinsX()+1)<<"MeV"<<endl;


//uncomment line below for debugging
// cout<<Form("bin # |  center (MeV) |   value  |  xs(b/sr)")<<endl;

//some variables 
double xscalc, bw, value ;


float transformation_factor = 2*TMath::Pi()*sin(telN*20*TMath::Pi()/180);
cout<<"transformation_factor = "<<transformation_factor<<endl;
//for each bin in the neutron energy measure histogram:
for(int b=1; b<= hEn->GetNbinsX(); b++){

//uncomment line below for debugging
//     //cout<<Form("bin %03d |  %02.01f |   %e  |  %02.04f", b,hSubtraction->GetBinCenter(b),hSubtraction->GetBinContent(b), xs->Eval(hSubtraction->GetBinCenter(b)))<<endl;

    //Get the corresponding XS: 
    xscalc = transformation_factor*xs->Eval(hEn->GetBinCenter(b));

    //get the bin-width from the neutrons' energy histo
    bw = hEn->GetBinWidth(b);

    //calculate: Φ(E) = Cte * #(E) /[dσ/dΩ(E)]
    //here I divided to binwidth to be comparible with the other results
    value = Cte*hEn->GetBinContent(b)/(xscalc*bw);

    //fill the value in hflux 
    hflux->SetBinContent(b,value);

}


//Axis labels and drawing 
 hflux->GetYaxis()->SetTitle("neutrons/#muC/cm^{2}/1-MeV");
 hflux->GetYaxis()->SetMaxDigits(2);
 hflux->Draw();

//plot PPACs reference
 PPACsFlux->Draw("same");
gPad->SetGridx();
gPad->SetGridy();

cout<<"==================================================================================="<<endl;
cout<<"|                                     SUMMARY                                     |"<<endl;
cout<<"-----------------------------------------------------------------------------------"<<endl;
cout<<"..................................................................................."<<endl;
std::scientific;
cout<<"mPol        :"<<mPol<<" g"<<endl;
cout<<"Cte  :"<<Cte<<endl;
cout<<"QPol        :"<<QPol<<" µC"<<endl;

//uncomment this to take a look at the XS:
//Checking XS
// TGraph *xs = XS_fecther(telN,0,0,"/home/e802/Analysis/pre_analysis/cross_secs/endf_B-VI.8_20deg.40deg.60deg.80deg.csv");
// TCanvas *XScv = new TCanvas("XScv", "cross section", 580,150,800,600);
// xs->GetXaxis()->SetTitle("En (MeV)");
// xs->GetXaxis()->SetRangeUser(1,40);
// xs->Draw();



cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;



}
