// passer du spectre proton au spectre initial en neutrons
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include <cmath>
#include "Math/InterpolationTypes.h"
#include "Math/Interpolator.h"
#include "Math/Polynomial.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TString.h"
#include "TObjArray.h"
#include "TLegend.h"
#include <TCanvas.h>
#include <TPostScript.h>
#include <TPaveText.h>

///////////////////////////////////////////////////////////////////////////////////////////
//  on cherche la section efficace doublement differentielle pour Ein et cos(theta)cm
Double_t sigma_diff(Double_t E,Double_t cos)

{

TFile *f=new TFile("XS/dsigma_np_online.root","readonly");

//printf("E=%f   cos=%f\n",E,cos);
TH2F *hh=(TH2F *)f->Get("hh");
Int_t nx=hh->GetXaxis()->GetNbins();
Double_t x[nx];
hh->GetXaxis()->GetLowEdge(x);
//for(Int_t i=0;i<nx;i++){printf("x=%f\n",x[i]);}

Int_t ny=hh->GetYaxis()->GetNbins();
Double_t Ener[ny];
hh->GetYaxis()->GetLowEdge(Ener);

Int_t j0=0;
for(Int_t j=1;j<ny;j++){
	if(E>=Ener[j-1] &&E <Ener[j]) j0=j;   // on cherche la ligne y qui correspond à l'energie E
	if(E>=Ener[ny-1]) j0=ny-1;
}

// on cree le TGraph    dsig/dOmega = f(cos(theta)) pour energie E
Double_t y[nx];
for(Int_t i=0;i<nx;i++){
	y[i]=hh->GetBinContent(i+1,j0+1);
}
TGraph *g=new TGraph(nx,x,y);

Double_t sig=g->Eval(cos);
//printf("E=%f   J0=%i     sig=%f\n",E,j0,sig);
f->Close();
return sig;
}

//////////////////////////////////////////////////////////
void Inverse4(Int_t itel=1, TString cible="CH2")
{
gROOT->Reset();
Int_t icolor[10]={1,2,3,4,6,16,46,9,32,21};

TLegend *leg0 = new TLegend(0.612,0.564,0.86,0.876);
TLegend *leg1 = new TLegend(0.612,0.695,0.86,0.876);

Int_t irun=0;
Double_t masse=0,Nat=0;
Double_t Omega=0.01821;
Double_t dist=464;
Double_t pi=TMath::ACos(-1.);
	
TString nomHISTO;

if(cible=="CH2") {
        masse=2.37E-02;
        Nat=masse/14.*2.*6.02e23;
        irun=407;
	nomHISTO="hp_Pol";
}

else if(cible=="C") {
        masse=4.08000E-02;
        Nat=masse/12.*6.02e23;
        irun=406;
	nomHISTO="hp_C";
}

//////////////////////////////////////////////////////////
//  		flux du PPAC
//
TString file_flux="run111.root";
TFile *aa = new TFile(file_flux,"readonly");
TGraph *gPPAC=(TGraph *)aa->Get("Nspectrum");
leg0->AddEntry(gPPAC,"PPAC Yield");

Double_t angle=float(itel)*20.;
TH1F *hExp;
TString file_experiment;
//file_experiment.Form("Experiment/run%03d_tel%i_protons.root",irun,itel);
file_experiment="Sous_388_370.root";
printf("file experiment=%s\n",file_experiment.Data());
TFile *bb = new TFile(file_experiment,"readonly");

//////////////////////////////////////////////////////////////////////////////:
//      Spectre proton MEDLEY
//
//hExp=(TH1F *)bb->Get(nomHISTO);if(!hExp) {printf("l histo %s n existe pas\n",nomHISTO.Data());exit(-1.);}

hExp=(TH1F *)bb->Get("hEprot20");if(!hExp) {printf("l histo %s n existe pas\n","hEprot20");exit(-1.);}

hExp->GetXaxis()->SetTitle("Energie (MeV)");
hExp->GetYaxis()->SetTitle("cps/MeV");	
hExp->SetLineColor(icolor[2]);
hExp->SetLineWidth(2);
leg1->AddEntry(hExp,"Proton spectrum");

/////////////////////////////////////////////////////////////////////////::

Double_t co=TMath::Cos(angle/180.*pi);
Double_t En[hExp->GetNbinsX()],Fn[hExp->GetNbinsX()],Xp[hExp->GetNbinsX()+1];
Double_t Enprim[hExp->GetNbinsX()];  // energie du neutron diffusé dans La

Double_t cos_phi[hExp->GetNbinsX()];  // angle du neutron dans le CM
Double_t cos_theta[hExp->GetNbinsX()];  // angle du neutron dans le Lab

for(Int_t i=0;i<hExp->GetNbinsX();i++){
	Xp[i]=hExp->GetBinLowEdge(i)/co/co;
}
Xp[hExp->GetNbinsX()]=hExp->GetXaxis()->GetXmax()/co/co;
TH1F *hFn=new TH1F("Flux_inverse","Flux_inverse",hExp->GetNbinsX(),Xp);
hFn->GetXaxis()->SetTitle("Energie (MeV)");
hFn->GetYaxis()->SetTitle("Flux (n/cm2/MeV/#muC");	
hFn->SetLineColor(icolor[1]);
hFn->SetLineWidth(2);

printf("nombre de bins = %i\n",hExp->GetNbinsX());

printf("%i   Xp[i]=%4.2f\n",0,Xp[0]);
printf("%i   Xp[i]=%4.2f\n",hFn->GetNbinsX(),Xp[hFn->GetNbinsX()]);


for(Int_t i=0;i<hExp->GetNbinsX();i++){
	Double_t Ep=hExp->GetBinCenter(i);				// energie du proton diffusé
	Double_t Np=hExp->GetBinContent(i)*hExp->GetBinWidth(i);	// nombre de protons dans le bin
	En[i]=Ep/co/co;							// Energie du neutron incident deduit du Ep et de l'angle de diffusion

	Double_t theta_lab=180/pi*TMath::ACos(co);
	Double_t dsig=sigma_diff(En[i],theta_lab)*1.e-3*1e-24;   	// dsig/dOmega dans labo pour proton à theta_lab passage en barn puis en cm2

	Fn[i]=Np/(Nat*dsig*Omega)/hFn->GetBinWidth(i); 
	hFn->SetBinContent(i,Fn[i]);

	//if(Np>0) {printf("%i  Ep=%4.2f En=%4.2f CosTheta=%6.5f theta_lab=%4.2f dsig=%5.4g Np=%4.3g    Fn=%f\n",i,Ep,En[i],co,theta_lab,dsig,Np,Fn[i]);}
}


TH1F *hYield=(TH1F*) hFn->Clone();
hYield->GetXaxis()->SetTitle("Energie (MeV)");
hYield->GetYaxis()->SetTitle("Yield (n/MeV/sr/#muC)");	
hYield->SetLineColor(icolor[1]);
hYield->SetLineWidth(2);
hYield->Scale(dist*dist);

cout<<"dist*dist = "<<dist*dist<<endl;

leg0->AddEntry(hYield,"Medley Yield");

TGraph *gFluxN=new TGraph(hExp->GetNbinsX(),En,Fn);
gFluxN->GetXaxis()->SetTitle("Energie (MeV)");
gFluxN->GetYaxis()->SetTitle("Flux (n/cm2/#muC)");


///////////////////////////////////////////////////////////////////////////
printf("Masse de la cible m=%4.3g g\n",masse);
printf("Angle solide Omega=%4.3g sr\n",Omega);
printf("Distance convertisseur - cible =%4.1f cm\n",dist);
printf("Telescope num %i    angle=%4.2f deg\n",itel,angle);

// printf("Integrale du Yield MEDLEY=%4.3g    \n",integrale_histo(hYield,1.));
// printf("Integrale du Yield PPAC=%4.3g      \n",integrale_graph(gPPAC,1.));
// printf("Integrale de Flux Medley=%4.3g     \n",integrale_histo(hFn,1.));
// printf("Integrale de Flux Medley=%4.3g     \n",integrale_graph(gFluxN,1.));

///////////////////////////////////////////////////////////////////////////
TCanvas *c1=new TCanvas("spec","",100,70,900,750);

leg1->SetTextSize(0.058);

c1->Divide(1,3);

///////////////////     spectre   protons  ///////////////////////////
c1->cd(1);
gStyle->SetOptStat(0);gStyle->SetOptTitle(0);
gPad->SetRightMargin(0.1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.15);

hExp->GetXaxis()->SetTitleSize(0.06);
hExp->GetYaxis()->SetTitleSize(0.06);
hExp->GetXaxis()->SetLabelSize(0.06);
hExp->GetYaxis()->SetLabelSize(0.06);
hExp->Draw();
leg1->Draw("same");
TPaveText *pt=new TPaveText(0.357,0.89,0.65,0.99,"brNDC");
TString text;
text.Form("%i deg",int(angle));
pt->AddText(text);
pt->Draw("same");

///////////////////     Flux neutron at MEDLEY target    ///////////////////////////
c1->cd(2);
gPad->SetRightMargin(0.1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.2);
gStyle->SetOptStat(0);

hFn->Draw();
hFn->GetXaxis()->SetTitleSize(0.06);
hFn->GetYaxis()->SetTitleSize(0.06);
hFn->GetXaxis()->SetLabelSize(0.06);
hFn->GetYaxis()->SetLabelSize(0.06);
gFluxN->Draw("*");


///////////////////     Neutron Yield    ///////////////////////////
c1->cd(3);

gPad->SetRightMargin(0.1);
gPad->SetLeftMargin(0.15);
gPad->SetBottomMargin(0.2);
gStyle->SetOptStat(0);

hYield->GetXaxis()->SetTitleSize(0.06);
hYield->GetYaxis()->SetTitleSize(0.06);
hYield->GetXaxis()->SetLabelSize(0.06);
hYield->GetYaxis()->SetLabelSize(0.06);
//hYield->SetMaximum(2.2e10);
//hYield->GetYaxis()->SetRangeUser(0);
hYield->Draw();
//gPPAC->Draw("*");
leg0->Draw("same");


///////////////////     Write    ///////////////////////////
TFile *gg=new TFile("toto.root","update");
hExp->Write();
gg->Close();


}




