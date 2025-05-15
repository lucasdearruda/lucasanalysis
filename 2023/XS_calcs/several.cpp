#include "plotDDX.cpp" // for the energy calibration


//IBM palette
const Int_t nColors = 5;
Int_t colors[nColors] = {
        TColor::GetColor("#648FFF"),  // Blue
        TColor::GetColor("#785EF0"),  // Purple
        TColor::GetColor("#DC267F"),  // Salmon
        TColor::GetColor("#FE6100"),  // Orange
        TColor::GetColor("#FFB000"),  // Gold
};

void several(){


int ang;
int E; 

char part[]={'p','d'};
char ch;

TGraph *mygr; 
TH1D *myh;

Float_t maxval = 0;
//for(int i=0;i<2;i++){
for(int i=0;i<2;i++){
    cout<< "plotting "<<part[i]<<endl;
    for(ang = 20; ang<=80; ang=ang+20){
        for(E = 20; E<=35;E=E+5){
            
            TCanvas *c = plotDDX(E,part[i],ang,1,2);

            myh = (TH1D*)gROOT->FindObject("hddx");


            plotMeAxs(Form("/mnt/medley/LucasAnalysis/talys_calcs/Fe56/best/%cddxE00%02d.000A0%02d.0.deg",part[i],E,ang));
            mygr = (TGraph*)gROOT->FindObject("gr");
            maxval = TMath::MaxElement(mygr->GetN(),mygr->GetY());
            
            for(int i=1;i<=myh->GetNbinsX();i++){myh->SetBinContent(i,myh->GetBinWidth(i)*myh->GetBinContent(i));}
            
            if(myh->GetMaximum() < maxval) myh->GetYaxis()->SetRangeUser(0,1.20*maxval);
            makeTL(myh, mygr,Form("^{nat}Fe(n,X%c) E_{NN} = %02d MeV, %02ddeg",part[i],E,ang), "Medley","Talys 2.0 (^{56}Fe)");
            c->SaveAs(Form("results/F/Fe_n_X%c_%02ddeg_%02dMeV.root",part[i],ang,E));
            c->SaveAs(Form("results/F/png/Fe_n_X%c_%02ddeg_%02dMeV.png",part[i],ang,E));

            // std::cout << "Press Enter to continue..." << std::endl;
            // std::cin.ignore();  // opcional pra limpar buffer
            // std::cin.get();     // espera enter
        }
    }
}
return;

}


Float_t angularXS(int E = 20, Float_t *xs_sigma =nullptr){

    TGraphErrors *xs_ang;
    TGraph *xs_talys, *xs_min, *xs_max; 

    int ang;
     // I fixed one energy to try, first. 
    
    //char part[]={'p','d'};
    
    char part[]={'p'};
    
    
    //TGraph *mygr; 
    TH1D *myh[4];
    Float_t sigmaxs;


    Float_t maxval = 0;
    Float_t xs_value;
    Int_t index; 




    TCanvas *histosCv = new TCanvas("histosCv","histosCv",800,600);
    TLegend *tl = new TLegend(0.4,0.53,0.88,0.86);
    cout<< "plotting "<<part[0]<<endl;
    
    xs_ang = new TGraphErrors();
    
    
    
    xs_talys = new TGraph();
    
    xs_min = new TGraph();
    xs_max = new TGraph();
    xs_min->SetName("xs_min");
    xs_max->SetName("xs_max");

    Int_t pn = 0;
    Float_t IntHisto;
    Float_t IntTalys;


    TGraph *talysgr; // Talys TGraph 
    for(ang = 20; ang<=80; ang=ang+20){

        TCanvas *c = plotDDX(E,part[0],ang,2);
        talysgr = plotMeAxs(Form("/mnt/medley/LucasAnalysis/talys_calcs/Fe56/best/%cddxE00%02d.000A0%02d.0.deg",part[0],E,ang));
        IntTalys = TrapezoidalIntegration(talysgr,0,40);

        cout<<"Integral Talys: "<<IntTalys<<endl;

        index = ang/20 - 1;
        cout<<"index = "<<index<<endl;

        myh[index] = (TH1D*)gROOT->FindObject("hddx");
        myh[index]->SetName(Form("hddx_%02d",ang));
        myh[index]->SetFillStyle(0);
        myh[index]->SetLineWidth(2);
        myh[index]->SetLineColor( colors[index]);
        tl->AddEntry(myh[index],Form("%02d deg",ang),"l");
        
    
        IntHisto = myh[index]->Integral("width");
        cout<<"Angle:"<<ang<<", Integral = "<<IntHisto<<endl;
        xs_ang->AddPoint(cos(TMath::Pi()*ang/180),IntHisto);
        xs_min->AddPoint(cos(TMath::Pi()*ang/180),IntHisto - sqrt(IntHisto));
        xs_max->AddPoint(cos(TMath::Pi()*ang/180),IntHisto + sqrt(IntHisto));
        xs_ang->SetPointError(pn,0.0,sqrt(IntHisto));


        xs_talys->AddPoint(cos(TMath::Pi()*ang/180),IntTalys);
        pn++;

        histosCv->cd();
        if(!index){
            myh[index]->Draw();
        }else{
            myh[index]->Draw("same");
        }

    }
    histosCv->cd();
    tl->Draw();    
    TCanvas *gg = new TCanvas("gg","gg",800,600);
    xs_ang->SetName("xs_ang");
    xs_ang->SetTitle(Form("^{nat}Fe(n,X%c) E_{NN} = %02d MeV",part[0],E));
    xs_ang->GetXaxis()->SetTitle("cos(#theta)");
    xs_ang->GetYaxis()->SetTitle("d#sigma/d#Omega (mb/sr)");
    xs_ang->SetMarkerStyle(20);
    xs_ang->SetMarkerColor(kRed);
    xs_ang->SetMarkerSize(1.5);
    xs_ang->SetLineColor(kRed);
    xs_ang->SetLineWidth(2);
    
    
    float maxval_exp = TMath::MaxElement(xs_ang->GetN(),xs_ang->GetY());
    //float maxval_talys = TMath::MaxElement(xs_ang->GetN(),xs_ang->GetY());      
    //if(myh->GetMaximum() < maxval) myh->GetYaxis()->SetRangeUser(0,1.20*maxval);
    xs_ang->GetYaxis()->SetRangeUser(0,1.2*maxval_exp);
    xs_ang->Draw();

    xs_min->SetLineColor(kBlue);
    xs_max->SetLineColor(kBlue);

    xs_min->Draw("same");
    xs_max->Draw("same");


    //*xs_sigma = TrapezoidalIntegration(xs_max,0,1) - TrapezoidalIntegration(xs_min,0,1);

    if (xs_sigma != nullptr) {
        *xs_sigma = TrapezoidalIntegration(xs_max,0,1) - TrapezoidalIntegration(xs_min,0,1);;
    }

    xs_value = 2*TrapezoidalIntegration(xs_ang,0,1);
    cout<<"Production XS: "<<xs_value<<" "<<sigmaxs<< " mb."<<endl;
    
    xs_talys->Draw("same PL");

    return xs_value;
    
}

void plotXSprod(){
    TCanvas *c = new TCanvas("c","c",800,600);
    
    TGraphErrors *xsprod = new TGraphErrors();
    xsprod->SetName("xsprod");
    xsprod->SetTitle("XS production");
    xsprod->GetXaxis()->SetTitle("E_{NN} (MeV)");
    xsprod->GetYaxis()->SetTitle("#sigma_{prod} (mb)");
    Float_t sigma_xs;
    Float_t E;
    Int_t npoint =0;
    for(E = 10; E<=35;E=E+5){
        xsprod->AddPoint(E,angularXS(E, &sigma_xs));
        xsprod->SetPointError(npoint++,0.0,sigma_xs);
        cout<<"E: "<<E<<endl;
    }

    for(int i=0;i<xsprod->GetN();i++){
        cout<<Form("E = %.1f, xs = %.1f ( %.1f)",xsprod->GetPointX(i),xsprod->GetPointY(i),xsprod->GetErrorY(i))<<endl;
    }
    

    
    TTree *tx = new TTree("xs","xs");
    //tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MPreeq/F:Compound/F:PreeqR/F:BKr/F:Str/F:KO/F:BrU/F");
    //##     E-out           xs           Direct     Preequilibrium Multiple_preeq    Compound
    tx->ReadFile("/mnt/medley/LucasAnalysis/talys_calcs/Fe56/best/pprod.tot", "Eprod/F:Total/F:Something");
    Float_t xtalys_pp, Etalys_pp;
    Long64_t nEntries = tx->GetEntries();
    tx->SetBranchAddress("Eprod", &Etalys_pp);
    tx->SetBranchAddress("Total", &xtalys_pp);


    TGraph *grT = new TGraph();
    for(Long64_t i = 0; i < nEntries; ++i) {
        tx->GetEntry(i);

        //cout<<Eprod<<" "<<total<<endl;
        grT->AddPoint(Etalys_pp, xtalys_pp);
        // Plotting code here
        // For example:
        // hist->Fill(Eplot, total);
    }
    grT->SetName("grT");
    grT->SetLineWidth(2);
    grT->SetLineColor(kRed);


    xsprod->Draw();
    grT->Draw("same");
    return;

}