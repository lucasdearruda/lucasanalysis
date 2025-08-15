//version 22.2025.08.14.001
Double_t mc;
Double_t M;
Double_t omega_tel = 0.010 ;//sr ---> after I will add fine adjustments to i
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Int_t A, Z;
Double_t L = 504.4; //cm
Double_t A_tgt = TMath::Pi()*pow(2.5/2,2); //cm^2
Double_t cm2_to_barn = 1e+24 ;//cm²

Float_t binWidth = 0.5; //MeV

using namespace std;




unsigned long long version(bool show = true) {
    string versionStr = "2025.08.08.001";

    // Remove os pontos
    versionStr.erase(remove(versionStr.begin(), versionStr.end(), '.'), versionStr.end());

    // Converte para número inteiro grande
    unsigned long long versionNum = stoull(versionStr);

    if (show) {
        cout << "// Set of functions developed by Lucas de Arruda //" << endl;
        cout << "version2025.08.07.001" << endl << endl;
        cout << "2025-08-07 :: Adding several functions:" << endl;
        cout << "  -> getAZ" << endl;
        cout << "  -> Attribute_Target" << endl;
        cout << "  -> getNflux" << endl;
        cout << "  -> getFluxHisto" << endl;
        cout << "  -> getCurrentTime" << endl;
        cout << "  -> thicknessFirstTel" << endl;
        cout << "  -> giveMeTheAngle2, from useful.h, version2025.07.22.001" << endl;
        cout << "  -> GetPfactor: to calculate the dead fraction of particles within the target" << endl;
        cout << "  -> splitPath: splits a full file path into directory and filename." << endl;
        cout << "  -> Fcorr: computes F(E) correction factor from input and lost particles, and optionally draws the graph." << endl;
        cout << "  -> defineParticle: sets Z and A based on particle symbol (p, d, t, h, a) and prints its name." << endl;
        cout << "  -> correctSpec: return the corrected spectrum after applying the TTC" << endl;
        cout << "  -> particleName: return the particlename from the particle char" << endl;

        cout << endl;

        // Imprimir variáveis
        cout << "mc = " << mc << " g" << endl;
        cout << "M = " << M << " g/mol" << endl;
        cout << "omega_tel = " << omega_tel << " sr" << endl;
        cout << "NA = " << NA << " mol^-1" << endl;
        cout << "Nc = " << Nc << " atoms/cm^3" << endl;
        cout << "A = " << A << endl;
        cout << "Z = " << Z << endl;
        cout << "binWidth = " << binWidth << " MeV" << endl;
        cout << "L = " << M << " cm" << endl;
        cout << "A_tgt = " << A_tgt << " cm^2" << endl;
        cout << "cm2_to_barn = " <<  cm2_to_barn<< " b/cm^2" << endl;
    }

    return versionNum;
}


void getAZ(char particle) {
    switch (particle) {
        case 'p': Z = 1; A = 1; break;  // próton
        case 'd': Z = 1; A = 2; break;  // deutério
        case 't': Z = 1; A = 3; break;  // trítio
        case 'h': Z = 2; A = 3; break;  // hélio-3
        case 'a': Z = 2; A = 4; break;  // alfa (hélio-4)
        default:  Z = 0; A = 0; break;  // inválido
    }
}


void Attribute_Target(string target_name = "MedleyCarbon"){
    //default values:
    mc = 0.0;//g 
    M = 1; //g/mol
    if(target_name == "Fe_thick_Medley" || target_name == "Fe_thick"){
        mc = 0.0851;//g 
        M = 55.845; //g/mol
        Nc = (mc/M) * NA;
    }else if(target_name == "Fe_thin_Medley" || target_name == "Fe_thin"){
        mc = 0.0182;//g 
        M = 55.845; //g/mol
        Nc = (mc/M) * NA;
    }else if(target_name == "MedleyCarbon" || target_name == "C"|| target_name == "CMedley"){
        mc = 0.0408;//g 
        M =  12.01070; //g/mol
        Nc = (mc/M) * NA;
    }else if(target_name == "Polyethylene" || target_name == "CH2"){
        mc = 0.0237;//g 
        M = 14.026; //g/mol
        Nc = (mc/M) * NA;
    }else if(target_name == ""){
        cerr<<"Error: Target name not specified."<<endl;
        cout<<">>------------------------------------------------------- "<<endl;
        cout<<"List of targets:"<<endl;
        cout<<" * Fe_thick_Medley / Fe_thick \n * Fe_thin_Medley / Fe_thin\n * MedleyCarbon / C / CMedley \n * Polyethylene /  CH2"<<endl;
        return;
    }

        cout<<">>=.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.== "<<endl;
        cout<<target_name<<endl;
        cout<<">>------------------------------------------------------- "<<endl;
        cout<<"m = "<<mc<<" g."<<endl;
        cout<<"M = "<<M<<" g/mol."<<endl;
        cout<<"Nc = "<<Nc<<" ."<<endl;
        cout<<">>-.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.-- "<<endl;
}

Float_t getNflux(
    Float_t Ea = 10,
    Float_t Eb = 12, 
    string fluxfilename = "/mnt/medley/LucasAnalysis/2022/nflux22.root", 
    bool verbose = false,
    bool drawIt = false
    ){

    if(verbose){
        cout<<"Loading nflux from: "<<fluxfilename<<endl;
    }
    
    TFile *ff = new TFile(fluxfilename.c_str(), "READ");
    TH1D *nflux = (TH1D*)ff->Get("nflux");
    
    if(drawIt){
        TCanvas * cflux = new TCanvas("cflux","cflux");
        nflux->Draw();
    }
    //nflux->Draw();

    Float_t nflux_En = nflux->Integral(nflux->FindBin(Ea),nflux->FindBin(Eb),"w")/(Eb-Ea);
    if(verbose) cout << "nflux_En ("<<(Ea+Eb)/2.0<<" MeV): " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;

    return nflux_En;
} 

TH1D *getFluxHisto(Float_t Ea = 10,Float_t Eb = 12, string fluxfilename = "/mnt/medley/LucasAnalysis/2022/nflux22.root", bool verbose = false, bool drawIt = false){

    if(verbose){
        cout<<"Loading nflux from: "<<fluxfilename<<endl;
    }
    
    TFile *ff = new TFile(fluxfilename.c_str(), "READ");
    TH1D *nflux = (TH1D*)ff->Get("nflux");
    TCanvas * cflux = new TCanvas("cflux","cflux");
    if(drawIt){
        nflux->Draw();
    }
    //nflux->Draw();

    Float_t nflux_En = nflux->Integral(nflux->FindBin(Ea),nflux->FindBin(Eb),"w")/(Eb-Ea);
    if(verbose) cout << "nflux_En ("<<(Ea+Eb)/2.0<<" MeV): " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;

    return nflux;
} 


std::string getCurrentTime() {
    char cur_time[128];
    time_t t = time(NULL);
    struct tm* ptm = localtime(&t);

    // Formata o tempo no formato desejado: "YYYY-MM-DD_HH:MM:SS"
    strftime(cur_time, sizeof(cur_time), "%Y-%m-%d_%H:%M:%S", ptm);

    return std::string(cur_time);
}


Int_t pCode(char particle = 'p'){

    if(particle == 'p'){
        return 1;
    }else if(particle == 'd'){
        return 2;
    }else if(particle == 't'){
        return 3;
    }else if(particle == 'h'){
        return 4;   
    }else if(particle == 'a'){
        return 5;
    }else{
        cerr<<"Error: Unknown particle code."<<endl;
        return 0;
    }   
    return 0;
}

Float_t thicknessFirstTel(int telN = 1){
    switch (telN)
    {
    case 1:
        return 53.4; 
        break;
    case 2:
        return 50; 
        break;
    case 3:
        return 61.2;
        break;
    case 4:
        return 50;
        break;
    case 5://ONLY CONSIDERING 'CLASSIC' TELESCOPES
        return 50;
        break;
    case 6:
        return 61.2;
        break;
    case 7:
        return 50;
        break;
    case 8:
        return 53.4;
        break;

    default:
        return 53.4;
        break;
    }

}


std::vector<Double_t> giveMeTheAngle2(
    Long64_t N = 1e3,
    Double_t angle_telescope = 20.0, 
    Double_t D_mm = 149.7, 
    Double_t target_radius_mm = 25.0 / 2, 
    Double_t siA_radius_mm = sqrt(450.0 / TMath::Pi()), 
    TVector3 target_n = TVector3(cos(TMath::Pi() / 4), 0., cos(TMath::Pi() / 4))
) {
    TRandom2 *rn = new TRandom2();
    std::time_t seed = std::time(nullptr);
    rn->SetSeed(seed);

    std::vector<Double_t> angles; 

    Double_t xt, yt, xd, yd;

    TVector3 Point_detector, Point_target, distance;
    TVector3 displacement(0.,0.,D_mm);

    if(angle_telescope<=90){
        displacement.RotateY(-(angle_telescope/180.)*TMath::Pi());
    }else{
        displacement.RotateY((angle_telescope/180.)*TMath::Pi());
    }
        

    for(Int_t i=0;i<N;i++){
        Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
        Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_target.X(),2) + pow(Point_target.Y(),2) >= pow(target_radius_mm,2)){
            Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
            Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
        Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_detector.X(),2) + pow(Point_detector.Y(),2) >= pow(siA_radius_mm,2)){
            Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
            Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_target.SetZ(0);
        Point_detector.SetZ(0);


        Point_target.RotateY(TMath::Pi()/4); //rotate point in target because the target has 45 deg inclination! 
        Point_detector.RotateY(-TMath::Pi()/9);//rotate the coordinate to match the detector position
    
        Point_detector += displacement;

        distance = Point_detector - Point_target;
        //here we have coordinates in target (xt,yt) and detector (xd,xy)! 
        //Now we will calculate the angle

        angles.push_back(distance.Angle(displacement));

    }

    return angles;

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Calculate corrected spectrum
//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.


#include "TVector3.h"
#include "TMatrixT.h"


#include <sys/stat.h>


TGraph * GetPfactor(
    string filename = "eloss_p_20.0deg_050.0um.root", 
    string path = "/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick_Medley/v7/",
    Float_t nbins = 800,
    bool drawit = false
    ){
    

    TFile *f = new TFile(Form("%s/%s",path.c_str(),filename.c_str()),"READ");
    TTree *t = (TTree*)f->Get("SIM");
    Float_t Pfactor;

    TH1D *hE = new TH1D("hE","hEnergy",nbins,0,40);
    TH1D *hEl = new TH1D("hEl","hEnergyLostParticles",nbins,0,40);
    
    
    if(drawit){
        TCanvas *cc = new TCanvas("cc","cc",1600,600);
        cc->Divide(2,1);
        cc->cd(1);

        t->Draw("E>>hE");
        hE->GetYaxis()->SetRangeUser(0,1.1*hE->GetMaximum());
        t->Draw("E>>hEl","detected<1", "same");
        cc->cd(2);

    }else{
        t->Draw("E>>hE","","goff");
        t->Draw("E>>hEl","detected<1", "goff");

    }
    

    t->Draw("E>>hE","","goff");
    hE->GetYaxis()->SetRangeUser(0,1.1*hE->GetMaximum());
    t->Draw("E>>hEl","detected<1", "goff");


    TGraph *gr = new TGraph();

    for (int i = 0; i < nbins; i++) {
        if(hE->GetBinContent(i) > hEl->GetBinContent(i)){
            Pfactor = hE->GetBinContent(i)/(hE->GetBinContent(i) - hEl->GetBinContent(i));
            gr->AddPoint(hE->GetBinCenter(i),Pfactor);
            //if(drawit) cout<<"Ei = "<<hE->GetBinCenter(i)<<"MeV, Pfactor = "<<Pfactor<<endl;
        }
    }

    gr->SetMarkerStyle(20);
    if(drawit){
        gr->Draw("ALP");
    }
    

return gr;
}



string ganil_folder= "/home/pi/ganil/"; 



void splitPath(const std::string& fullPath, std::string& directory, std::string& filename) {
    size_t found = fullPath.find_last_of("/\\"); // Encontra última barra
    directory = fullPath.substr(0, found + 1);  // Inclui a barra final
    filename = fullPath.substr(found + 1);      // Parte após a barra
}

TGraph *Fcorr(TTree *S, Int_t nbins = 400, bool drawit = true){

TH1D *hprod = new TH1D("hprod","hprod",nbins,0,40);
TH1D *hlost = new TH1D("hlost","hlost",nbins,0,40);

S->Draw("E>>hprod","","goff");
S->Draw("E>>hlost","!transmitted","goff");

TGraph *Ffunc = new TGraph();
Double_t F; 
for(int i = 1; i <= nbins; i++){
    F = hprod->GetBinContent(i)/(hprod->GetBinContent(i) - hlost->GetBinContent(i));
    //F = F/hprod->GetBinWidth(i);
    Ffunc->AddPoint(hprod->GetBinCenter(i),F);
}
Ffunc->GetXaxis()->SetTitle("Energy (MeV)");
Ffunc->GetYaxis()->SetTitle("F factor ");

if(drawit){
    TCanvas *corrFcv = new TCanvas("corrFactorCv","correction Factor -- lost particles",150,150,800,678);
    corrFcv->SetLeftMargin(0.14);

    Ffunc->Draw("ALP");
    gPad->SetGridx();
    gPad->SetGridy();
}

return Ffunc;

}


void defineParticle(char particle, Int_t* Zfis, Int_t* Afis){


    switch(particle){
        case 'p':
            *Zfis =1;
            *Afis =1;
            cout<<"|   particle    |   PROTON"<<endl;
            break;
        case 'd':
            *Zfis =1;
            *Afis =2;
            cout<<"|   particle    |   DEUTERON"<<endl;
            break;	
        case 't':
            *Zfis =1;
            *Afis =3;
            cout<<"|   particle    |   TRITON"<<endl;
            break;
        case 'h':
            *Zfis =2;
            *Afis =3;
            cout<<"|   particle    |   He-3"<<endl;
            break;
        case 'a':
            *Zfis =2;
            *Afis =4;
            cout<<"|   particle    |   ALPHA"<<endl;
            break;
        default: 
            *Zfis =1;
            *Afis =1;	
            cout<<"PROTON!!"<<endl;
    }


}


TH1D* correctSpec(TH1D *expSpec = nullptr,
                        bool correctDead = true,
                        char particle = 'p',
                        Float_t ang = 20,
                        string target_mat = "MedleyCarbon",
                        Double_t th = 75,
                        bool saveCanvas = false,
                        string pSpecFileName = "myspec.root"
                        ){

TStopwatch timer;


string pElossFileName = "/home/pi/ganil/kalscripts/eloss/results/UniformZ/"+target_mat+"/v7/"+Form("eloss_%c_%.1fdeg_%05.1fum.root",particle,(ang < 100) ? ang : 180 - ang,th); 

//defining the particle
Int_t Afis, Zfis;
defineParticle(particle, &Zfis, &Afis);



TChain *ElTree = new TChain("SIM");

ElTree->Add(pElossFileName.c_str());

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//defining the response function:

//First of all, the number of bins will depend on the spectrum I want to correct:
Int_t nbins = expSpec->GetNbinsX();
cout<<"Definig the number of bins for the response function equal to NbinsX for input spec = "<<nbins<<endl;


//Now, I will create the response function
TH2D *h= new TH2D("h","h",nbins, expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1),nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));


//Here I will draw the response function
ElTree->Draw("E:Erem>>h","","goff");
gStyle->SetOptStat("e");
h->SetTitle(Form("protons leaving %s target, %4.1f#mum",target_mat.c_str(),th));
h->GetXaxis()->SetTitle("Measured proton energy (MeV)");
h->GetYaxis()->SetTitle("Initial proton energy (MeV)");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//creating the projections over Y axis
TH1D *projY[nbins+1];
for (int b = 1; b <=nbins; b++)
{
    projY[b] =  h->ProjectionY(Form("hproj%d",b),b,b);
}

///////////////////////////////////////////////////////////////////
//get new energy
//First: create the random generator 
TRandom2 *rand = new TRandom2();
std::time_t seed = std::time(nullptr);
rand->SetSeed(seed);

//Create a canvas to temporarly store the corrected spectrum
TCanvas *corrCv = new TCanvas("corrCv","corrCv",50,50,1598,678);
corrCv->SetLeftMargin(0.14);
corrCv->Divide(2,1);



int cvID = 1;
corrCv->cd(cvID);

gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
gStyle->SetTextSize(0.04);
gStyle->SetTitleYSize(0.04);

gPad->SetGridx();
gPad->SetGridy();
// Configurações do canvas
corrCv->cd(cvID)->SetRightMargin(0.14);
corrCv->cd(cvID)->SetTopMargin(0.08);
corrCv->cd(cvID)->SetLeftMargin(0.14);
corrCv->cd(cvID)->SetBottomMargin(0.14);

h->Draw("colz");
h->GetYaxis()->SetTitleSize(0.06);
h->GetYaxis()->SetLabelSize(0.06);
h->GetZaxis()->SetLabelSize(0.06);
h->GetXaxis()->SetTitleSize(0.06);
h->GetXaxis()->SetLabelSize(0.06);

//create the correct spectrum
TH1D *corrected_exp = new TH1D("corrected_exp","corrected_exp",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));

Int_t counts; 
for (Int_t bin = 1; bin <= nbins; bin++)
{
  counts = expSpec->GetBinContent(bin);
  for (Int_t i = 0; i < counts; i++)
  {
    corrected_exp->Fill(projY[bin]->GetRandom(rand));
    //corrected_exp->Fill(expSpec->GetBinCenter(bin));
  }
  
}
corrected_exp->SetFillColor(kRed);
corrected_exp->SetFillStyle(3005);
corrected_exp->SetLineColor(kRed);
corrected_exp->SetLineWidth(2);
expSpec->SetLineColor(kBlack);
expSpec->SetLineWidth(2);

cvID =2;
corrCv->cd(cvID);


// Configurações do canvas
corrCv->cd(cvID)->SetRightMargin(0.03);
corrCv->cd(cvID)->SetTopMargin(0.08);
corrCv->cd(cvID)->SetLeftMargin(0.14);
corrCv->cd(cvID)->SetBottomMargin(0.14);

expSpec->SetTitle("");
expSpec->Draw();

gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
gStyle->SetTextSize(0.04);
gStyle->SetTitleYSize(0.04);
expSpec->GetXaxis()->SetTitle("Energy (MeV)");
expSpec->GetXaxis()->SetTitleSize(0.06);
expSpec->GetXaxis()->SetLabelSize(0.06);

expSpec->GetYaxis()->SetTitle("counts");
expSpec->GetYaxis()->SetMaxDigits(2);
expSpec->GetYaxis()->SetTitleSize(0.06);
expSpec->GetYaxis()->SetLabelSize(0.06);


gPad->SetGridx();
gPad->SetGridy();
corrected_exp->Draw("same");

// Find the maximum counts between corrected_exp and expSpec
Double_t maxCounts = TMath::Max(corrected_exp->GetMaximum(), expSpec->GetMaximum());
std::cout << "The maximum counts between corrected_exp and expSpec is: " << maxCounts << std::endl;

expSpec->GetYaxis()->SetRangeUser(0,maxCounts*1.15);


//regarding the directory    
std::string directory = "figures";

// Verificar se o diretório existe usando gSystem->AccessPathName
if (gSystem->AccessPathName(directory.c_str())) { // Diretório NÃO existe
    if (gSystem->mkdir(directory.c_str(), true) != 0) { // Tentar criar o diretório
        std::cerr << "Error creating directory: " << directory << std::endl;
        return nullptr;
    }
    std::cout << "Directory created successfully: " << directory << std::endl;
} else {
    std::cout << "Directory already exists: " << directory << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//generate correction factor for lost particles
//Create a canvas to temporarly store the corrected spectrum

TH1D *c2_exp;
if(correctDead){
    TCanvas *corrFcv = new TCanvas("corrFactorCv","correction Factor -- lost particles",150,150,800,678);
    corrFcv->SetLeftMargin(0.14);
    cout<<"\n.\n.\n.nbins = "<<nbins<<endl;

    TH1D *hprod = new TH1D("hprod","hprod",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));
    TH1D *hlost = new TH1D("hlost","hlost",nbins,expSpec->GetBinLowEdge(1),expSpec->GetBinLowEdge(nbins+1));

    std::string directory, filename;

    splitPath(pElossFileName, directory, filename);

    cout<<"Debugging splitPath:"<<endl;
    cout<<"directory = "<<directory<<endl;
    cout<<"filename = "<<filename<<endl;

    TGraph *Ffunc = GetPfactor(filename,directory,800, false);//= Fcorr(ElTree,nbins,false);
    Ffunc->GetXaxis()->SetTitle("Energy (MeV)");
    Ffunc->GetYaxis()->SetTitle("F factor ");
    Ffunc->Draw("ALP");

    gPad->SetGridx();
    gPad->SetGridy();


    c2_exp = (TH1D*)corrected_exp->Clone();
    c2_exp->SetNameTitle("c2_exp","c2_exp");

    for(int i = 1; i <= c2_exp->GetNbinsX(); i++){
        c2_exp->SetBinContent(i,corrected_exp->GetBinContent(i)*Ffunc->Eval(corrected_exp->GetBinCenter(i)));
    }
    corrCv->cd(2);
    c2_exp->SetLineColor(kBlue);
    c2_exp->SetFillColor(kCyan);
    c2_exp->SetFillStyle(3004);
    c2_exp->SetLineWidth(2);
    c2_exp->Draw("same");
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(expSpec,"Experimental Spectrum","lpf");
    leg->AddEntry(corrected_exp,"Corrected Spectrum","lpf");
    leg->AddEntry(c2_exp,"Corrected Spectrum with F factor","lpf");
    leg->Draw();

}else{ // if correctDead is false...
    corrCv->cd(2);
    TLegend *leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(expSpec,"Experimental Spectrum","lpf");
    leg->AddEntry(corrected_exp,"Corrected Spectrum","lpf");
    leg->Draw();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

if(saveCanvas){
    corrCv->SaveAs(Form("%s/%s_TTC.root",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.png",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    corrCv->SaveAs(Form("%s/%s_TTC.pdf",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()));
    cout << "Saving canvas as: " << Form("%s/%s_TTC.root (.pdf and .png)",directory.c_str(),pSpecFileName.substr(0, pSpecFileName.find_last_of('.')).c_str()) << endl;
    
}


timer.Print();
delete ElTree;
if(correctDead)return c2_exp;
else return corrected_exp;

}


string particleName(char particle = 'p'){

    string name;
    switch(particle){
        case 'p':
            name = "proton";
            break;
        case 'd':
            name = "deuteron";
            break;	
        case 't':
            name = "triton";
            break;
        case 'h':
            name = "He-3";
            break;
        case 'a':
            name = "alpha";
            break;
        default: 
            name = "proton";	
    }

    return name;

}