//version 2025.08.07.001
Double_t mc;
Double_t M;
Double_t omega_tel = 0.018 ;//sr ---> after I will add fine adjustments to i
Double_t NA = 6.02214076e23; //mol^-1
Double_t Nc = (mc/M) * NA; //atoms/cm^3
Int_t A, Z;

Float_t binWidth = 0.5; //MeV

using namespace std;




unsigned long long version(bool show = true) {
    string versionStr = "2025.08.07.001";

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
    }

    return versionNum;
}


void getAZ(char particle, Int_t &Z, Int_t &A) {
    switch (particle) {
        case 'p': Z = 1; A = 1; break;  // próton
        case 'd': Z = 1; A = 2; break;  // deutério
        case 't': Z = 1; A = 3; break;  // trítio
        case 'h': Z = 2; A = 3; break;  // hélio-3
        case 'a': Z = 2; A = 4; break;  // alfa (hélio-4)
        default:  Z = 0; A = 0; break;  // inválido
    }
}


void Attribute_Target(string target_name = "Fe_thick_Medley"){
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

Float_t getNflux(Float_t Ea = 10,Float_t Eb = 12, string fluxfilename = "/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", bool verbose = false, bool drawIt = false){

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

    return nflux_En;
} 

TH1D *getFluxHisto(Float_t Ea = 10,Float_t Eb = 12, string fluxfilename = "/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", bool verbose = false, bool drawIt = false){

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
    case 5:
        return 20;
        break;
    case 6:
        return 20;
        break;
    case 7:
        return 20;
        break;
    case 8:
        return 20;
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
