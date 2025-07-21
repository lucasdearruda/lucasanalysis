// #include <iostream>
// #include <fstream>
// #include <string>
// #include "TChain.h"
// #include "TTree.h"
// #include "TFile.h"
// #include "TStopwatch.h"
// #include "KVMaterial.h" // Ensure this header is included

//#include "/mnt/medley/LucasAnalysis/useful.h"

//function to tranform ToF neutrons' histo to Energy 
TH1D *ScaleXhisto(TH1D *toBeScaled, double scaleFactor = 1.0, string name = "scaled",bool verbose = false){//verbose put some info regarding the function working 

	Int_t nbins = toBeScaled->GetNbinsX();
	if(verbose)cout<<"There are "<<nbins<<" bins. The scale factor provided is "<<scaleFactor<<"."<<endl; 
	double xb[nbins+1];
	for(Int_t i=0;i<nbins+1;i++){
		if(verbose){
			cout<<"bin "<<i+1<<": "<<toBeScaled->GetBinLowEdge(i+1)<<" -->"<<scaleFactor*toBeScaled->GetBinLowEdge(i+1)<<endl;
		}
		xb[i] = scaleFactor*toBeScaled->GetBinLowEdge(i+1);
	}

	TH1D *hscaled = new TH1D(name.c_str(),name.c_str(),nbins,xb);
	for(Int_t i=1;i<nbins+1;i++){
		hscaled->SetBinContent(i,toBeScaled->GetBinContent(i));
	}

	return hscaled;

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

TH1D *newE(
    TH1D *hh = nullptr, 
    string nametitle = "h", 
    bool considerangle = false,
    Double_t tSi1 = 53.4, 
    Int_t Z = 1,
    Int_t A = 1
){
    TStopwatch timer;
    timer.Start();
    
    if(!hh){
        cerr<<"Error: The input histogram is null."<<endl;
        return nullptr;
    }
    
    KVMaterial *det1 = new KVMaterial("Si");
    det1->SetThickness(tSi1 * KVUnits::um);
    cout << "Thickness of Si1: " << tSi1 << " µm" << endl;

    Int_t nbins = hh->GetNbinsX();
    Int_t ncounts = hh->GetEntries();
    TH1D *h = new TH1D(nametitle.c_str(), nametitle.c_str(), nbins, hh->GetBinLowEdge(0),  hh->GetBinLowEdge(nbins+1));


    // TRandom2 *rn = new TRandom2();
    // std::time_t seed = std::time(nullptr);
    // rn->SetSeed(seed);
    
    
    vector<Double_t> angles;
    if(considerangle) vector<Double_t> angles = giveMeTheAngle2(ncounts, 20);
    
    
    
    Double_t si1, newE,angle_factor;

    angle_factor = 1;
    Long64_t countsOff = 0;
    double maxE = det1->GetEIncOfMaxDeltaE(Z, A);
    for(Int_t i=0;i<ncounts;i++){
        si1 = hh->GetRandom();
        
        if(considerangle){
            angle_factor = 1 /cos(angles[i]);
            det1->SetThickness(tSi1 * angle_factor * KVUnits::um);    
        }
        
        if(si1 > maxE){
            countsOff++;
            cerr <<countsOff<< ", Warning: si1 = " << si1 << " MeV > " << maxE<<" MeV!"<<endl;
            continue; // Skip this iteration if si1 is too high
        }
        newE = si1 + det1->GetEResFromDeltaE(Z, A, si1)/KVUnits::MeV;
        h->Fill(newE);
    }
    cout<<"Efficiency:  "<<countsOff<<"/"<<ncounts<<" = "<<(1.0 - (Double_t)countsOff/ncounts)*100<<"%"<<endl;

    timer.Print();
    return h;

}
