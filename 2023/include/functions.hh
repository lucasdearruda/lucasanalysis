#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <string>
#include <cmath>
#include <iostream>
#include <TTree.h> // TTree é necessário porque você usa a função GetGflash
#include <time.h> //for the time


using namespace std;

double En(double tof, double L = 4.6472){//L in meters, 2023 reference
	//return neutron energy in MeV, given its time of fligt.

	//if(tof<=0){ cout<<"ERROR. Tof<=0 is not allowed. Check the other steps. I will return -1 for that."<<endl; return -1;}
	if(tof<=0){ return 0;}
	double v = L/tof;
	double c = 0.299792458; //m/ns
	double gamma = 1/sqrt(1-(v/c)*(v/c));
	if(gamma!=gamma) return 0; //case where gamma = nan;
	double mc2 = 939.56542052; //MeV
	return (gamma -1)*mc2; //MeV
}


double ToFparticleNumber(double Eparticle, int particle, double L = 4.6472)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'
// 1: proton
// 2: deuteron
// 3: triton
// 4: helium-3
// 5: alpha
// 6: neutron -- not used of course (but the map follows it in the function)


	double c= 0.299792458;
	double m;

	switch(particle){//data from codata
        	case 6://neutrons
        		m = 939.56542052;
        		break;
		case 1://protons
			m = 938.27208816;
			break;
		case 2://deuterons
			m = 1875.61294257;
			break;
		case 3://tritons 
			m = 2808.92113298;
			break;
		case 4: //helium3
			m =  2808.39160743;
			break;
		case 5://alphas
			m = 3727.3794066;
			break;
		default:
			m = 938.27208816;//proton
	}

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


Float_t GetGflash(TTree *tx, const std::string& branchname = "Medley_1_dE2_ToF", float guess = 500, bool closecanvas = false) {
//Float_t GetGflash(TTree *tx, string branchname="Medley_1_dE2_ToF", float guess = 500, bool closecanvas = false){


	TCanvas *Ctof = new TCanvas("ctof",Form("Time of flight"),50,50,600,600);

    TH1D *htof;
    htof = new TH1D("htof","htof",500,100,600);

    tx->Draw(Form("%s>>htof",branchname.c_str()));

    double tof_peak[2];
    TSpectrum *spec = new TSpectrum();
    Int_t npeaks = spec->Search(htof, 6,"",0.1);//originally 0.3
    Double_t *xpeaks = spec->GetPositionX();
    cout<<"ToF plot: "<<npeaks<<" peak found.\n";
    for(int t =0;t<npeaks;t++){
        cout<<"position peak"<<t<<": "<<xpeaks[t]<<" ns.\n";
    }

    if (npeaks >= 2) {
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = xpeaks[npeaks -1] - 10;
        double fitRangeMax = xpeaks[npeaks -1] + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[npeaks -1], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }else{ //if did not find the peak
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = guess - 10;
        double fitRangeMax = guess  + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[0], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }


	if(closecanvas) Ctof->Close();
	
	return tof_peak[0];
}

double provideTgama(int telN, int year = 2023){//just to save some lines in my script 

	double tgamma;

	if(year == 2023){
		switch(telN){
			case 1:
				tgamma = 15.742;
				break;
			case 2:
				tgamma = 15.700;
				break;
			case 3:
				tgamma = 15.635;
				break;
			case 4:
				tgamma = 15.554;
				break;
			case 5:
				tgamma = 15.465;
				break;
			case 6:
				tgamma = 15.381;
				break;
			case 7:
				tgamma = 15.311;
				break;
			case 8:
				tgamma = 15.264;
				break;
			default:
				tgamma = 15.742;
		}
	}else{
		tgamma = 15;
		cout<<"****************************************************************"<<endl;
		cout<<"	year "<<year<<" nor available!"<<endl;
		cout<<"****************************************************************"<<endl;
	}
	return tgamma;
}

std::string getCurrentTime() {
    char cur_time[128];
    time_t t = time(NULL);
    struct tm* ptm = localtime(&t);

    // Formata o tempo no formato desejado: "YYYY-MM-DD_HH:MM:SS"
    strftime(cur_time, sizeof(cur_time), "%Y-%m-%d_%H:%M:%S", ptm);

    return std::string(cur_time);
}


#endif // FUNCTIONS_HH
