#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <string>
#include <cmath>
#include <optional>
#include <iostream>
#include <TTree.h> // TTree é necessário porque você usa a função GetGflash
#include <time.h> //for the time


using namespace std;


// Definição do struct para armazenar os dados
struct MedleyData {
    Int_t RunN;
    std::string MedleyConfig;
    std::string Target;
    Int_t RunTime;
    Float_t RunCharge;
};

// Função para obter os dados de uma run específica
std::optional<MedleyData> GetRunData(const std::string& filepath, Int_t requestedRunN) {
    // Criar o TTree e ler o arquivo
    TTree* InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(filepath.c_str(), "RunN/I:MedleyConfig/C:Target/C:RunTime/I:RunCharge/F");

    // Declarar variáveis para leitura dos dados
    Int_t RunN;
    char MedleyConfig[256];
    char Target[256];
    Int_t RunTime;
    Float_t RunCharge;

    // Associar as variáveis às branches do TTree
    InfoTree->SetBranchAddress("RunN", &RunN);
    InfoTree->SetBranchAddress("MedleyConfig", MedleyConfig);
    InfoTree->SetBranchAddress("Target", Target);
    InfoTree->SetBranchAddress("RunTime", &RunTime);
    InfoTree->SetBranchAddress("RunCharge", &RunCharge);

    // Loop sobre as entradas no TTree
    Long64_t nEntries = InfoTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        InfoTree->GetEntry(i);

        // Verificar se o RunN atual é o solicitado
        if (RunN == requestedRunN) {
            MedleyData entry;
            entry.RunN = RunN;
            entry.MedleyConfig = MedleyConfig;
            entry.Target = Target;
            entry.RunTime = RunTime;
            entry.RunCharge = RunCharge;

            delete InfoTree; // Liberar memória
            return entry; // Retorna os dados da run encontrada
        }
    }

    // Caso não encontre a run, liberar memória e retornar vazio
    delete InfoTree;
    return std::nullopt;
}

double En(double tof, double L = 5.044){//L in meters, 2022 reference
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


double ToFparticleNumber(double Eparticle, int particle, double L = 5.044)
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
    //Int_t npeaks = spec->Search(htof, 6,"",0.1);//originally 0.3
	Int_t npeaks = spec->Search(htof, 4,"",0.01);//originally 0.3
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

double provideTgama(int telN, int year = 2022, bool normal_config = true){//just to save some lines in my script 

	double tgamma = 15.0;
	if(year == 2023){

		switch(telN){
			case 1:
				tgamma = 15.742;
				if(!normal_config) tgamma = 15.282;
				break;
			case 2:
				tgamma = 15.700;
				if(!normal_config) tgamma = 15.324;
				break;
			case 3:
				tgamma = 15.635;
				if(!normal_config) tgamma = 15.388;
				break;
			case 4:
				tgamma = 15.554;
				if(!normal_config) tgamma = 15.467;
				break;
			case 5:
				tgamma = 15.465;
				if(!normal_config) tgamma = 15.553;
				break;
			case 6:
				tgamma = 15.381;
				if(!normal_config) tgamma = 15.632;
				break;
			case 7:
				tgamma = 15.311;
				if(!normal_config) tgamma = 15.698;
				break;
			case 8:
				tgamma = 15.264;
				if(!normal_config) tgamma = 15.742;
				break;
			default:
				tgamma = 15.742;
				if(!normal_config) tgamma = 15.282;
		}



	}else if(year == 2022){
		switch(telN){
			case 1:
				tgamma = 17.179;
				if(!normal_config) tgamma = 16.510;
				break;
			case 2:
				tgamma = 17.038;
				if(!normal_config) tgamma = 16.615;
				break;
			case 3:
				tgamma = 16.967;
				if(!normal_config) tgamma = 16.689;
				break;
			case 4:
				tgamma = 16.880;
				if(!normal_config) tgamma = 16.785;
				break;
			case 5:
				tgamma = 16.785;
				if(!normal_config) tgamma = 16.880;
				break;
			case 6:
				tgamma = 16.689;
				if(!normal_config) tgamma = 16.967;
				break;
			case 7:
				tgamma = 16.615;
				if(!normal_config) tgamma = 17.038;
				break;
			case 8:
				tgamma = 16.615;
				if(!normal_config) tgamma = 17.179;
				break;
			default:
				tgamma = 17.179;
				if(!normal_config) tgamma = 15.282;
		}
	}else{
		cout<<"****************************************************************"<<endl;
		cout<<"	year "<<year<<" nor available! Tgamma_std = "<<tgamma<<" ns."<<endl;
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

Double_t getDistanceTel(int telN = 1){
	switch (telN)
	{
	case 1:
		return 218.4*1e-3; //meters
		break;
	case 2:
		return 160.9*1e-3; //meters
		break;
	case 3:
		return 160.9*1e-3; //meters
		break;
	case 4:
		return 160.9*1e-3; //meters
		break;
	case 5:
		return 173.9*1e-3; //meters
		break;
	case 6:
		return 173.9*1e-3; //meters
		break;
	case 7:
		return 170.9*1e-3; //meters
		break;
	case 8:
		return 207.9*1e-3; //meters
		break;
	default:
		return 218.4*1e-3; //meters
		break;
	}

}


double ToFneutron(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'

	if(Eparticle<=0){ cout<<Form("ToFproton ERROR. Eparticle(%.3f MeV)<=0 is not allowed. Check the other steps.",Eparticle)<<endl; return 0;}
	double c= 0.299792458;
	double m = 939.56542052;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


#endif // FUNCTIONS_HH
