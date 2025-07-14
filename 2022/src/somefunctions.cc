#include "../include/somefunctions.hh"
#include <iostream>
#include <cmath>
#include <ctime>
#include <TSpectrum.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1D.h>

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


TChain* loadRuns(int nruni =63 ,int nrunf =63, double *Qtot = nullptr, Int_t *TotalTime = nullptr, string infoPath = "/mnt/medley/LucasAnalysis/2022/runlist2022.csv") {

	// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . 
	//charge loading: 
	TTree* InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(infoPath.c_str(), "RunN/I:MedleyConfig/C:Target/C:RunTime/I:RunCharge/F");

    // Declarar variáveis para leitura dos dados
    
    Float_t RunCharge;
 	Int_t RunN;
	Float_t TotalChargeFaraday = 0;
    Int_t runtime = 0;

	InfoTree->SetBranchAddress("RunN", &RunN);
	InfoTree->SetBranchAddress("RunCharge", &RunCharge);
    InfoTree->SetBranchAddress("RunTime", &runtime);
    if(TotalTime) *TotalTime = 0;

	// -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . 

    if(!nruni){
        cout<<"No run number given, using default run 111"<<endl;
        nruni = 111;
        nrunf = 111;
    }

    char name[1000];
    if (!tx)
        tx = new TChain("RD");

    vector<Int_t> entries;

    for (int i = nruni; i <= nrunf; i++)
    {
        Bool_t fExist = true;
        for (Int_t j = 0; j <= 999 && fExist; j++)
        {
            ifstream mfile;
            sprintf(name, "%s/r%04d_%03dr.root",pathRuns.c_str(), i, j);
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
                cout << "Entries: " << entries.back() / 1000 << "k " << endl;

				//charge:
				for(int t=0;t<InfoTree->GetEntries();t++)
				{
				    InfoTree->GetEntry(t);
					if(RunN == i){
						cout << "RunN: " << RunN << " ChargeFaraday: " << RunCharge <<" µC."<< endl;
						TotalChargeFaraday += RunCharge;
						cout<<"."<<endl;
                        cout<< "RunTime: " << runtime << " seconds."<< endl;
                        if(TotalTime) *TotalTime += runtime;
					}
            	}	
			}
        }
    }
	
	if (Qtot)
        *Qtot = TotalChargeFaraday;
    return tx;

}


void loadCuts(int telN=1){
    cout<<"\n\n --- --- --- LOADING FOR TELESCOPE "<<telN<<" --- --- ---\n"<<endl;
         
    //NPT CUTS dEdE (upper bands) =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //protons npt
        cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading NPT PROTONS CUT: %s/p_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/p_npt_%d.C", cuts_path.c_str(),telN));
        cut_protons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_npt_%d",telN));
        if(cut_protons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_protons[telN]->SetName(Form("p_npt_%d",telN));
        cut_protons[telN]->SetLineWidth(2);
        cut_protons[telN]->SetLineColor(kRed);
    
        //deuterons npt
        cout << Form("\n->[tel %d] loading NPT DEUTERONS CUT: %s/d_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/d_npt_%d.C", cuts_path.c_str(),telN));
        cut_deuterons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_npt_%d",telN));
        if(cut_deuterons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_deuterons[telN]->SetName(Form("d_npt_%d",telN));
        cut_deuterons[telN]->SetLineWidth(2);
        cut_deuterons[telN]->SetLineColor(kGreen);

        //tritons npt
        cout << Form("\n->[tel %d] loading NPT TRITONS CUT: %s/t_npt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/t_npt_%d.C", cuts_path.c_str(),telN));
        cut_tritons[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_npt_%d",telN));
        if(cut_tritons[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_tritons[telN]->SetName(Form("t_npt_%d",telN));
        cut_tritons[telN]->SetLineWidth(2);
        cut_tritons[telN]->SetLineColor(kBlack);
    
    //NPT CUTS dEdE (upper bands) =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //helium 
        cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH (he,alpha) CUTS FOR dE1:dE2 (UPPER BRANCH)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading NPT HE3 CUT: %s/he3_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/he3_%d.C", cuts_path.c_str(),telN));
        cut_he[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("he3_%d",telN));
        if(cut_he[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_he[telN]->SetName(Form("he3_%d",telN));
        cut_he[telN]->SetLineWidth(2);
        cut_he[telN]->SetLineColor(kRed);
    
        //alphas 
        cout << Form("\n->[tel %d] loading NPT ALPHAS CUT: %s/alphas%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/alphas%d.C", cuts_path.c_str(),telN));
        cut_alphas[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("alphas%d",telN));
        if(cut_alphas[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_alphas[telN]->SetName(Form("alphas%d",telN));
        cut_alphas[telN]->SetLineWidth(2);
        cut_alphas[telN]->SetLineColor(kGreen);

    //PT CUTS in dEdE (bottom band).=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //protons
        cout << Form("\n- - -   LOADING PUNCH-THROUGH CUTS FOR dE1:dE2 (BOTTOM BRANCH)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading PT PROTONS CUT: %s/p_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/p_pt_dE_%d.C", cuts_path.c_str(),telN));
        cut_protonsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_dE_%d",telN));
        if(cut_protonsPT[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_protonsPT[telN]->SetName(Form("p_pt_dE_%d",telN));
        cut_protonsPT[telN]->SetLineWidth(2);
        cut_protonsPT[telN]->SetLineColor(kRed);
        //deuterons:
        cout << Form("\n->[tel %d] loading PT DEUTERONS CUT: %s/d_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/d_pt_dE_%d.C", cuts_path.c_str(),telN));
        cut_deuteronsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_dE_%d",telN));
        if(cut_deuteronsPT[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_deuteronsPT[telN]->SetName(Form("d_pt_dE_%d",telN));
        cut_deuteronsPT[telN]->SetLineWidth(2);
        cut_deuteronsPT[telN]->SetLineColor(kGreen);
        //tritons:
        cout << Form("\n->[tel %d] loading PT TRITONSCUT: %s/t_pt_dE_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/t_pt_dE_%d.C", cuts_path.c_str(),telN));
        cut_tritonsPT[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_dE_%d",telN));
        if(cut_tritonsPT[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_tritonsPT[telN]->SetName(Form("t_pt_dE_%d",telN));
        cut_tritonsPT[telN]->SetLineWidth(2);
        cut_tritonsPT[telN]->SetLineColor(kBlack);

    // //PT csi CUTS =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //protons
        cout << Form("\n- - -   LOADING PUNCH-THROUGH CUTS FOR dE2:Eres (CsI)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading PT PROTONS (CSI) CUT: %s/p_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/p_pt_%d.C", cuts_path.c_str(),telN));
        cut_protonsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("p_pt_%d",telN));
        if(cut_protonsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_protonsCsI[telN]->SetName(Form("p_pt_%d",telN));
        cut_protonsCsI[telN]->SetLineWidth(2);
        cut_protonsCsI[telN]->SetLineColor(kRed);
        //deuterons:
        cout << Form("\n->[tel %d] loading PT DEUTERONS CUT: %s/d_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/d_pt_%d.C", cuts_path.c_str(),telN));
        cut_deuteronsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("d_pt_%d",telN));
        if(cut_deuteronsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_deuteronsCsI[telN]->SetName(Form("d_pt_%d",telN));
        cut_deuteronsCsI[telN]->SetLineWidth(2);
        cut_deuteronsCsI[telN]->SetLineColor(kGreen);
        //tritons:
        cout << Form("\n->[tel %d] loading PT TRITONS CUT: %s/d_pt_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/t_pt_%d.C", cuts_path.c_str(),telN));
        cut_tritonsCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("t_pt_%d",telN));
        if(cut_tritonsCsI[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_tritonsCsI[telN]->SetName(Form("t_pt_%d",telN));
        cut_tritonsCsI[telN]->SetLineWidth(2);
        cut_tritonsCsI[telN]->SetLineColor(kBlack);

    // //NPT csi CUTS =.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.=.
        //protons
        cout << Form("\n- - -   LOADING NON-PUNCH-THROUGH CUTS FOR dE2:Eres (CsI)   - - -\n")<< "... " ;
        cout << Form("\n->[tel %d] loading NPT PROTONS (CSI) CUT: %s/pnptcsi_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/pnptcsi_%d.C", cuts_path.c_str(),telN));
        cut_protonsNPTCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("pnptcsi_%d",telN));
        if(cut_protonsNPTCsI[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_protonsNPTCsI[telN]->SetName(Form("pnptcsi_%d",telN));
        cut_protonsNPTCsI[telN]->SetLineWidth(2);
        cut_protonsNPTCsI[telN]->SetLineColor(kRed);
        //deuterons:
        cout << Form("\n->[tel %d] loading NPT DEUTERONS CUT: %s/dnptcsi_%d.C",telN, cuts_path.c_str(),telN)<< "... " ;
        gROOT->ProcessLine(Form(".L  %s/dnptcsi_%d.C", cuts_path.c_str(),telN));
        cut_deuteronsNPTCsI[telN] = (TCutG *)gROOT->GetListOfSpecials()->FindObject(Form("dnptcsi_%d",telN));
        if(cut_deuteronsNPTCsI[telN]!=NULL) cout << "  --> OK."<< endl;
        cut_deuteronsNPTCsI[telN]->SetName(Form("dnptcsi_%d",telN));
        cut_deuteronsNPTCsI[telN]->SetLineWidth(2);
        cut_deuteronsNPTCsI[telN]->SetLineColor(kGreen);

}


void setParticlesAliases(int telN=1 ){
   
    loadCuts(telN);
    tx->SetAlias("protons", Form("(p_npt_%d && pnptcsi_%d)||(p_pt_dE_%d && p_pt_%d)",telN,telN,telN,telN));
    tx->SetAlias("deuterons", Form("!protons && ((d_npt_%d && dnptcsi_%d)||(d_pt_dE_%d && d_pt_%d))",telN,telN,telN,telN));
    tx->SetAlias("tritons", Form("!protons &&!deuterons ((t_npt_%d && Eres<0.5)||(t_pt_dE_%d && t_pt_%d))",telN,telN,telN));
    tx->SetAlias("he3", Form("he3_%d",telN));
    tx->SetAlias("alphas", Form("alphas%d",telN));
    
    cout<<"Aliases set for telescope "<<telN<<":\n";
    cout<<"protons = (p_npt_"<<telN<<" && pnptcsi_"<<telN<<")||(p_pt_dE_"<<telN<<" && p_pt_"<<telN<<")"<<endl;
    cout<<"deuterons = !protons && ((d_npt_"<<telN<<" && dnptcsi_"<<telN<<")||(d_pt_dE_"<<telN<<" && d_pt_"<<telN<<"))"<<endl;
    cout<<"tritons = !protons &&!deuterons ((t_npt_"<<telN<<" && Eres<0.5)||(t_pt_dE_"<<telN<<" && t_pt_"<<telN<<"))"<<endl;
    cout<<"he3 = he3_"<<telN<<endl;
    cout<<"alphas = alphas"<<telN<<endl;
    
    if(cut_protons[telN] != nullptr){
        cut_protons[telN]->SetVarY("dE1");
        cut_protons[telN]->SetVarX("dE2");
        cout<<"Cut for protons set to dE1:dE2"<<endl;
    }
    
    if(cut_protonsPT[telN] != nullptr){    
        cut_protonsPT[telN]->SetVarY("dE1");
        cut_protonsPT[telN]->SetVarX("dE2");
        cout<<"Cut for protons PT set to dE1:dE2"<<endl;
    }

    if(cut_protonsCsI[telN]!= nullptr){
        cut_protonsCsI[telN]->SetVarY("dE2");
        cut_protonsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for protons CsI set to dE2:Eres"<<endl;
    }

    if(cut_protonsNPTCsI[telN]!=nullptr){
        cut_protonsNPTCsI[telN]->SetVarY("dE2");
        cut_protonsNPTCsI[telN]->SetVarX("Eres");    
        cout<<"Cut for protons NPT CsI set to dE2:Eres"<<endl;
    }

    if(cut_deuterons[telN] != nullptr){
        cut_deuterons[telN]->SetVarY("dE1");
        cut_deuterons[telN]->SetVarX("dE2");
        cout<<"Cut for deuterons set to dE1:dE2"<<endl;
    }
    
    if(cut_deuteronsPT[telN] != nullptr){
        cut_deuteronsPT[telN]->SetVarY("dE1");
        cut_deuteronsPT[telN]->SetVarX("dE2");
        cout<<"Cut for deuterons PT set to dE1:dE2"<<endl;
    }

    if(cut_deuteronsCsI[telN] != nullptr){
        cut_deuteronsCsI[telN]->SetVarY("dE2");
        cut_deuteronsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for deuterons CsI set to dE2:Eres"<<endl;
    }
    if(cut_deuteronsNPTCsI[telN] != nullptr){
        cut_deuteronsNPTCsI[telN]->SetVarY("dE2");
        cut_deuteronsNPTCsI[telN]->SetVarX("Eres");
        cout<<"Cut for deuterons NPT CsI set to dE2:Eres"<<endl;
    }

    if(cut_tritons[telN] != nullptr){
        cut_tritons[telN]->SetVarY("dE1");
        cut_tritons[telN]->SetVarX("dE2");
        cout<<"Cut for tritons set to dE1:dE2"<<endl;
    }
    if(cut_tritonsPT[telN] != nullptr){
        cut_tritonsPT[telN]->SetVarY("dE1");
        cut_tritonsPT[telN]->SetVarX("dE2");
        cout<<"Cut for tritons PT set to dE1:dE2"<<endl;
    }
    if(cut_tritonsCsI[telN] != nullptr){
        cut_tritonsCsI[telN]->SetVarY("dE2");
        cut_tritonsCsI[telN]->SetVarX("Eres");
        cout<<"Cut for tritons CsI set to dE2:Eres"<<endl;  
    }
    if(cut_he[telN] != nullptr){
        cut_he[telN]->SetVarY("dE1");
        cut_he[telN]->SetVarX("dE2");
        cout<<"Cut for he3 set to dE1:dE2"<<endl;
    }
    if(cut_alphas[telN] != nullptr){
        cut_alphas[telN]->SetVarY("dE1");
        cut_alphas[telN]->SetVarX("dE2");
        cout<<"Cut for alphas set to dE1:dE2"<<endl;
    }
    
    return;
}

void setEnergyAliases(int telN=1 ){
    cout <<"g's values: g1 = "<<g1<<", g2 = "<<g2<<", g3 = "<<g3<<endl;
    tx->SetAlias("dE1", Form("Medley_%d_dE1*%f",telN,g1));
    tx->SetAlias("dE2", Form("Medley_%d_dE2*%f",telN,g2));
    tx->SetAlias("Eres", Form("Medley_%d_Eres*%f",telN,g3));
    tx->SetAlias("E", "dE1 + dE2 + Eres");
	cout<<" - - -"<<endl;        
    cout<<"Aliases set for telescope "<<telN<<":\n";
    cout<<"dE1 = Medley_"<<telN<<"_dE1 * "<<g1<<endl;
    cout<<"dE2 = Medley_"<<telN<<"_dE2 * "<<g2<<endl;
    cout<<"Eres = Medley_"<<telN<<"_Eres * "<<g3<<endl;
    cout<<"E = dE1 + dE2 + Eres"<<endl;
}
void setGuesses(Double_t guess1=1, Double_t guess2 =1, Double_t guess3 = 1){
    g1 = guess1;
    g2 = guess2;
    g3 = guess3;
    cout<<"\n.\n.\n.Guesses successfully set to:"<<endl;
    cout<<g1<<endl;
    cout<<g2<<endl;
    cout<<g3<<endl;
    cout<<" - - -"<<endl;
    return;
}


Float_t GetGflash(TTree *tx, const std::string& branchname = "Medley_1_dE2_ToF", float guess = 400, bool closecanvas = false) {



	TCanvas *Ctof = new TCanvas("ctof",Form("Time of flight"),50,50,600,600);

    TH1D *htof;
    htof = new TH1D("htof","htof",500,100,600);

    tx->Draw(Form("%s>>htof",branchname.c_str()));

    double tof_peak[2];
    TSpectrum *spec = new TSpectrum();
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


