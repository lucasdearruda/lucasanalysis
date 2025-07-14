#ifndef SOMEFUNCTIONS_HH
#define SOMEFUNCTIONS_HH

#include <string>
#include <optional>
#include <TTree.h>



//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//CUTS VARIABLES FOR PID:

//cuts for npt particles
TCutG *cut_protons[5], *cut_deuterons[5], *cut_tritons[5]; 
TCutG *cut_he[5], *cut_alphas[5];

//cuts for pt particles
TCutG *cut_protonsPT[5], *cut_deuteronsPT[5], *cut_tritonsPT[5];
TCutG *cut_protonsCsI[5], *cut_deuteronsCsI[5], *cut_tritonsCsI[5]; 

//cut npt csi 
TCutG *cut_protonsNPTCsI[5], *cut_deuteronsNPTCsI[5];

string cuts_path =   "/mnt/medley/LucasAnalysis/2022/PID";

string pathRuns = "/mnt/medley/runswithMM_2022";


//TFile *g[6];

//TTree *TR[6];

//std::vector<TGraph*> myplots;
TChain *tx = NULL;

Long64_t nEntries[6];
string filename[6];
string pathCal = "/mnt/medley/kaliveda_results/Telescopes_dE_E_results2023";

Double_t g1, g2, g3;

//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.


struct MedleyData {
    Int_t RunN;
    std::string MedleyConfig;
    std::string Target;
    Int_t RunTime;
    Float_t RunCharge;
};

// Funções declaradas:
std::optional<MedleyData> GetRunData(const std::string& filepath, Int_t requestedRunN);
double En(double tof, double L = 5.044);
double ToFparticleNumber(double Eparticle, int particle, double L = 5.044);
double provideTgama(int telN, int year = 2022, bool normal_config = true);
std::string getCurrentTime();
Double_t getDistanceTel(int telN = 1);
double ToFneutron(double Eparticle, double L = 0.1497);

//TChain* loadRuns(int nruni =63 ,int nrunf =63);
TChain* loadRuns(int nruni =63 ,int nrunf =63, double *Qtot = nullptr, string infoPath = "/mnt/medley/LucasAnalysis/2022/runlist2022.csv");
void loadCuts(int telN=1);
void setParticlesAliases(int telN=1 );
void setEnergyAliases(int telN=1 );
void setGuesses(Double_t guess1=1, Double_t guess2 =1, Double_t guess3 = 1);
//Float_t GetGflash(TTree *tx, const std::string& branchname = "Medley_1_dE2_ToF", float guess = 400, bool closecanvas = false);

#endif // SOMEFUNCTIONS_HH
