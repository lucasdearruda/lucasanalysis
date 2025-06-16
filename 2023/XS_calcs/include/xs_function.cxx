#include "headers/xs_functions.hh"
#include "bin_width.hh"
#include "/mnt/medley/LucasAnalysis/2023/Etot/newE.hh"  
//#include "/mnt/medley/LucasAnalysis/useful.h"
#include "/mnt/medley/LucasAnalysis/2023/applyTTC/include/correctSpec2.cxx"
// Function to get the thickness of Si1 based on the angle
Double_t GetSi1Thickness(Double_t angle = 20){

    if(abs(angle - 20.0f) < 1e-3f) {
        return 53.4; // 25 µm for 20 degrees
    } else if (abs(angle - 40.0f) < 1e-3f) {
        return 50; // 25 µm for 40 degrees
    } else if (abs(angle - 60.0f) < 1e-3f) {
        return 61.6; // 25 µm for 60 degrees
    } else if (abs(angle - 80.0f) < 1e-3f) {
        return 50; // 25 µm for 80 degrees
    } else if (abs(angle - 100.0f) < 1e-3f) {
        return 53.4; // 25 µm for 100 degrees
    } else if (abs(angle - 120.0f) < 1e-3f) {
        return 50; // 25 µm for 120 degrees
    } else if (abs(angle - 140.0f) < 1e-3f) {
        return 61.6; // 25 µm for 140 degrees
    } else if (abs(angle - 160.0f) < 1e-3f) {
        return 50; // 25 µm for other angles
    }else{
        cerr << "Error: Angle not recognized. Returning default thickness." << endl;
        return 53.4; // Default thickness for unrecognized angles
    }

}


//Function to convert a double value to a string in the format "XpY" (e.g., 2.3 -> "2p3")
std::string DoubleToXpY(double value) {
    int intPart = static_cast<int>(value);
    int fracPart = static_cast<int>(std::round(std::abs((value - intPart) * 10))); // one digit

    std::ostringstream oss;
    oss << intPart << "p" << fracPart;
    return oss.str();
}

//// Function to get the atomic number (Z) and mass number (A) for a given particle
std::pair<int, int> get_ZA(char particle) {
    switch (particle) {
        case 'p': return {1, 1};  // proton
        case 'd': return {1, 2};  // deuteron
        case 't': return {1, 3};  // triton
        case 'h': return {2, 3};  // helium-3 (He3)
        case 'a': return {2, 4};  // alpha (He4)
        default:  return {0, 0};  // unknown
    }
}
// Definição das variáveis globais
Double_t mc = 0;  // inicialize com valor padrão ou deixe para setar em runtime
Double_t M = 0;

Double_t L = 464.72; //cm
Double_t cm2_to_mbarn = 1e+27 ;//cm²


const Double_t NA = 6.02214076e23; // mol^-1
Double_t Nc = 0;                   // precisa ser calculada, depois que mc e M forem definidos
const Double_t omega_tel = 0.040;  // sr

void Attribute_Target(string target_name = "Fe"){
    if(target_name == "Fe"){
        mc = 0.0851;//g 
        M = 55.845; //g/mol
        Nc = (mc/M) * NA;
        cout<<">>=.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.==.== "<<endl;
        cout<<"m = "<<mc<<" g."<<endl;
        cout<<"M = "<<M<<" g/mol."<<endl;
        cout<<"Nc = "<<Nc<<" ."<<endl;
        cout<<">>-.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.--.-- "<<endl;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to calculate the error in the integral of a histogram
Float_t histIntegralError(TH1D * histo){
    //Function to calculate the error in the integral of a histogram
    Float_t sum = 0;
    Int_t nbins = histo->GetNbinsX();
    for (Int_t i = 1; i <= nbins; i++) {
        if (histo->GetBinContent(i) > 0) {
            Float_t error = histo->GetBinError(i);
            sum += error * error * histo->GetBinWidth(i) * histo->GetBinWidth(i);
        }
    }

    return sqrt(sum);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//return the eq. number for a given particle 
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
//-------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//return the telescope number
Float_t omega_telescope(Float_t angle = 20){

    if (fabs(angle - 20.0f) < 1e-3f) {
    // safe comparison for float
        return 0.01800;
    }else if (fabs(angle - 40.0f) < 1e-3f) {
        return 0.01785;
    }else if (fabs(angle - 60.0f) < 1e-3f) {
        return 0.01781;
    }else if (fabs(angle - 80.0f) < 1e-3f) {
        return 0.01768;   
    }else if (fabs(angle - 100.0f) < 1e-3f) {
        return 0.01772;
    }else if (fabs(angle - 120.0f) < 1e-3f) {
        return 0.01805;
    }else if (fabs(angle - 140.0f) < 1e-3f) {
        return 0.01815;
    }else if (fabs(angle - 100.0f) < 1e-3f) {
        return 0.01814;
    }
    return 0.018;
}
//-------------------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to calculate the total Charge of a selection of runs 
//returns the total charge in µC
//  -> 'runs' is the vector of runs to be summed
//  -> 'verbose' is a boolean to print the runs and their charges
//  -> 'infofile' is the path to the runlist file containing the charge information
Float_t totalCharge(vector<Int_t> runs, bool verbose = true, string infofile =  "/mnt/medley/LucasAnalysis/2023/runlist.csv") {
    
    Float_t TotalCharge = 0;
    TTree *InfoTree = new TTree("InfoTree", "InfoTree");
    InfoTree->ReadFile(infofile.c_str(), "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F");

    if(verbose) {
        cout << "---------------------------------------------------" << endl;
        cout << "Total Charge Calculation" << endl;
        cout << "---------------------------------------------------" << endl;
        cout << "Info file: " << infofile << endl;
        cout << "Total Charge will be calculated in µC." << endl;
        cout << "Number of runs: " << runs.size() << endl;
        if (runs.empty()) {
            cout << "No runs provided. Returning 0." << endl;
            return 0;
        }
        cout << "Calculating total charge for runs: ";  
        for (const auto& run : runs) {
            cout << run << " ";
        }
        cout << endl;
    }

    Int_t RunN;
    Float_t ChargeFaraday;

    InfoTree->SetBranchAddress("RunN", &RunN);
    InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

 // Coloque os runs em um set para busca rápida
    std::set<Int_t> runSet(runs.begin(), runs.end());

    int nEntries = InfoTree->GetEntries();
    for(int i = 0; i < nEntries; i++) {
        InfoTree->GetEntry(i);
        if(runSet.find(RunN) != runSet.end()) {
            if(verbose) {
                cout << "RunN: " << RunN << ", ChargeFaraday: " << ChargeFaraday << " µC" << endl;
            }
            // Adiciona a carga Faraday ao total
            TotalCharge += ChargeFaraday;
        }
    }

    if(verbose) cout<<"---------------------------------------------------"<<endl;
    return TotalCharge;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to load several files into a TChain
// 'chain' is the TChain to which the files will be added
// 'runIntervals' is a vector containing pairs of run numbers, where each pair represents a range of runs to be loaded
// 'basePath' is the base path where the files are located

Long64_t LoadRunsToChain(TChain* chain, const std::vector<int>& runIntervals, const std::string& dir, bool verbose = true) {
    Long64_t totalEntries = 0;
    for (size_t i = 0; i + 1 < runIntervals.size(); i += 2) {
        int startRun = runIntervals[i];
        int endRun = runIntervals[i + 1];
        for (int run = startRun; run <= endRun; ++run) {
            std::string filename = Form("%s/%03d.root", dir.c_str(), run);

            if (verbose)
                std::cout << Form("%03d : ", run) << filename << std::endl;

            std::ifstream file(filename);
            if (!file) {
                if (verbose)
                    std::cout << "File does not exist: " << filename << std::endl;
                continue;
            }
            file.close();

            if (verbose)
                std::cout << "Adding " << filename << std::endl;

            chain->Add(filename.c_str());
        }
    }

    if (verbose) {
        totalEntries = chain->GetEntries();
        std::cout << "Total entries: " << totalEntries / 1000 << "k" << std::endl;
    }
    return totalEntries;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to transform from histo to TGraph
TGraph * histoToGraph(TH1D *hh){
    TGraph *gr = new TGraph();
    gr->SetName(hh->GetName());
    for(int b = 1; b<= hh->GetNbinsX();b++){
        gr->AddPoint(hh->GetBinCenter(b),hh->GetBinContent(b));
    }
    return gr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Function to get the differential cross-section for a given energy range and angles
// Ea: lower energy limit in MeV for the neutron energy
// Eb: upper energy limit in MeV for the neutron energy
// angles: vector of angles in degrees for which the differential cross-section will be calculated
// particles: vector of particles for which the differential cross-section will be calculated
// runsForward: vector of run numbers for the forward angle
// runsBackward: vector of run numbers for the backward angle
    //Disclaimer: ---------------------------------------------------------------------------------------------------------------
    // runsForward and runsBackward are vectors containing intvs of runs:
    //  std::vector<int> runsForward = {12, 12, 15, 15, 18, 21, 44, 44, 56, 59}; -->12, 15, 18 to 21, 44 to 59
    //to input this into totalCharge, we need to iterate over the vector two by two:
    // std::vector<int> runs = {12, 15, 18, 19, 20, 21, 44, 56, 57, 58, 59};
// Match_CORR: boolean to activate of not the matching correction
// TTC_CORR: boolean to activate of not the thick targe correction
// target: string with the name of the target; it is forwarded to the function Attribute_Target
// DrawFlux: boolean to draw the flux in the canvas
std::vector<std::vector<TH1D*>> GetDDX(
    float Ea = 25,
    float Eb = 26,
    const std::vector<float>& angles =  {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0},
    const std::vector<char>& particles = {'p', 'd', 't', 'h', 'a'},
    const std::vector<int>& runsForward = {},
    const std::vector<int>& runsBackward = {},
    bool Match_CORR = true,
    bool TTC_CORR = true,
    string target = "Fe",
    bool DrawFlux = false
){
    Attribute_Target(target);
    //First thing, calculate the total charge for the runs forward and backward: 
    //the vectors runsForward(Backward) contains intvs, so they whould be read two by two
    //reading charges:: 
    Int_t intvs = runsForward.size()/2; //Should always be pair
        if(runsForward.size() % 2 != 0 || runsBackward.size() % 2 != 0) {
        std::cerr << "Error: The number of runs in runsForward and runsBackward must be even." << std::endl;
        return {};
    } 
    
    
    //Disclaimer: ---------------------------------------------------------------------------------------------------------------
    // runsForward and runsBackward are vectors containing intvs of runs:
    //  std::vector<int> runsForward = {12, 12, 15, 15, 18, 21, 44, 44, 56, 59}; -->12, 15, 18 to 21, 44 to 59
    //to input this into totalCharge, we need to iterate over the vector two by two:
    // std::vector<int> runs = {12, 15, 18, 19, 20, 21, 44, 56, 57, 58, 59};
    Float_t TotalChargeForward = 0;
    Float_t TotalChargeBackward = 0;
    

    // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --
    //Sanity checking..
    // Check if runsForward has an odd number of elements

    // Make mutable copies to adjust if needed
    std::vector<int> runsForwardAdjusted = runsForward;
    std::vector<int> runsBackwardAdjusted = runsBackward;

    if (runsForwardAdjusted.size() % 2 != 0) {
        std::cerr << "Warning: runsForward has an odd number of elements. The last element will be ignored." << std::endl;
        runsForwardAdjusted.pop_back();  // remove last element
    }

    if (runsBackwardAdjusted.size() % 2 != 0) {
        std::cerr << "Warning: runsBackward has an odd number of elements. The last element will be ignored." << std::endl;
        runsBackwardAdjusted.pop_back();
    }
    // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --  // --

        std::vector<int> runs;

        for (size_t i = 0; i < runsForward.size(); i += 2) {
            int start = runsForward[i];
            int end = runsForward[i + 1];

            for (int run = start; run <= end; ++run) {
                runs.push_back(run);
            }
        }
    TotalChargeForward = totalCharge(runs, true);

    runs.clear();

    for (size_t i = 0; i < runsBackward.size(); i += 2) {
        int start = runsBackward[i];
        int end = runsBackward[i + 1];

        for (int run = start; run <= end; ++run) {
            runs.push_back(run);
        }
    }
    TotalChargeBackward = totalCharge(runs, true);
    //Now we have the total charge for the forward and backward runs, so we will print them for the user
    std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    std::cout << ".." << std::endl;
    std::cout << "Total Charge Forward: " << TotalChargeForward << " µC" << std::endl;
    std::cout << "Total Charge Backward: " << TotalChargeBackward << " µC" << std::endl;

    //Charge thing: ---------------------------------------------------------------------------------------------------------------


    TChain *tx_Forward = new TChain("M");
    TChain *tx_Backward = new TChain("M");

    Long64_t EntriesForward = 0;
    Long64_t EntriesBackward = 0;
    //Load the files into the TChains

    string processed_runs = "/mnt/medley/LucasAnalysis/2023/reducedv61";
    string name;
    
    EntriesForward =    LoadRunsToChain(tx_Forward, runsForward, processed_runs, true);
    EntriesBackward =    LoadRunsToChain(tx_Backward, runsBackward, processed_runs, true);
        


    std::cout << "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
    std::cout << ".." << std::endl;
    std::cout << "Total Charge Forward: " << TotalChargeForward << " µC" << std::endl;
    std::cout << "Total Charge Backward: " << TotalChargeBackward << " µC" << std::endl;
    std::cout << "Entries Forward: " << EntriesForward / 1000 << "k" << std::endl;
    std::cout << "Entries Backward: " << EntriesBackward / 1000 << "k" << std::endl;
    // ---------------------------------------------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Open nflux -- from Medley 
    //__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

    TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/nflux_direct/nflux.root", "READ");
    TH1D *nflux = (TH1D*)ff->Get("nflux");
    TCanvas * cflux = new TCanvas("cflux","cflux");
    nflux->Draw();

    TGraph *fluxGr = histoToGraph(nflux);
    // fluxGr->SetLineColor(kRed);
    // fluxGr->SetLineWidth(2);
    // fluxGr->Draw("same");

    Double_t nflux_En = fluxGr->Eval((Ea+Eb)/2.0);
    cout << "nflux_En ("<<(Ea+Eb)/2.0<<" MeV): " << nflux_En << " n / sr / µC / 1-MeV."<<endl;// --> "<<nflux_En*(Eb-Ea)<< " n / sr / µC ."<<endl;

    //__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

    std::vector<std::vector<TH1D*>> histograms;
    
    //histograms to be corrected
    std::vector<std::vector<TH1D*>> histograms_Ecorr;

    Double_t maxE = 40;
    Double_t minE = 0;
    Double_t binsdE = (maxE-minE)/bindE;
    
    
    //useful only if Match_CORR is true:
        Double_t binsFactor = 60; //factor to increase the number of bins

    for (size_t i = 0; i < particles.size(); ++i) {
        std::vector<TH1D*> hist_per_angle;
        std::vector<TH1D*> hist_per_angle_Ecorr;

        for (size_t j = 0; j < angles.size(); ++j) {
            char particle = particles[i];
            float angle = angles[j];

            std::string hname = Form("%c_ang_%.0f_deg_%.1f_%.1fMeV", particle, angle, Ea, Eb);
            std::string htitle = Form("%c_ang_%.0f_deg_%.1f_%.1fMeV", particle, angle, Ea, Eb);

            TH1D* h = new TH1D(hname.c_str(), htitle.c_str(), binsdE, minE,maxE);
            h->Sumw2();

            //To be Corrected
              std::string hname_Ecorr = Form("%c_ang_%.0f_deg_%.1f_%.1fMeV_MatchC", particle, angle, Ea, Eb);
              std::string htitle_Ecorr = Form("%c_ang_%.0f_deg_%.1f_%.1fMeV_MatchC", particle, angle, Ea, Eb);
            
              TH1D* h_Ecorr = new TH1D(hname_Ecorr.c_str(), htitle_Ecorr.c_str(), binsdE, minE, maxE);
              h_Ecorr->Sumw2();
              
              
              if(Match_CORR){//in this case we need more bins
                
                h_Ecorr = new TH1D(hname_Ecorr.c_str(), htitle_Ecorr.c_str(), binsFactor*binsdE, minE, maxE);
                h_Ecorr->Sumw2();
              }
              
            

            std::string cut = Form("ENN > %.3f && ENN < %.3f && PID == %d && ang == %.1f",
                                Ea, Eb, pCode(particle), angle);

            if (angle > 80.0) {
                tx_Backward->Draw(Form("E >> %s", hname.c_str()), cut.c_str(), "goff");
                
                //Ecorr

                if(Match_CORR){
                    tx_Backward->Draw(Form("si1 >> %s", hname_Ecorr.c_str()), cut.c_str(), "goff");
                }else{
                    tx_Backward->Draw(Form("E >> %s", hname_Ecorr.c_str()), cut.c_str(), "goff");
                }
            } else {
                tx_Forward->Draw(Form("E >> %s", hname.c_str()), cut.c_str(), "goff");

                //Ecorr
                if(Match_CORR){
                    tx_Forward->Draw(Form("si1 >> %s", hname_Ecorr.c_str()), cut.c_str(), "goff");
                }else{
                    tx_Forward->Draw(Form("E >> %s", hname_Ecorr.c_str()), cut.c_str(), "goff");
                } 
            }

            hist_per_angle.push_back(h);
            hist_per_angle_Ecorr.push_back(h_Ecorr);

        }

        histograms.push_back(hist_per_angle);
        histograms_Ecorr.push_back(hist_per_angle_Ecorr);
    }

    //Any transformations we may want for the histograms 

    Double_t binWidthNeutron = Eb-Ea; 


    //Double_t FactorF = pow(L,2)*cm2_to_mbarn/(Nc*omega_tel*TotalChargeForward*nflux_En);
    //Double_t FactorB = pow(L,2)*cm2_to_mbarn/(Nc*omega_tel*TotalChargeBackward*nflux_En);
    
    Double_t FactorF = pow(L,2)*cm2_to_mbarn/(Nc*omega_tel*TotalChargeForward*nflux_En*binWidthNeutron*bindE);
    Double_t FactorB = pow(L,2)*cm2_to_mbarn/(Nc*omega_tel*TotalChargeBackward*nflux_En*binWidthNeutron*bindE);

    cout<<"FactorF: "<<FactorF<<endl;
    cout<<"FactorB: "<<FactorB<<endl;


    //apply the factor bi by bin to transform counts to DDX
    for (Int_t i = 0; i < histograms.size(); ++i) {
        for (Int_t j = 0; j < histograms[i].size(); ++j) {

            float angle = angles[j];
            TH1D* h = histograms[i][j];
            TH1D* h_Ecorr = histograms_Ecorr[i][j];
            if(Match_CORR){
                ///needs to add the particle:
                std::pair<int, int> ZA  = get_ZA(particles[i]);
                int Z = ZA.first;
                int A = ZA.second;
                cout<<"Z: "<<Z<<", A: "<<A<<endl;
                
                
                //TCanvas *cc = new TCanvas("tempCanvas", "tempCanvas", 1400, 600);
                //cc->Divide(2,1);
                //cc->cd(1);
                //h_Ecorr->Draw("hist");
                ///wait for the user to close the canvas
                //cc->WaitPrimitive();
                //cc->cd(2);

                h_Ecorr = newE(h_Ecorr, h_Ecorr->GetTitle(), false, GetSi1Thickness(angle), Z, A);
                //wait user press any key to continue
                //h_Ecorr->Draw("hist");
                
                h_Ecorr->Rebin(binsFactor);
                // h_Ecorr->SetLineColor(kRed);
                // h_Ecorr->SetLineWidth(2);
                // h_Ecorr->SetDrawOption("hist");
                // h->SetLineColor(kBlack);
                // h->Draw("same hist");
                histograms_Ecorr[i][j] = h_Ecorr;
                
            }
            

            //TTC:
            if(TTC_CORR){
                
                //apply the correction to the histograms
                //h_Ecorr = correctSpec(h, true, Form("/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick/v7/eloss_%c_%02.1fdeg_025.0um.root",particles[i], angle), "Fe_thick_Medley", 25, particles[i], false);
               Double_t angle_ttc;
                //mapping the angles to angle_ttc, 100 deg -> 80 deg, 120 deg -> 60 deg, 140 deg -> 40 deg, 160 deg -> 20 deg
                if(angle > 80.0){
                    angle_ttc = -angle + 180.0;
                }else{
                    angle_ttc = angle;
                }
                cout<<"angle_ttc: "<<angle_ttc<<endl;
                string histoname = h_Ecorr->GetName();
                string histotitle = h_Ecorr->GetTitle();
                
                // // DEBUGGING...
                // Double_t sum_of_errors, sum_of_errors_after;
                // sum_of_errors = 0;
                // sum_of_errors_after = 0;
                // //calculate the sum of errors of all bins in h_Ecorr
                // for(Int_t b = 1; b <= h_Ecorr->GetNbinsX(); b++){
                //     sum_of_errors += h_Ecorr->GetBinError(b);
                // }
                

                h_Ecorr = correctSpec(h_Ecorr, true, Form("/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick/v7/eloss_%c_%02.1fdeg_025.0um.root",particles[i], angle_ttc), "Fe_thick_Medley", 25, particles[i], false);
                histograms_Ecorr[i][j] = h_Ecorr;
                
                // // DEBUGGING...
                // for(Int_t b = 1; b <= h_Ecorr->GetNbinsX(); b++){
                //     sum_of_errors_after += h_Ecorr->GetBinError(b);
                // }

                // cout<<"Sum of errors before correction: "<<sum_of_errors<<endl;
                // cout<<"Sum of errors after correction: "<<sum_of_errors_after<<endl;
                // //wait for the user to press ENTER
                // std::cout << "Pressione ENTER para continuar...";
                // std::cin.get();
                
                h_Ecorr->SetName(Form("%s",histoname.c_str()));
                h_Ecorr->SetTitle(Form("%s",histotitle.c_str()));

                histoname = h->GetName();
                histotitle = h->GetTitle();

                h = correctSpec(h, true, Form("/home/pi/ganil/kalscripts/eloss/results/UniformZ/Fe_thick/v7/eloss_%c_%02.1fdeg_025.0um.root",particles[i], angle_ttc), "Fe_thick_Medley", 25, particles[i], false);
                histograms[i][j] = h;

                h->SetName(Form("%s",histoname.c_str()));
                h->SetTitle(Form("%s",histotitle.c_str()));
            }

            h_Ecorr->Sumw2();
            h->Sumw2();
            for(Int_t b = 1; b<= h->GetNbinsX();b++){
                if (angle < 90.0) {
                    h->SetBinContent(b,h->GetBinContent(b)*FactorF);
                    h_Ecorr->SetBinContent(b,h_Ecorr->GetBinContent(b)*FactorF);

                    h->SetBinError(b,h->GetBinError(b)*FactorF);
                    h_Ecorr->SetBinError(b,h_Ecorr->GetBinError(b)*FactorF);
                }else{
                    h->SetBinContent(b,h->GetBinContent(b)*FactorB);
                    h_Ecorr->SetBinContent(b,h_Ecorr->GetBinContent(b)*FactorB);

                    h->SetBinError(b,h->GetBinError(b)*FactorB);
                    h_Ecorr->SetBinError(b,h_Ecorr->GetBinError(b)*FactorB);
                }

            }
        }
    }


    if(Match_CORR) return histograms_Ecorr;
    return histograms;

}

