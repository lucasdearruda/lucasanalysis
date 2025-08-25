//prod_DXX.cpp, version 2025-08-25.0
// I added some TNameds to the output file, with the runs used
// and some TNameds with the input parameters, and target info

//prod_DXX.cpp, version 2025-08-21.0
//I am revising it
// - first thing: add the string for the output file name inside the function


//prod_DXX.cpp, version 2025-08-08.0
// Script for evaluation the ddx for Fe in 2023 runs
// This script accepts multiple runs


#include <iostream>     // for std::cout, std::cerr, std::endl
#include <fstream>      // for std::ifstream
#include <string>       // for std::string


//for the time:
#include <time.h>
//for the vectors and pairs:
#include <vector>
#include <utility> 
#include "TVector3.h"

using namespace std;


//This should be included manually, but it is here for convenience:
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/runLoader.cxx" // for GetRunData
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/functions.cxx" // for my functions 
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/newE.cxx" // for the matching correction 

//////////
//The 'global' variables (like binWidth, mc, M, ...) are defined in functions.cxx for simplicity

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here we will procude DDX for Fe in 2023 runs:: this new version is more organized and save everything in a ROOT file 
// Different than before, many settings will be defined within the script, not passed as arguments
void prod_DDX(
    bool matchingCORR = true,
    bool TTC = true,
    bool plotIt = true, 
    bool saveIt = true, 
    bool pause_each = false,
    bool verbose = true
    ) {
    
    string cur_time = getCurrentTime();
    clock_t tStart = clock();
    string outputFileName = "prod_DDX_thick"  ; 

    //create the output file name based on the parameters:
    if(matchingCORR) {
        outputFileName += "_MC";
    }
    if(TTC) {
        outputFileName += "_TTC";
    }
    outputFileName += ".root";

    //For iron: 
    string mytarget= "Fe_thick_Medley";
    Attribute_Target(mytarget);
        //Attribute_Target("Fe_thin");
    //for carbon:
    //Attribute_Target("MedleyCarbon");
    version();

    //cout<< "If not compiled, remember loading the .CXX files:" << endl;
    //cout<< ".L /mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/runLoader.cxx" << endl;
    //cout<< ".L /mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/functions.cxx" << endl;
    cout<<endl<<endl<<endl<<endl;

    //: : : Defining what we want to get: 
        //std::vector<float> angles =  {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0};

    //std::vector<float> angles =  {20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0};
    std::vector<float> angles =  {20.0};
    
    //std::vector<char> particles = {'p', 'd', 't', 'h', 'a'};
    std::vector<char> particles = { 'p', 'a'};

    
    cout << "Output file name: " << outputFileName << endl;

    //: : : Here we can define using loops os stuff like that::

    //I defined the pairs because we know that the relationship bet time and energi is non-linear
    //and we measured neutrons's energy (calc. from ToF)...
    std::vector<std::pair<float, float>> energy_bins = {
        //{25, 26}//,
        {3.5, 4.5},
        {4.5, 5.5},
        {5.5, 6.5},
        {6.5, 7.5},
        {7.5, 8.5},
        {8.5, 9.5},
        {9.5, 10.5},
        {10.5, 11.5},
        {11.5, 12.5},
        {12.5, 13.5},
        {13.5, 14.5},
        {14.5, 15.5},
        {15.5, 16.5},
        {16.5, 17.6},
        {17.6, 18.8},
        {18.8, 20.1},
        {20.1, 21.6},
        {21.6, 23.2},
        {23.2, 25},
        {25, 27},
        {27, 29.3},
        {29.3, 32},
        {32, 35},
        {35, 38.4},
        {38.4, 40}
    };


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Getting runs and charges
//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

    cout<< "Loading runs and charges..." << endl;  
    
    
    //for Fe-thick: _ _ _ _ _ _ _ _ _ 
    std::vector<int> runsForward = {397,405};
                        //std::vector<int> runsBackward = {383,384}; PROBLEMATIC RUNS (120° very weird)!!!! 
    std::vector<int> runsBackward = {410,412};

    //for Fe-thin: _ _ _ _ _ _ _ _ _ 
    // std::vector<int> runsForward = {389,394};
    //                     //std::vector<int> runsBackward = {383,384}; PROBLEMATIC RUNS!!!! 
    // std::vector<int> runsBackward = {385,387};



    //for carbon
    //std::vector<int> runsForward = {35,39};
    //std::vector<int> runsBackward = {383,384};

    TChain* Fr = nullptr; 
    TChain* Br = nullptr;
    cout<< "Starting prod_DDX..." << endl;

    //quick benchmark
    Float_t charge_fr = 0;
    Float_t charge_br = 0;

    for(size_t i = 0; i + 1 < runsForward.size(); i += 2) {//two by two 
        Fr = loadRuns(runsForward[i], runsForward[i+1], Fr, &charge_fr);
    }

    for(size_t i = 0; i + 1 < runsBackward.size(); i += 2) {//two by two 
        Br = loadRuns(runsBackward[i], runsBackward[i+1], Br, &charge_br);
    }


    //TChain* Fr = loadRuns(397, 405, nullptr, &charge_fr); //ForwardRuns
    //TChain* Br = loadRuns(383, 384, nullptr, &charge_br); //BackwardRuns

    cout<< " _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _ - _" << endl;
    cout << "Total charge for forward runs: " << charge_fr << " µC" << endl;
    cout << "Total charge for backward runs: " << charge_br << " µC" << endl;


//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

// Defining some auxiliary variables and the visualization canvas:

Float_t Ea, Eb; // Ea and Eb are defined based on the binwidth and the energy:: 
Int_t Nbins = 200;
Float_t Factor = 0;
        
string targetmat = "Fe_thick_Medley"; // Target material
TCutG * cut_ee[4];

if(targetmat == "MedleyCarbon" || targetmat == "C" || targetmat == "CMedley") {
    //load the carbon masks...
// I want to exclude    ee1.C ee2.C
//and consider only ee4.C, for this case we have to remove angles for tel3 as well 

    gROOT->ProcessLine(".L  /mnt/medley/LucasAnalysis/2023/XS_calcs/eeCcuts/ee1.C");
    gROOT->ProcessLine(".L  /mnt/medley/LucasAnalysis/2023/XS_calcs/eeCcuts/ee2.C");
    gROOT->ProcessLine(".L  /mnt/medley/LucasAnalysis/2023/XS_calcs/eeCcuts/ee3.C");
    cut_ee[0] = (TCutG *)gROOT->GetListOfSpecials()->FindObject("ee1");
    cut_ee[1] = (TCutG *)gROOT->GetListOfSpecials()->FindObject("ee2");
    //cut_ee[2] = (TCutG *)gROOT->GetListOfSpecials()->FindObject("ee3"); THERE IS NO ee3 FOR CARBON
    cut_ee[3] = (TCutG *)gROOT->GetListOfSpecials()->FindObject("ee4");
}

std::vector<std::vector<std::vector<TH1D*>>> Hddx; // vector for TH1D for each particle and angle 

TCanvas *c;
if(plotIt) c =  new TCanvas("c","c",800,600);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Process everything, looping over energies, angles and particles
//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.

cout<<"..L = " << L << " cm"<<endl;;
cout<<"..A_tgt = " << A_tgt <<" cm²"<<endl;
cout<< "..Nc = " << Nc << " atoms/cm³"<<endl;
cout<<"..omega_tel = " << omega_tel << " sr"<<endl;
cout<<"...charge_fr = " << charge_fr << " µC."<<endl;
cout<<"...charge_br = " << charge_br << " µC."<<endl;


for (const auto& bin : energy_bins) {
    Float_t Ea = bin.first;
    Float_t Eb = bin.second;
    Float_t En = (Ea + Eb) / 2.0; // Energy in MeV
    Float_t flux = getNflux(Ea, Eb); // Get the flux for the given energy range, in n/sr/µC/MeV - from functions.cxx
    if(verbose)cout << "Processing energy: " << En << " MeV" << endl;
    if(verbose)cout << " ---> Ea: " << Ea << " MeV, Eb: " << Eb << " MeV" << endl;
    if(verbose)cout << "   ---> neutron flux: " << flux << " n/sr/µC/MeV" << endl;

    std::vector<std::vector<TH1D*>> angle_vec; // angle vector


    //Now we can loop over angles and particles
    for(float angle : angles){

        
        std::vector<TH1D*> particle_vec; //vector for particles
      
        if(angle<100){
            Factor = 1e3*L*L*cm2_to_barn/(flux*Nc*omega_tel*charge_fr); // *TMath::Cos(45*TMath::DegToRad())
        }else{
            Factor = 1e3*L*L*cm2_to_barn/(flux*Nc*omega_tel*charge_br);
        }
        if(verbose)cout<< " ** Factor = " << Factor << endl;
        
        
        for(char particle : particles){
            
            getAZ(particle); //update Z and A based on the particle
            
            cout <<"Processing angle: " << angle << " deg, particle: " << particle << endl;
            
            // Here we call the function that does the actual work       

            //Create the name and title for the given histo.. 
            TString hname = Form("P%c_E%.1f_A%.1fdeg", particle, En, angle);
            hname.ReplaceAll(".", "p");

            TString htitle = Form("E_{NN}=%.1f-%.1f MeV, ang=%.0fdeg, %s", Ea, Eb, angle, particleName(particle).c_str());

            // Create the histogram to be saved (has to be created out of the if): 
            TH1D* hist = new TH1D(hname.Data(), htitle.Data(), Nbins, 0, 40); 
            //TH1D* hist = nullptr;
            
            //conditions for the cuts:
            string conditions = Form("ENN>%f && ENN<%f && PID==%d && ang == %f", Ea, Eb, pCode(particle), angle);
            cout << "  - Conditions: " << conditions << endl;


            if(
                (targetmat == "MedleyCarbon" || targetmat == "C" || targetmat == "CMedley")
                && (angle == 20 || angle == 40 ||  angle == 80 || angle == 100 || angle == 120 || angle == 160)
                ) {
                //if we are using the carbon target, we have to apply the cuts:
                if(angle<=80)
                    conditions = Form("ENN>%f && ENN<%f && PID==%d && ang == %f && ee%d", Ea, Eb, pCode(particle), angle,(int)angle/20);
                else
                    conditions = Form("ENN>%f && ENN<%f && PID==%d && ang == %f && ee%d", Ea, Eb, pCode(particle), angle, (int)(180 - angle)/20);
                cout << "  - [Carbon] -->> Corrected conditions: " << conditions << endl;
            }

            //if p d or t && match correction ON: 
            if(matchingCORR && (particle == 'p' || particle == 'd' || particle == 't')) { //if we need to apply the matching correction:

                TH1D* histSi1 = new TH1D(Form("hsi1_%c_%.1f",particle,angle), Form("hsi1_%c_%.1f_title",particle,angle), 20*Nbins, 0, 50); //create the histogram for si1

                if(angle<=80){
                    //Fr->Draw("si1>>hsi1", Form("ENN>%f && ENN<%f && PID==%d && ang == %f", Ea, Eb, pCode(particle), angle));
                    Fr->Draw(Form("si1>>hsi1_%c_%.1f",particle,angle), conditions.c_str(),"goff");
                }else{
                    //Br->Draw("si1>>hsi1", Form("ENN>%f && ENN<%f && PID==%d && ang == %f", Ea, Eb, pCode(particle), angle));
                    Br->Draw(Form("si1>>hsi1_%c_%.1f",particle,angle), conditions.c_str(),"goff");
                }
                
                //hsi1 constructed, now we can apply the newE function to the histogram
                
                hname = string(hname.Data()) + "_MC";
                htitle = string(htitle.Data()) + ", MC";
                // cout <<"\n-- -- \n -- -- \n -- -- -- -- \n -- -- -- -- -- \n --- -- -- -- -- -- \n";
                // cout<<"Info :: particle = "<< particle << ", angle= " << angle << ", A = " << A << ", Z = " << Z << endl;
                // cout<< "Conditions:"<<conditions<<endl;
                // cout << "Thickness first tel:"<<thicknessFirstTel(int(angle / 20))<<endl;

                //now we attribute

                Float_t thSi1 = thicknessFirstTel(int(angle / 20));//µm
                hist = newE(histSi1, hname.Data(), htitle.Data(), thSi1, Z, A); // Apply the newE function to the histogram
                hist->Rebin(20); // Rebin the histogram to 20 bins
                //hist->SetNameTitle(hname.Data(), htitle.Data()); // Set the name and title of the histogram

                //gDirectory->Remove(histSi1); // Remove the histogram from the directory
                //delete histSi1; // Delete the histogram to free memory

            }else{
                
                if(angle<=80){
                    //Fr->Draw(Form("E>>%s", hname.Data()), Form("ENN>%f && ENN<%f && PID==%d && ang == %f", Ea, Eb, pCode(particle), angle));
                    Fr->Draw(Form("E>>%s", hname.Data()), conditions.c_str());
                }else{
                    //Br->Draw(Form("E>>%s", hname.Data()), Form("ENN>%f && ENN<%f && PID==%d && ang == %f", Ea, Eb, pCode(particle), angle));
                    Br->Draw(Form("E>>%s", hname.Data()), conditions.c_str());
                }
                
            }


            //required transformations on hist:: 
            if(TTC){
                if(verbose)cout<< "Applying TTC correction to " << hist->GetName() <<endl;
                hist = correctSpec(hist, true,particle, angle, targetmat); //from functions.cxx
                hname = string(hname.Data()) + "_TTC";
                htitle = string(htitle.Data()) + ", TTC";
                hist->SetNameTitle(hname.Data(), htitle.Data());
            }


            //visualization hitso-by-histo:
            if(pause_each) {
                c->cd();
                hist->Draw();
                c->Update();
                cout << "Press Enter to continue..." << endl;
                cin.get(); // Wait for user input
            }

            //////////////////////////////////
            //the last part is converting the histograms of particle counts to DDX:
            //I defined before the variables L, A_tgt and cm2_to_barn in functions.cxx
            //
            // dσ/dΩdE (En,Ep) = (1/binW_part)) * Np(Ep)L²/( φ(En) * Nc *ΔΩ_tel * Q * cos ζ) * [b/cm²],
            //
            // ζ = target angle
            // [φ(En)] = n/sr/µC/MeV
            //
            
            for(int bn = 1;bn<=hist->GetNbinsX();bn++){
                Double_t binContent = hist->GetBinContent(bn);
                Double_t binWidth = hist->GetBinWidth(bn);
                if(binContent > 0) {
                    hist->SetBinContent(bn, binContent * Factor / binWidth);
                } else {
                    hist->SetBinContent(bn, 0); // Avoid division by zero
                }
            }
            

            particle_vec.push_back(hist);
        }
        angle_vec.push_back(particle_vec);
    }

    Hddx.push_back(angle_vec);

}
 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Create and FILL file output...
//__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.


TFile *fOut = new TFile(outputFileName.c_str(), "RECREATE");
TNamed *processing_info= new TNamed("File produced by prod_DDX.cpp:",Form("version2025-08-25.0 in %s",cur_time.c_str()));


std::string forwardStr = makeRunsString(runsForward);
std::string backwardStr = makeRunsString(runsBackward);

TNamed *runs_infoF = new TNamed("Forward runs:", forwardStr.c_str());
TNamed *runs_infoB = new TNamed("Backward runs:", backwardStr.c_str());


// Loop for saving the histograms::0
for (const auto& angle_vec : Hddx) {
    for (const auto& particle_vec : angle_vec) {
        for (const auto& hist : particle_vec) {
            if(hist) hist->Write();
        }
    }
}

 // --- Create a TNamed for each flag ---
    TNamed* flag_matchingCORR = new TNamed("matchingCORR", matchingCORR ? "true" : "false");
    TNamed* flag_TTC         = new TNamed("TTC", TTC ? "true" : "false");
    TNamed* flag_plotIt      = new TNamed("plotIt", plotIt ? "true" : "false");
    TNamed* flag_saveIt      = new TNamed("saveIt", saveIt ? "true" : "false");
    TNamed* flag_pause_each  = new TNamed("pause_each", pause_each ? "true" : "false");
    TNamed* flag_verbose     = new TNamed("verbose", verbose ? "true" : "false");
    TNamed* target_named     = new TNamed("Target", mytarget.c_str());



flag_pause_each->Write();
flag_matchingCORR->Write();
flag_TTC->Write();
flag_plotIt->Write();
flag_saveIt->Write();
flag_verbose->Write();
target_named->Write();

runs_infoF->Write();
runs_infoB->Write();
processing_info->Write();
fOut->Write();
fOut->Close();
cout << "Output file " << outputFileName << " created." << endl;

cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
cout<<"prod_DXX.cpp, version 1.2025-08-25.0"<<endl;

    


}