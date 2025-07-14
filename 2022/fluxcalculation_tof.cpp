#include "src/somefunctions.cc"
#include "include/somefunctions.hh"
#include <any>

//////
// This script is in its vertions version1.2025-07-14.0

void fluxcalculation_tof() {
    //this function will calculate the flux for the 2022 data using the time of flight method.

    std::cout << "Flux calculation function called." << std::endl;

    TFile *fcorr = new TFile("pcorr.root", "READ");
    TGraph *pcorr = (TGraph*)fcorr->Get("pcorr");

    Int_t runA = 35; // Starting run number
    Int_t runB = 38; // Ending run number
    cout<<"Loading runs from " << runA << " to " << runB << endl;
    double Qtot = 0; // Total charge variable
    Int_t Ttot = 0; // Total charge variable
    TChain *tx = loadRuns(runA, runB, &Qtot,&Ttot); // Load runs from 63 to 63

    cout<<"The total charge for runs " << runA << " to " << runB << " is: " << Qtot << " µC." << endl;
    
    loadCuts(1); // Load cuts for telescope 1
    setGuesses(0.000402, 0.002275, 0.0016942782); // Set default guesses for g1, g2, g3
    setEnergyAliases(1); // Set energy aliases for telescope 1
    setParticlesAliases(1); // Set particle aliases for telescope 1
    
    TCutG *cutEE;
    cout<<" Reading cut for EEBand from specialSelections/eeband.C..."<<endl;
    gROOT->ProcessLine(".L  specialSelections/eeband.C");
    cutEE = (TCutG *)gROOT->GetListOfSpecials()->FindObject("eeband");
    if(cutEE!=NULL) 
        cout << "  --> OK."<< endl;
    else 
        cout<< "  --> ERROR: EEBand cut not found." << endl;

    tx->SetAlias("itof", "Medley_1_dE2_ToF"); // Set alias for total energy
    cout<<"Defining alias itof = Medley_1_dE2_ToF"<<endl;
    Float_t gflash = GetGflash(tx, "Medley_1_dE2_ToF", 400); // Get Gflash for the total energy
    cout << "Gflash value: " << gflash << endl;

    //using aliases:

    tx->SetAlias("tof", Form("%f - Medley_1_dE2_ToF + %f",gflash,provideTgama(1))); // Set alias for total energy
    cout<<"Defining alias tof = "<<gflash<<" - Medley_1_dE2_ToF + "<<provideTgama(1)<<" ns"<<endl;
    
    double dT1 = getDistanceTel(1);// in meters
    
    cout<<".\nDistance for telescope 1: "<<dT1<<" m.\n."<<endl;
    tx->SetAlias("tofn", Form("tof - ToFparticleNumber(E,1,%f)", dT1)); // Set alias for time of flight for the neutron
    cout<<"Defining alias tofn = tof - ToFparticleNumber(E,1," << dT1 << ")" << endl;
    tx->SetAlias("En", "En(tofn)"); // Set alias for time of flight minus particle number for second particle
    cout << "Defining alias En = En(tofn)" << endl;
    
    //Create a new TTree with the events I am interested in:
    const char * outputFileName  = Form("runs_%03d-%03d.root",runA, runB);
    cout<<"--> OUTPUT FILENAME: "<<outputFileName<<"."<<endl;

    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
    TTree *newtree = new TTree("M","Medley's EE data TTree");

    //Int_t particle; 
    
    Float_t energy; 
    
    Float_t si1; 
    Float_t si2; 
    Float_t csi; 
    Float_t rawtof; // raw ToF
    Float_t tof_measured;
    Float_t tof_correction; 
    Float_t tofnn; 
    Float_t tofpart; 
    Float_t tofnn_woc;// without correction 
    
    Float_t ENN;
    Float_t ENNwoc; 

   // TBranch * pidBranch  = newtree->Branch("PID", &particle, "PID/I");
    TBranch * EtotBranch  = newtree->Branch("E", &energy, "energy/F");
    
    TBranch * rawtof_Branch = newtree->Branch("rawtof", &rawtof, "rawtof/F");
    TBranch * tof_measured_Branch = newtree->Branch("tof_measured", &tof_measured, "tof_measured/F");
    TBranch * tof_correction_Branch  = newtree->Branch("tof_correction", &tof_correction, "tof_correction/F");
    
    
    TBranch * tofnn_Branch  = newtree->Branch("tofnn", &tofnn, "tofnn/F");
    TBranch * tofpart_Branch  = newtree->Branch("tofpart", &tofpart, "tofpart/F");
    TBranch * tofnn_woc_Branch  = newtree->Branch("tofnn_woc", &tofnn_woc, "tofnn_woc/F");
    
    TBranch * EnBranch  = newtree->Branch("ENN", &ENN, "ENN/F");
    TBranch * EnWocBranch  = newtree->Branch("ENNwoc", &ENNwoc, "ENNwoc/F");
    
    //debbugging part - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    TBranch * branch_si1  = newtree->Branch("si1", &si1, "si1/F");
    TBranch * branch_si2  = newtree->Branch("si2", &si2, "si2/F");
    TBranch * branch_csi  = newtree->Branch("csi", &csi, "csi/F");
    
    //information
    string cur_time = getCurrentTime();
    cout<<"Current time: "<<cur_time<<endl;
    TNamed *processing_info= new TNamed("processed by 2022/fluxcalculation_tof.cpp:version1.2025-07-14.0",cur_time);
    TNamed *charge_in_uC = NULL;
    TNamed *MedleyTarget= NULL;
    TNamed *duration_of_the_run_in_sec= NULL;
    TNamed *configuration_telescopes= NULL;


    std::optional<MedleyData> rundata =  GetRunData("/mnt/medley/LucasAnalysis/2022/runlist2022.csv",runA);
    cout<<"Run data for run "<<runA<<": "<<endl;
    if(rundata.has_value()){
        cout<<"RunN: "<<rundata->RunN<<endl;
        cout<<"MedleyConfig: "<<rundata->MedleyConfig<<endl;
        cout<<"Target: "<<rundata->Target<<endl;
        cout<<"RunTime: "<<Ttot<<endl;
        cout<<"RunCharge: "<<Qtot<<endl;
        MedleyTarget = new TNamed("MedleyTarget", rundata->Target);
        duration_of_the_run_in_sec = new TNamed("duration_of_the_run_in_sec", std::to_string(Ttot));
        configuration_telescopes = new TNamed("configuration_telescopes", rundata->MedleyConfig);
        charge_in_uC = new TNamed("charge_in_uC", std::to_string(Qtot));
        charge_in_uC->Write();
        MedleyTarget->Write();
        duration_of_the_run_in_sec->Write();
    }
    processing_info->Write();
    
    

    Float_t itof = 0;
    

    int telN = 1 ;
    Long64_t nentries = tx->GetEntries();
    cout << "Total entries in the chain: " << nentries << endl;

    double tgamma = provideTgama(telN); // Get the gamma value for the telescope
    cout << "Tgamma for telescope " << telN << ": " << tgamma << " ns." << endl;

    

    TH2D *hh = new TH2D("hh", "hh", 500, 0, 500, 400, 0, 40);
    //for (Long64_t i = 0; i < nentries/50; i++) {

    int p=0;
    for (Long64_t i = 0; i < 2000000; i++) {

        tx->GetEntry(i);
        itof = tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue();
        
        if(itof > gflash || itof ==0) continue; //if the event is faster than gamma... go to next event
      
                //cout<<i<<", "<<p++<<": itof = "<<itof<<", "<<gflash<<endl;

        tof_measured = gflash - tx->GetLeaf(Form("Medley_%d_dE2_ToF",telN))->GetValue() + tgamma; // calculate the time of flight

                //cout<<"   tof : "<<tof<<" ns . . . tofNN = ";
        
        si1 = tx->GetLeaf(Form("Medley_%d_dE1",telN))->GetValue()*g1;
        si2 = tx->GetLeaf(Form("Medley_%d_dE2",telN))->GetValue()*g2;
        csi = tx->GetLeaf(Form("Medley_%d_Eres",telN))->GetValue()*g3;
        energy = si1 + si2 + csi;
        tofpart = ToFparticleNumber(energy, 1, dT1);
        tofnn_woc = tof_measured - tofpart; // calculate the time of flight for the neutron
    
        tof_correction = pcorr->Eval(energy);

        tofnn = tofnn_woc + tof_correction; // calculate the time of flight for the neutron
                // cout<<tofnn_woc<<" ns / "<<tofnn<<" ns*"<<endl;
        
                // cout<<" . . .: "<<si1<<", "<<si2<<", "<< csi<<" =>"<<energy<<endl;

        //if it has suitable energies:
        if( si1 >= 0.02 && si2 >= 0.02 ){
            //and its a proton:
            if( 
                (cut_protons[telN]->IsInside(si2,si1) && cut_protonsNPTCsI[telN]->IsInside(csi,si2))
                ||
                (cut_protonsPT[telN]->IsInside(si2,si1) && cut_protonsCsI[telN]->IsInside(csi,si2)) 
            ){
                //cout<<"Proton found: dE1 = "<<si1<<", si2 = "<<si2<<", Eres = "<<csi<<", E = "<<energy<<", tof = "<<tof_measured<<", tofn = "<<tofnn<<endl;
                //from the elastic scattering band:
                if(
                    cutEE->IsInside(tofnn_woc,energy)
                ){  

                    ENN = En(tofnn);
                    ENNwoc = En(tofnn_woc);

                    hh->Fill(tofnn, energy);
                    
                }
            
            }

        }


        //if (i/nentries) {
            double progress = 100.0 * i / nentries;
            std::cout << "\rProgress: " << std::fixed << std::setprecision(2)
                    << progress << "%" << std::flush;
        //}


    }
    hh->Draw("colz");
    //pcorr->Draw();




}