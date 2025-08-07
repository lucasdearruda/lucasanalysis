//
//
// Description:
// Tgis file contains the function loadRuns() which loads the runs from the defined folder, providing also the integrated chage



// Here I will load the runs in the TChain tx
// and return the TChain tx, and the integrated charge in µC
//If there were runs previouslu loaded, we will just update the TChain with the new runs
TChain* loadRuns(Int_t runa = 28, 
                Int_t runb = 30,
                TChain* tx = nullptr,
                Float_t *charge = nullptr, 
                string folder = "/mnt/medley/LucasAnalysis/2023/reducedv7"){



    
    if (tx == nullptr) {
        tx = new TChain("M");
    } else {
        cout << "Updating existing TChain." << endl;
    }

    
    string name;

    for(Int_t i = runa; i <= runb; i++) {
        name = Form("%s/%03d.root", folder.c_str(), i);
        cout << name << endl;

        ifstream mfile(name);
        if (mfile) {
            mfile.close();
            cout << "Adding " << name << endl;


            Int_t before = tx->GetEntries();
            tx->Add(name.c_str());
            Int_t after = tx->GetEntries();
            Int_t added = after - before;

            cout << "Entries = +" << added / 1000 << "k " << endl;
            cout << " . . .  =  " << after / 1000 << "k " << endl;
        } else {
            cout << "File not found: " << name << endl;
            cout << "Skipping..." << endl;
        }
    
    
}

cout<<"---------------------------------------------------"<<endl;
TTree *InfoTree = new TTree("InfoTree", "InfoTree");
InfoTree->ReadFile("/mnt/medley/LucasAnalysis/2023/runlist.csv", "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F");

Int_t RunN;
Float_t ChargeFaraday;
Float_t TotalChargeFaraday = 0;

InfoTree->SetBranchAddress("RunN", &RunN);
InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

for(int i=0;i<InfoTree->GetEntries();i++)
{
    InfoTree->GetEntry(i);
    //cout<< i<<", nrun = "<<RunN<< endl;
    if (RunN<= runb && RunN>=runa)
    {
        //InfoTree->GetEntry(i);
        cout << "RunN: " << RunN << " ChargeFaraday: " << ChargeFaraday << endl;
        TotalChargeFaraday += ChargeFaraday;
    }
}
cout<<"---------------------------------------------------"<<endl;
cout<< "Total Charge Faraday: " << TotalChargeFaraday <<" µC"<< endl;
if (charge != nullptr) {
    *charge += TotalChargeFaraday;
}else{
    cout<< "No charge pointer provided, not updating charge." << endl;
}

return tx;

}