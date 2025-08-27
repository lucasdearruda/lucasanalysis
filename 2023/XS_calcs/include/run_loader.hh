//
//
// Description:
// Tgis file contains the function loadRuns() which loads the runs from the defined folder, providing also the integrated chage


TChain* loadRuns(Int_t runa = 28, Int_t runb = 30,Float_t *charge = nullptr, string folder = "/mnt/medley/LucasAnalysis/2023/reducedv61"){



Int_t Ta;
TChain *tx = new TChain("M");

string name;
for(Int_t i=runa;i<=runb;i++)
{
    Bool_t fExist = true;
    name = Form("%s/%03d.root",folder.c_str(),i);
    cout << name << endl;
    ifstream mfile;
    mfile.open(name);
    if(mfile)
    {
        mfile.close();
    }
    else
        fExist=false;
    if(fExist)
    {	  
        cout << "Adding " << name << endl;
        tx->Add(name.c_str());
        Ta = (Int_t) (tx->GetEntries());
        cout << "Entries " << Ta/1000 << "k "  << endl;
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
    *charge = TotalChargeFaraday;
}

return tx;

}