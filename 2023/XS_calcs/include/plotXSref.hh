

TGraphErrors *myXS(string namefile = "/mnt/medley/LucasAnalysis/exfor/exfor_files/datasetsC.root", Float_t Ang = 20, string ref = "12838.002", Float_t EN = 27.4){


TChain *t = new TChain("DataTree");
t->Add(namefile.c_str());
t->Print();

// TLeaf* leaf = t->GetLeaf("reference");
//  cout << "Type of 'Reference': " << leaf->GetTypeName() << endl;
//  cout << "Value of 'Reference': " << leaf->GetValue() << endl;

Long64_t nentries = data->GetEntries();
for(Long64_t jentry=0; jentry<nentries; jentry++){
    t->GetEntry(jentry);
    //cout << t->GetLeaf("Ang")->GetValue()<< " "<<t->GetLeaf("EN")->GetValue()<<endl; 
    if(t->GetLeaf("Ang")->GetValue() == Ang && t->GetLeaf("EN")->GetValue() == EN && std::strcmp(t->GetLeaf("reference")->GetValue(), ref.c_str()) == 0){
        cout<<"Found the entry: "<<jentry<<endl;
        cout << t->GetLeaf("Ang")->GetValue()<< " "<<t->GetLeaf("EN")->GetValue()<<endl<<endl; 
        //break;
    }
}


return nullptr;

}