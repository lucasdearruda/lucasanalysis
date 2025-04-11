void plotMeAxs(string namefile = "", bool same = true){

    if(namefile == ""){
        cerr<<"Error: No file name provided."<<endl;
        return;
    }else{


        TTree *tx = new TTree("xs","xs");
        //tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MPreeq/F:Compound/F:PreeqR/F:BKr/F:Str/F:KO/F:BrU/F");
        //##     E-out           xs           Direct     Preequilibrium Multiple_preeq    Compound
        tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MultPreeq/F:Compound/F");
        Float_t Eprod, total;
        Long64_t nEntries = tx->GetEntries();
        tx->SetBranchAddress("Eprod", &Eprod);
        tx->SetBranchAddress("Total", &total);


        TGraph *gr = new TGraph();
        for(Long64_t i = 0; i < nEntries; ++i) {
            tx->GetEntry(i);

            cout<<Eprod<<" "<<total<<endl;
            gr->AddPoint(Eprod, total);
            // Plotting code here
            // For example:
            // hist->Fill(Eplot, total);
        }
        if(same){
            gr->Draw("same");
        }else{
            gr->Draw();
        }

    }


}