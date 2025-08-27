

TGraphErrors *plotXS(Double_t ENa = 25, Double_t ENb = 26, Double_t setAng = 20, string filename = "datasets.root", bool plot = false){

TChain *txs = new TChain("DataTree");
txs->Add(filename.c_str());

Float_t EN, E, XS, XSerr, Ang;
txs->SetBranchAddress("EN", &EN);
txs->SetBranchAddress("E", &E);
txs->SetBranchAddress("Ang", &Ang);
txs->SetBranchAddress("XS", &XS);
txs->SetBranchAddress("XSerr", &XSerr);

Long64_t nEntries = txs->GetEntries();

TGraphErrors *graph = new TGraphErrors();
for (size_t i = 0; i < nEntries; i++)
{
    txs->GetEntry(i);
    if (Ang == setAng && EN >= ENa && EN <= ENb)
    {
        graph->SetPoint(graph->GetN(), E, XS);
        graph->SetPointError(graph->GetN() - 1, 0, XSerr);
    }
}

graph->SetMarkerStyle(8);
graph->SetMarkerColor(kBlack);
graph->SetLineColor(kBlack);
graph->SetLineWidth(2);

if (plot)
{
    TCanvas *xs_cv = new TCanvas("XS", "XS");
    graph->Draw("ALP");
}


return graph;

}