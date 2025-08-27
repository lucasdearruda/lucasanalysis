#include <TTree.h>
#include <TGraphErrors.h>
#include <string>
#include <cmath> // para fabs

TGraphErrors* readExforGraphForAngle(std::string filename, float targetAngle) {
    TTree* txs = new TTree("txs", "txs");
    txs->ReadFile(filename.c_str(), "EN/F:E/F:Ang/F:XS/F:XSerr/F");

    float EN, E, Ang, XS, XSerr;
    txs->SetBranchAddress("EN", &EN);
    txs->SetBranchAddress("E", &E);
    txs->SetBranchAddress("Ang", &Ang);
    txs->SetBranchAddress("XS", &XS);
    txs->SetBranchAddress("XSerr", &XSerr);

    TGraphErrors* graph = new TGraphErrors();

    int nEntries = txs->GetEntries();
    for (int i = 0; i < nEntries; ++i) {
        txs->GetEntry(i);
        if (std::fabs(Ang - targetAngle) < 0.5) {
            int n = graph->GetN();
            graph->SetPoint(n, E, XS);
            graph->SetPointError(n, 0, XSerr);
        }
    }

    graph->SetTitle(Form("EXFOR - %.0f°", targetAngle));
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.3);
    graph->SetMarkerColor(kBlue+1);
    graph->SetLineColor(kBlue+1);
    graph->GetXaxis()->SetTitle("E_{p} (MeV)");
    graph->GetYaxis()->SetTitle("d#sigma/d#Omega dE (mb/sr/MeV)");

    return graph;
}
