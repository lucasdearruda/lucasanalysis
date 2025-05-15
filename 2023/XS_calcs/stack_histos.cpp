void stack_histos() {
    const int N = 8;
    TH1D* h[N];

    // Criar histos de exemplo
    for (int i = 0; i < N; ++i) {
        TString name = Form("h_%d", i);
        h[i] = new TH1D(name, name, 100, 0, 60);
        for (int j = 0; j < 1000; ++j)
            h[i]->Fill(gRandom->Gaus(10 + i*5, 2));
        h[i]->SetFillColor(i + 2);
    }

    TCanvas* c1 = new TCanvas("c1", "Stacked Histos", 800, 600);
    THStack* hs = new THStack("hs", "Stacked Histograms");

    for (int i = 0; i < N; ++i) {
        hs->Add(h[i]);
    }

    hs->Draw("hist");
    hs->GetXaxis()->SetTitle("Energy (MeV)");
    hs->GetYaxis()->SetTitle("Counts");
}