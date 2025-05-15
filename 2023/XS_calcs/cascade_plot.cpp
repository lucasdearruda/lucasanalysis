void cascade_plot() {
    const int N = 8;
    TH1D* h[N];
    double yShift = 15.0;  // Espaço entre histogramas
    double angle[N] = {15, 30, 45, 60, 75, 90, 115, 135};

    TCanvas* c = new TCanvas("c", "Cascade", 1000, 600);
    c->SetRightMargin(0.05);

    double maxY = 0;

    for (int i = 0; i < N; ++i) {
        TString name = Form("h_%d", i);
        h[i] = new TH1D(name, Form("Angle %.0f deg", angle[i]), 100, 0, 60);
        for (int j = 0; j < 1000; ++j)
            h[i]->Fill(gRandom->Gaus(20 + i*2, 3));

        h[i]->SetLineColor(i + 2);
        h[i]->SetFillColorAlpha(i + 2, 0.6);
        h[i]->SetStats(0);
        h[i]->Scale(1.0);  // opcional: normaliza

        // Encontra o maior valor
        double locMax = h[i]->GetMaximum();
        if (locMax > maxY) maxY = locMax;
    }

    // Cria hist base invisível para o eixo
    TH2F* frame = new TH2F("frame", "Cascade Plot;Energy (MeV);Shifted Counts", 
                           100, 0, 60, 100, 0, maxY + yShift * N);
    frame->Draw();

    // Desenha cada histograma com deslocamento em Y
    for (int i = 0; i < N; ++i) {
        h[i]->DrawCopy("same hist");
        h[i]->GetXaxis()->SetLabelSize(0);
        h[i]->GetYaxis()->SetLabelSize(0);

        // Aplica deslocamento vertical (em y)
        for (int bin = 1; bin <= h[i]->GetNbinsX(); ++bin) {
            double old = h[i]->GetBinContent(bin);
            h[i]->SetBinContent(bin, old + i * yShift);
        }

        h[i]->Draw("same hist");
    }
}