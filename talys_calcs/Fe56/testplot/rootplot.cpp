{
    TTree *tt = new TTree("tt", "dados");
    tt->ReadFile("pddxE0067.000A015.0.deg", "E/F:xs/F:dir/F:preeq/F:mult/F:compound/F");

    TGraph *g_dir     = new TGraph();
    TGraph *g_preeq   = new TGraph();
    TGraph *g_mult    = new TGraph();
    TGraph *g_comp    = new TGraph();
    TGraph *g_xs      = new TGraph();

    Float_t E, xs, dir, preeq, mult, compound;
    tt->SetBranchAddress("E", &E);
    tt->SetBranchAddress("xs", &xs);
    tt->SetBranchAddress("dir", &dir);
    tt->SetBranchAddress("preeq", &preeq);
    tt->SetBranchAddress("mult", &mult);
    tt->SetBranchAddress("compound", &compound);

    Long64_t n = tt->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tt->GetEntry(i);
        g_dir->SetPoint(i, E, dir);
        g_preeq->SetPoint(i, E, preeq);
        g_mult->SetPoint(i, E, mult);
        g_comp->SetPoint(i, E, compound);
        g_xs->SetPoint(i, E, xs);
    }

    g_dir->SetLineColor(kRed);
    g_preeq->SetLineColor(kBlue);
    g_mult->SetLineColor(kMagenta);
    g_comp->SetLineColor(kGreen+2);
    g_xs->SetLineColor(kBlack);

    g_dir->SetLineWidth(2);
    g_preeq->SetLineWidth(2);
    g_mult->SetLineWidth(2);
    g_comp->SetLineWidth(2);
    g_xs->SetLineWidth(3);

    g_dir->SetTitle("Variáveis vs Energia;Energia (MeV);Valor");

    g_dir->Draw("AL");
    g_preeq->Draw("L SAME");
    g_mult->Draw("L SAME");
    g_comp->Draw("L SAME");
    g_xs->Draw("L SAME");

    auto leg = new TLegend(0.65, 0.6, 0.88, 0.88);
    leg->AddEntry(g_dir, "dir", "l");
    leg->AddEntry(g_preeq, "preeq", "l");
    leg->AddEntry(g_mult, "mult", "l");
    leg->AddEntry(g_comp, "compound", "l");
    leg->AddEntry(g_xs, "xs (total)", "l");
    leg->Draw();
}

