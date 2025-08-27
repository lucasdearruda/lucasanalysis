#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TGClient.h>
#include <TApplication.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TLatex.h>
#include <TString.h>
#include <TVirtualPad.h>

class FeXSViewer : public TGMainFrame {
private:
    TFile *f;
    int nEnergyBins;
    int currentBin;

    TCanvas *ecanvas;
    TVirtualPad *subPad[9];

    TGTextButton *btnNext;
    TGTextButton *btnPrev;

public:
    FeXSViewer(const char *filename);
    virtual ~FeXSViewer();

    void DrawBin(int bin);
    void OnNext();
    void OnPrev();
};

// Implementações aqui (sem ClassImp)

FeXSViewer::FeXSViewer(const char *filename) : TGMainFrame(gClient->GetRoot(), 1200, 900), currentBin(0) {
    f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        printf("Erro ao abrir o arquivo %s\n", filename);
        return;
    }

    nEnergyBins = 0;
    while(true) {
        TString testName = Form("h_20deg_En_%.1f_%.1f_MeV_bin%d", 25.0 + nEnergyBins * 1.0, 25.0 + (nEnergyBins+1) * 1.0, nEnergyBins);
        if (f->Get(testName)) nEnergyBins++;
        else break;
    }
    if (nEnergyBins == 0) {
        printf("Nenhum histograma encontrado no arquivo.\n");
        return;
    }
    printf("Encontrados %d bins de energia.\n", nEnergyBins);

    ecanvas = new TCanvas("ecanvas", "Fe(n,Xp) Cross Sections", 1200, 900);
    ecanvas->Divide(3,3);

    for (int i = 0; i < 9; i++) {
        subPad[i] = ecanvas->cd(i+1);
    }

    btnPrev = new TGTextButton(this, "&Anterior");
    btnNext = new TGTextButton(this, "&Próximo");
    AddFrame(btnPrev, new TGLayoutHints(kLHintsLeft | kLHintsTop, 5,5,5,5));
    AddFrame(btnNext, new TGLayoutHints(kLHintsRight | kLHintsTop, 5,5,5,5));

    btnPrev->Connect("Clicked()", "FeXSViewer", this, "OnPrev()");
    btnNext->Connect("Clicked()", "FeXSViewer", this, "OnNext()");

    SetWindowName("Visualizador Fe(n,Xp) XS");
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    DrawBin(currentBin);
}

FeXSViewer::~FeXSViewer() {
    if (f) f->Close();
    delete ecanvas;
}

void FeXSViewer::DrawBin(int bin) {
    if (bin < 0) bin = 0;
    if (bin >= nEnergyBins) bin = nEnergyBins - 1;
    currentBin = bin;

    Float_t Ea_zero = 25.0 + bin * 1.0;
    Float_t Eb_zero = 25.0 + (bin + 1) * 1.0;
    TString tag = Form("%.1f_%.1f_MeV_bin%d", Ea_zero, Eb_zero, bin);

    for (int i = 0; i < 9; i++) {
        subPad[i]->Clear();
        subPad[i]->SetGrid();
    }

    for (int i = 0; i < 8; i++) {
        int angle = 20 + 20*i;
        TString hname = Form("h_%ddeg_En_%s", angle, tag.Data());
        TH1D *h = (TH1D*)f->Get(hname);
        if (!h) continue;

        subPad[i]->cd();
        h->SetLineColor(kBlue+2);
        h->SetLineWidth(2);
        h->Draw("hist");
    }

    subPad[8]->cd();
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.1, 0.5, Form("Bin de energia: %d\nE_{a}=%.1f MeV\nE_{b}=%.1f MeV", bin, Ea_zero, Eb_zero));

    ecanvas->Update();
}

void FeXSViewer::OnNext() {
    if (currentBin < nEnergyBins - 1) {
        DrawBin(currentBin + 1);
    }
}

void FeXSViewer::OnPrev() {
    if (currentBin > 0) {
        DrawBin(currentBin - 1);
    }
}

int main(int argc, char **argv) {
    TApplication app("app", &argc, argv);

    if (argc < 2) {
        printf("Uso: %s arquivo.root\n", argv[0]);
        return 1;
    }

    FeXSViewer *viewer = new FeXSViewer(argv[1]);

    app.Run();

    delete viewer;
    return 0;
}
