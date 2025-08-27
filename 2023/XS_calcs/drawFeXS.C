#include <TFile.h>
#include <TKey.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TList.h>
#include <vector>
#include <string>
#include <iostream>

std::vector<std::string> binTags;
int currentBin = 0;
TCanvas* cView = nullptr;
TFile* f = nullptr;

void drawBin(int binIndex) {
    if (!f || binIndex < 0 || binIndex >= (int)binTags.size()) return;

    currentBin = binIndex;
    TString tag = binTags[binIndex].c_str();

    TPaveText* status = (TPaveText*)cView->FindObject("statusBox");
    if (status) {
        status->Clear();
        status->AddText(Form("Energia bin %d - %s", binIndex, tag.Data()));
        status->Draw();
    }

    for (int j = 0; j < 8; ++j) {
        TString hname = Form("h_%ddeg_En_%s", 20 + 20 * j, tag.Data());
        TH1D* h = (TH1D*)f->Get(hname);
        if (!h) continue;
        cView->cd(j + 1);
        gPad->SetGrid();
        h->SetLineColor(kBlue + 2);
        h->SetLineWidth(2);
        h->Draw("hist");
    }
    cView->cd(0);
    cView->Update();
}

void nextBin() {
    drawBin((currentBin + 1) % binTags.size());
}

void prevBin() {
    drawBin((currentBin - 1 + binTags.size()) % binTags.size());
}

void drawFeXS(const char* filename = "Fe_pXS.root") {
    f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }

    binTags.clear();
    TIter next(f->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TString name = key->GetName();
        if (name.BeginsWith("h_20deg_En_")) {
            TString tag = name;
            tag.ReplaceAll("h_20deg_En_", "");
            binTags.push_back(tag.Data());
        }
    }

    if (binTags.empty()) {
        std::cerr << "Nenhum histograma encontrado." << std::endl;
        return;
    }

    std::cout << "Bins energéticos encontrados: " << binTags.size() << std::endl;
    std::cout << "Use nextBin() e prevBin() para navegar." << std::endl;

    cView = new TCanvas("cView", "Fe(n,Xp) viewer", 1200, 800);
    cView->Divide(3, 3);

    cView->cd(0);
    TPaveText* status = new TPaveText(0.3, 0.94, 0.7, 0.99, "NDC");
    status->SetName("statusBox");
    status->SetFillColor(0);
    status->SetTextFont(42);
    status->AddText("Fe(n,Xp) Navegação Interativa");
    status->Draw();

    drawBin(0);
}
