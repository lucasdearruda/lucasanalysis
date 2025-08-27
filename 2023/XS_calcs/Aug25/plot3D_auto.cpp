#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TKey.h>
#include <TClass.h>
#include <TString.h>
#include <iostream>
#include <map>
#include <vector>


#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLegend.h>


#include <TH1D.h>
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TStyle.h>

void plot3D_auto() {
    // Ajuste aqui se quiser trocar partícula
    char particle = 'p'; // 'p' -> Pp (prótons), 'a' -> Pa (alfas)
    
    // Abre arquivo ROOT
    TFile *f = TFile::Open("prod_DDX_thick_MC.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erro ao abrir arquivo ROOT!" << std::endl;
        return;
    }

    // Mapeia energias -> lista de histogramas por ângulo
    std::map<TString, std::vector<TString>> energyMap;

    // Percorre todos os objetos do arquivo
    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*) next())) {
        TObject *obj = key->ReadObj();
        if (!obj->InheritsFrom(TH1D::Class())) continue;

        TString hname = obj->GetName();

        // Checa se é do tipo desejado (ex.: "Pp_E..._A..._MC")
        if (!hname.BeginsWith(Form("P%c_", particle))) continue;

        // Extrai parte da energia (ex.: de "Pp_E4p0_A20p0deg_MC" vira "E4p0")
        Ssiz_t posE = hname.Index("_E");
        Ssiz_t posA = hname.Index("_A");
        if (posE == kNPOS || posA == kNPOS) continue;

        TString energy = hname(posE+1, posA-posE-1); // pega "E4p0"

        // Adiciona ao mapa
        energyMap[energy].push_back(hname);
    }

    // Agora percorremos cada energia e montamos o TH2D
    for (auto &entry : energyMap) {
        TString energy = entry.first;
        auto &histNames = entry.second;

        if (histNames.empty()) continue;

        // Usa o primeiro histograma como referência para os bins de energia
        TH1D *h_ref = (TH1D*) f->Get(histNames[0]);
        if (!h_ref) continue;

        int nBinsX = h_ref->GetNbinsX();
        double xMin = h_ref->GetXaxis()->GetXmin();
        double xMax = h_ref->GetXaxis()->GetXmax();

        // Cria o TH2D com eixo Y = ângulo (180 bins de 1° cada)
        TH2D *h2 = new TH2D(Form("h2_%s", energy.Data()),
                            Form("Fe(n,X%c), %s;Energy (MeV);Angle (deg);d^{2}#sigma/d#Omega dE", particle, energy.Data()),
                            nBinsX, xMin, xMax,
                            180, 0, 180);

        // Preenche com todos os histogramas dessa energia
        for (auto &hname : histNames) {
            TH1D *h = (TH1D*) f->Get(hname);
            if (!h) continue;

            // Extrai ângulo do nome
            Ssiz_t posA = hname.Index("_A");
            Ssiz_t posDeg = hname.Index("deg");
            if (posA == kNPOS || posDeg == kNPOS) continue;
            TString angleStr = hname(posA+2, posDeg - (posA+2)); // ex: "20p0"
            double angle = angleStr.Atof(); // converte p0 -> float

            int iy = h2->GetYaxis()->FindBin(angle); // acha o bin certo no eixo Y

            for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
                double value = h->GetBinContent(ix);
                double x = h->GetXaxis()->GetBinCenter(ix);
                h2->SetBinContent(h2->GetXaxis()->FindBin(x), iy, value);
            }
        }

        // Desenha
        TCanvas *c = new TCanvas(Form("c_%s", energy.Data()), "3D plot", 900, 700);
        gStyle->SetOptStat(0);
        h2->SetStats(0);
        h2->Draw("lego2z"); // opções: "lego2", "surf1", "colz"
    }
}


void plotGraphs_auto() {
    char particle = 'p'; // 'p' -> Pp (prótons), 'a' -> Pa (alfas)
    
    TFile *f = TFile::Open("prod_DDX_thick_MC.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erro ao abrir arquivo ROOT!" << std::endl;
        return;
    }
  // Agrupa histogramas por energia
    std::map<TString, std::vector<TString>> energyMap;

    TIter next(f->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*) next())) {
        TObject *obj = key->ReadObj();
        if (!obj->InheritsFrom(TH1D::Class())) continue;

        TString hname = obj->GetName();
        if (!hname.BeginsWith(Form("P%c_", particle))) continue;

        // Extrai energia
        Ssiz_t posE = hname.Index("_E");
        Ssiz_t posA = hname.Index("_A");
        if (posE == kNPOS || posA == kNPOS) continue;

        TString energy = hname(posE+1, posA-posE-1);
        energyMap[energy].push_back(hname);
    }

    // Para cada energia, desenha vários TPolyLine3D com um eixo TH3D de fundo
    for (auto &entry : energyMap) {
        TString energy = entry.first;
        auto &histNames = entry.second;

        if (histNames.empty()) continue;

        // Primeiro, varrer para achar ranges globais
        double xmin = 1e9, xmax = -1e9;
        double zmin = 1e9, zmax = -1e9;
        double ymin = 1e9, ymax = -1e9;

        for (auto &hname : histNames) {
            TH1D *h = (TH1D*) f->Get(hname);
            if (!h) continue;

            // Extrai ângulo
            Ssiz_t posA = hname.Index("_A");
            Ssiz_t posDeg = hname.Index("deg");
            if (posA == kNPOS || posDeg == kNPOS) continue;
            TString angleStr = hname(posA+2, posDeg - (posA+2));
            double angle = angleStr.Atof();

            ymin = std::min(ymin, angle);
            ymax = std::max(ymax, angle);

            xmin = std::min(xmin, h->GetXaxis()->GetXmin());
            xmax = std::max(xmax, h->GetXaxis()->GetXmax());

            for (int ix = 1; ix <= h->GetNbinsX(); ix++) {
                double z = h->GetBinContent(ix);
                zmin = std::min(zmin, z);
                zmax = std::max(zmax, z);
            }
        }

        // Cria um TH3D vazio só para os eixos
        TH3D *frame = new TH3D("frame", 
            Form("Fe(n,X%c), %s;Energy (MeV);Angle (deg);d^{2}#sigma/d#Omega dE", 
                 particle, energy.Data()),
            10, xmin, xmax,
            10, ymin-5, ymax+5,
            10, zmin, zmax*1.1);

        TCanvas *c = new TCanvas(Form("c_graph3d_%s", energy.Data()), "3D lines", 900, 700);
        gStyle->SetOptStat(0);

        frame->SetStats(0);
        frame->Draw(); // só os eixos

        // Agora desenha as linhas em cima
        int colorIndex = 2;

        for (auto &hname : histNames) {
            TH1D *h = (TH1D*) f->Get(hname);
            if (!h) continue;

            // Extrai ângulo
            Ssiz_t posA = hname.Index("_A");
            Ssiz_t posDeg = hname.Index("deg");
            TString angleStr = hname(posA+2, posDeg - (posA+2));
            double angle = angleStr.Atof();

            int nBins = h->GetNbinsX();
            std::vector<double> x(nBins), y(nBins), z(nBins);

            for (int ix = 1; ix <= nBins; ix++) {
                x[ix-1] = h->GetXaxis()->GetBinCenter(ix);
                y[ix-1] = angle;
                z[ix-1] = h->GetBinContent(ix);
            }

            TPolyLine3D *pline = new TPolyLine3D(nBins, x.data(), y.data(), z.data());
            pline->SetLineColor(colorIndex);
            pline->SetLineWidth(2);
            pline->Draw("SAME");

            colorIndex++;
            if (colorIndex > 9) colorIndex = 2;
        }

        c->Update();
    }
}
