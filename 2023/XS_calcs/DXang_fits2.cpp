/////////////////////////////////
//Code pour faire le fit des angles DX à partir du fichier Fe_pXS.root
///////

#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TLatex.h"
#include "TChain.h"
#include <string>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TKey.h"
#include "TGraph.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TMath.h>
//#include "/mnt/medley/LucasAnalysis/useful.h" //version6.10.2024.0002

#include "TStopwatch.h"
#include <time.h>
TGraph* fitToGraph(TF1* f, int npoints=1000, double xmin=-1, double xmax=1) {
    TGraph* gr = new TGraph(npoints);
    double step = (xmax - xmin) / (npoints-1);
    for (int i = 0; i < npoints; ++i) {
        double x = xmin + i*step;
        double y = f->Eval(x);
        gr->SetPoint(i, x, y);
    }
    return gr;
}

// Define uma função com combinação de P0, P2, P4, P6
Double_t legendreFitFunc(Double_t *x, Double_t *par) {
    Double_t c = x[0];
    return par[0]*1.0 +                           // P0
           par[1]*c +                             // P1
           par[2]*1.0/2*(3*c*c - 1) +               // P2
           par[3]*1.0/2*(5*c*c*c - 3*c) +           // P3
           par[4]*1.0/8*(35*c*c*c*c - 30*c*c + 3) + // P4
           par[5]*1.0/8*(63*c*c*c*c*c - 70*c*c*c + 15*c);// + // P5
          // par[6]*0.0625*(231*c*c*c*c*c*c - 315*c*c*c*c + 105*c*c - 5); // P6
}

void DXang_fits2(){


    TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_pXS.root", "READ");


    // Arquivo de saída para salvar os parâmetros do fit
    std::ofstream out("fit_results_Legendre.txt");
    if (!out.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    out << "Graph_Name\tParameter_Index\tValue\tError\n";



    // Vector para guardar os TF1 após o fit
    std::vector<TF1*> fitsVec;
    std::vector<TGraph*> graphsVec; // também guarda os graphs (opcional)


    TIter next(ff->GetListOfKeys());
    TKey *key;

    TCanvas *cc = new TCanvas("cc", "Fits de DXang", 800, 600);
    while ((key = (TKey *)next())) {
        std::string name = key->GetName();

        // Filtra apenas os TGraph que terminam com "_cos"
        if (name.find("_cos") == std::string::npos) continue;

        TObject *obj = key->ReadObj();
        TGraph *graph = dynamic_cast<TGraph *>(obj);
        if (!graph) continue;

        //TF1 *fitFunc = new TF1(("fit_" + name).c_str(), "pol4", -1.0, 1.0);
        TF1 *fitFunc = new TF1(("fit_" + name).c_str(), legendreFitFunc, -1.0, 1.0,6);
        fitFunc->SetParNames("a0", "a1", "a2", "a3", "a4", "a5");


        cc->Clear();
        graph->SetMarkerStyle(20);      // Estilo 20 = círculo cheio
        graph->SetMarkerSize(1.2);      // Tamanho do marcador (ajuste conforme necessário)
        graph->SetLineWidth(2);         // Espessura da linha do marcador (não afeta linha entre pontos se não tiver)
        graph->Draw("AP");   // "AP" = Axis and Points (eixos e pontos)


        graph->Fit(fitFunc, "R");  // R: range fit, Q: quiet
        fitsVec.push_back(fitFunc);
        graphsVec.push_back(graph);  // opcional, para desenhar os pontos depois



        cc->SaveAs(("fits/"+name + "_fit.png").c_str()); // Salva o gráfico com o fit
        // Salva os parâmetros
        //for (int i = 0; i <= 4; ++i) {
        out << name << endl;
        for (int i = 0; i <= 5; ++i) {
            double value = fitFunc->GetParameter(i);
            double error = fitFunc->GetParError(i);
            //out << name << "\t" << i << "\t" << value << "\t" << error << "\n";
            out << "a" << i << "\t" << value << "\t" << error << "\n";
        }

        //delete fitFunc;
    }

    // Cria o canvas e o multigraph
    TCanvas *c0 = new TCanvas("c1", "multigraph fits pol4", 1200, 700);
    TMultiGraph *mg = new TMultiGraph();
    // Dentro do seu loop for:
    for (size_t i = 0; i < fitsVec.size(); ++i) {
        fitsVec[i]->SetLineWidth(2);
        fitsVec[i]->SetLineColor(kBlue + (i % 7));
        
        TGraph* fitGraph = fitToGraph(fitsVec[i]);  // usa o TF1 armazenado no vetor
        
        //fitGraph->SetLineColor(kBlue + (i % 7));
        fitGraph->SetLineWidth(2);
        fitGraph->SetLineColor(kRed);
        fitGraph->SetNameTitle(Form("%zu",i+1),Form("%zu",i+1));
        mg->Add(fitGraph, "L");

    }



    mg->SetTitle("Fits Pol4; cos(#theta); Cross section (mb/sr)");
    mg->Draw("a fb l3d");  // Desenha axes + as linhas

    c0->Update();

//Extra section 
    
    TGraph2D* g2d_points = new TGraph2D();
    TGraph2D* g2d_fits   = new TGraph2D();

    int pointIndex = 0;
    for (size_t i = 0; i < graphsVec.size(); ++i) {
        TGraph* g = graphsVec[i];
        TF1* f = fitsVec[i];

        int npoints = g->GetN();
        for (int j = 0; j < npoints; ++j) {
            double x, y;
            g->GetPoint(j, x, y);
            g2d_points->SetPoint(pointIndex, x, i, y);  // z = valor
            ++pointIndex;
        }

        // Adiciona curva suavizada do fit (100 pontos)
        const int nfit = 100;
        for (int j = 0; j < nfit; ++j) {
            double x = -1.0 + j * (2.0 / (nfit - 1));
            double y = f->Eval(x);
            g2d_fits->SetPoint(i * nfit + j, x, i, y);  // mesmo índice i
        }
    }

    TCanvas *c3d = new TCanvas("c3d", "3D Fit + Data", 1000, 700);
    g2d_points->SetMarkerStyle(20);
    g2d_points->SetMarkerSize(0.7);
    g2d_points->SetMarkerColor(kBlack);
    g2d_points->GetXaxis()->SetTitle("cos(#theta)");
    g2d_points->GetYaxis()->SetTitle("bin number");
    g2d_points->GetZaxis()->SetTitle("Cross section (mb/sr)");
    g2d_points->Draw("p0");  // pontos no espaço 3D

    std::vector<TPolyLine3D*> curves3D;

    for (size_t i = 0; i < fitsVec.size(); ++i) {
        TF1* f = fitsVec[i];

        const int nfit = 100;
        TPolyLine3D* line = new TPolyLine3D(nfit);
        
        for (int j = 0; j < nfit; ++j) {
            double x = -1.0 + j * (2.0 / (nfit - 1));  // cos(theta)
            double y = f->Eval(x);                     // sigma(theta)
            double z = i;                              // índice (ou energia)
            line->SetPoint(j, x, z, y);                // ordem: x, y, z
        }

        line->SetLineColor(kRed);
        line->SetLineWidth(2);
        curves3D.push_back(line);
    }
    for (auto line : curves3D) {
        line->Draw();  // Desenha curva ajustada como linha 3D
    }

    // Limpar memoria depois se quiser
    for (auto f : fitsVec) delete f;
    for (auto g : graphsVec) delete g;

    ff->Close();
    out.close();
    //ff->Close();

    std::cout << "Fits concluídos e resultados salvos em 'fit_results_pol6.txt'.\n";
}