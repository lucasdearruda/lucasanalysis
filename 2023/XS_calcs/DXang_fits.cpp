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

// Define uma função com combinação de P0, P2, P4, P6
Double_t legendreFitFunc(Double_t *x, Double_t *par) {
    Double_t c = x[0];
    return par[0]*1.0 +                           // P0
           par[1]*c +                             // P1
           par[2]*0.5*(3*c*c - 1) +               // P2
           par[3]*0.5*(5*c*c*c - 3*c) +           // P3
           par[4]*0.125*(35*c*c*c*c - 30*c*c + 3); //+  // P4
          // par[5]*0.125*(63*c*c*c*c*c - 70*c*c*c + 15*c) + // P5
          // par[6]*0.0625*(231*c*c*c*c*c*c - 315*c*c*c*c + 105*c*c - 5); // P6
}

void DXang_fits(){


    TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_pXS.root", "READ");


    // Arquivo de saída para salvar os parâmetros do fit
    std::ofstream out("fit_results_pol4.txt");
    if (!out.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    out << "Graph_Name\tParameter_Index\tValue\tError\n";

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


        // Cria e executa o fit pol6
        // TF1 *fitFunc = new TF1("fitFunc", "pol6", -1.0, 1.0);  // Range genérico em cos(theta)

        //Legendre
        TF1 *fitFunc = new TF1("f_leg", legendreFitFunc, -1, 1, 4); // 4 parâmetros: P0, P2, P4, P6


        cc->Clear();
        graph->SetMarkerStyle(20);      // Estilo 20 = círculo cheio
        graph->SetMarkerSize(1.2);      // Tamanho do marcador (ajuste conforme necessário)
        graph->SetLineWidth(2);         // Espessura da linha do marcador (não afeta linha entre pontos se não tiver)
        graph->Draw("AP");   // "AP" = Axis and Points (eixos e pontos)

        graph->Draw();
        //graph->Fit(fitFunc, "Q"); // "Q" = Quiet mode (sem prints)
        //graph->Fit(fitFunc, "R"); // 
        graph->Fit("pol4", "R"); // 
        cc->SaveAs(("fits/"+name + "_fit.png").c_str()); // Salva o gráfico com o fit
        // Salva os parâmetros
        for (int i = 0; i <= 4; ++i) {
            double value = fitFunc->GetParameter(i);
            double error = fitFunc->GetParError(i);
            out << name << "\t" << i << "\t" << value << "\t" << error << "\n";
        }

        delete fitFunc;
    }

    out.close();
    ff->Close();

    std::cout << "Fits concluídos e resultados salvos em 'fit_results_pol6.txt'.\n";
}