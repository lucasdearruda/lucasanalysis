#include "TF1.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "/mnt/medley/LucasAnalysis/useful.h" //version2025.02.26.002

//
//
//  Script adapted to deal with TGraphErrors' results, procuded by 'prodduce_pXSv2.cpp'.
//
//

//=========================================
// Function to convert TF1 into a TGraph
//=========================================
TGraph* fitToGraph(TF1* f, int npoints = 1000, double xmin = -1.0, double xmax = 1.0) {
    TGraph* gr = new TGraph(npoints);
    double step = (xmax - xmin) / (npoints - 1);
    for (int i = 0; i < npoints; ++i) {
        double x = xmin + i * step;
        gr->SetPoint(i, x, f->Eval(x));
    }
    return gr;
}


//=========================================
// Function to extract neutron energy value
// from a graph name string.
//
// It searches for the energy substring located 
// between the last underscore ('_') before 
// "_MeV_cos" and the "_MeV_cos" itself.
//
// Converts notation like "4p5" to "4.5" and returns 
// the energy as a Float_t.
//
// Returns -1 if the expected pattern is not found.
//=========================================
Float_t EfromName(const std::string& name) {
    size_t posMeV = name.find("_MeV_cos");
    if (posMeV == std::string::npos) return -1; // erro se não achar

    // procura o '_' antes do trecho '_MeV_cos'
    size_t posUnder = name.rfind('_', posMeV - 1);
    if (posUnder == std::string::npos) return -1;

    std::string energy_str = name.substr(posUnder + 1, posMeV - posUnder - 1);

    // troca 'p' por '.'
    for (auto& c : energy_str) {
        if (c == 'p') c = '.';
    }

    return atof(energy_str.c_str());
}



//======================================================
// Legendre fit function with P0 to P5 terms (up to l=5)
//======================================================
Double_t legendreFitFunc(Double_t *x, Double_t *par) {
    double c = x[0];
    return par[0]*1.0 +
           par[1]*c +
           par[2]*0.5*(3*c*c - 1) +
           par[3]*0.5*(5*c*c*c - 3*c)+
           par[4]*0.125*(35*c*c*c*c - 30*c*c + 3);//+
           //par[5]*0.125*(63*c*c*c*c*c - 70*c*c*c + 15*c);
    // return par[0]*(1.0 +
    // par[1]*c +
    // 0*par[2]*0.5*(3*c*c - 1) +
    // 0*par[3]*0.5*(5*c*c*c - 3*c)) ;//+
    //par[4]*0.125*(35*c*c*c*c - 30*c*c + 3);// +
    //par[5]*0.125*(63*c*c*c*c*c - 70*c*c*c + 15*c);
}
Int_t N_param = 5;
//======================================================
// Main function to read graphs, fit with Legendre, and draw
//======================================================
void DXang_fits_wErrors() {
    // Open input ROOT file
    TFile *ff = new TFile("/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_pXS.root", "READ");

    // Output text file for fit parameters
    std::ofstream out("fit_results_Legendre.txt", std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    //
    //Temporary output 
    //
    std::ofstream out_integral("integral_results_Legendre.txt", std::ios::app);
    if (!out_integral.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    out << "#"<<N_param<<" parameters:"<<endl;
    out_integral << "#"<<N_param<<" parameters: [En Int(mb)]"<<endl;
    
    // Store fitted functions and graphs
    std::vector<TF1*> fitsVec;
    std::vector<TGraph*> graphsVec;
    std::vector<Float_t> neutron_En;

    // Iterate over keys in the file
    TIter next(ff->GetListOfKeys());
    
    TKey *key;
    TCanvas *cc = new TCanvas("cc", "Fits de DXang", 800, 600);
    TGraphErrors *totXS = new TGraphErrors();
    
    float parameters[] = {7,1.3,4.5,8.5,0.2}; // Initialize parameters for the fit

    while ((key = (TKey*)next())) {
        std::string name = key->GetName();
        if (name.find("_cos") == std::string::npos) continue; //if didnt find "_cos" in the name, skip it

        TGraphErrors *graph = dynamic_cast<TGraphErrors*>(key->ReadObj());
        if (!graph) continue;
        neutron_En.push_back(EfromName(name)); // Store neutron energy
        cout<< "\n. .\n. .\nProcessing graph: " << name << " with energy: " << neutron_En.back() << " MeV" << endl;
        cout<< ". .\n. ."<<endl;

        TF1 *fitFunc = new TF1(("fit_" + name).c_str(), legendreFitFunc, -1.0, 1.0, N_param);
        fitFunc->SetParNames("a0", "a1", "a2", "a3", "a4", "a5");
        //fitFunc->SetParameters(7,1.3,4.5,8.5,0.2);
        fitFunc->SetParameters(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
        // Style and draw
        cc->Clear();
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetLineWidth(2);
        graph->Draw("AP");
        // Fit and store
        graph->Fit(fitFunc, "RF");
        fitsVec.push_back(fitFunc);
        graphsVec.push_back(graph);
        
        
        cout<<"press some key to continue..."<<endl;
        cin.get(); // Wait for user input to continue
        //Create folder for each parnum
        
        char folderName[100];
        sprintf(folderName, "fits/par%d", N_param);  // Cria o nome da pasta com N_param

        // Cria a pasta se não existir
        gSystem->mkdir(folderName, true);

        // Monta o caminho completo do arquivo
        char FigfileName[200];
        sprintf(FigfileName, "%s/%s_fit.png", folderName, name.c_str());

        // Save canvas
        cc->SaveAs(FigfileName); // Salva o gráfico com o fit

        // Save parameters to text file
        out << name << "\n";
        for (int i = 0; i < N_param; ++i) {
            out << "a" << i << "\t" 
                << fitFunc->GetParameter(i) << "\t" 
                << fitFunc->GetParError(i) << "\n";

            parameters[i] = fitFunc->GetParameter(i); // Update parameters for next fit
        }
        // Float_t Leg_int = fitFunc->Integral(-1, 1);
        // totXS->AddPoint(points++, Leg_int);
         //out << "Integral\t" << Leg_int << "\n\n"; // Save integral of the it
        // out_integral << neutron_En[neutron_En.size()-1] <<"\t"<< Leg_int <<endl; // Save integral of the fit
        
    }

    

    // //========================
    // //Integral of the fits -- total XS
    // //========================
    // TCanvas *cXS = new TCanvas("cXS", "Integral of the fits", 1000, 700);
    // totXS->SetName("totalXSgraph");
    // totXS->SetTitle("Total XS from Legendre fits;cos(#theta);Cross section (mb/)");
    // totXS->Draw();
    

    // //========================
    // // Multigraph of all fits
    // //========================
    // TCanvas *c0 = new TCanvas("c1", "multigraph fits", 1200, 700);
    // TMultiGraph *mg = new TMultiGraph();
    
    // Int_t points =0;
    
    // for (size_t i = 0; i < fitsVec.size(); ++i) {
    //     TGraph* fitGraph = fitToGraph(fitsVec[i]);
    //     Float_t Leg_int = (Float_t)TrapezoidalIntegration(fitGraph, -1, 1); // Integral using trapezoidal rule
    //     //Float_t Leg_int = (Float_t)TrapezoidalIntegration(fitGraph, -1, 1, 1000); // Integral using trapezoidal rule
    //     totXS->AddPoint(points++, Leg_int);
        
        
    //     //out << "Integral\t" << Leg_int << "\n\n"; // Save integral of the fit
    //     out_integral << neutron_En[i] <<"\t"<< Leg_int <<endl; // Save integral of the fit

    //     fitGraph->SetLineWidth(2);
    //     fitGraph->SetLineColor(kRed);
    //     fitGraph->SetNameTitle(Form("%.1f", neutron_En[i]), Form("%.1f", neutron_En[i]));
    //     mg->Add(fitGraph, "L");
    //     //out_integral << neutron_En[i] << "\t" << Leg_int << endl; // Save integral of the fit
    //     cout << "Integral for " << neutron_En[i] << " MeV: " << 2*TMath::Pi()*Leg_int << endl;
    // }

    // mg->SetTitle("Fits with Legendre Polynomials;cos(#theta);Cross section (mb/sr)");
    // mg->Draw("a fb l3d");
    // c0->Update();

    // //==============================
    // // 3D Plot: Raw Data + Fit Curves
    // //==============================
    // TGraph2D* g2d_points = new TGraph2D();

    // // Fill data points into 3D graph
    // int pointIndex = 0;
    // for (size_t i = 0; i < graphsVec.size(); ++i) {
    //     TGraph* g = graphsVec[i];
    //     for (int j = 0; j < g->GetN(); ++j) {
    //         double x, y;
    //         g->GetPoint(j, x, y);
    //         g2d_points->SetPoint(pointIndex++, x, neutron_En[i], y);  // x = cos(theta), y = binIndex, z = cross section
    //     }
    // }

    // // Fit curves in 3D using PolyLine3D
    // std::vector<TPolyLine3D*> curves3D;
    // for (size_t i = 0; i < fitsVec.size(); ++i) {
    //     TF1* f = fitsVec[i];
    //     const int nfit = 100;
    //     TPolyLine3D* line = new TPolyLine3D(nfit);

    //     for (int j = 0; j < nfit; ++j) {
    //         double x = -1.0 + j * (2.0 / (nfit - 1));
    //         double y = f->Eval(x);
    //         double z = i;  // bin index
    //         line->SetPoint(j, x, neutron_En[i], y);  // order: x, y, z
    //     }

    //     line->SetLineColor(kRed);
    //     line->SetLineWidth(2);
    //     curves3D.push_back(line);
    // }

    // // Draw in 3D
    // TCanvas *c3d = new TCanvas("c3d", "3D Fit + Data", 1000, 700);
    // g2d_points->SetMarkerStyle(20);
    // g2d_points->SetMarkerSize(0.7);
    // g2d_points->SetMarkerColor(kBlack);
    // g2d_points->GetXaxis()->SetTitle("cos(#theta)");
    // g2d_points->GetYaxis()->SetTitle("E_{NN} (MeV)");
    // g2d_points->GetZaxis()->SetTitle("Cross section (mb/sr)");
    // g2d_points->Draw("p0");

    // for (auto line : curves3D)
    //     line->Draw();

    // //========================
    // // Clean up
    // //========================
    // for (auto f : fitsVec) delete f;
    // for (auto g : graphsVec) delete g;
    // ff->Close();
    // out.close();
    // out_integral.close();
    // cc->Close();

    // std::cout << "Fits concluídos e resultados salvos em 'fit_results_Legendre.txt'.\n";
}
