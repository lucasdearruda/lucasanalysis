#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>

void PlotXs() {
    // File path — adjust if needed
    std::string filename = "XS_output.csv";

    // Open file
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    // Prepare angles (matching your columns after the first)
    std::vector<double> angles = {20, 40, 60, 80, 100, 120, 140, 160};

    TGraph2D *graph = new TGraph2D();
    TGraph *sGr[angles.size()]; //single Graphs

    //declaring the single graphs
    for(int i=0;i<angles.size();i++){
        sGr[i] = new TGraph();
    }

    std::string line;
    int pointIndex = 0;
    //index for single plots
    std::vector<int> grIndex(angles.size(), 0);  // Um índice para cada sGr[i]
    // it the same as this: std::vector<int> grIndex = {0, 0, 0, 0, 0, 0, 0, 0}; 


    
    while (std::getline(infile, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;

        // Split line by comma
        while (std::getline(ss, cell, ',')) {
            cells.push_back(cell);
        }

        // Check if enough columns
        if (cells.size() < angles.size() + 1) {
            std::cerr << "Warning: line with insufficient columns:\n" << line << std::endl;
            continue;
        }

        // First column is X
        double X = std::stod(cells[0]);

        // Loop over angle columns
        for (size_t i = 0; i < angles.size(); ++i) {
            double val = std::stod(cells[i + 1]);
            if (X != 0) {
                sGr[i]->SetNameTitle(Form("sGr%d", (int)angles[i]), Form("sGr%d", (int)angles[i]));
                sGr[i]->SetPoint(grIndex[i]++, X, val);  // Usa o índice próprio do gráfico
            }
            graph->SetPoint(pointIndex++, X, angles[i], val);
        }
    }

    infile.close();

    // Draw graph
    TCanvas *c = new TCanvas("c", "TGraph2D Example", 800, 600);
    c->SetRightMargin(0.15);
    graph->SetTitle(";E_{NN} (MeV);#theta_{LAB} (deg);d#sigma/d#Omega (mb/sr)");
    
    //graph->Draw("P0 COLZ");
    graph->Draw("surf1Z");
    graph->GetYaxis()->SetTitleOffset(3.62);
    graph->GetXaxis()->SetTitleOffset(3.62);
    gPad->Update(); 
    c->SaveAs("graph2D.png");

    // Save graph to ROOT file
    TFile out("output.root", "RECREATE");
    graph->Write("graph2D");

    for(int i=0;i<angles.size();i++){
        sGr[i]->Write();
    }

    out.Close();
}
