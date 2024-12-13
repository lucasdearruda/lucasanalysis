#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TColor.h"

void readRGBFileAndSetPalette(const std::string& filename) {
    std::vector<Int_t> colors;
    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        float r, g, b;
        if (iss >> r >> g >> b) {
            // Convert normalized RGB values (0 to 1) to integer values (0 to 255)
            Int_t r_int = static_cast<Int_t>(r * 255);
            Int_t g_int = static_cast<Int_t>(g * 255);
            Int_t b_int = static_cast<Int_t>(b * 255);
            
            Int_t color = TColor::GetColor(r_int, g_int, b_int);
            colors.push_back(color);
            
            // Uncomment this for debug output
            // std::cout << "Read RGB (normalized): (" << r << ", " << g << ", " << b << ") -> "
            //          << "RGB (integer): (" << r_int << ", " << g_int << ", " << b_int << ") -> "
            //          << "Color Index: " << color << std::endl;
        } else {
            std::cerr << "Error reading RGB values from line: " << line << std::endl;
        }
    }

    // Set the palette in gStyle directly
    int nColors = colors.size();
    gStyle->SetPalette(nColors, &colors[0]);
}


void example() {

    readRGBFileAndSetPalette("colors.txt");

    TCanvas* c1 = new TCanvas("c1", "Palette Example", 800, 200);
    // Create a 2D histogram
    TH2F *h2 = new TH2F("h2", "Example 2D Histogram", 100, -4, 4, 100, -4, 4);

    // Fill 2D histogram with random data
    for (int i = 0; i < 10000; ++i) {
        h2->Fill(gRandom->Gaus(0, 1), gRandom->Gaus(0, 1));
    }

    // Create a canvas and draw the 2D histogram
    TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
    h2->Draw("COLZ");  // Use COLZ option to show the color palette
}
