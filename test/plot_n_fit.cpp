#include <iostream>
#include <fstream>
#include <vector>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

void plot_n_fit() {
    // Open the file
    std::ifstream infile("points.txt");
    if (!infile.is_open()) {
        std::cerr << "Error: Cannot open points.txt!" << std::endl;
        return;
    }

    // Read the data into vectors
    std::vector<double> x_vals, y_vals;
    double x, y;
    while (infile >> x >> y) {
        x_vals.push_back(x);
        y_vals.push_back(y);
    }
    infile.close();

    // Create a TGraph from the data
    int n_points = x_vals.size();
    TGraph *graph = new TGraph(n_points, &x_vals[0], &y_vals[0]);
    graph->SetTitle("Fitted Data;X;Y");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.6);

    // Fit a 4th-degree polynomial in the range [-3, 5]
    TF1 *fitFunc = new TF1("fitFunc", "pol4", -3, 5);
    graph->Fit(fitFunc, "R");

    // Find the maximum value of the fit function
    double max_x = fitFunc->GetMaximumX(-3, 5);
    double max_y = fitFunc->Eval(max_x);

    // Print the result
    std::cout << "Maximum value of the fit: Y = " << max_y << " at X = " << max_x << std::endl;

    // Draw the graph and the fit
    TCanvas *c1 = new TCanvas("c1", "Fit Plot", 800, 600);
    graph->Draw("AP");
    fitFunc->Draw("same");

    // Save the plot as an image
    c1->SaveAs("fit_plot.png");
}

