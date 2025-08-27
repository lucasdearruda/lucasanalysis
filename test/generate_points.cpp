#include <iostream>
#include <fstream>
#include <TRandom3.h>
#include <TF1.h>

void generate_points() {
    // Open the output file
    std::ofstream out("points.txt");

    // Create the Landau function with parameters
    TF1* fb2 = new TF1("fb2", "TMath::Landau(x,[0],[1],0)", -5, 10);
    fb2->SetParameters(0.2, 1.3);  // Set the parameters for the Landau distribution

    // Random number generator for jitter
    TRandom3 r;

    // Generate 100 points with jitter
    for (int i = 0; i < 100; ++i) {
        // Sample an x-value randomly in the range [-5, 10]
        double x = fb2->GetRandom();

        // Evaluate the Landau function at x
        double y = fb2->Eval(x);

        // Add random Gaussian jitter to y
        y += r.Gaus(0, 0.01);  // Mean 0, stddev 0.1 for jitter

        // Write the point to the file
        out << x << " " << y << "\n";
    }

    // Close the file
    out.close();
    
    std::cout << "Points generated and saved to points.txt" << std::endl;
}

