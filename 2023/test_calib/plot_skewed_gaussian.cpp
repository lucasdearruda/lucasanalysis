#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "Math/DistFunc.h" // for normal_cdf

// Skewed Gaussian function definition
double skewed_gaussian(double* x, double* par) {
    double A  = par[0]; // amplitude
    double mu = par[1]; // mean
    double sigma = par[2]; // standard deviation
    double alpha = par[3]; // skewness

    // Normalization factor
    double norm = A * 2.0 / (sigma * sqrt(2.0 * TMath::Pi()));
    double arg1 = (x[0] - mu) / sigma;
    double gauss = exp(-0.5 * arg1 * arg1); // Gaussian part
    double cdf = ROOT::Math::normal_cdf(alpha * arg1, 1.0, 0.0); // CDF part

    // Return the skewed Gaussian
    return norm * gauss * cdf;
}

// Function to generate data, plot histogram, and fit the skewed Gaussian
void gfit(Int_t Npoints = 1e4) {
    // Initialize TRandom2 for random number generation
    TRandom2 rand(0); // No seed, or you can set a specific value

    // Retrieve the skewed Gaussian function
    TF1* fskew = (TF1*)gROOT->GetFunction("fskew");
    if (!fskew) {
        std::cerr << "Function 'fskew' not found. Run plot_skewed_gaussian() first." << std::endl;
        return;
    }

    // Create histogram to store generated data
    TH1D* h = new TH1D("h", "Generated Data from Skewed Gaussian; x; Counts", 200, -10, 10);

    // Generate Npoints according to the skewed Gaussian distribution
    for (Int_t i = 0; i < Npoints; ++i) {
        double val = fskew->GetRandom(); // Generate a random value based on the skewed Gaussian
        h->Fill(val); // Fill histogram with the generated value
    }

    // Create a canvas and draw the histogram
    TCanvas* c2 = new TCanvas("c2", "Generated Data & Fit", 800, 600);
    h->SetLineColor(kBlack);
    h->SetMarkerStyle(20);
    h->Draw();

    // Fit the histogram with the skewed Gaussian function
    h->Fit(fskew, "R"); // "R" = use the function's range for fitting

    // Print the fit results in the terminal
    std::cout << "\n--- Fit Results ---" << std::endl;
    for (int i = 0; i < 4; ++i) {
        std::cout << fskew->GetParName(i) << " = " 
                  << fskew->GetParameter(i) << " ± " 
                  << fskew->GetParError(i) << std::endl;
    }

    // Create a TPaveText to display the parameters on the plot
    TPaveText* pave = new TPaveText(0.15, 0.6, 0.4, 0.80, "NDC");
    pave->SetBorderSize(0);
    pave->SetFillColor(0);
    pave->SetTextAlign(12); // Centered text

    // Add fit parameter text
    pave->AddText(Form("A = %.3f #pm %.3f", fskew->GetParameter(0), fskew->GetParError(0)));
    pave->AddText(Form("mu = %.3f #pm %.3f", fskew->GetParameter(1), fskew->GetParError(1)));
    pave->AddText(Form("sigma = %.3f #pm %.3f", fskew->GetParameter(2), fskew->GetParError(2)));
    pave->AddText(Form("alpha = %.3f #pm %.3f", fskew->GetParameter(3), fskew->GetParError(3)));

    pave->Draw();
}

// Function to plot the skewed Gaussian and call gfit()
void plot_skewed_gaussian() {
    TCanvas* c1 = new TCanvas("c1", "Skewed Gaussian using normal_cdf", 800, 600);

    // Define the skewed Gaussian function
    TF1* fskew = new TF1("fskew", skewed_gaussian, -10, 10, 4);
    fskew->SetParameters(1.0, 2.3, 1.0, 3.0); // A, mu, sigma, alpha
    fskew->SetParNames("Amplitude", "Mean", "Sigma", "Alpha");
    fskew->SetLineColor(kRed+1);

    // Increase the number of points for a smooth curve
    fskew->SetNpx(1000); 

    // Set plot title and axis labels
    fskew->SetTitle("Skewed Gaussian using normal_cdf; x; f(x)");

    // Draw the skewed Gaussian function
    fskew->Draw();

    
    // Generate data, plot the histogram, and fit the skewed Gaussian
    gfit();
}
