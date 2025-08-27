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
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/include/bin_width.hh" 
//
//
//  Script adapted to deal with TGraphErrors' results, procuded by 'prodduce_pXSv2.cpp'.
//
//
#include <string>
#include <utility> // for std::pair
#include <regex>


std::pair<double, double> ExtractEnergyRange(const std::string& name) {
    std::regex re("ENN_([0-9]+)p([0-9]+)_([0-9]+)p([0-9]+)");
    std::smatch match;

    if (std::regex_search(name, match, re) && match.size() == 5) {
        std::string low_str = match[1].str() + "." + match[2].str();
        std::string high_str = match[3].str() + "." + match[4].str();
        double low = std::stod(low_str);
        double high = std::stod(high_str);
        return {low, high};
    }

    // Return default (or throw exception if needed)
    return {0.0, 0.0};
}

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
    return par[0]*(1.0 +
           par[1]*c +
           par[2]*0.5*(3*c*c - 1) +
           par[3]*0.5*(5*c*c*c - 3*c)+
           par[4]*0.125*(35*c*c*c*c - 30*c*c + 3));//+
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
void GetXSfromDDX_night(string filename = "Fe_a_DDX.root", bool step_by_step = false) {
    // Open input ROOT file
    
    TFile* fout = new TFile(Form("XSfrom%s_results2.root",filename.substr(0, filename.rfind(".root")).c_str()), "RECREATE"); // Output ROOT file for graphs and fits
    cout<<"Writting results to: "<<fout->GetName()<<endl;
    TFile *ff = new TFile(filename.c_str(), "READ");
    cout<< "Opening file: " << filename << endl;

    // Output text file for fit parameters
    std::ofstream out("fit_results_Legendre.txt", std::ios::app);
    if (!out.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    //
    //Temporary output 
    //
    std::ofstream out_integral("integral_results_Legendre2.txt", std::ios::app);
    if (!out_integral.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }


    //
    //Temporary output -- numerical
    //
    std::ofstream out_integral_num("integral_results_TRAPEZ2.txt", std::ios::app);
    if (!out_integral.is_open()) {
        std::cerr << "Erro ao criar o arquivo de saída." << std::endl;
        return;
    }

    out << "#"<<N_param<<" parameters:"<<endl;
    out_integral << "#"<<N_param<<" parameters: [En Int(mb)]"<<endl;
    out_integral_num << "E, Int, errInt"<<endl;
    
    // Iterate over keys in the file
    TIter next(ff->GetListOfKeys());
    
    TKey *key;
    //TCanvas *cc = new TCanvas("cc", "Fits de DXang", 800, 600);

    TCanvas *cc = new TCanvas("cc", "Angles", 1700, 600);
    cc->Divide(3,1); // Divide canvas into 3 pads for better visualization
    TGraphErrors *totXS = new TGraphErrors();
    Float_t integral_function = 0;
    Float_t parameters[] = {5.20957, 3.59675, 1.81913, -0.0940362, 0.836435}; // Initialize parameters for the fit

    // Loop through all keys in the file
    // This will read all TGraphErrors objects in the file
    // and perform the fitting and graphing operations.
        while ((key = (TKey*)next())) {
            if (TString(key->GetClassName()) == "TGraphErrors") {
                TGraphErrors* g1 = (TGraphErrors*)key->ReadObj();
                TString name1 = g1->GetName();
                // You can also print its name:
                std::cout << "Loaded: " << name1 << std::endl;
                

                // Only process gAngle graphs (skip gCosAngle for now)
                if (!name1.BeginsWith("gAngle")) continue;

                std::cout << "Processing pair: " << name1 << std::endl;

                // Get the corresponding gCosAngle graph (next key)
                TKey* nextKey = (TKey*)next();
                if (!nextKey || TString(nextKey->GetClassName()) != "TGraphErrors") continue;

                TGraphErrors* g2 = (TGraphErrors*)nextKey->ReadObj(); // gCosAngle
                string name2 = g2->GetName();
                std::cout << "               : " << name2 << std::endl;
                ///////// Pad 1 - plot gAngle and integrate //////////
                        cc->cd(1);
                        g1->SetTitle(Form("gAngle-%s",name2.c_str()));
                        g1->SetMarkerStyle(20);
                        g1->Draw("APL");
                        gPad->SetGrid();


                ///////// Pad 2 - plot gCosAngle //////////////////////
                cc->cd(2);
                g2->SetTitle("gCosAngle - for Fitting");
                g2->SetMarkerStyle(21);
                g2->SetMarkerColor(kBlue);
                g2->Draw("AP");
                gPad->SetGrid();

                TF1 *fitFunc = new TF1(("fit_" + name2).c_str(), legendreFitFunc, -1.0, 1.0, N_param);
                fitFunc->SetParNames("a0", "a1", "a2", "a3", "a4", "a5");
                fitFunc->SetParameters(parameters[0], parameters[1], parameters[2], parameters[3], parameters[4]);
                //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Fumilli");
                g2->Fit(fitFunc, "RFU");
                parameters[0] = fitFunc->GetParameter(0);
                parameters[1] = fitFunc->GetParameter(1);
                parameters[2] = fitFunc->GetParameter(2);
                parameters[3] = fitFunc->GetParameter(3);
                parameters[4] = fitFunc->GetParameter(4);
                ///////// Pad 3 - major/minor variation //////////////
                cc->cd(3);
                TGraph* major = new TGraph();
                TGraph* minor = new TGraph();

                for (int i = 0; i < g2->GetN(); i++) {
                    double x, y;
                    g2->GetPoint(i, x, y);
                    double ey = g2->GetErrorY(i);
                    major->SetPoint(i, x, y + ey);
                    minor->SetPoint(i, x, y - ey);
                }
                major->SetLineColor(kBlue);
                major->SetLineWidth(2);
                major->SetNameTitle("Major Error", "Major Error");
                major->Draw("AP");
                major->SetMarkerStyle(20);
                major->SetMarkerColor(kBlue);
                major->SetMarkerSize(1.2);

                minor->SetLineColor(kBlack);
                minor->SetLineWidth(2);
                minor->SetNameTitle("Minor Error", "Minor Error");
                minor->Draw("P SAME");
                minor->SetMarkerStyle(20);
                minor->SetMarkerColor(kBlack);
                minor->SetMarkerSize(1.2);

                major->Draw("AP");
                minor->Draw("P SAME");

                TF1 *major_fitFunc = new TF1(("major_fit_" + name2).c_str(), legendreFitFunc, -1.0, 1.0, N_param);
                major_fitFunc->SetParameters(fitFunc->GetParameters());
                major->Fit(major_fitFunc, "RF");
                TF1 *minor_fitFunc = new TF1(("minor_fit_" + name2).c_str(), legendreFitFunc, -1.0, 1.0, N_param);
                minor_fitFunc->SetParameters(fitFunc->GetParameters());
                minor->Fit(minor_fitFunc, "RF");


                ///////////////////////////
                       //adding to the output
                fout->cd();
                cc->Write((name2 + "_canvas").c_str());
                cc->SaveAs(("fits/XSfromDDXs/"+name2 + "_canvas.png").c_str());
                g2->Write((name2 + "_data").c_str());
                fitFunc->Write((name2 + "_fit").c_str());
                major->Write((name2 + "_Major").c_str());
                major_fitFunc->Write((name2 + "_major_fit").c_str());
                minor->Write((name2 + "_Minor").c_str());
                minor_fitFunc->Write((name2 + "_minor_fit").c_str());
        
                //////////////////////////////////////////////////////////
                // Calculate the integral of the major and minor fits

                auto range = ExtractEnergyRange(name2);
                Float_t neutron_En = (range.first + range.second)*0.5;

                cout<<"Energies ::  "<<range.first<<" MeV - "<<range.second<<" MeV"<<endl;
                cout<<"Neutron Energy :: "<<neutron_En<<" MeV"<<endl;
                cout<<"/ :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: / :: /" << endl;
                cout << "§§\nCalculating integrals. . .  " << name2 << ":" << endl;
                cout<<"_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ "<<endl;
                Float_t major_integral =  2 * TMath::Pi() * major_fitFunc->Integral(-1, 1);
                Float_t minor_integral =  2 * TMath::Pi() * minor_fitFunc->Integral(-1, 1);
                // Calculate the error in the integral
                Float_t integral_error =   TMath::Abs(major_integral - minor_integral) / 2.0; // Average of the two errors
                // Print the results
                cout << "Major Integral: " << major_integral << ", Minor Integral: " << minor_integral 
                    << ", Integral Error: " << integral_error << endl;

                //cout<< endl <<endl<<"Summarizing integral results for " << name << ":" << endl;
                Float_t integral_function =  2 * TMath::Pi() * fitFunc->Integral(-1, 1);
                cout << "Integral of the fit: " << integral_function << " \u00B1 " << integral_error << endl;
                // Now, we can save the integral results to the output file
                out_integral << neutron_En << "\t" << integral_function << "\t" << integral_error << std::endl;
                //and add the point to the total XS graph

                //numerical trapezoidal integration
                Float_t major_integralT =  2 * TMath::Pi() * TrapezoidalIntegration(major,-1,1);
                Float_t minor_integralT =  2 * TMath::Pi() * TrapezoidalIntegration(minor,-1,1);
                // Calculate the error in the integral
                Float_t integral_errorT =   TMath::Abs(major_integralT - minor_integralT) / 2.;
                Float_t integralT =   2 * TMath::Pi() * TrapezoidalIntegration(g2,-1,1);
                cout<< "Numerical Integral of the fit: " << integralT << " \u00B1 " << integral_errorT << endl;
                out_integral_num << neutron_En << "\t" << integralT << "\t" << integral_errorT << std::endl;

                totXS->AddPoint(neutron_En, integral_function);
                totXS->SetPointError(totXS->GetN()-1, neutron_En, integral_error); // Set error for the point

                        
                gPad->SetGrid();
                if(step_by_step){
                    cout<<"press some key to continue..."<<endl;
                    cin.get(); // Wait for user input to continue
                }
                        
            
            
            }
        }
    fout->Close();
    delete fout;
}
