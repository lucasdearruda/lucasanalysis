
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <TGraph.h>
#include <TTree.h>

using namespace std;

void plotMeAxs_gpt(string namefile = "", string path = "/mnt/medley/LucasAnalysis/talys_calcs", string subpath = "Fe0/best_31.5", bool same = true, string outfile = "output.out") {

    // Check if the namefile is provided
    if(namefile == "") {
        cerr << "Error: No file name provided." << endl;
        return;
    }

    // Construct the full file path using path, subpath, and namefile
    string fullPath = path + "/" + subpath + "/" + namefile;
    string fullPath_out = path + "/" + subpath + "/" + outfile;

    // Open the output file for reading isotope abundances
    cout<< "Opening file: " << fullPath_out << endl;
    ifstream infile(fullPath_out);
    if (!infile.is_open()) {
        cerr << "Error: Could not open output file " << outfile << endl;
        return;
    }
 // Read the isotope abundances section
string line;
bool readingAbundances = false;
map<string, float> isotopeAbundances;

// Skipping to the "Isotope Abundance" section in the output file
while (getline(infile, line)) {
    //cout << "Line: " << line << endl;  // Debugging output
    
    if (line.find("Isotope Abundance") != string::npos) {
        readingAbundances = true;  // Start reading isotope abundances
        continue;  // Skip the "Isotope Abundance" header line
    }

    if (readingAbundances) {
        // Skip any empty lines after the header
        if (line.empty()) {
            continue;  // Skip the empty line
        }

        // Read isotope and abundance from the file
        stringstream ss(line);
        string isotope;
        float abundance;

        ss >> isotope >> abundance;
        cout << "Isotope: " << isotope << ", Abundance: " << abundance << endl;  // Debugging output
        if (!isotope.empty()) {
            isotopeAbundances[isotope] = abundance;
        }

        // Stop reading abundances when reaching the next section (after isotope abundances)
        // For example, when we encounter "TALYS-2.1" in the next line or any other header
        if (line.find("TALYS-2.1") != string::npos) {
            break;  // Stop reading after this point
        }
    }
}


    // Create a graph for the total weighted sum of all isotopes
    TGraph *grTotal = new TGraph();
    bool firstGraph = true;

    // Loop through the isotopes in the map and process the corresponding files
    for (const auto& entry : isotopeAbundances) {
        string isotope = entry.first;
        float abundance = entry.second;

        // Construct the isotope-specific file name
        string isotopeFile = fullPath + "." + isotope.substr(0, 3);  // e.g., "dddxE0031.500A020.0.deg.054"
        
        // Open the isotope file
        TTree *tx = new TTree("xs", "xs");
        tx->ReadFile(isotopeFile.c_str(), "Eprod/F:Total/F");

        Float_t Eprod, total;
        Long64_t nEntries = tx->GetEntries();
        tx->SetBranchAddress("Eprod", &Eprod);
        tx->SetBranchAddress("Total", &total);

        // Create a graph for the individual isotope plot
        TGraph *grIndividual = new TGraph();

        // Loop through the entries in the tree and accumulate the weighted sum
        for (Long64_t i = 0; i < nEntries; ++i) {
            tx->GetEntry(i);

            // Add the weighted data to the individual graph
            grIndividual->AddPoint(Eprod, total);

            // Add the weighted data to the total weighted sum
            grTotal->AddPoint(Eprod, total * abundance);
        }

        // Plot the individual isotope data
        if (same && !firstGraph) {
            grIndividual->Draw("same");
        } else {
            grIndividual->Draw();
            firstGraph = false;
        }
    }

    // Plot the total weighted sum
    if (same) {
        grTotal->Draw("same");
    } else {
        grTotal->Draw();
    }
}
