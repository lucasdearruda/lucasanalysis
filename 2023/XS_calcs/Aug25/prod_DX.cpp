//prod_DX.cpp, version 2025-08-25.0
// I added some TNameds to the output file, with the runs used
// and some TNameds with the input parameters, and target info

#include <iostream>     // for std::cout, std::cerr, std::endl
#include <fstream>      // for std::ifstream
#include <string>       // for std::string

//for the time:
#include <time.h>
//for the vectors and pairs:
#include <vector>
#include <utility> 
#include "TVector3.h"

using namespace std;

//////////
//The 'global' variables (like binWidth, mc, M, ...) are defined in functions.cxx for simplicity
#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Aug25/src/functions.cxx" // for my functions 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here we will procude DX for Fe in 2023 runs, from DDX file processed by prod_DDX.cpp::

void prod_DX(
    string inputFileName = "prod_DDX_thick_MC.root",
    bool plotIt = true, 
    bool saveIt = true, 
    bool pause_each = false,
    bool verbose = true
    ) {
    
    string cur_time = getCurrentTime();
    clock_t tStart = clock();
    string outputFileName = "myDX_thick.root"  ; 


    TFile *fIn = TFile::Open(inputFileName.c_str());
    if(!fIn || fIn->IsZombie()) {
        cerr << "Error opening input file: " << inputFileName << endl;
        return;
    }

   // Map: "particle_energy" -> TGraph
    // map<string, TGraph*> graphs;

    // TIter next(fIn->GetListOfKeys()); 
    // TKey *key; 
    // while ((key = (TKey*)next())){
    //     TObject *obj = key->ReadObj(); 
    //     TH1D *h = dynamic_cast<TH1D*>(obj); 
    //     if(h) {
    //         string name = h->GetName(); // ex: "h_Proton_E100_theta30"

            
    //         // separar por "_"
    //         vector<string> tokens;
    //         stringstream ss(name);
    //         string token;
    //         while(getline(ss, token, '_')) tokens.push_back(token);

    //         // Variáveis para armazenar info
    //         string particle;
    //         string e_str, a_str;
    //         bool MC_flag = false;
    //         bool TTC_flag = false;

    //         for (auto &t : tokens) {
    //             if (t[0] == 'P') {
    //                 particle = t.substr(1);
    //             } else if (t[0] == 'E') {
    //                 e_str = t.substr(1); // remove 'E'
    //                 replace(e_str.begin(), e_str.end(), 'p', '.');
    //             } else if (t[0] == 'A') {
    //                 a_str = t.substr(1);
    //                 auto deg_pos = a_str.find("deg");
    //                 if (deg_pos != string::npos) a_str = a_str.substr(0, deg_pos);
    //                 replace(a_str.begin(), a_str.end(), 'p', '.');
    //             } else if (t == "MC" || t == "_MC") {
    //                 MC_flag = true;
    //             } else if (t == "TTC" || t == "_TTC") {
    //                 TTC_flag = true;
    //             }
    //         }

    //         // Converter para números
    //         double energy = stod(e_str);
    //         double theta  = stod(a_str);

    //         cout << h->GetName() 
    //             << ", particle = " << particle << endl
    //             << "  E = " << energy << " MeV"<< endl
    //             << "  theta = " << theta << " deg"<< endl
    //             << "  MC = " << MC_flag<< endl
    //             << "  TTC = " << TTC_flag<< endl
    //             << endl<< endl<< endl;

    //     }
    // }

        
    TIter next(fIn->GetListOfKeys()); 
    TKey *key; 
    map<pair<char,double>, TGraph*> graphs;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj(); 
        TH1D *h = dynamic_cast<TH1D*>(obj); 
        if(!h) continue;

        string name = h->GetName();
        vector<string> tokens;
        stringstream ss(name);
        string token;
        while(getline(ss, token, '_')) tokens.push_back(token);

        char particle = 0;
        string e_str, a_str;
        bool MC_flag = false, TTC_flag = false;

        for(auto &t : tokens) {
            if(t[0] == 'P') particle = t[1];
            else if(t[0] == 'E') { e_str = t.substr(1); replace(e_str.begin(), e_str.end(), 'p', '.'); }
            else if(t[0] == 'A') { 
                a_str = t.substr(1); 
                auto deg_pos = a_str.find("deg"); 
                if(deg_pos != string::npos) a_str = a_str.substr(0, deg_pos);
                replace(a_str.begin(), a_str.end(), 'p', '.');
            }
            else if(t=="MC" || t=="_MC") MC_flag = true;
            else if(t=="TTC" || t=="_TTC") TTC_flag = true;
        }

        double energy = stod(e_str);
        double theta  = stod(a_str);

        cout    << ", particle = " << particle << endl
                << "  E = " << energy << " MeV"<< endl
                << "  theta = " << theta << " deg"<< endl
                << "  MC = " << MC_flag<< endl
                << "  TTC = " << TTC_flag<< endl
                << endl<< endl<< endl;

        // criar/recuperar TGraph
        auto key_graph = make_pair(particle, energy);
        if(graphs.find(key_graph) == graphs.end()) {
            graphs[key_graph] = new TGraph();
            graphs[key_graph]->SetName(Form("g_%c_E%.1f", particle, energy));
            graphs[key_graph]->SetTitle(Form("Particle %c, E=%.1f MeV", particle, energy));
        }

        double integral = h->Integral("width");
        double cos_theta = cos(theta * M_PI / 180.0);
        TGraph *g = graphs[key_graph];
        int n = g->GetN();
        g->SetPoint(n, cos_theta, integral);
    }

    TFile *fOut = new TFile(outputFileName.c_str(), "RECREATE");
    if(!fOut || fOut->IsZombie()) {
        cerr << "Error creating output file: " << outputFileName << endl;
        return;
    }
    // salvar todos os TGraph do map
    for (auto &pair : graphs) {
        TGraph *g = pair.second;
        g->Write(); // salva no ROOT file
    }

    fOut->Close();
    cout << "All graphs saved in " << outputFileName << endl;
    //cout << "Output file " << outputFileName << " created." << endl;

    cout <<"\nTotal execution time: "<< double(clock() - tStart) / (double)CLOCKS_PER_SEC<<" s."<<endl;
    cout<<"prod_DX.cpp, version 1.2025-08-25.0"<<endl;

    


}