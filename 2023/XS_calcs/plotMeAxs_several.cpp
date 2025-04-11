#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <regex>  // Para usar expressões regulares
using namespace std;


void parseIso(map<string, float>& isotopeAbundances, string outfile = "path/to/out") {
    ifstream infile(outfile);
    if (!infile.is_open()) {
        cerr << "Erro: Não foi possível abrir o arquivo de saída " << outfile << endl;
        return;
    }

    string line;
    bool readingAbundances = false;
    bool skippingEmptyLine = false;

    while (getline(infile, line)) {
        // Verifica se estamos na seção "Isotope Abundance"
        if (line.find("Isotope Abundance") != string::npos) {
            readingAbundances = true;
            continue;  // Pula a linha "Isotope Abundance"
        }

        if (readingAbundances) {
            // Se a linha for vazia, apenas pula
            if (line.empty()) {
                if (skippingEmptyLine) {
                    break;  // Se já pulou uma linha vazia, agora é hora de parar
                } else {
                    skippingEmptyLine = true;
                    continue;  // Pula a primeira linha vazia
                }
            }

            // Caso contrário, a linha contém dados de isótopos
            stringstream ss(line);
            string isotope;
            float abundance;

            ss >> isotope >> abundance;
            if (!isotope.empty()) {
                // Remove o sufixo do isótopo (se houver) e mantém o número
                string isotopeNumber = isotope.substr(0, isotope.find_first_not_of("0123456789"));
                
                    isotopeAbundances[isotopeNumber] = abundance;
                    cout << "Isotope: " << isotopeNumber << ", Abundance: " << abundance << endl;
                
            }
        }
    }

    infile.close();
}



void plotMeAxsSeveral(char particle = 'p', float angle = 20.,  bool same = true, string path = "/mnt/medley/LucasAnalysis/talys_calcs/", string outfile = "output.out"){

    if(path == ""){
        cerr<<"Error: No PATH name provided."<<endl;
        return;
    }

        string outputfile = path + outfile;
        cout<<"Output file: "<<outputfile<<endl;
        ifstream file(outputfile);
        if(!file){
            cerr<<"Error: File not found."<<endl;
            return;
        }
        
        map<string, float> isoptopeAbundances;
        parseIso(isoptopeAbundances, outputfile);
/*

        TTree *tx = new TTree("xs","xs");
        //tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MPreeq/F:Compound/F:PreeqR/F:BKr/F:Str/F:KO/F:BrU/F");
        //##     E-out           xs           Direct     Preequilibrium Multiple_preeq    Compound
        tx->ReadFile(fullPath.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MultPreeq/F:Compound/F");
        Float_t Eprod, total;
        Long64_t nEntries = tx->GetEntries();
        tx->SetBranchAddress("Eprod", &Eprod);
        tx->SetBranchAddress("Total", &total);


        TGraph *gr = new TGraph();
        for(Long64_t i = 0; i < nEntries; ++i) {
            tx->GetEntry(i);

            cout<<Eprod<<" "<<total<<endl;
            gr->AddPoint(Eprod, total);
            // Plotting code here
            // For example:
            // hist->Fill(Eplot, total);
        }
        if(same){
            gr->Draw("same");
        }else{
            gr->Draw();
        }

    }
*/

}