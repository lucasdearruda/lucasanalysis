#include <TFile.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>  // mkdir
#include <sys/types.h>

// Função para extrair energia de nêutrons do nome do histograma


float EfromName(const std::string& name) {
    size_t posEn = name.find("En_");
    if (posEn == std::string::npos) return -1;

    size_t posMid = name.find('_', posEn + 3);
    size_t posMeV = name.find("_MeV", posMid);
    if (posMid == std::string::npos || posMeV == std::string::npos) return -1;

    std::string low_str = name.substr(posEn + 3, posMid - (posEn + 3));
    std::string high_str = name.substr(posMid + 1, posMeV - (posMid + 1));

    // substitui 'p' por '.'
    for (auto& c : low_str) if (c == 'p') c = '.';
    for (auto& c : high_str) if (c == 'p') c = '.';

    float low = atof(low_str.c_str());
    float high = atof(high_str.c_str());

    return (low + high) / 2.0f;
}

// Função para ler arquivo TALYS e retornar TGraph
int AngleFromName(const std::string& name) {
    size_t posDeg = name.find("deg");
    if (posDeg == std::string::npos) return -1;

    // Procura o underscore antes do "deg"
    size_t posUnder = name.rfind('_', posDeg);
    if (posUnder == std::string::npos) return -1;

    // Extrai substring entre o underscore e "deg"
    std::string angle_str = name.substr(posUnder + 1, posDeg - posUnder - 1);

    return atoi(angle_str.c_str());
}


std::vector<TGraph*> plotMeAxs(const std::string& namefile = "", bool same = true, bool convertToCounts = false, Float_t binWidth = 1.0) {
    std::vector<TGraph*> graphs;

    if (namefile == "") {
        std::cerr << "Error: No file name provided." << std::endl;
        return graphs;  // vetor vazio
    }

    TTree* tx = new TTree("xs", "xs");
    // Lê as colunas do arquivo
    tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MultPreeq/F:Compound/F");

    Float_t Eprod, total, direct, preeq, multpreeq, compound;
    Long64_t nEntries = tx->GetEntries();

    tx->SetBranchAddress("Eprod", &Eprod);
    tx->SetBranchAddress("Total", &total);
    tx->SetBranchAddress("Direct", &direct);
    tx->SetBranchAddress("Preeq", &preeq);
    tx->SetBranchAddress("MultPreeq", &multpreeq);
    tx->SetBranchAddress("Compound", &compound);

    // Cria os 5 TGraphs
    TGraph* grTotal     = new TGraph();
    TGraph* grDirect    = new TGraph();
    TGraph* grPreeq     = new TGraph();
    TGraph* grMultPreeq = new TGraph();
    TGraph* grCompound  = new TGraph();

    grTotal->SetName("Total");
    grDirect->SetName("Direct");
    grPreeq->SetName("Preeq");
    grMultPreeq->SetName("MultPreeq");
    grCompound->SetName("Compound");

    grTotal->SetLineWidth(2);
    grDirect->SetLineWidth(2);
    grPreeq->SetLineWidth(2);
    grMultPreeq->SetLineWidth(2);
    grCompound->SetLineWidth(2);

    for (Long64_t i = 0; i < nEntries; ++i) {
        tx->GetEntry(i);

        float norm = convertToCounts ? (1.0 / binWidth) : 1.0;

        grTotal->AddPoint(Eprod, total * norm);
        grDirect->AddPoint(Eprod, direct * norm);
        grPreeq->AddPoint(Eprod, preeq * norm);
        grMultPreeq->AddPoint(Eprod, multpreeq * norm);
        grCompound->AddPoint(Eprod, compound * norm);
    }

    // Desenha os gráficos
    if (same) {
        grTotal->Draw();
        grDirect->Draw("same");
        grPreeq->Draw("same");
        grMultPreeq->Draw("same");
        grCompound->Draw("same");
    } else {
        grTotal->Draw();
        grDirect->Draw();
        grPreeq->Draw();
        grMultPreeq->Draw();
        grCompound->Draw();
    }

    graphs.push_back(grTotal);
    graphs.push_back(grDirect);
    graphs.push_back(grPreeq);
    graphs.push_back(grMultPreeq);
    graphs.push_back(grCompound);

    return graphs;
}

TGraph* plotMeAxsOld(string namefile = "", bool same = true, bool convertToCounts = false,       Float_t binWidth = 1.0){
    if(namefile == ""){
        cerr<<"Error: No file name provided."<<endl;
        return nullptr;
    }else{


        TTree *tx = new TTree("xs","xs");
        //tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MPreeq/F:Compound/F:PreeqR/F:BKr/F:Str/F:KO/F:BrU/F");
        //##     E-out           xs           Direct     Preequilibrium Multiple_preeq    Compound
        tx->ReadFile(namefile.c_str(), "Eprod/F:Total/F:Direct/F:Preeq/F:MultPreeq/F:Compound/F");
        Float_t Eprod, total;
        Long64_t nEntries = tx->GetEntries();
        tx->SetBranchAddress("Eprod", &Eprod);
        tx->SetBranchAddress("Total", &total);


        TGraph *gr = new TGraph();
 

        for(Long64_t i = 0; i < nEntries; ++i) {
            
            tx->GetEntry(i);
            
            if(convertToCounts){   
                // Convert to counts
                //cout<<"Converting to counts: " << total << " / " << binWidth <<" = "<<total/binWidth<< endl;
                gr->AddPoint(Eprod, total/binWidth);
            }else{
                gr->AddPoint(Eprod, total);
            }

            
            // Plotting code here
            // For example:
            // hist->Fill(Eplot, total);
        }
        gr->SetName("gr");
        gr->SetLineWidth(2);
        if(same){
            gr->Draw("same");
        }else{
            gr->Draw();
        }
        return gr;
    }


}
void saveCanvasIfRequested(TCanvas* c, float targetEnergy, int targetAngle, bool saveFig) {
    if (!saveFig) return;

    // Criar pasta results/Comp se não existir
    struct stat info;
    if (stat("results/Comp", &info) != 0) {  // pasta não existe
        int status = mkdir("results/Comp", 0777);
        if (status != 0) {
            std::cerr << "Erro ao criar diretório results/Comp" << std::endl;
            return;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        std::cerr << "results/Comp existe mas não é diretório" << std::endl;
        return;
    }

    // Montar nome do arquivo
    TString filename = Form("results/Comp/Fe_En%.1fdeg_%d.png", targetEnergy, targetAngle);
    // Aumentar resolução PNG (canvas maior)
    //int w = 1600, h = 1200;
    //c->SetCanvasSize(w, h);

    TString baseName = Form("results/Comp/Fe_En%.1fdeg_%d", targetEnergy, targetAngle);

    c->SaveAs(baseName + ".png");   // PNG alta res (~300 dpi)
    c->SaveAs(baseName + ".pdf");   // PDF vetorial
    c->SaveAs(baseName + ".root");  // ROOT file
}
void plot_DDX_wRaw(bool saveFig = true,  float targetEnergy = 10.5, int targetAngle = 40) {
    // // Parâmetros
    // float targetEnergy = 10.5;
    // int targetAngle = 40;

    // Caminhos
    std::string rawPath = "Fe_pXS_3MeVbin_3to40MeV_rawE.root";
    std::string newEPath = "Fe_pXS_3MeVbin_3to40MeV.root";
    std::string talysBase = "/mnt/medley/LucasAnalysis/talys_calcs/Fe56/best_9June/";
    

    
    std::string talysFullPath = talysBase + Form("pddxE%08.3fA%05.1f.deg", targetEnergy, (float)targetAngle);

    // Abrir arquivos ROOT
    TFile* f_raw = TFile::Open(rawPath.c_str());
    TFile* f_proc = TFile::Open(newEPath.c_str());
    if (!f_raw || !f_proc) {
        std::cerr << "Erro ao abrir arquivos ROOT!" << std::endl;
        return;
    }

    // Procurar histograma correspondente
    TH1D* h_raw = nullptr;
    TH1D* h_proc = nullptr;

    TIter next_raw(f_raw->GetListOfKeys());
    TObject* obj_raw;
    while ((obj_raw = next_raw())) {
        std::string name = obj_raw->GetName();
        cout<< name<<endl;
        if (name.find("Fe_thick") != std::string::npos) {
            float e = EfromName(name);
            int a = AngleFromName(name);
            if (a == targetAngle && fabs(e - targetEnergy) < 0.6) {
                h_raw = (TH1D*)f_raw->Get(name.c_str());
                break;
            }
        }
    }

    TIter next_proc(f_proc->GetListOfKeys());
    TObject* obj_proc;
    while ((obj_proc = next_proc())) {
        std::string name = obj_proc->GetName();
        if (name.find("Fe_thick") != std::string::npos) {
            float e = EfromName(name);
            int a = AngleFromName(name);
            if (a == targetAngle && fabs(e - targetEnergy) < 0.6) {
                h_proc = (TH1D*)f_proc->Get(name.c_str());
                break;
            }
        }
    }

    if (!h_raw || !h_proc) {
        std::cerr << "Histograma não encontrado!" << std::endl;
        return;
    }

    // Ler TALYS
    // TGraph* g_talys = plotMeAxs(talysFullPath, false, false, 0.1);
    // if (!g_talys) return;

    // // Estilos
    // h_raw->SetLineColor(kRed);
    // h_proc->SetLineColor(kBlue);
    // g_talys->SetLineColor(kGreen+2);
    // g_talys->SetLineWidth(2);

    // // Plot
    // TCanvas* c = new TCanvas("c", "Comparação ROOT vs TALYS", 800, 600);
    // h_raw->SetTitle(Form("Fe(n,Xp) - %d deg - En=%.1f MeV", targetAngle, targetEnergy));
    // h_raw->GetXaxis()->SetTitle("Ep (MeV)");
    // h_raw->GetYaxis()->SetTitle("Counts (a.u.)");

    // h_raw->Draw("HIST E1");
    // h_proc->Draw("HIST SAME E1");
    // g_talys->Draw("L SAME");

    // TLegend* leg = new TLegend(0.6,0.7,0.88,0.88);
    // leg->AddEntry(h_raw, "Raw data", "l");
    // leg->AddEntry(h_proc, "Processed data", "l");
    // leg->AddEntry(g_talys, "TALYS", "l");
    // leg->Draw();

    // c->Update();


    // Ler TALYS (vários TGraphs)
   std::vector<TGraph*> g_talys_vec = plotMeAxs(talysFullPath, false, false, 0.1);
if (g_talys_vec.empty()) return;

// Estilos - cores e estilos definidos
const int colors[] = {kGreen+2, kRed, kBlack, kBlue, kMagenta};
const int styles[] = {1, 3, 2, 3, 2};  // 1=solid, 2=dashed, 3=dotted

for (size_t i = 0; i < g_talys_vec.size(); ++i) {
    g_talys_vec[i]->SetLineColor(colors[i]);
    g_talys_vec[i]->SetLineStyle(styles[i]);
    g_talys_vec[i]->SetLineWidth(2);
}

// Estilo dos histogramas
h_raw->SetLineColor(kRed);
h_proc->SetLineColor(kBlue);

// Ajustar eixo Y para acomodar tudo
double maxY = h_raw->GetMaximum();
maxY = std::max(maxY, h_proc->GetMaximum());
for (auto& g : g_talys_vec) {
    double gxmax = g->GetHistogram() ? g->GetHistogram()->GetMaximum() : 0;
    maxY = std::max(maxY, gxmax);
}
maxY *= 1.1; // margem de 20%

h_raw->SetMaximum(maxY);

// Labels e títulos maiores
h_raw->GetXaxis()->SetTitleSize(0.05);
h_raw->GetXaxis()->SetLabelSize(0.045);
h_raw->GetYaxis()->SetTitleSize(0.05);
h_raw->GetYaxis()->SetLabelSize(0.045);

h_raw->SetTitle(Form("Fe(n,Xp) - %d deg - En=%.1f MeV", targetAngle, targetEnergy));
h_raw->GetXaxis()->SetTitle("Ep (MeV)");
h_raw->GetYaxis()->SetTitle("d^{2}#sigma / d#Omega dE (mb/sr/MeV)");

// Criar canvas e desenhar
TCanvas* c = new TCanvas("c", "comparison ROOT vs TALYS", 800, 600);
c->SetLeftMargin(0.15); 
c->SetRightMargin(0.05); 
h_raw->Draw("HIST E1");
h_proc->Draw("HIST SAME E1");

gPad->SetGrid();

g_talys_vec[0]->Draw("L SAME"); // Total - verde sólido
for (size_t i = 1; i < g_talys_vec.size(); ++i) {
    g_talys_vec[i]->Draw("L SAME");
}

// Legenda com textos maiores
TLegend* leg = new TLegend(0.45, 0.44, 0.93, 0.88);
leg->SetTextSize(0.04);
leg->SetFillColor(kWhite);         // fundo branco
//leg->SetFillStyle(3001);           // fundo com transparência leve
leg->SetBorderSize(1);             // borda visível
leg->AddEntry(h_raw, "Raw data", "l");
leg->AddEntry(h_proc, "Processed data", "l");
leg->AddEntry(g_talys_vec[0], "TALYS Total", "l");
leg->AddEntry(g_talys_vec[1], "TALYS Direct", "l");
leg->AddEntry(g_talys_vec[2], "TALYS Preeq", "l");
leg->AddEntry(g_talys_vec[3], "TALYS MultPreeq", "l");
leg->AddEntry(g_talys_vec[4], "TALYS Compound", "l");
leg->Draw();


c->Update();

saveCanvasIfRequested(c, targetEnergy, targetAngle, saveFig);

}
