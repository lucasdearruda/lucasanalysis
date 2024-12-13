// #include "TCanvas.h"
// #include "TH1.h"
// #include "TStyle.h"
// #include "TROOT.h"

// Função com parâmetros para títulos dos eixos e range do eixo X
void LoadVisual(TCanvas* NNCanvas, TH1* hneutrons_unc,
                const char* xTitle = "E (MeV)", const char* yTitle = "n  sr^{-1} 1-MeV^{-1} #muC^{-1}",
                double Xi = 0.0, double Xf = 45.0) {
    
    // Configurações de estilo globais
    gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
    gStyle->SetTextSize(0.04);
    gStyle->SetTitleYSize(0.04);
    
    // Configurações do canvas
    NNCanvas->cd()->SetRightMargin(0.03);
    NNCanvas->cd()->SetTopMargin(0.08);
    NNCanvas->cd()->SetLeftMargin(0.14);
    NNCanvas->cd()->SetBottomMargin(0.12);
    NNCanvas->SetWindowSize(790, 652);

    // Configurações do histograma
    hneutrons_unc->GetXaxis()->SetRangeUser(Xi, Xf);        // Define o intervalo do eixo X
    hneutrons_unc->GetYaxis()->SetTitle(yTitle);            // Define o título do eixo Y
    hneutrons_unc->GetXaxis()->SetTitle(xTitle);            // Define o título do eixo X
    hneutrons_unc->SetTitle("");                            // Título do histograma vazio

    // Ajustar tamanhos das fontes dos eixos
    hneutrons_unc->GetXaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetYaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetXaxis()->SetTitleSize(0.06);
    hneutrons_unc->GetYaxis()->SetTitleSize(0.06);
    
    // Atualizar o canvas com as novas configurações
    NNCanvas->Modified();
    NNCanvas->Update();
}
// Função com parâmetros para títulos dos eixos e range dos eixos X e Y
void LoadVisual_th2(TCanvas* NNCanvas, TH2* hneutrons_unc,
                    const char* xTitle = "X axis title", const char* yTitle = "Y axis title",
                    double Xi = 0.0, double Xf = 45.0, double Yi = 0.0, double Yf = 45.0) {

    // Configurações de estilo globais
    gStyle->SetOptStat(0); // Desativar a caixa de estatísticas
    gStyle->SetTextSize(0.04);
    gStyle->SetTitleYSize(0.04);
    
    // Configurações do canvas
    NNCanvas->cd()->SetRightMargin(0.13);
    NNCanvas->cd()->SetTopMargin(0.08);
    NNCanvas->cd()->SetLeftMargin(0.14);
    NNCanvas->cd()->SetBottomMargin(0.14);
    NNCanvas->SetWindowSize(875, 554);
    //hneutrons_unc->GetZaxis()->SetMaxDigits(3);
    auto palette = (TPaletteAxis*)hneutrons_unc->GetListOfFunctions()->FindObject("palette");
    palette->SetX1NDC(0.866762);
    palette->SetY1NDC(0.14285714);
    palette->SetX2NDC(0.89828080);
    palette->SetY2NDC(0.922269);

    // Configurações do histograma
    hneutrons_unc->GetXaxis()->SetRangeUser(Xi, Xf);        // Define o intervalo do eixo X
    hneutrons_unc->GetYaxis()->SetRangeUser(Yi, Yf);        // Define o intervalo do eixo Y
    hneutrons_unc->GetYaxis()->SetTitle(yTitle);            // Define o título do eixo Y
    hneutrons_unc->GetXaxis()->SetTitle(xTitle);            // Define o título do eixo X
    hneutrons_unc->SetTitle("");                            // Título do histograma vazio
    
    

    // Ajustar tamanhos das fontes dos eixos
    hneutrons_unc->GetXaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetYaxis()->SetLabelSize(0.06);
    hneutrons_unc->GetXaxis()->SetTitleSize(0.06);
    hneutrons_unc->GetYaxis()->SetTitleSize(0.06);
    hneutrons_unc->GetZaxis()->SetLabelSize(0.06);
    
    // Atualizar o canvas com as novas configurações
    NNCanvas->Modified();
    NNCanvas->Update();

}