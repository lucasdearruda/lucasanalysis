#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_ddx_mRuns.cpp" 
using namespace std;

std::vector<TH1D*> Fe_allAngles2(Float_t Ea = 25, Float_t Eb = 26) {
    //TH1::AddDirectory(kFALSE);
    std::vector<TH1D*> h(8);

    vector<Int_t> runs_a = {396};//, 47, 29, 396, 396, 396, 396, 396};
    vector<Int_t> runs_b = {405};// 49, 34, 405, 405, 405, 405, 405};
    
    vector<Int_t> runs_aR = {410};//, 396, 396, 396, 396, 396, 396, 396};
    vector<Int_t> runs_bR = {412};//, 405, 405, 405, 405, 405, 405, 405};
    

    h[0] = Fe_ddx_mRuns(Ea, Eb, 1, 20.0, 'p', runs_a, runs_b);
    h[1] = Fe_ddx_mRuns(Ea, Eb, 1, 40.0, 'p', runs_a, runs_b);
    h[2] = Fe_ddx_mRuns(Ea, Eb, 1, 60.0, 'p', runs_a, runs_b);
    h[3] = Fe_ddx_mRuns(Ea, Eb, 1, 80.0, 'p', runs_a, runs_b);
    h[4] = Fe_ddx_mRuns(Ea, Eb, 1, 100.0, 'p', runs_aR, runs_bR);  
    h[5] = Fe_ddx_mRuns(Ea, Eb, 1, 120.0, 'p', runs_aR, runs_bR);
    h[6] = Fe_ddx_mRuns(Ea, Eb, 1, 140.0, 'p', runs_aR, runs_bR);
    h[7] = Fe_ddx_mRuns(Ea, Eb, 1, 160.0, 'p', runs_aR, runs_bR);

     TCanvas *cv = new TCanvas("cv","cv");
    h[0]->Draw();
    for(int i=1;i<8;i++){
        h[i]->Draw("same");
    }
    delete cv;
    return h;
}

void histosPlot(std::vector<TH1D*> h){
    TCanvas *cv = new TCanvas("cv","cv");
    //h[0]->Draw();

    TCanvas *cc = new TCanvas("cc","cc");

    TH2D *h2 = new TH2D("h2","h2",h[0]->GetNbinsX(),h[0]->GetBinLowEdge(0),h[0]->GetBinLowEdge(h[0]->GetNbinsX()),180,0,180);
    Int_t angles[8] = {20,40,60,80,100,120,140,160};
    for(int i=0;i<8;i++){
        for(int b= 1;b<=h[i]->GetNbinsX();b++){
            h2->SetBinContent(b,angles[i],h[i]->GetBinContent(b));
        }
    }
    h2->Draw("lego2");
    h2->GetXaxis()->SetTitle("E_{p} (MeV)");
    h2->GetYaxis()->SetTitle("Angle (deg)");
    h2->GetZaxis()->SetTitle("mb/sr#dot1-MeV");
    
}