#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_allAngles2.cpp" 

#include <vector>
#include <iostream>
#include "TStopwatch.h"

void XS_angularDX(){
    //TH1::AddDirectory(kFALSE);
    TStopwatch timer;
    TCanvas *cc = new TCanvas("cc","cc");
    
    std::vector<std::vector<TH1D*>> allHistos;


//    std::vector<TH1D*> h = Fe_allAngles2(25,26);

    Float_t Ezero = 5;
    for (int i = Ezero; i < 40; i++) {
        std::vector<TH1D*> h = Fe_allAngles2(i, i + 1);
        allHistos.push_back(h);
        //gROOT->GetListOfCanvases()->Delete();
    }


    TGraph2D * XSgr = new TGraph2D();
    Float_t intg;
    // Exemplo de acesso:
    for (size_t i = 0; i < allHistos.size(); ++i) {
        std::cout << "Pair " << (i + Ezero) << "," << (i + Ezero+1) << ": "
                  << allHistos[i].size() << " histograms" << std::endl;
        
        for(int j=0;j<allHistos[i].size();j++){
            intg = allHistos[i][j]->Integral("width");
            XSgr->AddPoint(i+Ezero+0.5,20+20*j,intg);
            std::cout << "Angle: " << 20 + 20 * j << " MeV: " << i + 2 << " Value: " << intg << std::endl;
        }
    }

    XSgr->Draw();


    timer.Print();
}