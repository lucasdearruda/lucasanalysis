#include "/mnt/medley/LucasAnalysis/2023/XS_calcs/Fe_allAngles2.cpp" 

#include <vector>
#include <iostream>
#include "TStopwatch.h"

void XS_angularDX2(){
    //TH1::AddDirectory(kFALSE);
    TStopwatch timer;
    TCanvas *cc = new TCanvas("cc","cc");

    std::vector<TH1D*> h;


    TGraph2D * XSgr = new TGraph2D();
    Float_t intg;

    Float_t Ezero = 5;
    for (int i = Ezero; i < 38; i=i+3) {
 
        h = Fe_allAngles2(i, i + 3);
        for(int j=0;j<8;j++){
            intg = h[j]->Integral("width");
            XSgr->AddPoint(20+20*j,i+Ezero+1.5,intg);
            std::cout << "Angle: " << 20 + 20 * j << " MeV: " << i + 2 << " Value: " << intg << std::endl;
        }


    }

    XSgr->Draw();


    timer.Print();
}