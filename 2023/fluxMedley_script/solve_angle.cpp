
#include "TF1.h"
#include "TH1D.h"
#include "TVirtualFitter.h"
#include "TPaletteAxis.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TLatex.h"
#include "TStopwatch.h"
#include <string>
#include <iostream>
#include <fstream>
#include <Math/RootFinderAlgorithms.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
using namespace ROOT::Math;
using namespace std;
#include "/home/e802/Analysis/pre_analysis/libs/useful.h"



const double mn = 939.5654133; //MeV/c²
const double mp = 938.2720813; //MeV/c²

double FindTheta(double u =1, int solution=+1)
{


TF1 *f = new TF1("f", "-[0] -[3]*([1]/[2])*( cos(x)*sqrt( pow([2]/[1],2)- pow(sin(x),2) ) -pow(sin(x),2) )", 0, 3.142);

// 	//defining parameters:
//     //[0]-> u
//     //[1]-> m_n.c^2 [N]
//     //[2]-> m_p.c^2 [P]
//     //[3]-> solution (+1 or -1)


 	f->SetParameter(0,u);

 	//mass of the neutron
     f->SetParameter(1,mn);

	//mass of the proton
     f->SetParameter(2,mp);
     
     f->SetParameter(3,solution);

     //f->Draw();
    //gPad->SetGridx();
    //gPad->SetGridy();

    RootFinder *k = new RootFinder();
    k->SetMethod(RootFinder::kBRENT);
    k->Solve(*f,0,3.142);

    double c = k->Root();
    //cout << c*180/Pi() << endl;

     return c;
}
using namespace TMath;


void solve_angle(){

vector<double> a_values = {20,40,60,80};
vector<double> u_values = {-0.8816,-0.5857,-0.2560,-0.0445};

vector <double> theta1;
vector <double> theta2;
cout<<" Formula D.70 for theta calculation based on mu (cos(β)). See https://www.oecd-nea.org/dbdata/data/manual-endf/endf102.pdf"<<endl;
for(int i=0;i<a_values.size();i++){

    theta1.push_back(FindTheta(u_values[i],1));
    theta2.push_back(FindTheta(u_values[i],-1));
    cout<<Form("%04.1f° [µ= %06.4f]: %04.3frad, %04.3frad; %04.3f°, %04.3f°",a_values[i],u_values[i],theta1[i],theta2[i],theta1[i]*180/TMath::Pi(),theta2[i]*180/TMath::Pi())<<endl;
}


cout<<" \n.\n.\nFormula D.71 (same reference):"<<endl;
vector <double> dcosb_dcostheta1;
vector <double> dcosb_dcostheta2;
for(int i=0;i<a_values.size();i++){
    dcosb_dcostheta1.push_back( 2*cos(theta1[i])*mn/mp + ( 1+ (2*pow(cos(theta1[i]),2) -1)*pow(mn/mp,2) )/( sqrt( 1 - pow(sin(theta1[i]),2)*pow(mn/mp,2) ) ) );
    dcosb_dcostheta2.push_back( 2*cos(theta2[i])*mn/mp + ( 1+ (2*pow(cos(theta2[i]),2) -1)*pow(mn/mp,2) )/( sqrt( 1 - pow(sin(theta2[i]),2)*pow(mn/mp,2) ) ) );
    cout<<Form("%04.1f° [µ= %06.4f]: %04.3f, %04.3f",a_values[i],u_values[i],dcosb_dcostheta1[i],dcosb_dcostheta2[i] )<<endl;
}




}