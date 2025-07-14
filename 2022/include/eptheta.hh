
//#include <iostream>
#include <Math/RootFinderAlgorithms.h>
#include <TF1.h>
#include <Math/RootFinder.h>
#include <Math/Functor.h>
using namespace ROOT::Math;
using namespace std;

double EptTheta(double TOF, double theta, int tel=1)
{
//these calculations just make sense for ELASTIC SCATTERING (because we assume Ep = En.cos^2(theta) )
//return the ENERGY of a elastic scattered proton

	//equation: TOF[event] = ToF[NN] + ToF[proton]


	//for elastic scattering:
	double M = 1.0; // mass of recoil nucleus (target) in u
	double mp = 1.007276467; //mass o proton (scattered) in u

	//E_particle = E_neutron*Factor (factor is the share of energy lost by the neutron)
	double Factor = 4*M*mp*pow(cos(TMath::Pi()*theta/180),2)/pow(M+mp,2);
	//cout<<Form("Factor = %f", Factor)<<endl;

	// ct = [ L/sqrt( 1 - ( m_n.c^2/(E_p/Factor+m_n.c^2) )^2 ) +  L1/sqrt( 1 - ( m_p.c^2/(E_p+m_p.c^2) )^2 ) ]

	TF1 *f = new TF1("f", "[0]/sqrt(1 - pow( [2]/( x/[4]  +[2] )   ,2) )+[1]/sqrt( 1 - pow( [3]/( x + [3])  ,2)   )-[5]*[6]", 0, 200);

	//defining parameters:
    //[0]-> L
    //[1]-> L1
    //[2]-> m_n.c^2 [N]
    //[3]-> m_p.c^2 [P]
    //[4]-> Factor
    //[5]-> c
    //[6]-> t


	//Distance between rotating converter and medley target (L):
	double L = 5.044;

	//Distance between medley target and telescope:
	double L1; //distance from the target to Si2


	switch(tel){ //which depends on the telescope we are using.
		case 1:
			L1=0.2184;
			break;
		case 2:
			L1=0.1609;
			break;
		case 3:
			L1=0.1609;
			break;
		case 4:
			L1=0.1609;
			break;
		case 5:
			L1=0.1739;
			break;
		case 6:
			L1=0.1739;
			break;
		case 7:
			L1=0.1709;
			break;
		case 8:
			L1=0.2079;
			break;
		default:
			L1=0.2184;
	}

	f->SetParameter(0,L);

	f->SetParameter(1,L1);

	//mass of the proton
    f->SetParameter(2,938.2720813);

	//mass of the neutron
    f->SetParameter(3,939.5654133);

	//angle convertion: degree -> radian
    f->SetParameter(4,Factor);

	//velocity of light
	f->SetParameter(5,0.299792458);

	//TOF
    f->SetParameter(6,TOF);

    //f->Draw();
    //gPad->SetGridx();
    //gPad->SetGridy();


     RootFinder *k = new RootFinder();
       k->SetMethod(RootFinder::kBRENT);
     k->Solve(*f,0,20000);

     double c = k->Root();
    //cout << c << endl;

    return c;
}

