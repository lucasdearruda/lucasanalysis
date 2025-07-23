//string listfiles(const char *ext="Rx.root",const char *dirname="/mnt/c/Users/lucas_det5oo0/neutron/root2rx/agosto2020" ){

// #include <iostream>     // std::cout
//#include <algorithm>    // std::sort
//#include <vector>       // std::vector


bool compareFunction (string a, string b) {return a<b;}
//compare any way you like, here I am using the default string comparison

string listfiles(const char *dirname=getenv("PWD"),const char *ext="" ){
//show the files and the return the name of the chosen one
vector <string> names;

TSystemDirectory dir(dirname, dirname);
TList *files = dir.GetListOfFiles();
		int i;
	if(files){
		i=1;
		cout << "  n  |" <<'\t'<< "Filename" << endl;
		TSystemFile *file;
		TString fname;
		TIter next(files);
			while((file=(TSystemFile*)next())){
				fname = file->GetName();
				if(!file->IsDirectory() && fname.EndsWith(ext)){
					//cout << Form("%5d", i) <<"|\t"<< fname.Data() << endl;
					names.push_back(fname.Data());
					i++;
				}
			}
		}
		int number=-1;

		bool b;
		b= false;
		while(b == false){
			b = true;
			for(int p=0;p<names.size()-1;p++){
				if(names[p]>names[p+1]){
					b = false;
					swap(names[p], names[p+1]);
				}
			}

		}
		for(int p =0; p<names.size();p++){
			cout << Form("%5d", p+1) <<"|\t"<< names[p] << endl;
		}
		while(number<1||number>=i){
			cout <<Form("Inform number of file to be opened (%d to %d): ",1, i-1)<< endl;
			cin>>number;
		}//cout <<names[number-1]<< endl;
		return names[number-1];
}


// 	bool b;
// 	b= false;
// 	while(b == false){
// 		b=true;
// 		for(int i=0; i<calib.size()-1;i++){
// 			if(calib[i].detec>calib[i+1].detec){
// 				b = false;
// 				swap(calib[i],calib[i+1]);
// 			}
// 		}
// 	}//organize detectors


void showfiles(const char *dirname=getenv("PWD"),const char *ext="" ){
//show the files names with the numbers, but doesnt return nothing
vector <string> names;

TSystemDirectory dir(dirname, dirname);
TList *files = dir.GetListOfFiles();
		int i;
	if(files){
		i=1;
		cout << "  n  |" <<'\t'<< "Filename" << endl;
		TSystemFile *file;
		TString fname;
		TIter next(files);
			while((file=(TSystemFile*)next())){
				fname = file->GetName();
				if(!file->IsDirectory() && fname.EndsWith(ext)){
					cout << Form("%5d", i) <<"|\t"<< fname.Data() << endl;
					names.push_back(fname.Data());
					i++;
				}
			}
		}
}

vector<string> namesfiles(const char *dirname=getenv("PWD"),const char *ext=""){
//return a vector with N names
//the inputs are:
//- name of the directory [dirname]
//- extension of wanted file [ext]
vector <string> names;
vector <string> chosen_names;


TSystemDirectory dir(dirname, dirname);
TList *files = dir.GetListOfFiles();
		int i;
	if(files){
		i=1;
		cout << "  n  |" <<'\t'<< "Filename" << endl;
		TSystemFile *file;
		TString fname;
		TIter next(files);
			while((file=(TSystemFile*)next())){
				fname = file->GetName();
				if(!file->IsDirectory() && fname.EndsWith(ext)){
					//cout << Form("%5d", i) <<"|\t"<< fname.Data() << endl;
					names.push_back(fname.Data());
					i++;
				}
			}
			sort(names.begin(),names.end(),compareFunction);
			for(int t=0;t<i-1;t++){
				cout << Form("%5d", t+1) <<"|\t"<<names[t]<< endl;
			}
		}
		int number=-1;

		// cout <<Form("Inform number of file to be opened (%d to %d): ",1, i-1)<< endl;
		// string line;
		// getline( cin, line );
		// istringstream is( line );
		// int n;
		// while( is >> n ) {
		// 	if(n>0 && n<i){
		// 		cout<<"File "<<n<<" selected: ";
		// 		chosen_names.push_back(names[n-1]);
		// 		cout<<chosen_names.back()<<endl;
		// 	}else{
		// 		cout<<"-> File "<<n<<" does not exist. Ignored."<<endl;
		// 	}
		// }

		//return chosen_names;
		return names;
}

double gauswbg(double * x, double * p){
		double arg = (x[0] - p[1])/p[2];

		return  p[0]*TMath::Exp(-0.5*arg*arg)  + p[3];
}



double FWHM(TH1* h1){
	int bin1 = h1->FindFirstBinAbove(h1->GetMaximum()/2);
   int bin2 = h1->FindLastBinAbove(h1->GetMaximum()/2);
   double fwhm = h1->GetBinCenter(bin2) - h1->GetBinCenter(bin1);
   return fwhm;
}

double erf(double * x, double * p){
		return  TMath::Erf((x[0]-p[0])/p[1])*p[2]   + p[3];
}
double erfc(double * x, double * p){
		return  TMath::Erfc((x[0]-p[0])/p[1])*p[2]  + p[3];
}


double linear(double x, double A, double B){
		return  A*x + B;
}

void sortAsc(double x[], double v[], int s){
	//void sortAsc(vector a, vector b, vector size (number of elements) )
	// sort double vector X ascending, doing the same trnasformations on vector V
	// example: let x[] = {3, 7, 4, 8, 2, 9, 1} and b[] = {0, 1, 2, 3, 4, 5, 6}.
	//calling sortAsc(x,b,7) will result:
	// x ={1,2,3,4,7,8,9} and b={6,4,0,2,1,3,5}
	bool b= false;
	double aux, aux1;
	while(b == false){
		b=true;
		for(int i=0; i<s-1;i++){
			if(x[i]>x[i+1]){
				b = false;
				aux = x[i]; aux1= v[i];
				x[i] = x[i+1]; v[i] = v[i+1];
				x[i+1] = aux; v[i+1]=aux1;
			}
		}
	}
}

struct cal
{
  char detec;//A or B or C
  int Tel;//1-8
  float A, B, sA, sB;
};




double nneg(double x){
	if(x<0){
		return 0;
	}else{
		return x;
	}

}

vector <cal> sortCal(vector <cal> calib){
	//sort calibration vector as:
	//Telescope 1> det 1 ; det 2; det3 > Telescope 2 > det 1 ; ....
	// calib is the vector
	bool b;
	b= false;
	while(b == false){
		b=true;
		for(int i=0; i<calib.size()-1;i++){
			if(calib[i].detec>calib[i+1].detec){
				b = false;
				swap(calib[i],calib[i+1]);
			}
		}
	}//organize detectors
	b= false;
	while(b == false){
		b=true;
		for(int i=0; i<calib.size()-1;i++){
			if(calib[i].Tel>calib[i+1].Tel){
				b = false;
				swap(calib[i],calib[i+1]);
			}
		}
	}//organize telescopes

	return calib;
}

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


double EpToFL(double TOF, double Ep)
{
//return the DISTANCE L in centimeters, of a proton with energy Ep and ToF tof.
	// ct = [ L/sqrt( 1 - ( m_n.c^2/(E_p/cos^2(theta)+m_n.c^2) )^2 )

	TF1 *f = new TF1("f", "x*[2]*sqrt(1 - pow( [0]/([1] +[0]),2) ) ", 0, 1e6);
	//defining parameters:
    //x -> L
    //[0]-> m_p.c^2 [N]
	//[1]-> Ep
    //[2]-> c

	//mass of the proton*c^2
    f->SetParameter(0,938.2720813);

	//Energy of the proton
    f->SetParameter(1,Ep);

	//velocity of light
	f->SetParameter(2,0.299792458);

	//TOF
    f->SetParameter(3,TOF);

    return 100*f->Eval(TOF);
}

double EptThetaM(double tof, char particle, double Ep)
{
//returns the energy of the neutron
//particle = 'p', 'd', 't', 'h', 'a' for each particle
//Ep = particle's energy
    TF1 *f = new TF1("f", "[0]/sqrt(1-pow([2]/(x+[2]),2))+[1]/sqrt(1-pow([3]/([4]+[3]),2))-[5]*[6]", 0, 200);

    //[0]-> L
    //[1]-> L'
    //[2]-> mn.c^2
    // x-> En (our incognita)
    //[3]-> mparticle.c^2
    //[4]-> Ep
    //[5]-> t
	//[6]-> c

    f->SetParameter(0,5);//L
    f->SetParameter(1,0.02184);//L'
	f->SetParameter(2,939.56542052);//Lmn.c2

	switch(particle){
		case 'p':
			f->SetParameter(3,938.2720882);
			break;
		case 'd':
			f->SetParameter(3,1875.612943);
			break;
		case 't':
			f->SetParameter(3,2808.921133);
			break;
		case 'h':
			f->SetParameter(3,2808.391607);
			break;
		case 'a':
			f->SetParameter(3,3727.379407);
			break;
		default:
			f->SetParameter(3,938.2720882);//proton

	}

    f->SetParameter(4,Ep);
    f->SetParameter(5,tof);
    f->SetParameter(6,0.299792458);//c in m/ns


    //f->Draw();
    //gPad->SetGridx();
    //gPad->SetGridy();


     RootFinder *k = new RootFinder();
       k->SetMethod(RootFinder::kBRENT);
     k->Solve(*f,0,10000);

     double c = k->Root();
    //cout << c << endl;

    return c;
}



double ToFp(double Eparticle, char particle, int tel=1)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, given the distance 'L' between target and the ToF detector

	double L; //distance from the target to Si2

	switch(tel){ //which depends on the telescope we are using.
		case 1:
			L=0.2184;
			break;
		case 2:
			L=0.1609;
			break;
		case 3:
			L=0.1609;
			break;
		case 4:
			L=0.1609;
			break;
		case 5:
			L=0.1739;
			break;
		case 6:
			L=0.1739;
			break;
		case 7:
			L=0.1709;
			break;
		case 8:
			L=0.2079;
			break;
		default:
			L=0.2184;
	}

	double c= 0.299792458;
	double m = 938.27208816;

	switch(particle){//data from codata
		case 'p':
			m = 938.27208816;
			break;
		case 'd':
			m = 1875.61294257;
			break;
		case 't':
			m = 2808.92113298;
			break;
		case 'h':
			m =  2808.39160743;
			break;
		case 'a':
			m = 3727.3794066;
			break;
		default:
			m = 938.27208816;//proton
	}

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}

double ToFn(double Eneutron, double L)
{
//returns the ToF for a neutron with energy Eneutron, given the distance 'L' (in meters) between Li target and Medley's target

	double c= 0.299792458;//m/ns
	double m = 939.56542052;

    return L/(c*sqrt(1 - pow((m/(Eneutron+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


double En(double tof, double L = 5.044){
	//return neutron energy in MeV, given its time of fligt.
	//double L = 5.044; //meters
	if(tof==0){ cout<<"ERROR. Tof==0 is not allowed. Check the other steps."<<endl; return -1;}
	double v = L/tof;
	double c = 0.299792458; //m/ns
	double gamma = 1/sqrt(1-pow(v/c,2));
	double mc2 = 939.56542052; //MeV
	return (gamma -1)*mc2; //MeV
}


string lastname(string name){
	int l = name.size();
	l--;
	string namefile;
	while(l>=0){
		if(name[l] == '/'){
			break;

		}else{
			namefile  = name[l] + namefile;
		}
		l--;
	}
	return namefile;
}

string lastnameNoExt(string name){
	// script that gives you the last name of file without extention:
	// example:
	// root [1] string namefile = "home/lucas/dearruda/ganil/filename.root"
	// (std::string) "filename"

	string namefile = lastname(name);

	int t=0;
	int l=namefile.size();
	//cout<<"size = "<<l<<endl;
	while(namefile[l]!='.'){
	//	cout<<Form("l = %d : %c",l,(char)namefile[l])<<endl;
		l--;
		t++;
		if(l==0){
			break;
		}
	}
	//cout<<"final l ="<<l<<endl;
	namefile = namefile.substr(0,namefile.size()-t);

	return namefile;
}

struct calpoint
{
  int ch, sigch;
  float E;
};

vector <calpoint>getAddPts(string filename, char D, int ID){
  char c;
  int tel, n, tempch, tempsigch, p,st;
  float tempf;
  p=-1;//index for the vector, starts at -1 because it is updated as  the new position is created
 st=0;
  vector<calpoint> addpts; //create a vector of calpoint struct

  cout<<"reading additional points..."<<endl;

  ifstream file(filename.c_str()); //read file
  string line; // declare the line string

  if (!file) // give error if the file does not exist
    cout << "Couldn't open file : " << filename << endl ;
    //return -1;
  else{   //if exists...
      bool control=true; //starts with true
      getline(file,line); //read the first line
	  			//cout<<"first line: "<<line<<endl;

	  if(line[0]=='=') control = false;

      while(control == false){ //while it is false
        getline(file,line);//keep reading
				//cout<<"another line: "<<line<<endl;
        if(line[0]=='=') control = true;
      }
      while(!file.eof()){
        //getline(file,line);//keep reading
		file>>c;
		file>>tel;
		file>>n;
		st++;
		if(st==10) break;
        //cout<<"c = "<<c<<" tel= "<<tel<<" n= "<<n<<endl;
        if(c == D && tel == ID){
            for(int j=0;j<n;j++){
				//cout<<" j = "<<j<<" n= "<<n<<endl;
                addpts.push_back(calpoint());// add element fo the vector
                p++;//update index
                file>>addpts[p].ch;
                file>>addpts[p].sigch;
                file>>addpts[p].E;
				cout<<"	...addpts["<<p<<"].ch = "<<addpts[p].ch<<" addpts["<<p<<"].sigch = "<<addpts[p].sigch<<" addpts["<<p<<"].E = "<<addpts[p].E<<endl;
            }
        }else{
			for(int j=0;j<n;j++){
				//cout<<"--> j = "<<j<<" n= "<<n<<endl;
                file>>tempch;
                file>>tempsigch;
                file>>tempf;
				//cout<<"--> TEMP.ch = "<<tempch<<" TEMP.sigch = "<<tempsigch<<" TEMP.E = "<<tempf<<endl;
            }
		}
      }
  }
    return addpts;
}

vector<cal> fetchCal(string path="/home/dearruda/ganil/completed_calibrations-3alpha/second_calibrations/"){
//fetch calibrations inside a directory.

cout<<"\nFetching calibration...\n\n"<<endl;//info to the user

vector <string> names; //string vector for the names it found
names = namesfiles(path.c_str(),".txt");

int N;
N = names.size(); // the size of the vector == number of names
cout<<"...Number of calibration files found = "<<N<<" files."<<endl; // show the information

string tp;
vector<cal> calibration; // create calibration vector
int ctrl, det, p;
float value;

p=-1;//p is the index for the struct vector; starts at -1 because it is updated as soon as a new struct is created
cout<<"\nFiles and identified telescopes:"<<endl;
for (int i=0;i<N;i++) {//i is the file number
  cout << names[i] <<"\t---> "<<names[i][5]<<".....if#"<<i<<endl; // which contains the telescope's type

  ifstream cfile(Form("%s/%s",path.c_str(),names[i].c_str()));//to open the file we want

  if (cfile.is_open()){ //checking whether the file is open
    //cout<<"INSIDE IF #:"<<i<<endl;
    ctrl = -1;
    while(getline(cfile, tp)){ //read data from file object and put it into string.
      if(tp[0] == '=') ctrl++; //ckech in which line we are
                               //cout<<">>CTRL "<<ctrl<<endl;
      if(ctrl == 1){//control that determines when we are out of the header area
        while(!cfile.eof()){
          cfile>>det;
          //cout<<"det#"<<det<<endl;
          calibration.push_back(cal());//for a new detector, create space on structure vector
          p++;//update index
          calibration[p].detec=names[i][5]; //assign detector type (character out of the name of the file) to the p-th struct
          calibration[p].Tel = det; //telescope

			if(calibration[p].detec == 'C'){
					//cfile>>value;
					//cfile>>value;//read peak location in chn

				cfile>>value;
				calibration[p].A = value;
				//cout<<"A="<<value<<endl;

				cfile>>value;
				calibration[p].B = value;
				//cout<<"sA="<<value<<endl;

			}else{// not C detector
				for(int t=1; t<=6;t++){
					cfile>>value;
					//cout<<"value"<<t<<"="<<value<<endl;
				}
				cfile>>value;
				calibration[p].A = value;
				//cout<<"A="<<value<<endl;

				cfile>>value;
				calibration[p].sA = value;
				//cout<<"sA="<<value<<endl;

				calibration[p].B = 0;
				calibration[p].sB = 0;//now it is zero; before was:
											//cfile>>value;
											//calibration[p].B = value;
											//cout<<"B="<<value<<endl;

											//cfile>>value;
											//calibration[p].sB = value;
											//cout<<"sB="<<value<<endl;
			}
        }
      }
    }
    cfile.close();
  } //close the file object.
}

// cout<<"\n--- Summary of the processed information: ---\n"<<endl;

// for(int i=0; i<=p; i++){
//   cout<<Form("calibration[%d]\n.detec=%c\n.Tel=%d\n",i,calibration[i].detec,calibration[i].Tel);
//   cout<<Form(".A=%f\t.sA=%f\n",calibration[i].A,calibration[i].sA);
//   //cout<<Form(".B=%f\t.sB=%f\n:\n",calibration[i].B,calibration[i].sB);
// }

return calibration;

}



// cout<<"\n--- Summary of the processed information: ---\n"<<endl;

// for(int i=0; i<=p; i++){
//   cout<<Form("calibration[%d]\n.detec=%c\n.Tel=%d\n",i,calibration[i].detec,calibration[i].Tel);
//   cout<<Form(".A=%f\t.sA=%f\n",calibration[i].A,calibration[i].sA);
//   //cout<<Form(".B=%f\t.sB=%f\n:\n",calibration[i].B,calibration[i].sB);
// }

// return calibration;

// }

double pol3(double a, double b,double c,double d, double x){
    //polinomial-4 edited
    return a + b*x + c*x*x + d*x*x*x;
}

double pol4(double a, double b,double c,double d,double e, double x){
    //polinomial-4 edited
    return a + b*x + c*pow(x,2) + d*pow(x,3) + e*pow(x,4);
}

double calibLin(double A, double B, double x){
	double cal = A*x+B;
	if(cal<0){
		return 0;
	}else{
		return cal;
	}
}

double TTcSec(TTree *tree, double Er, char* branch1, char* branch2){
//return the interpolated value of cross section for a provided E value.
//You have to provide a TTree with branches "E"(eV) and <branch> (b/sr).
// return -1 if you are out of range.

Er = Er*1e6;

Double_t E, s;
tree->SetBranchAddress(branch1, &E);
tree->SetBranchAddress(branch2, &s);

//cout<<"branch1 = "<<branch1<<endl;
//cout<<"branch2 = "<<branch2<<endl;

Int_t Nevents = tree->GetEntries();

Double_t Emin, Emax;
Emin = tree->GetMinimum(branch1);
Emax = tree->GetMaximum(branch1);

if(Er<Emin || Er>Emax){
	cout<<Form("E_required (%.2f MeV) is out of range (%.2f to  %.2f) MeV", Er*1e-6, Emin*1e-6,Emax*1e-6)<<endl;
	return -1;
}

Int_t i=0;
E = -1; //initial condition
double Ea, Eb, sa, sb, alpha, beta;
while(E<Er){
	tree->GetEntry(i);

	if(i==0){
		Eb = E;
		sb = s;
	}else{
		Ea = Eb;
		sa = sb;
		Eb = E;
		sb = s;
	}

	//cout<<Form("E = %f s= %f\n", E, s);
	if(i == Nevents){
		return s;
	}
	i++;
}

if(Er == E)return s;

//cout<<Form("Eb = %f, Ea = %f, sb= %f, sa = %f\n", Eb, Ea, sb, sa);
//s(Er) = alpha*Er + beta

alpha = (sb-sa)/(Eb-Ea);
beta = sb - alpha*Eb;

return alpha*Er + beta;

}

Int_t IntegralC(TH1D *h, double a, double b){
	//return the number of counts between 2 bins
	return h->Integral(h->GetXaxis()->FindBin(a),h->GetXaxis()->FindBin(b));
}

double IntegralH(TH1D *h, double a, double b){//make the integral on a histogram
	int nbina, nbinb;
	double integral=0;

	nbina= h->GetXaxis()->FindBin(a);//Findbin considers the bin his low edge but not the up edge
	nbinb= h->GetXaxis()->FindBin(b);

	if(nbina == nbinb){//they are in the same bin
		integral  = (b-a)*h->GetBinContent(nbina);
	}else{//separated bins

		integral = (h->GetBinLowEdge(nbina+1) - a)*h->GetBinContent(nbina)+(b-h->GetBinLowEdge(nbinb))*h->GetBinContent(nbinb);

		if(nbinb-nbina>1){
			for(int t=nbina+1; t<nbinb;t++){
				integral += h->GetBinContent(t)*h->GetBinWidth(t);
			}
		}

	}
	return integral;
}

TH1D *RebinHisto(TH1D *histo, int nbins, double xi, double xf, string name = "newhisto"){ //make the benchmark of it!
//receives histo and returns a newhisto(nbins, xi, xf). Counts ouside [xi,xf] will be not considered.

	TH1D * nhisto = new TH1D(name.c_str(), name.c_str(), nbins, xi, xf);

	double le, ue, intg,bw;
	for(int i=1;i<=nbins;i++){//looping over the new histo

	le = nhisto->GetBinLowEdge(i);
	ue = nhisto->GetBinLowEdge(i+1);
	intg= IntegralH(histo,nhisto->GetBinLowEdge(i),nhisto->GetBinLowEdge(i+1));
	bw= nhisto->GetBinWidth(i);
	//cout<<Form("le: %f\t ue: %f\t int: %f\t bw: %f",le, ue, intg, bw)<<endl;
	nhisto->SetBinContent(i,IntegralH(histo,le,ue)/bw);

	}
	return nhisto;
}

double interpol(TTree *tree, double xr, char* branch1, char* branch2){
//return the interpolated value of a MONOTONIC TTree.
//You provide the ttree with the branches 'branch1' (x) and 'branch2' (y(x)).
//it will return y(x) :D

Double_t y, x;
tree->SetBranchAddress(branch1, &x);
tree->SetBranchAddress(branch2, &y);

//cout<<"branch1 = "<<branch1<<endl;
//cout<<"branch2 = "<<branch2<<endl;

Int_t Nevents = tree->GetEntries();

Double_t xmin, xmax;
xmin = tree->GetMinimum(branch1);
xmax = tree->GetMaximum(branch1);

if(xr<xmin || xr>xmax){
	cout<<Form("x_required (%.2f) is out of range (%.2f to  %.2f)", xr,xmin,xmax)<<endl;
	return -1;
}

Int_t i=0;
double ya, yb, xa, xb, alpha, beta;

tree->GetEntry(i++); //read first ev and assign to x and y;

xa = x; ya=y; //record this value as event A

while(i<Nevents){//keep doing this until the end of the  ttree

	tree->GetEntry(i++); // read next event

	xb=x;//attributes to ev2
	yb=y;

	if(xa>xb){ //verify the order of events to compare correctly
		swap(xa,xb);//if the order is wrong, order them!
		swap(ya,yb);
	}
	if(xr<=xb && xr>=xa && y!=0){//compare: if xr belongs to [xa,xb]
		alpha = (yb-ya)/(xb-xa);
		beta = yb - alpha*xb;

		return alpha*xr + beta;
	}
	//the new ev1 (ev1') is the old ev2
	xa = x;
	ya = y;
}
//if ure out the loop, you reached the last event
return y; //return the last y read

}


double xdmin(TF1 *f, double x0, double y0, int step = 1e3) {//produced by chatgpt
  // Define a function to calculate the distance between two points
  auto distance = [&](double x, double y) {
    return std::sqrt(std::pow(x - x0, 2) + std::pow(y - y0, 2));
  };

  double xMin = f->GetXmin();
  double xMax = f->GetXmax();
  double xStep = (xMax - xMin) / step;  // set the step size for x

  double closestDist = std::numeric_limits<double>::max();
  double closestX = 0.0;

  for (double x = xMin; x <= xMax; x += xStep) {
    double y = f->Eval(x);  // evaluate the function at x

    double d = distance(x, y);  // calculate the distance

    if (d < closestDist) {  // check if this point is the closest
      closestDist = d;
      closestX = x;
    }
  }

  return closestX;
}

double ToF(double Eparticle, char particle, double L = 1)
{
//returns the ToF (ns) for the a particle 'particle' with energy Eparticle (MeV), given the distance 'L' (m)

	double c= 0.299792458; // m/ns
	double m = 938.27208816; // proton mass by standard

	switch(particle){//data from codata
		case 'p':
			m = 938.27208816;
			break;
		case 'd':
			m = 1875.61294257;
			break;
		case 't':
			m = 2808.92113298;
			break;
		case 'h':
			m =  2808.39160743;
			break;
		case 'a':
			m = 3727.3794066;
			break;
		default:
			m = 938.27208816;//proton
	}

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}

double timeRun(TTree *tree, bool plots = true,Int_t *counts=0, char* nameplot = Form("%s","avg"), char* branch1 = Form("%s","Micromegas_TOFTS"), char* branch2 = Form("%s","Micromegas")){
//Function created in 04-05-2023
//Last update: 24-05-2023
//	-> added the possibility of title for the plots, by informing 'nameplot'
//	-> added possibility of returning the number of counts above threshold
//	-> added the line to state useful info:
///		 "There are 123784 (99.1 %) events in Micromegas above the 5000 threshold.""
//	previous update: 12-05-2023
//  	-> add option to suppress plots
//typically processes about 15k events/second.
//Returns the proper run time in seconds, based on the MM's timestamp.
//It receives the TTree with 'Micromegas' and <branch1> branches
//and plot their events accordingly to timestamp. Then it evaluates the number of events for second
//and plot it, removing from the total time span the amount of time when the avg counts are smaller than
//a givern threshold (thcounts)

bool first = true;

cout<<"\n>>>>> Computing run time..."<<endl;
if(tree->GetListOfBranches()->FindObject(branch1)){
	 cout<<Form("branch %s found. :)", branch1)<<endl;
}else{
	cout<<Form("branch %s NOT found.", branch1)<<endl; return 0;
}

if(tree->GetListOfBranches()->FindObject(branch2)){
	 cout<<Form("branch %s found. :)", branch2)<<endl;
}else{
	cout<<Form("branch %s NOT found.", branch2)<<endl; return 0;
}

ULong64_t ts;
Int_t evts = tree->GetEntries();
ULong64_t primeiro;
vector<ULong64_t> tme;
vector<float> signal;

float sg;

UShort_t sig;
///cout<<"comparison between "<<branch2<<" and Micro: "<<strcmp(branch2, "Micromegas")<<endl;



tree->SetBranchAddress(branch1,&ts);

if(strcmp((char*)branch2, "Micromegas")){
	tree->SetBranchAddress(branch2,&sig);
}else{
	tree->SetBranchAddress(branch2,&sg);
}

//tree->Draw(Form("%s>>h",branch1)); it is working

TGraph *MM = new TGraph();// create a plot to analyse it
Int_t ev_above_threshold=0;
int threshold = 5000;
int p=1;//define p as 1, because the p=0 will be defined within the if
for(Int_t i=0; i<evts; i++){

	tree->GetEntry(i);

	if(ts>0 && !first){ // if it is NOT the first event and ts>0...
		tme.push_back(10*(ts - primeiro)); //evaluate time

		if(strcmp((char*)branch2, "Micromegas")){
			signal.push_back((float)sig);		  //evaluate signal
			if(plots) MM->SetPoint(p,tme.at(p),(float)sig);	  //setpoint
			if(sig>threshold) ev_above_threshold++;
		}else{
			signal.push_back(sg);		  //evaluate signal
			if(plots) MM->SetPoint(p,tme.at(p),sg);	  //setpoint
			if(sg>threshold) ev_above_threshold++;
		}

		//cout<<"time = "<<time[p]<<" signal = "<<sg<<endl;

		p++; //update counter

	}else if(ts>0){ //otherwise if it is the first event
		primeiro = ts;	 // primeiro evento
		first = false; //is not the first anymore (so will enter in the previous loop next time)
		tme.push_back(0);	// time for this is zero

		if(strcmp((char*)branch2, "Micromegas")){
			signal.push_back((float)sig); // signal is its own signal
			if(plots) MM->SetPoint(0,0.0,(float)sig);	 // setpoint
		}else{
			signal.push_back(sg); // signal is its own signal
			if(plots) MM->SetPoint(0,0.0,sg);	 // setpoint
		}
		//cout<<"time = "<<time[0]<<" ";
	}

	//cout<<i<<" timestamp = "<<ts<<endl;
}
cout<<"There are "<<p<<" MM events."<<endl;
p--;// set the counter to the last event.

cout<<Form("There are %d (%03.1f %%) events in %s above the %d threshold.",ev_above_threshold,100.0*ev_above_threshold/(p+1),branch2,threshold)<<endl;

*counts = ev_above_threshold;

//define the window to average:
//Int_t window = 1e9; //1 second expressed in ns
Int_t window = 1e9; //1 second expressed in ns

//express the last one in windows unit (seconds)
int final_time = ceil(1.0*tme.back()/window); // get the closest high integer

cout << "\ntme.back() = " <<tme.back()<<" ns"<<endl;
cout << "final_time = "<<final_time <<" s." << endl;


p=0; // p is the position at the time vector


TH1I *avg = new TH1I("avg", "avg", final_time,0,final_time);
//TH1I *avg = new TH1I("avg", "avg", ceil(final_time/window),0,ceil(final_time/window));

for(Int_t g=0; g<tme.size();g++){
	avg->Fill(round(tme.at(g)*1e-9));
}



//now I will count the number of seconds in wich the counts  are below a threshold

int thcounts = 30;
int rem=0; //threshold counts/sec and 'rem' which is the number of bins to remove
for(int b = 1; b<=avg->GetNbinsX();b++){
	if(avg->GetBinContent(b)<thcounts) rem++;
}
cout << "time to remove = " << rem<<'.'<<endl;

cout << "\n.\n.\n. RUNTIME = " << final_time - rem <<" s.";
if(plots){
	new TCanvas(Form("%s_1-%s",branch1,nameplot),Form("%s_1-%s",branch1,nameplot));
	avg->SetTitle(Form("avg th=%d - %s", threshold,nameplot));
	avg->GetXaxis()->SetTitle("time (s)");
	avg->GetYaxis()->SetTitle("avg. counts/s");
	avg->Draw();
	gPad->SetGridx();
	gPad->SetGridy();
	new TCanvas(Form("%s_2-%s",branch1,nameplot),Form("%s_2-%s",branch1,nameplot));
	MM->SetTitle(Form("Events %s",nameplot));
	MM->GetXaxis()->SetTitle("Event#");
	MM->GetYaxis()->SetTitle("channel");
	MM->Draw("AP");
}
return final_time - rem;
}

Int_t MM(TTree * tree, bool plot = true, char* branch1 = Form("%s","Micromegas"), double threshold = 5000){
//Function created in 12-05-2023
//Count how many events of a branch are above a certain threshold.
//It allows you to delete the plot branch1:Entry$ (by setting plot = false).

	TCanvas *ctemp = new TCanvas("ctemp","ctemp");

	Int_t counts;
	counts = tree->Draw(Form("%s:Entry$", branch1), Form("%s>%f", branch1,  threshold));
	cout<<"counts = "<<counts<<endl;
	if(!plot) delete ctemp;

	return counts;
}

double timeRun2(TTree *tree, bool plots = true,Int_t *counts=0, char* nameplot = Form("%s","avg"), char* branch1 = Form("%s","Micromegas_TOFTS"), char* branch2 = Form("%s","Micromegas")){
//Function created in 30-06-2023
//Returns the proper run time in seconds, based on the MM's timestamp.
//It receives the TTree with <branch1>='Micromegas' and <branch2>=Micromegas branches



cout<<"\n>>>>> Computing run time..."<<endl;
if(tree->GetListOfBranches()->FindObject(branch1)){
	 cout<<Form("branch %s found. :)", branch1)<<endl;
}else{
	cout<<Form("branch %s NOT found.", branch1)<<endl; return 0;
}

if(tree->GetListOfBranches()->FindObject(branch2)){
	 cout<<Form("branch %s found. :)", branch2)<<endl;
}else{
	cout<<Form("branch %s NOT found.", branch2)<<endl; return 0;
}

ULong64_t ts;
Int_t evts = tree->GetEntries();
double primeiro;
vector<double> tme;
vector<float> signal;

float sg;

UShort_t sig;

tree->SetBranchAddress(branch1,&ts);

if(strcmp((char*)branch2, "Micromegas")){
	tree->SetBranchAddress(branch2,&sig);
}else{
	tree->SetBranchAddress(branch2,&sg);
}

//tree->Draw(Form("%s>>h",branch1)); it is working

TGraph *MM = new TGraph();// create a plot to analyse it
Int_t ev_above_threshold=0;
int threshold = 5000;
int p=1;//define p as 1, because the p=0 will be defined within the if

bool first = true;

for(Int_t i=0; i<evts; i++){

	tree->GetEntry(i);

	if(ts>0 && !first){ // if it is NOT the first event and ts>0...
		tme.push_back(10*(ts - primeiro)*1e-9); //evaluate time in seconds

		if(strcmp((char*)branch2, "Micromegas")){
			signal.push_back((float)sig);		  //evaluate signal
			if(plots) MM->SetPoint(p,tme.at(p),(float)sig);	  //setpoint
			if(sig>threshold) ev_above_threshold++;
		}else{
			signal.push_back(sg);		  //evaluate signal
			if(plots) MM->SetPoint(p,tme.at(p),sg);	  //setpoint
			if(sg>threshold) ev_above_threshold++;
		}

		//cout<<"time = "<<time[p]<<" signal = "<<sg<<endl;

		p++; //update counter

	}else if(ts>0){ //otherwise if it is the first event
		primeiro = ts;	 // primeiro evento
		first = false; //is not the first anymore (so will enter in the previous loop next time)
		tme.push_back(0);	// time for this is zero

		if(strcmp((char*)branch2, "Micromegas")){
			signal.push_back((float)sig); // signal is its own signal
			if(plots) MM->SetPoint(0,0.0,(float)sig);	 // setpoint
		}else{
			signal.push_back(sg); // signal is its own signal
			if(plots) MM->SetPoint(0,0.0,sg);	 // setpoint
		}
		//cout<<"time = "<<time[0]<<" ";
	}

	//cout<<i<<" timestamp = "<<ts<<endl;
}
cout<<"There are "<<p<<" MM events."<<endl;
p--;// set the counter to the last event.

cout<<Form("There are %d (%03.1f %%) events in %s above the %d threshold.",ev_above_threshold,100.0*ev_above_threshold/(p+1),branch2,threshold)<<endl;

*counts = ev_above_threshold;

//define the window to average:
//Int_t window = 1e9; //1 second expressed in ns
double window = 15; //in seconds

//express the last one in windows unit (seconds)
int final_time = ceil(1.0*tme.back()/window); // get the closest high integer

cout << "\ntme.back() = " <<tme.back()<<" s"<<endl;
cout << "final_time = "<<final_time <<" [windows]." << endl;


p=0; // p is the position at the time vector


TH1I *avg = new TH1I("avg", "avg", final_time,0,final_time);
//TH1I *avg = new TH1I("avg", "avg", ceil(final_time/window),0,ceil(final_time/window));

for(Int_t g=0; g<tme.size();g++){
	avg->Fill(round(tme.at(g)/window));
}



//now I will count the number of seconds in wich the counts  are below a threshold

int thcounts = 30;
int rem=0; //threshold counts/sec and 'rem' which is the number of bins to remove
for(int b = 1; b<=avg->GetNbinsX();b++){
	if(avg->GetBinContent(b)<thcounts) rem++;
}
cout << "seconds to remove = " << rem*window;

cout << "\n.\n.\n. RUNTIME = " << (final_time - rem)*window <<" s.";
if(plots){
	new TCanvas(Form("%s_1-%s",branch1,nameplot),Form("%s_1-%s",branch1,nameplot));
	avg->SetTitle(Form("avg th=%d - %s", threshold,nameplot));
	avg->GetXaxis()->SetTitle("time (s)");
	avg->GetYaxis()->SetTitle("avg. counts/s");
	avg->Draw();
	new TCanvas(Form("%s_2-%s",branch1,nameplot),Form("%s_2-%s",branch1,nameplot));
	MM->SetTitle(Form("Events %s",nameplot));
	MM->GetXaxis()->SetTitle("time (s)");
	MM->GetYaxis()->SetTitle("channel");
	MM->Draw("AP");
}
return (final_time - rem)*window ;
}


double tgamma(int telescope){

//Compute the time the gamma flash needs to travel from the target to the telescope


//Distance between rotating converter and medley target (L):
	double L = 5.044;

	//                 Lf              ______:::::__;-
	//                   ______:::::'''''     __;--
	//   ______:::::'''''            X __;--         L1
	//  -------------------------------'   y
	//                L
	//
	// Lf^2 = L^2 + L1^2 - L1.L.cos(X), but X = pi/2 - y
	// Lf = sqrt[L^2 + L1^2 - L1.L.sin(y)], but X = pi/2 - y




	//Distance between medley target and telescope:
	double L1, angle, Lf; //distance from the target to Si2


	switch(telescope){ //which depends on the telescope we are using.
		case 1:
			L1=0.2184;
			angle = Pi()/9;
			break;
		case 2:
			L1=0.1609;
			angle = 2*Pi()/9;
			break;
		case 3:
			L1=0.1609;
			angle = 3*Pi()/9;
			break;
		case 4:
			L1=0.1609;
			angle = 4*Pi()/9;
			break;
		case 5:
			L1=0.1739;
			angle = 5*Pi()/9;
			break;
		case 6:
			L1=0.1739;
			angle = 6*Pi()/9;
			break;
		case 7:
			L1=0.1709;
			angle = 7*Pi()/9;
			break;
		case 8:
			L1=0.2079;
			angle = 8*Pi()/9;
			break;
		default:
			angle = 1*Pi()/9;
			L1=0.2184;
	}
	Lf = sqrt(L1*L1 + L*L + L1*L*cos(angle));


return Lf/0.299792458;
}


TH1D * hsub( TH1D* histoA, TH1D* histoB, double factor, string name="Subtracted_histo"){

	// Check for nullptr
    if (!histoA || !histoB) {
        cerr << "Error: One or both input histograms are nullptr." << endl;
        return nullptr;
    }

    // Check for negative factor
    if (factor < 0) {
        cerr << "Warning: Negative factor provided. Using it as an addition, not subtraction." << endl;
    }

	TH1D * hsub;
	hsub = (TH1D*)histoA->Clone();
	hsub->SetName(name.c_str());
	hsub->SetNameTitle(name.c_str(),name.c_str());

	hsub->Add(histoB,-factor);
	//we need to correct negative bins
	for(int i=1; i<=hsub->GetNbinsX();i++){
		if(hsub->GetBinContent(i)<0) hsub->SetBinContent(i,0);
	}

	return hsub;
}


vector<double> ReadDistances() {
    // Create a vector to store the distances
    std::vector<double> d;

    // Open the input file
    std::ifstream inputFile("/home/dearruda/ganil/runs_logbook22/defs/distancesTelescopes.txt");

    // Check if the file is open successfully
    if (!inputFile.is_open()) {
        std::cerr << "Error: Unable to open the file defs/distanceTelescopes.txt" << std::endl;
        return d;  // Return an empty vector in case of an error
    }

    // Read distances from the file and store them in the vector, skipping lines starting with '#'
    double distance;
    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') {
            // Skip empty lines and lines starting with '#'
            continue;
        }
        // Try to convert the line to a double and store it in the vector
        std::istringstream iss(line);
        if (iss >> distance) {
            d.push_back(distance);
            //std::cout << "Read distance: " << distance << std::endl;
        }
    }

    // Close the input file
    inputFile.close();

    return d; // Return the vector containing distances
}

// void normalize_3vec(double x[3]) {
//     double magnitude = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

//     if (magnitude != 0) {
//         x[0] /= magnitude;
//         x[1] /= magnitude;
//         x[2] /= magnitude;
//     } else {
//         // Handle the case where the magnitude is zero to avoid division by zero.
//         x[0] = 0.0;
//         x[1] = 0.0;
//         x[2] = 0.0;
//     }
// }


TVector3 coordTransform(TVector3 D,TVector3 sb){
//Traform coordinates from B to A; where B is a coordinate system with origin at vector D (represented in A), whose z-axis
//is parallel to it.
//inputs:
// D = {xd,yd,zd} (represented in A)
// sb = {sx,sy,sz} (represented in B)
//
//you have to #include "TVector3.h"


//defining a versors:
TVector3 i_a(1.,0.,0.);
TVector3 j_a(0.,1.,0.);
TVector3 k_a(0.,0.,1.);


//angle between Z axis in A and D
Double_t angle = k_a.Angle(D);

//this defines de axis of rotation to make Z reach D direction
TVector3 axisR = k_a.Cross(D.Unit());

sb.Rotate(angle,axisR);

return sb + D;

}


//Created 22/03/2024
int IntegralCut(TCutG *cut, TH2D * histo){

	double x, y, counts;
	int sum=0;
		for(int i=1; i<=histo->GetNbinsY();i++){
			for(int j=1; j<=histo->GetNbinsX();j++){
				x = histo->GetXaxis()->GetBinCenter(j);
				y = histo->GetYaxis()->GetBinCenter(i);
				counts = histo->GetBinContent(j,i);
				if(cut->IsInside(x,y)) sum += counts;
			}
		}

return sum;

}
