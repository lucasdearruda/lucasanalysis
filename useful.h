/*
//Set of functions developed by Lucas de Arruda//
version2025.07.22.001

2025-07-22:: Updating GetGflash function to include sigma and threshold parameters
2025-05-27:: creation of function NumIntegrationTF1
2025-01-26.2::  creation of function GiveMetheCharge
2025-01-26.1::  beamlivetime update. The method tx->GetMaximum() was causing a segmentation fault 
2025-01-24::  giveMeTheAngle2 added
2025-01-24:: definParticle added. Taken originally from kalscripts/eloss/includes/functions.hh
2024-10-18::  new version of giveMeTheAngle added
*/

double NumIntegrationTF1(TF1 *f, double x1, double x2, int npts = 1000) {

	double integral = 0.0;

	double xa,xb;

	//calculate step size
	double step = (x2 - x1) / (npts - 1);

	for (int i = 0; i < npts-1; ++i) {
		xa = x1 + i * step;
		xb = xa + step;
		integral += (f->Eval(xa) + f->Eval(xb)) * step / 2.0; // Trapezoidal rule
	}

	return integral;
}


Float_t giveMeTheCharge(
	Int_t runa = 31,
	Int_t runb = 32,
	string namefile = "/mnt/medley/LucasAnalysis/2023/runlist.csv", 
	string fileformat = "RunN/I:Or/C:Target/C:Time_s/I:Time_h/F:TimeEval_s/I:TimeEval_h/F:ChargeIntegrator/F:ChargeFaraday/F"){

cout<<"---------------------------------------------------"<<endl;
TTree *InfoTree = new TTree("InfoTree", "InfoTree");
InfoTree->ReadFile(namefile.c_str(), fileformat.c_str());

Int_t RunN;
Float_t ChargeFaraday;
Float_t TotalChargeFaraday = 0;

InfoTree->SetBranchAddress("RunN", &RunN);
InfoTree->SetBranchAddress("ChargeFaraday", &ChargeFaraday);

for(int i=0;i<InfoTree->GetEntries();i++)
{
    InfoTree->GetEntry(i);
    //cout<< i<<", nrun = "<<RunN<< endl;
    if (RunN<= runb && RunN>=runa)
    {
        //InfoTree->GetEntry(i);
        cout << "RunN: " << RunN << " ChargeFaraday: " << ChargeFaraday << endl;
        TotalChargeFaraday += ChargeFaraday;
    }
}
cout<<"---------------------------------------------------"<<endl;
cout<< "Total Charge Faraday: " << TotalChargeFaraday <<" µC"<< endl;


return TotalChargeFaraday;

}


void defineParticle(char particle, Int_t* Zfis, Int_t* Afis){


    switch(particle){
        case 'p':
            *Zfis =1;
            *Afis =1;
            cout<<"|   particle    |   PROTON"<<endl;
            break;
        case 'd':
            *Zfis =1;
            *Afis =2;
            cout<<"|   particle    |   DEUTERON"<<endl;
            break;	
        case 't':
            *Zfis =1;
            *Afis =3;
            cout<<"|   particle    |   TRITON"<<endl;
            break;
        case 'h':
            *Zfis =2;
            *Afis =3;
            cout<<"|   particle    |   He-3"<<endl;
            break;
        case 'a':
            *Zfis =2;
            *Afis =4;
            cout<<"|   particle    |   ALPHA"<<endl;
            break;
        default: 
            *Zfis =1;
            *Afis =1;	
            cout<<"PROTON!!"<<endl;
    }


}

std::vector<Double_t> giveMeTheAngle(
    Long64_t N = 1e3,
    Double_t D_mm = 149.7, 
	TVector3 rotation_target = TVector3(0., TMath::Pi() / 4, 0.), //x y z
	TVector3 rotation_detector =  TVector3(0., TMath::Pi() / 9, 0.), //x y z
    Double_t target_radius_mm = 25.0 / 2, 
	Double_t target_thickness_um  = 25.,
    Double_t siA_radius_mm = sqrt(450.0 / TMath::Pi()), 
	Double_t siA_thickness_um = 0.,
	Int_t wn = 2 //which_normal: gives the angle with respect to what.
				//0 : with respect to Z axis
				//1: with respect to TARGET's normal (leaving target's face (beam's direction positive))
				//2 : with respect to Detectors's normal (Si's front face towards it's backface)
) {
    TRandom2 *rn = new TRandom2();
    std::time_t seed = std::time(nullptr);
    rn->SetSeed(seed);

    std::vector<Double_t> angles; 

    Double_t xt, yt, xd, yd;

    TVector3 Point_detector, Point_target, distance;
    TVector3 displacement(0.,0.,D_mm);
    TVector3 normal(0.,0.,1.0);
	TVector3 normalSi = normal;

	normalSi.RotateX(rotation_detector.X()); 
	normalSi.RotateY(rotation_detector.Y()); 
	normalSi.RotateZ(rotation_detector.Z()); 

	normal.RotateX(rotation_target.X()); 
	normal.RotateY(rotation_target.Y()); 
	normal.RotateZ(rotation_target.Z()); 

	/*cout<<" . . . . . . . . . normal: "<<endl;
	normal.Print();
	cout<<" . . . . . . . . . normal Si: "<<endl;
	normalSi.Print();
	cout<<"wn ="<<  wn <<endl;*/
    for(Int_t i=0;i<N;i++){
        Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
        Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_target.X(),2) + pow(Point_target.Y(),2) >= pow(target_radius_mm,2)){
            Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
            Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
        Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_detector.X(),2) + pow(Point_detector.Y(),2) >= pow(siA_radius_mm,2)){
            Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
            Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_target.SetZ(target_thickness_um*(rn->Rndm()-0.5)/1000);
        Point_detector.SetZ(siA_thickness_um*(rn->Rndm()-0.5)/1000);

		//Apply proper rotations in the target and normal


        Point_target.RotateX(rotation_target.X()); 
		Point_target.RotateY(rotation_target.Y()); 
		Point_target.RotateZ(rotation_target.Z()); 
        
		Point_detector += displacement;
		Point_detector.RotateX(rotation_detector.X());//rotate the coordinate to match the detector position
		Point_detector.RotateY(rotation_detector.Y());//rotate the coordinate to match the detector position
    	Point_detector.RotateZ(rotation_detector.Z());//rotate the coordinate to match the detector position

        distance = Point_detector - Point_target;
        //here we have coordinates in target (xt,yt) and detector (xd,xy)! 
        //Now we will calculate the angle
		if(wn == 0){
			angles.push_back(distance.Angle(normal));
		}else if(wn == 1){
			angles.push_back(distance.Angle(displacement));
		}else if(wn == 2){
			angles.push_back(distance.Angle(normalSi));
		}
			

    }

    return angles;

}



std::vector<Double_t> giveMeTheAngle2(
    Long64_t N = 1e3,
    Double_t angle_telescope = 20.0, 
    Double_t D_mm = 149.7, 
    Double_t target_radius_mm = 25.0 / 2, 
    Double_t siA_radius_mm = sqrt(450.0 / TMath::Pi()), 
    TVector3 target_n = TVector3(cos(TMath::Pi() / 4), 0., cos(TMath::Pi() / 4))
) {
    TRandom2 *rn = new TRandom2();
    std::time_t seed = std::time(nullptr);
    rn->SetSeed(seed);

    std::vector<Double_t> angles; 

    Double_t xt, yt, xd, yd;

    TVector3 Point_detector, Point_target, distance;
    TVector3 displacement(0.,0.,D_mm);

    if(angle_telescope<=90){
        displacement.RotateY(-(angle_telescope/180.)*TMath::Pi());
    }else{
        displacement.RotateY((angle_telescope/180.)*TMath::Pi());
    }
        

    for(Int_t i=0;i<N;i++){
        Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
        Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_target.X(),2) + pow(Point_target.Y(),2) >= pow(target_radius_mm,2)){
            Point_target.SetX(target_radius_mm*2*(rn->Rndm()-0.5));
            Point_target.SetY(target_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
        Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        while(pow(Point_detector.X(),2) + pow(Point_detector.Y(),2) >= pow(siA_radius_mm,2)){
            Point_detector.SetX(siA_radius_mm*2*(rn->Rndm()-0.5));
            Point_detector.SetY(siA_radius_mm*2*(rn->Rndm()-0.5));
        }

        Point_target.SetZ(0);
        Point_detector.SetZ(0);


        Point_target.RotateY(TMath::Pi()/4); //rotate point in target because the target has 45 deg inclination! 
        Point_detector.RotateY(-TMath::Pi()/9);//rotate the coordinate to match the detector position
    
        Point_detector += displacement;

        distance = Point_detector - Point_target;
        //here we have coordinates in target (xt,yt) and detector (xd,xy)! 
        //Now we will calculate the angle

        angles.push_back(distance.Angle(displacement));

    }

    return angles;

}

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

bool compareFunction (string a, string b) {return a<b;}
//compare any way you like, here I am using the default string comparison

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



double NintGauss(double A=1, double sigma=1, double Nsigma=9, int npts = 1e5){ //numerical integration of gaussian
// it returns the numerical integration of a gaussian in (Nsigma)*sigma range (centered in mu).
// A - height of gaussian
// sigma - sigma
// Nsigma - number of sigmas  to integrate (it will go from -Nsigmas/s to Nsigmas/2)
// npts - number of points

double xi = -(Nsigma/2)*sigma, xf =  +(Nsigma/2)*sigma;
double step = (xf-xi)/npts;


double sum = 0;
double ya, yb, xa, xb;
for(int i=0; i<npts;i++){
	xa = xi + step*i;
	xb = xa + step;
	ya = A *  TMath::Exp( -0.5 * TMath::Power( xa/sigma ,2) );
	yb = A *  TMath::Exp( -0.5 * TMath::Power( xb/sigma ,2) );

	sum += (ya+yb)*step/2;
}

return sum;
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


double En(double tof, double L = 4.6472){//L in meters, 2023 reference
	//return neutron energy in MeV, given its time of fligt.

	//if(tof<=0){ cout<<"ERROR. Tof<=0 is not allowed. Check the other steps. I will return -1 for that."<<endl; return -1;}
	if(tof<=0){ return 0;}
	double v = L/tof;
	double c = 0.299792458; //m/ns
	double gamma = 1/sqrt(1-(v/c)*(v/c));
	if(gamma!=gamma) return 0; //case where gamma = nan;
	double mc2 = 939.56542052; //MeV
	return (gamma -1)*mc2; //MeV
}




double ToFparticle(double Eparticle, char particle, double L = 4.6472)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'


	double c= 0.299792458;
	double m;

	switch(particle){//data from codata
        	case 'n':
        		m = 939.56542052;
        		break;
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
			m = 939.56542052;//neutron
	}

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}



double ToFneutron(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'

	if(Eparticle<=0){ cout<<Form("ToFproton ERROR. Eparticle(%.3f MeV)<=0 is not allowed. Check the other steps.",Eparticle)<<endl; return 0;}
	double c= 0.299792458;
	double m = 939.56542052;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


double ToFproton(double Eparticle, double L = 0.1497)
{//updated 25/11/2024
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'

	if(Eparticle<=0){ 
	//	cout<<Form("ToFproton ERROR. Eparticle(%.3f MeV)<=0 is not allowed. Check the other steps.",Eparticle)<<endl; return 0;
		return 0;
	}
	double c= 0.299792458;
	double m = 938.27208816;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}
double ToFdeuteron(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'

	if(Eparticle<=0){ cout<<"ERROR. Eparticle<=0 is not allowed. Check the other steps."<<endl; return 0;}
	double c= 0.299792458;
	double m = 1875.61294257;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}

double ToFtriton(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'


	if(Eparticle<=0){ cout<<"ERROR. Eparticle<=0 is not allowed. Check the other steps."<<endl; return 0;}

	double c= 0.299792458;
	double m = 2808.92113298;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}
double ToFhelium(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'


	if(Eparticle<=0){ cout<<"ERROR. Eparticle<=0 is not allowed. Check the other steps."<<endl; return 0;}

	double c= 0.299792458;
	double m = 2808.39160743;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
}


double ToFalpha(double Eparticle, double L = 0.1497)
{
//returns the ToF for the a particle 'particle' with energy Eparticle, that crosses distance 'L'


	if(Eparticle<=0){ cout<<"ERROR. Eparticle<=0 is not allowed. Check the other steps."<<endl; return 0;}

	double c= 0.299792458;
	double m = 3727.3794066;

    return L/(c*sqrt(1 - pow((m/(Eparticle+m)),2))); // comes from E = (γ-1)mc^2 solved for t given that v = L/t.
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


double timeRun2(TTree *tree, bool plots = true, char* branch1 = Form("%s","Medley_1_SI_DE2TS"), char* branch2 = Form("%s","Medley_1_dE2"),Int_t *counts=0, char* nameplot = Form("%s","avg")){
//Function created in 30-06-2023, adapted to net version in 07/10/2023
//Returns the proper run time in seconds, based on the timestamp.
//It receives the TTree with <branch1>='Medley_1_SI_DE2TS' and <branch2>='Medley_1_dE2' branches

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
Float_t sig;

tree->SetBranchAddress(branch1,&ts);

if(strcmp((char*)branch2, "Medley_1_SI_DE2")){
	tree->SetBranchAddress(branch2,&sig);
}else{
	tree->SetBranchAddress(branch2,&sg);
}

//tree->Draw(Form("%s>>h",branch1)); it is working

TGraph *MM = new TGraph();// create a plot to analyse it
Int_t ev_above_threshold=0;
float threshold = 0.5;
int p=1;//define p as 1, because the p=0 will be defined within the if

bool first = true;

for(Int_t i=0; i<evts; i++){

	tree->GetEntry(i);

	if(ts>0 && !first){ // if it is NOT the first event and ts>0...
		tme.push_back(10*(ts - primeiro)*1e-9); //evaluate time in seconds

		if(strcmp((char*)branch2, "Medley_1_dE2")){
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

		if(strcmp((char*)branch2, "Medley_1_dE2")){
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
cout<<"There are "<<p<<" events."<<endl;
p--;// set the counter to the last event.

cout<<Form("There are %d (%03.1f %%) events in %s above the %f threshold.",ev_above_threshold,100.0*ev_above_threshold/(p+1),branch2,threshold)<<endl;

//*counts = ev_above_threshold;

//define the window to average:
//Int_t window = 1e9; //1 second expressed in ns
double window = 60; //in seconds

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
	avg->SetTitle(Form("avg th=%f - %s", threshold,nameplot));
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

//function to get the first non zero entry in given branch 
Long64_t GetFirstNonZero(TTree *G, const char* branchName = "Medley_1_SI_DE2TS"){
    ULong64_t value, nE;
    value = 0;
    G->SetBranchAddress(branchName,&value);
    nE=G->GetEntries();
    for(Int_t i=0; i<nE;i++){
    	G->GetEntry(i);
	if(value){
	//cout<<"evt "<<i<<endl;
	break;

	}

	}
    return value;
}


Long64_t GetLasttNonZero(TTree *G, const char* branchName = "Medley_1_SI_DE2TS"){
    ULong64_t value, nE;
    value = 0;
    G->SetBranchAddress(branchName,&value);
    nE=G->GetEntries();
    for(Int_t i=nE-1; i>=0;i--){
    	G->GetEntry(i);
	if(value){
	//cout<<"evt "<<i<<endl;
	break;

	}

	}
    return value;
}


//function to get the beam live time 
TH1D *beamlivetime(TTree * tx, string branchname="Medley_1_SI_DE2TS",float *time_sec=NULL, float *avc = NULL, float *totaltime = NULL,float threshold = 20,string histoname="hc", bool drawit =false){

//Float_t firstvalue = GetFirstNonZero(tx)/1.e8;
Long64_t firstvalue = GetFirstNonZero(tx,branchname.c_str())/1.e8;
Long64_t lastvalue = GetLasttNonZero(tx,branchname.c_str())/1.e8;

Int_t nbinstime = round(lastvalue - firstvalue) +1;
if(totaltime) *totaltime = lastvalue - firstvalue;

cout << ">>>> 'beamlivetime' function report: (considering TS branches)-------"<<endl;
cout<<".\n. first value found: "<<firstvalue<<" seconds."<<endl;
cout<<". last value found: "<<lastvalue<<" seconds."<<endl;
cout<<". bins TS histo: "<<nbinstime<<" seconds.\n."<<endl;
cout << ">>>>-----------------------------------------------------------------"<<endl;


  TH1D *hc=new TH1D(histoname.c_str(),Form("counts over time: %s",histoname.c_str()),nbinstime,0,nbinstime);

    //cout<<"Maximum TT = "<<tx->GetMaximum
if(drawit){
    hc->GetXaxis()->SetTitle("time (s)");
    tx->Draw(Form("%s/1.e8 - %lld>>%s",branchname.c_str(),firstvalue,histoname.c_str()),"");	// Divide by 1.e8 to get it in s, because time stamp happens every 10 ns;
    TF1 *f1 = new TF1("f1","pol0",0,nbinstime);
    f1->SetParameter(0,threshold);
    f1->SetLineColor(kRed);
    f1->SetLineWidth(2);
    f1->Draw("same");
}else{
    tx->Draw(Form("%s/1.e8 - %lld>>%s",branchname.c_str(),firstvalue,histoname.c_str()),"", "goff");	// Divide by 1.e8 to get it in s, because time stamp happens every 10 ns;
}
  
   Int_t nbins=hc->GetNbinsX();
  //cout << "nbins= " << nbins << " limits: " << h1->GetXaxis()->GetXmin() << " " << h1->GetXaxis()->GetXmax() << endl;
  Int_t realbeamtime=0;

  Float_t beamthreshold=threshold; // For C experiment: beamthreshold=40.; for Cr experiment: beamthreshold=10.;
  Float_t cc=0.;
  Float_t CR = 0.;
  for(Int_t i=0;i<nbins;i++){
   cc=hc->GetBinContent(i+1);

   if(cc>=beamthreshold)
   {
 //  cout << "bin " << i << " cc= " << cc << endl;
     realbeamtime++;
     CR += cc;
   }

  }

  CR = 1.*CR/realbeamtime; //in secs
  cout << "\n\n>> Livetime summary -------------------------------------------------"<<endl;

  cout << ". Real beam time: " << realbeamtime <<  "(s)  = " << realbeamtime/3600. << " (h). Avg CR = " << CR<< "/s" << endl;
  if(time_sec)*time_sec = realbeamtime;
  if(avc)*avc = CR;
    cout << ".\n---------------------------------------------------------------------"<<endl;

return hc;



}




//function for getting gamma flash
//Float_t GetGflash(TTree *tx, string branchname="Medley_1_dE2_ToF", float guess = 500, bool closecanvas = false){
Float_t GetGflash(TTree *tx, string branchname="Medley_1_dE2_ToF", float guess = 500, bool closecanvas = false, Float_t sigma= 4,Float_t threshold= 0.1){


	TCanvas *Ctof = new TCanvas("ctof",Form("Time of flight"),50,50,600,600);

    TH1D *htof;
    htof = new TH1D("htof","htof",500,100,600);

    tx->Draw(Form("%s>>htof",branchname.c_str()));

    double tof_peak[2];
    TSpectrum *spec = new TSpectrum();
    Int_t npeaks = spec->Search(htof, sigma,"",threshold);//originally 0.3
    Double_t *xpeaks = spec->GetPositionX();
    cout<<"ToF plot: "<<npeaks<<" peak found.\n";
    for(int t =0;t<npeaks;t++){
        cout<<"position peak"<<t<<": "<<xpeaks[t]<<" ns.\n";
    }

    if (npeaks >= 2) {
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = xpeaks[npeaks -1] - 10;
        double fitRangeMax = xpeaks[npeaks -1] + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[npeaks -1], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }else{ //if did not find the peak
        TF1 *gaussian = new TF1("gaussian_%02d", "gaus");
        double fitRangeMin = guess - 10;
        double fitRangeMax = guess  + 10;
        gaussian->SetParameters(htof->GetMaximum(), xpeaks[0], 5.0);
        htof->Fit(gaussian, "Q", "", fitRangeMin, fitRangeMax);
        tof_peak[0] = gaussian->GetParameter(1);
        tof_peak[1] = gaussian->GetParError(1);
    }


	if(closecanvas) Ctof->Close();
	
	return tof_peak[0];
}


TGraphErrors* invertGraph(TGraphErrors *g){

    TGraphErrors *gr;
    gr = new TGraphErrors();
    int np= g->GetN();
    float Xlp = g->GetPointX(np-1);
    int l=0;
    for(int p=np-1;p>=0;p--){
        gr->SetPoint(l++,2*Xlp-g->GetPointX(p),g->GetPointY(p));
    }

    return gr;
}


//summing Tgraphs...
TGraphErrors* sumTgraphs(TGraphErrors *a,TGraphErrors *b){

    TGraphErrors *gr;
    gr = new TGraphErrors();
    int npa, npb;
    npa= a->GetN();
    npb= b->GetN();
    for(int p=0;p<npa+npb;p++){
        if(p<npa){
            gr->SetPoint(p,a->GetPointX(p),a->GetPointY(p));
        }else{
            gr->SetPoint(p,b->GetPointX(p-npa),b->GetPointY(p-npa));
        }

    }

    return gr;
}


//trapezoidal integration of a TGraph
double TrapezoidalIntegration(TGraph* graph, double xmin, double xmax) {
	//function written partially by chatgpt
    int npoints = graph->GetN();

    // Ensure the graph has points
    if (npoints < 2) {
        return 0.0;
    }

    // Sort x values to ensure they are in ascending order
    graph->Sort();
    double integral = 0.0;

	//First I will count how many points are inside the range:
	int first_index = 0;
	int last_index = 0;
	int pointsinside =0;

	for(int i=0; i<graph->GetN();i++){
		if(xmin<graph->GetPointX(i) && graph->GetPointX(i)<xmax){ // if the point in inside
			if(pointsinside == 0) first_index = i; //if this is the first point inside, 'i' is the first index
			pointsinside++; // so the number of points inside gets bigger
			if(pointsinside>0) last_index=i; //and if there are points inside, 'i' is the last index..
		}
	}

	//cout<<"there are "<<pointsinside<<" points inside.";

	if(pointsinside == 0 ){
 		return 0.5*(graph->Eval(xmin)+graph->Eval(xmax))*(xmax-xmin); //in case there are no points inside
	}

    // Trapezoidal rule: evaluate integral of points inside
    for (int i = first_index; i < last_index; ++i) {
        double x1, x2, y1, y2;
        graph->GetPoint(i, x1, y1);
        graph->GetPoint(i + 1, x2, y2);
        integral += 0.5 * (x2 - x1) * (y1 + y2);
    }

	//as I used strictly smaller in L583, we can do that:
	double x1, x2, y1, y2;
	//beginning to first point:
	x1 = xmin;
	y1 = graph->Eval(xmin);
	graph->GetPoint(first_index, x2, y2);
	integral += 0.5 * (x2 - x1) * (y1 + y2);

	graph->GetPoint(first_index, x1, y1);
	x2 = xmax;
	y2 = graph->Eval(xmax);

	integral += 0.5 * (x2 - x1) * (y1 + y2);

    return integral;
}

//It works for TGraph and TF1 objects
double SimplifiedTrapezoidalIntegration(TObject* obj, double xmin, double xmax, int Nstep = 1e4) {
    if (!obj) {
        std::cerr << "Erro: Objeto nulo.\n";
        return 0;
    }

    double step = (xmax - xmin) / Nstep;
    double sum = 0;

    std::string className = obj->IsA()->GetName();

    if (auto func = dynamic_cast<TF1*>(obj)) {
        for (int i = 0; i < Nstep; ++i) {
            double x1 = xmin + i * step;
            double x2 = x1 + step;
            double y1 = func->Eval(x1);
            double y2 = func->Eval(x2);
            sum += 0.5 * step * (y1 + y2);
        }
    } else if (auto graph = dynamic_cast<TGraph*>(obj)) {
        for (int i = 0; i < Nstep; ++i) {
            double x1 = xmin + i * step;
            double x2 = x1 + step;
            double y1 = graph->Eval(x1);
            double y2 = graph->Eval(x2);
            sum += 0.5 * step * (y1 + y2);
        }
    } else {
        std::cerr << "Erro: Tipo '" << className << "' não suportado para integração.\n";
        return 0;
    }

    return sum;
}


//function to sort vector based on another vector
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

//function to sort vector<Int_t> based on another vector<Int_t>

void sortVecAsc(std::vector<int>& x, std::vector<int>& v) {
    // Check for size mismatch
    if (x.size() != v.size()) {
        std::cerr << "Error: vectors with different lengths!!" << std::endl;
        return;
    }

    // Bubble sort algorithm
    bool b = false;
    int s = x.size();

    while (!b) {
        b = true;
        for (int i = 0; i < s - 1; i++) {
            if (x[i] > x[i + 1]) {
                b = false;
                std::swap(x[i], x[i + 1]);
                std::swap(v[i], v[i + 1]);
            }
        }
    }
}
double provideTgama(int telN, int year = 2023){//just to save some lines in my script 

	double tgamma;

	if(year == 2023){
		switch(telN){
			case 1:
				tgamma = 15.742;
				break;
			case 2:
				tgamma = 15.700;
				break;
			case 3:
				tgamma = 15.635;
				break;
			case 4:
				tgamma = 15.554;
				break;
			case 5:
				tgamma = 15.465;
				break;
			case 6:
				tgamma = 15.381;
				break;
			case 7:
				tgamma = 15.311;
				break;
			case 8:
				tgamma = 15.264;
				break;
			default:
				tgamma = 15.742;
		}
	}else{
		tgamma = 15;
		cout<<"****************************************************************"<<endl;
		cout<<"	year "<<year<<" nor available!"<<endl;
		cout<<"****************************************************************"<<endl;
	}
	return tgamma;
}


//function to tranform ToF neutrons histo to Energy . v: 17/10/2024
// - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
//if density is true, returns density histo: En/ΔE
//Inputs:
// - htof: tof histogram 
// - density: bool to chose between density or counts histo
// - verbose: define if you want many info in the screen or not
// - L: distance to calculate the ToF
TH1D *NeutronFromToF(TH1D *tofHisto, bool density = false, bool verbose = false, double L = 4.6472 ){//verbose put some info regarding the function working 
	
	//First, I need to see if the ToF is not too small to be tranformed: 
	
	int firstBin =1; //we start in the first bin in tofHisto
	int nBinsTof =tofHisto->GetNbinsX();
	if(verbose)cout<<"Tof input:\n Nbins= "<<nBinsTof<<", from "<<tofHisto->GetBinLowEdge(1)<<" to "<< tofHisto->GetBinLowEdge(nBinsTof+1)<<" ns."<<endl;
	
	while(!En(tofHisto->GetBinLowEdge(firstBin),L)){ //while the low edge of the bin is a tof too small (== En is zero and then '!En' is true)
		if(verbose)cout<<"Not calculable tof. Bin #"<<firstBin<<" => Low edge = "<<tofHisto->GetBinLowEdge(firstBin)<<" ns."<<endl;
		firstBin++;//we go to next bin
		
		if(firstBin==nBinsTof+1){//if we reach the overflow bin
			cout<<"INCOMPATIBLE TOF RANGE!"<<endl;
			return NULL; //the whole thing is a bad placed region, so forget it 
		}
	}


	if(verbose)cout<<"FIRST CALCULABLE BIN =  bin#"<<firstBin<<" => Low edge = "<<tofHisto->GetBinLowEdge(firstBin)<<" ns."<<endl;
	//get the number of bins 
	int nbinsNeutron = nBinsTof+1 - firstBin;// ok 

	if(verbose){
		cout<<"Nbinsneutron ="<<nbinsNeutron<<". path considered = "<<L<<" meters."<<endl;
		cout<<"Proceeding to bin vector construction... "<<endl;
	}

	double lowEdge;
  	double v[nbinsNeutron+1]; //has to include overflow bin, this is the vector for the bin edges!  
	
    for(int t=nBinsTof+1, k=0;t>=firstBin;t--,k++){//t follows the bin number, from <firstBin> to Nbins+1 (overflow bin)
		lowEdge = tofHisto->GetBinLowEdge(t);
		if(verbose) cout<<"TOF bin #"<<t<<","<<lowEdge<<" ns (low edge) --> "<<En(lowEdge,L)<<" MeV."<<endl;
		v[k] = En(lowEdge,L);
	}

    TH1D *hEnn  = new TH1D("hEnn","hEnn",nbinsNeutron,v);
	double energy;
	Int_t counts;
	int binEn;
	
    for(int t=nBinsTof, k=1;t>=firstBin;t--,k++){
	    
		counts = tofHisto->GetBinContent(t);
	
		hEnn->SetBinContent(k,counts);
		if(verbose){
 			if(density) 
				cout<<"Fill[Tof bin # "<<t<<", center = "<<tofHisto->GetBinCenter(t)<<" ns] = "<<hEnn->GetBinCenter(k)<<" MeV, "<<tofHisto->GetBinContent(t)<<" counts. ==> "<<counts<<" corr_counts"<<endl;
			else
				cout<<"Fill[Tof bin # "<<t<<", center = "<<tofHisto->GetBinCenter(t)<<" ns] = "<<hEnn->GetBinCenter(k)<<" MeV, "<<tofHisto->GetBinContent(t)<<" counts."<<endl;
		}
		
	
	}

	if(density){
		for(int i=1;i<=hEnn->GetNbinsX();i++){
			hEnn->SetBinContent(i,hEnn->GetBinContent(i)/hEnn->GetBinWidth(i));
		}
	}
	
	return hEnn;


}



//function to integrate histogram
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



//function to rebin histogram (does not conserve the entries, but the area! I need to recheck it )
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


//function to fetch the cross sections from a csv file:
TGraph * XS_fecther(int telN = 1,bool drawXS = false,bool verbose = false, string xsfile = "/home/e802/Analysis/pre_analysis/cross_secs/H1_MT2_angXS.csv", string yaxis = "d#sigma/d#Omega (b/sr)"){
	 
  if(verbose) cout<<"\n\nLoading cross sections from file "<<lastname(xsfile).c_str()<<"...\n"<<endl;

  TTree *Cs = new TTree("CS", "CSTTree");
  Cs->ReadFile(xsfile.c_str(), "En/D:a/D:b/D:c/D:d/D");//:e/D:f/D:g/D:h/D");
  float Emin, Emax;
  Emin = Cs->GetMinimum("En");
  Emax = Cs->GetMaximum("En");
  Cs->Print();
  if(verbose) cout<<"Emin = "<<Emin*1e-6<<" MeV\tEmax = "<<Emax*1e-6<<" MeV"<<endl;

  int tel = telN;
   //__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__.__. Drawing option
  const char * branch2;

  switch(tel){//data from codata
      case 1:
          branch2 = "a";
          break;//   double D, Nat, NA, m_target, omega_telescope, barn_cm2,omega_target, r_target, a_target;
      case 2:
          branch2 = "b";
          break;
      case 3:
          branch2 = "c";
          break;
      case 4:
          branch2 = "d";
          break;
      case 5:
          branch2 = "e";
          break;
      case 6:
          branch2 = "f";
          break;
      case 7:
          branch2 = "g";
          break;
      case 8:
          branch2 = "h";
          break;
      default:
          branch2 = "a";
  }


//---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.---ooOOOoo---.
//plotting EE XS

      double er, s;
      if(drawXS){
		TCanvas *XS_canvas = new TCanvas("XS","XS",230,240,800,800);
      	XS_canvas->cd()->SetLeftMargin(0.16);
	  }

      TGraph *sgraph = new TGraph(Cs->GetEntries());
      Cs->SetBranchAddress("En", &er);
      Cs->SetBranchAddress(branch2, &s);

      for(Int_t i=0; i<Cs->GetEntries();i++){
          Cs->GetEntry(i);
          sgraph->SetPoint(i,er*1e-6,s);
      }
      sgraph->GetXaxis()->SetTitle("E (MeV)");
      sgraph->GetXaxis()->SetRangeUser(2,40);
      sgraph->GetYaxis()->SetTitle(yaxis.c_str());
      sgraph->SetTitle(Form("ee XS for #theta = %d^{o}", telN*20));
      sgraph->SetLineColor(kRed);
      sgraph->SetLineWidth(2.0);

      if(verbose)cout<<Form("XS value for 02 MeV = %.4f b/sr\nXS value for 40 MeV = %.4f b/sr", sgraph->Eval(2),sgraph->Eval(40))<<endl;

      if(drawXS){
		sgraph->Draw();
		gPad->SetGridx();
		gPad->SetGridy();
		gPad->SetLogy();
	  }

	  return sgraph;
}


//cuntcion to scale a TGraphErrors
void ScaleTGraphErrors(TGraphErrors *tg, double factor){
	for(int t=0;t<tg->GetN();t++)
	{
		tg->SetPointY(t, tg->GetPointY(t)*factor);
		tg->SetPointError(t, tg->GetErrorX(t),tg->GetErrorY(t)*factor);
	}
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


//function to tranform ToF neutrons' histo to Energy 
TH1D *ScaleXhisto(TH1D *toBeScaled, double scaleFactor = 1.0, string name = "scaled",bool verbose = false){//verbose put some info regarding the function working 

	Int_t nbins = toBeScaled->GetNbinsX();
	if(verbose)cout<<"There are "<<nbins<<" bins. The scale factor provided is "<<scaleFactor<<"."<<endl; 
	double xb[nbins+1];
	for(Int_t i=0;i<nbins+1;i++){
		if(verbose){
			cout<<"bin "<<i+1<<": "<<toBeScaled->GetBinLowEdge(i+1)<<" -->"<<scaleFactor*toBeScaled->GetBinLowEdge(i+1)<<endl;
		}
		xb[i] = scaleFactor*toBeScaled->GetBinLowEdge(i+1);
	}

	TH1D *hscaled = new TH1D(name.c_str(),name.c_str(),nbins,xb);
	for(Int_t i=1;i<nbins+1;i++){
		hscaled->SetBinContent(i,toBeScaled->GetBinContent(i));
	}

	return hscaled;

}



double Jacobian(double E, double L = 4.6472) {
    double c = 0.299792458; // m/ns
    //double mn = 939.56542052; // MeV/c^2
	double M = 939.56542052; // MeV

    return L*pow(M, 2)/(c*pow(E*(2*M+E), 3./2));
}




//possible jacobians:: - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
/*
Double_t dtdEn(Double_t En, Double_t L = 4.6472){//EN in MeV and L in m
	Double_t M =  939.56542052;
	Double_t c =  0.299792458;
	return (L/c) * pow(M,2)/pow((En+M),3)* pow(1 - pow(M/(En+M),2),-3.0/2);

}
*/

Double_t dtdEn(Double_t En, Double_t L = 4.6472){//EN in MeV and L in m
//Obtained by Wolfram, but I think its correct!
	Double_t M =  939.56542052;
	Double_t c =  0.299792458;
	return (L/c)* pow(M,2)/pow(En*(En+2*M),3.0/2);

}


Double_t dEndt(Double_t t, Double_t L = 4.6472){//t in ns in MeV and L in m
	Double_t M =  939.56542052;
	Double_t c =  0.299792458;
	return M*pow(L/c,2)/(pow(t,3)) *pow(1 - pow(L/(c*t),2),-3.0/2);
}
// - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - -  - - - - 
//function to tranform ToF neutrons' histo to Energy;
//if density is true, returns density histo: En/ΔE
//Inputs:
// - htof: tof histogram (counts)
// - density: bool to chose between density or counts histo
// - Emin: minimum energy of resulting histo
// - Emax: maximum energy of resulting histo
// - nbins: number of bins you want in the result
// - verbose: print info about the transformation
// - L: distance to calculate the ToF
TH1D *NeutronFromToFJacobian(TH1D *htof, bool density = true,Double_t Emin = 1,Double_t Emax = 50, Int_t nbins = 0, bool verbose = false, double L = 4.6472){//verbose put some info regarding the function working 


	TGraph *gtof = new TGraph();
	// Lets obtain the Tof density plot 
	for (int i = 1; i <= htof->GetNbinsX(); ++i){
		gtof->SetPoint(i-1,htof->GetBinCenter(i),htof->GetBinContent(i)/htof->GetBinWidth(i));
		//gtof->SetPoint(i-1,htof->GetBinCenter(i),htof->GetBinContent(i));
	}
	//test
	//gtof->Draw();


	//verify tof limits
	Int_t firstpoint = 0;
	while(En(gtof->GetPointX(firstpoint++),L)<=0){
		if(firstpoint == (gtof->GetN() -1) ){
			cerr<<"INVALID TOF DATA !"<<endl;
			return nullptr;
		}
	} 
	Double_t Emin_from_tof = ceil(En(htof->GetBinCenter(htof->GetNbinsX()),L));
	Double_t Emax_from_tof = floor(En(htof->GetBinCenter(firstpoint),L));

	if(Emax > Emax_from_tof){
		if(verbose) cout<<"Redefining Emax from "<<Emax<<" MeV to "<<Emax_from_tof<<" MeV."<<endl;
		Emax = Emax_from_tof;
	}
	if(Emin < Emin_from_tof){
		if(verbose) cout<<"Redefining Emin from "<<Emin<<" MeV to "<<Emin_from_tof<<" MeV."<<endl;
		Emin = Emin_from_tof;
	}

	if(verbose){
	cout<<"first point: "<<firstpoint<<" -->"<< gtof->GetPointX(firstpoint) << "ns ("<< En(gtof->GetPointX(firstpoint),L)<<" MeV)"<<endl; 
	cout<<"last point: "<<(gtof->GetN() -1)<<" -->"<< gtof->GetPointX(gtof->GetN() -1) << "ns ("<< En(gtof->GetPointX(gtof->GetN() -1),L)<<" MeV)"<<endl; 
	}


	// // Create the Energy histogram
	if(!nbins) nbins = htof->GetNbinsX(); // if nbins == 0, then it will assume the same number of bins from the ToF plot 

	TH1D* hEn = new TH1D("hEn", "Energy Histogram",nbins , Emin, Emax);

	// // Loop over each bin of the ToF histogram
	for (int i = 1; i <= nbins; ++i) {
		Double_t energy = hEn->GetBinCenter(i);
		Double_t tof = ToFparticle(energy,'n', L );
		Double_t correction = dtdEn(energy, L);
		hEn->SetBinContent(i,gtof->Eval(tof)*correction);
		
		if(verbose)cout<<"En = "<<energy<<" MeV ["<<tof<<" ns] -> correction = "<<correction<<endl;
	}

	//if density == false, multiply each binning by ΔE, to return the result in counts:
	if(!density){
		for(int i=1;i<hEn->GetNbinsX();i++)hEn->SetBinContent(i,hEn->GetBinContent(i)*hEn->GetBinWidth(i));
	}

 	return hEn;

}



//funtcion to scale a TGraphErrors
void ScaleTGraph(TGraph *tg, double factor){
	for(int t=0;t<tg->GetN();t++)
	{
		tg->SetPointY(t, tg->GetPointY(t)*factor);
	}
}
