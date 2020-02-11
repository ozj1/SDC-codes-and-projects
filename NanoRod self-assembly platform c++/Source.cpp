#include <iostream> // in order to using cin, cout, endl, ...
#include <sstream>   // using for getin data func--intel compiler in cluster
#include <fstream>   // for getting ofstream out("out.txt"); as a cout as a trick for avoiding high amount for dimensions because randNumx[plb][plb] will give us stackoverflow error
#include <string>   // for getting ofstream out("out.txt");
#include <stdlib.h>  // for using rand_s
#include <time.h>   // to use srand, time
#include <math.h>   // because of POW sqrt, abs, ...
#include <cmath>
#include<vector>    // this code use vector library which is said to be the most reliable way for pulling doubles from .txt
#include<cstdlib>
#include <iomanip>
#include <omp.h>   //parallel
#include <stdio.h>

#define x1 NanoRod{numrod,nanolength1,nanoCoreRad,coat_vlu} 
#define x2 x1, x1 
#define x4 x2, x2 
#define x8 x4, x4 
#define x16 x8, x8 
#define x32 x16, x16 
#define x64 x32, x32
#define x128 x64, x64 
#define x256 x128, x128 
#define x512 x256, x256 // numrod needs to be initialized by this (manually)


using namespace std;

//INTEL COMPILER
namespace patch
{
	template < typename T > std::string to_string(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

const char* generatefilename(int x, string prename) {
	string filename = "_";
	filename += patch::to_string(x);
	filename += prename;
	filename += "_";

	filename += ".txt";
	return filename.c_str();
}




// Global Parameters
const int numrod = 512;
const double padLgth_vlu = 150.; //hight of pad
const double coat_vlu = 25.; //coating length of nano rods
const double nanoCoreRad = 145.;
const double partRad_valu = coat_vlu + nanoCoreRad;
const double nanolength1 = 2379.;


const double plb = 2650., phb = 19800., plby = 6000., phby = 10800., lx = 22500., ly = 50000., z0 = 240.;

const double zetaau = -22.e-3, pH = 4.9;
const double zicoef = 0.02;

const int numRegionX = 10;  // for assigning regions to nanorods ( calculating dist and overlap for the neighbors.. not for all)
const int numRegionY = 20;


const double e0 = 1.60217657e-19, kT = 4.0453001636e-21, Na = 6.022e+23, eps = 80.1, eps0 = 8.854187817e-12;
const double Aau = 45.5e-20, Asio2 = 7.5e-20, Aw = 3.7e-20;
const double  pi = 3.1415926535897, zseg = 630.695*e0;
const double kBoltzmann = 1.38064852e-23;
const double Temp = 293.;

double maxlenght = nanolength1;
double cutoff1 = 5000., cutoff2 = 7000.;
double cutoffInSS = 1000., cutoffInSE = 900., cutoffInEE = 600.;

const int engSize = 100000;
const int timeStep = 100000;
const int totOrderSize = 100000;

//order parameter constants
const int xmeshn = 500; const int ymeshn = 1; //Kristen: no need for considering it in y directin 
const double cutfx = (lx / xmeshn) * 30;
const double cutfy = (ly / ymeshn) * 30;
double oriordprm[xmeshn][ymeshn][timeStep];
double oriordprmAvg[xmeshn][ymeshn];
double counter = 0.0;

/// parameters for printing
int printEng = 10000 * numrod;
int stepEng = 100 * numrod;
int printSphereCoor = 10000 * numrod;
int printCyl = 10000 * numrod;
int printOrder = 10000 * numrod;
int enoughStepOrder = 20000 * numrod;
int calOrder = 1000 * numrod;
int printCoor = 10000 * numrod;

int calTotOrder = 100 * numrod;
int printTotOrder = 10000 * numrod;

//RECIVING FILE

struct Info {
	double xvalue[numrod];
	double yvalue[numrod];
	double zvalue[numrod];
	double tethavalue[numrod];
	double length[numrod];
	double coreradf[numrod];
	double coatThickf[numrod];
	int numrod[numrod];
	int totalinfo;

};


Info getdata() {
	Info fromprevious;
	vector<double> values;
	double in = 0.0;

	ifstream infile;
	infile.open("_191Coordinates_.txt");//remember to save a text file in exact location that "out" file exist  
	while (infile >> in) {
		values.push_back(in);

	}
	int jj = 0;
	for (jj = 0; jj < values.size(); jj++) {


		//cout << values[j] << endl;
		if (jj % 7 == 0) {
			fromprevious.xvalue[jj / 7] = values[jj];
		}
		else if (jj % 7 == 1) {
			fromprevious.yvalue[jj / 7] = values[jj];
		}
		else if (jj % 7 == 2) {
			fromprevious.tethavalue[jj / 7] = values[jj];// it rotates periodiacally from x-y-...
		}
		else if (jj % 7 == 3) {
			fromprevious.length[jj / 7] = values[jj];// it rotates periodiacally from x-y-...
		}
		else if (jj % 7 == 4) {
			fromprevious.coreradf[jj / 7] = values[jj];// it rotates periodiacally from x-y-...
		}
		else if (jj % 7 == 5) {
			fromprevious.coatThickf[jj / 7] = values[jj];// it rotates periodiacally from x-y-...
		}
		else if (jj % 7 == 6) {
			fromprevious.numrod[jj / 7] = values[jj];// it rotates periodiacally from x-y-...
		}

	}
	fromprevious.totalinfo = jj;

	return fromprevious;
}




//usefull global functions
double sign(double x)
{
	if (x > 0) return 1;
	if (x < 0) return -1;
	return 0;
}

double fRand(double fMin, double fMax)
{

	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}


double Dist(double center1X, double center1Y, double center2X, double center2Y) {
	double dist = 0.;
	if ((abs(center1X - center2X) <= (lx / 2.)) && (abs(center1Y - center2Y) <= (ly / 2.))) {
		dist = sqrt((center2Y - center1Y)*(center2Y - center1Y) + (center2X - center1X)*(center2X - center1X));
	}
	else if ((abs(center1X - center2X) > (lx / 2.)) && (abs(center1Y - center2Y) <= (ly / 2.))) {
		dist = sqrt((center2Y - center1Y)*(center2Y - center1Y) + (lx - abs((center2X - center1X)))*(lx - abs((center2X - center1X))));
	}
	else if ((abs(center1X - center2X) <= (lx / 2.)) && (abs(center1Y - center2Y) > (ly / 2.))) {
		dist = sqrt((ly - abs((center2Y - center1Y)))*(ly - abs((center2Y - center1Y))) + (center2X - center1X)*(center2X - center1X));
	}
	else if ((abs(center1X - center2X) > (lx / 2.)) && (abs(center1Y - center2Y) > (ly / 2.))) {

		dist = sqrt((ly - abs((center2Y - center1Y)))*(ly - abs((center2Y - center1Y))) + (lx - abs((center2X - center1X)))*(lx - abs((center2X - center1X))));
	}

	return dist;
}




bool overlap(double dist, double rad1, double rad2) {
	bool value;
	if (dist > (rad1 + rad2)) {// come back to it
		value = false;

	}
	else {
		value = true;
	}

	return value;
}

class Energy {
public:






	double energySpherePadVDW(double x_pr, double z_pr, double r_pr, double plb_pr, double phb_pr, double coat_pr, double padLenght_pr) {

		double x = x_pr * 1.0e-9, z = z_pr * 1.0e-9, R = r_pr * 1.0e-9, plb = plb_pr * 1.0e-9, phb = phb_pr * 1.0e-9, coat = coat_pr * 1.0e-9, padLgth = padLenght_pr * 1.0e-9;
		double Agw, Ags;
		Ags = (pow(Aau, 1. / 2) - pow(Aw, 1. / 2)) * (pow(Asio2, 1. / 2) - pow(Aw, 1. / 2));
		Agw = pow((pow(Aau, 1. / 2) - pow(Aw, 1. / 2)), 2.);
		return (-1. / (12))*(Agw *Evdw4(x, z, R, plb, phb, padLgth) + Ags * Evdw4(x, z, coat + R, plb, phb, padLgth) - Ags * Evdw4(x, z, R, plb, phb, padLgth));
	}

	double energySphereSpherVDW(double dist, double R, double coat)
	{//it includes all core shell, shell shell, core core interactions

		double Agw, Ags, Asw;
		Ags = (pow(Aau, 1. / 2) - pow(Aw, 1. / 2)) * (pow(Asio2, 1. / 2) - pow(Aw, 1. / 2));
		Agw = pow((pow(Aau, 1. / 2) - pow(Aw, 1. / 2)), 2.);
		Asw = pow((pow(Asio2, 1. / 2) - pow(Aw, 1. / 2)), 2.);
		//no problem if R and coat are not in nm scale because x and y in Hr function are ratios

		//encourted a strange problem don't write 1/12 will make zero in results, use 1./12 !!!!
		return (-1. / 12)* (Agw*Hr((dist - 2.*R) / (2.*R), 1.) + 2. * Ags*Hr((dist - 2.*R - coat) / (2.*R + 2.*coat), R / (R + coat)) - 2. * Ags*Hr((dist - 2.*R) / (2.*R), 1.) + Asw * Hr((dist - 2.*R - 2.*coat) / (2.*R + 2.*coat), 1.) - Asw * Hr((dist - 2.*R) / (2.*R), 1.));
	}


	double ESt_SR_SR(double dist_prm)
	{
		double dist = dist_prm * 1.0e-9;
		double pH = 4.9, eps = 80.1, eps0 = 8.854187817e-12, e0 = 1.60217657e-19, kT = 4.0453001636e-21, R = 145.e-9, ds = 25.e-9;
		double Rt = R + ds;
		return eps * eps0*pow((kT / e0), 2.)* pow(yy(dist), 2)*(pow(Rt, 2.) / dist)*log(1 + exp(-kappa(pH) * (dist - 2.*R - 2.*ds)));

	}


	double ESt_SR_ER(double dist_pr)
	{
		double dist = dist_pr * 1.0e-9;
		double  partRad = partRad_valu * 1.e-9, zi = zicoef * zseg;

		return (exp(kappa(pH)*partRad) / (1. + kappa(pH)*partRad)) * (zi*zseg / (4. * pi*eps*eps0))*exp(-kappa(pH)*dist) / dist;
	}

	double ESt_ER_ER(double r)
	{
		r = r * 1.0e-9;
		double pH = 4.9, eps = 80.1, eps0 = 8.854187817e-12, e0 = 1.60217657e-19, kT = 4.0453001636e-21, pi = 3.1415926535897, zseg = 630.695*e0, zi = zicoef * zseg;

		return (pow(zi, 2.) / 1.)*(exp(-kappa(pH) * r) / (4. * pi*eps*eps0)) / r;
	}

	

private:
	//edw1 infinite plane...edw2 semi infinte plane..edw3 infinite -2* semi infinte(slit)..edw4 infinte pad.. edw5 modified sphere with pad
	double Evdw1(double x, double z, double r) {

		return (2. * r *x) / (-pow(r, 2.) + pow(x, 2.)) + (2. * r *z) / (-pow(r, 2.) + pow(z, 2.)) - (2. * r *x *z *(-2. * pow(r, 2.) + pow(x, 2.) + pow(z, 2.))) / ((-pow(r, 2.) + pow(x, 2.))*(-pow(r, 2.) + pow(z, 2.))* sqrt(-pow(r, 2.) + pow(x, 2.) + pow(z, 2.))) + log(pow(((-r + abs(x)) / (r + abs(x))), sign(x)) *pow(((-pow(r, 2.) + pow(x, 2.) + r * z + sqrt(-pow(r, 2.) + pow(x, 2.) + pow(z, 2.))*abs(x)) / (-pow(r, 2.) + pow(x, 2.) - r * z + sqrt(-pow(r, 2.) + pow(x, 2.) + pow(z, 2.))* abs(x))), sign(x))*pow(abs((-r + z) / (r + z)), (1 - sign(x))));

	}


	double Evdw2(double z, double r) {

		return (4. * r* z) / (-pow(r, 2) + pow(z, 2)) + 2. * log((-r + z) / (r + z));
	}


	double Evdw3(double x, double z, double r, double plb, double phb) {

		return Evdw2(z, r) - Evdw1((phb - plb) - x, z, r) - Evdw1(x, z, r);
	}


	double Evdw4(double x, double z, double r, double plb, double phb, double padLgth) {


		return Evdw3(x - plb, z, r, plb, phb) - Evdw3(x - plb, z + padLgth, r, plb, phb);
	}



	//edw1 infinite plane...edw2 semi infinte plane..edw3 infinite -2* semi infinte(slit)..edw4 infinte pad.. edw5 modified sphere with pad


	double Hr(double x, double y)  //Hr is an internal function for calculation of vdw among spheres
	{

		return y / (pow(x, 2.) + x * y + x) + y / (pow(x, 2.) + x * y + x + y) + 2.*log(((pow(x, 2.) + x * y + x) / (pow(x, 2.) + x * y + x + y)));
	}


	//electro static section calculation
	double kappa(double pH)
	{
		double c;
		c = pow(10., -pH)* 1.e+3;
		return sqrt((2.* pow(e0, 2.)*Na*c) / (eps* eps0* kT));
	}

	double ys(double zetaau)
	{
		return (zetaau*e0) / (kT);
	}

	double yy(double r)
	{
		double  R_SI = nanoCoreRad * 1.e-9, coat_SI = coat_vlu * 1.0e-9;
		return 4.*exp(0.5 *kappa(pH)* (r - 2.*R_SI - 2.*coat_SI)) * atanh(exp(-1.*kappa(pH)*(r - 2.*R_SI - 2.* coat_SI) / 2.)* tanh(ys(zetaau) / 4.));

	}







};



class NanoRod : public Energy {
public:
	NanoRod(int rodNum, double length, double coreRad, double coatThick)
		: rodNum{ rodNum }, length{ length }, coreRadf{ coreRad }, coatThickf{ coatThick }, theta{ 0 }, center{ 0,0,0 }\
		, numSphere{ int(floor(length / 2. / (coreRadf + coatThickf)) + 1) }, engWthPad{ 0 }, engWthOthr{ 0 }, engTot{ engWthPad + engWthOthr }\
		, spherRad{ length / 2 / numSphere }, coatThick{ coatThickf }, coreRad{ spherRad - coatThick }\
		, endCenter1X{ 0 }, endCenter1Y{ 0 }, endCenter2X{ 0 }, endCenter2Y{ 0 }\
		, regionx{ regionX(center[0], numRegionX) }, regiony{ regionX(center[0], numRegionX) }{

	} //constructor


	class NanoSphere {


	public:
		NanoSphere()
			:centerSphere{ 0,0,0 }, engWthPad{ 0 }, engWthOthr{ 0 }, engTot{ engWthOthr + engWthPad } {

		}





		double centerSphere[3];

		double engWthPad;
		double engWthOthr;
		double engTot;









	};

	NanoSphere *nanosphere;// for nanorods
	NanoSphere *nanospheres;// for nanorodtry




	int rodNum;
	double length;
	double coreRadf;
	double coatThickf;
	double theta;
	double center[3];



	int numSphere;


	double engWthPad;
	double engWthOthr;
	double engTot;
	double spherRad;
	double coatThick;
	double coreRad;


	double endCenter1X, endCenter1Y, endCenter2X, endCenter2Y;


	int regionx;
	int regiony;

	void PBCCheck() {
		if ((center[0] <= lx - length / 2.) && (center[0] >= length / 2.) && (center[1] <= ly - length / 2.) && (center[1] >= length / 2.)) {
			//do nothing
		}
		else {
			if (center[0] < lx) {
				center[0] = lx + center[0];
			}
			else if (center[0] > lx) {
				center[0] = -lx + center[0];
			}
			if (center[1] < ly) {
				center[1] = ly + center[1];
			}
			else if (center[1] > ly) {
				center[1] = -ly + center[1];
			}


		}

	}

	void updatepoints() {//for calculating electro energy

		endCenter1X = center[0] + length * cos(theta) / 2.;
		endCenter1Y = center[1] + length * sin(theta) / 2.;
		endCenter2X = center[0] - length * cos(theta) / 2.;
		endCenter2Y = center[1] - length * sin(theta) / 2.;

		///if pbc
		if ((center[0] <= lx - length / 2.) && (center[0] >= length / 2.) && (center[1] <= ly - length / 2.) && (center[1] >= length / 2.)) {
			//do nothing
		}
		else {
			if (center[0] > lx - length / 2.) {
				if (endCenter1X > lx) {
					endCenter1X = -lx + endCenter1X;
				}
				if (endCenter2X > lx) {
					endCenter2X = -lx + endCenter2X;
				}

			}

			if (center[0] < length / 2.) {
				if (endCenter1X < 0) {
					endCenter1X = lx + endCenter1X;
				}
				if (endCenter2X < 0) {
					endCenter2X = lx + endCenter2X;
				}
			}


			if (center[1] > ly - length / 2.) {
				if (endCenter1Y > ly) {
					endCenter1Y = -ly + endCenter1Y;
				}
				if (endCenter2Y > ly) {
					endCenter2Y = -ly + endCenter2Y;
				}
			}

			if (center[1] < length / 2.) {
				if (endCenter1Y < 0) {
					endCenter1Y = ly + endCenter1Y;
				}
				if (endCenter2Y < 0) {
					endCenter2Y = ly + endCenter2Y;
				}

			}


		}

	}

	void updateEnergy() {
		engTot = engWthOthr + engWthPad;
	}

	int regionX(double centerx, int numRegionX) {
		double result = (centerx* numRegionX / lx);
		return int(result);
	}

	int regionY(double centery, int numRegionY) {
		double result = (centery* numRegionY / ly);
		return int(result);

	}







	double findingBestFitRadOutside() {


		double remain = fmod(length, 2 * (coreRadf + coatThickf));
		double fitRad = (coreRadf + coatThickf) + (remain / numSphere);
		return fitRad;
	}





public:
	struct XY
	{
		double x;
		double y;

	};

	XY findingCenters(int iNum) {

		//double diang = sqrt(length*length + (coreRad + coatThick)*(coreRad + coatThick));  not using because of "spherocylendir model"

		double centerSphereX, centerSphereY;
		if (numSphere % 2 == 1) {

			centerSphereX = center[0] - 2 * (numSphere / 2 - iNum)*spherRad*cos(theta);
			centerSphereY = center[1] - 2 * (numSphere / 2 - iNum)*spherRad*sin(theta);
			//PBC
			if ((center[0] <= lx - length / 2.) && (center[0] >= length / 2.) && (center[1] <= ly - length / 2.) && (center[1] >= length / 2.)) {
				//do nothing
			}
			else {
				if (center[0] > lx - length / 2.) {
					if (centerSphereX > lx) {
						centerSphereX = -lx + centerSphereX;
					}
				}

				if (center[0] < length / 2.) {
					if (centerSphereX < 0) {
						centerSphereX = lx + centerSphereX;
					}
				}


				if (center[1] > ly - length / 2.) {
					if (centerSphereY > ly) {
						centerSphereY = -ly + centerSphereY;
					}
				}

				if (center[1] < length / 2.) {
					if (centerSphereY < 0) {
						centerSphereY = ly + centerSphereY;
					}

				}


			}

		}

		else {

			centerSphereX = center[0] - (2 * (numSphere / 2 - iNum) - 1)*spherRad*cos(theta);
			centerSphereY = center[1] - (2 * (numSphere / 2 - iNum) - 1)*spherRad*sin(theta);

			//PBC
			if ((center[0] <= lx - length / 2.) && (center[0] >= length / 2.) && (center[1] <= ly - length / 2.) && (center[1] >= length / 2.)) {
				//do nothing
			}
			else {
				if (center[0] > lx - length / 2.) {
					if (centerSphereX > lx) {
						centerSphereX = -lx + centerSphereX;
					}
				}

				if (center[0] < length / 2.) {
					if (centerSphereX < 0) {
						centerSphereX = lx + centerSphereX;
					}
				}


				if (center[1] > ly - length / 2.) {
					if (centerSphereY > ly) {
						centerSphereY = -ly + centerSphereY;
					}
				}

				if (center[1] < length / 2.) {
					if (centerSphereY < 0) {
						centerSphereY = ly + centerSphereY;
					}

				}


			}


		}
		XY  result = { centerSphereX,centerSphereY };

		return result;
	}








private:



};






int main(int argc, char* argv[]) {



	srand(time(NULL));

	NanoRod nanorod[numrod] = { x512 };// need to be adjusted manually.. number of arrays elements..for initialization




	ofstream out;

	double totOrderAvg[totOrderSize];
	double totOrder[totOrderSize];
	double energyTot[engSize];
	// initializing the system.. nanorods


	for (int j = 0; j < numrod; j++) {
		nanorod[j].nanosphere = new NanoRod::NanoSphere[nanorod[j].numSphere];
	}

	Info fromPrevious;
	fromPrevious = getdata();



	for (int j = 0; j < numrod; j++) {

		nanorod[j].center[0] = fromPrevious.xvalue[j];
		nanorod[j].center[1] = fromPrevious.yvalue[j];
		nanorod[j].theta = fromPrevious.tethavalue[j];
		nanorod[j].length = fromPrevious.length[j];
		nanorod[j].coreRadf = fromPrevious.coreradf[j];
		nanorod[j].coatThickf = fromPrevious.coatThickf[j];
		nanorod[j].rodNum = fromPrevious.numrod[j];

		for (int k = 0; k < nanorod[j].numSphere; k++) {// finds the sphere positions for each 

			nanorod[j].nanosphere[k].centerSphere[0] = nanorod[j].findingCenters(k).x;
			nanorod[j].nanosphere[k].centerSphere[1] = nanorod[j].findingCenters(k).y;
		}

	}






	//Assigning Energy

//#pragma omp parallel for private(k,i,l,dist1,dist2,dist3,dist4,dist5,dist6,dist7)
	for (int j = 0; j < numrod; j++) {


		if ((nanorod[j].center[0] < phb + cutoff2) && (nanorod[j].center[0] > plb - cutoff2)) {

			for (int k = 0; k < nanorod[j].numSphere; k++) {
				nanorod[j].engWthPad += nanorod[j].energySpherePadVDW(nanorod[j].nanosphere[k].centerSphere[0], z0, nanorod[j].coreRad, plb, phb, nanorod[j].coatThick, padLgth_vlu);
			}
		}

		for (int i = 0; i < numrod; i++) {

			if (i != j) {

				double rodDist = Dist(nanorod[j].center[0], nanorod[j].center[1], nanorod[i].center[0], nanorod[i].center[1]);

				if (rodDist > cutoff1) {
					//do nothing
				}
				else {
					nanorod[j].updatepoints();
					nanorod[i].updatepoints();
					double dist1 = Dist(nanorod[j].endCenter1X, nanorod[j].endCenter1Y, nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
					double dist2 = Dist(nanorod[j].endCenter2X, nanorod[j].endCenter2Y, nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
					double dist3 = Dist(nanorod[j].endCenter1X, nanorod[j].endCenter1Y, nanorod[i].endCenter2X, nanorod[i].endCenter2Y);
					double dist4 = Dist(nanorod[j].endCenter2X, nanorod[j].endCenter2Y, nanorod[i].endCenter2X, nanorod[i].endCenter2Y);

					if (dist1 < cutoffInEE) {
						nanorod[j].engWthOthr += nanorod[j].ESt_ER_ER(dist1);
					}
					if (dist2 < cutoffInEE) {
						nanorod[j].engWthOthr += nanorod[j].ESt_ER_ER(dist2);
					}
					if (dist3 < cutoffInEE) {
						nanorod[j].engWthOthr += nanorod[j].ESt_ER_ER(dist3);
					}
					if (dist4 < cutoffInEE) {
						nanorod[j].engWthOthr += nanorod[j].ESt_ER_ER(dist4);
					}

					
					for (int k = 0; k < nanorod[j].numSphere; k++) {
						double dist5 = Dist(nanorod[j].nanosphere[k].centerSphere[0], nanorod[j].nanosphere[k].centerSphere[1], nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
						double dist6 = Dist(nanorod[j].nanosphere[k].centerSphere[0], nanorod[j].nanosphere[k].centerSphere[1], nanorod[i].endCenter2X, nanorod[i].endCenter2Y);

						if (dist5 < cutoffInSE) {
							nanorod[j].engWthOthr += nanorod[j].ESt_SR_ER(dist5);
						}
						if (dist6 < cutoffInSE) {
							nanorod[j].engWthOthr += nanorod[j].ESt_SR_ER(dist6);
						}

						for (int l = 0; l < nanorod[i].numSphere; l++) {
							double dist7 = Dist(nanorod[j].nanosphere[k].centerSphere[0], nanorod[j].nanosphere[k].centerSphere[1], nanorod[i].nanosphere[l].centerSphere[0], nanorod[i].nanosphere[l].centerSphere[1]);
							//if (dist7 < 2 * (nanorod[j].coreRadf + nanorod[j].coatThickf)) { // Testing the overlap
							//	cout << dist7 << "\n" << nanorod[j].center[0] << "," << nanorod[j].center[1] << "\n" << nanorod[i].center[0] << "," << nanorod[i].center[1] << "\n\n\n";
							//}
							if (dist7 <= 385) {// testing
								dist7 = 385;
							}

							if (dist7 < cutoffInSE) {

								nanorod[j].engWthOthr += nanorod[j].ESt_SR_ER(dist7);
							}

						}




					}



				}
			}
		}
		nanorod[j].updateEnergy();


	}



	NanoRod nanorodtry{ numrod + 1,nanolength1,nanoCoreRad,coat_vlu };

	nanorodtry.nanospheres = new NanoRod::NanoSphere[nanorodtry.numSphere];


	bool acceptance = false;
	double acceptanceRate = 0;
	int total = 0;

	for (long long j = 0; j < 100000000 * long(numrod); j++) {
		total++;
		
			int min = 1, max = numrod;
			int randi = rand() % (max - min + 1) + min - 1;


			double radnum = (double)fRand(0., 1.);



			if (0. <= radnum && radnum < (1. / 3)) {
				double randdeltax = (double)fRand(-500., 500.);
				nanorodtry.center[0] = nanorod[randi].center[0] + randdeltax;
				nanorodtry.center[1] = nanorod[randi].center[1];
				nanorodtry.theta = nanorod[randi].theta;
			}

			if (1. / 3 <= radnum && radnum < (2. / 3)) {
				double randdeltay = (double)fRand(-500., 500.);
				nanorodtry.center[0] = nanorod[randi].center[0];
				nanorodtry.center[1] = nanorod[randi].center[1] + randdeltay;
				nanorodtry.theta = nanorod[randi].theta;

			}

			if (2. / 3 <= radnum && radnum <= 1) {
				double  randdeltaTetha = (double)fRand(-pi / 10, pi / 10);
				nanorodtry.center[0] = nanorod[randi].center[0];
				nanorodtry.center[1] = nanorod[randi].center[1];
				nanorodtry.theta = nanorod[randi].theta + randdeltaTetha;

			}
			nanorodtry.PBCCheck(); // center of nanorods need to be in window all the time

			for (int k = 0; k < nanorodtry.numSphere; k++) {// finds the sphere positions for each 

				nanorodtry.nanospheres[k].centerSphere[0] = nanorodtry.findingCenters(k).x;
				nanorodtry.nanospheres[k].centerSphere[1] = nanorodtry.findingCenters(k).y;

			}




			// Need to modified for 2nd version of cluster ready
			bool overlap_p = false; //overlap is assumed false unless otherwise stated
			for (int i = 0; i < numrod; i++) { // for rod j and rod i which was previously accepted
				if (i != randi) {
					double rodDist = (Dist(nanorodtry.center[0], nanorodtry.center[1], nanorod[i].center[0], nanorod[i].center[1]));
					double safeDist = (nanorodtry.length + nanorod[i].length) / 2.0;
					//bool xPBC = ((nanorodtry.center[0] - nanorod[i].center[0]) > lx - (nanorod[i].length + nanorodtry.length) / 2.);
					//bool yPBC = ((nanorodtry.center[1] - nanorod[i].center[1]) > ly - (nanorod[i].length + nanorodtry.length) / 2.);
					if ((rodDist < safeDist)) {// || xPBC || yPBC) { // if they are neighbors

						for (int k = 0; k < nanorodtry.numSphere; k++) { // for all spheres belong to rod try



							for (int l = 0; l < nanorod[i].numSphere; l++) { // for all spheres belong to rod i




								double dist = Dist(nanorodtry.nanospheres[k].centerSphere[0], nanorodtry.nanospheres[k].centerSphere[1], nanorod[i].nanosphere[l].centerSphere[0], nanorod[i].nanosphere[l].centerSphere[1]);
								if (overlap(dist, nanorodtry.coreRadf + nanorodtry.coatThickf, nanorod[i].coreRadf + nanorod[i].coatThickf) == true) {
									overlap_p = true;



									break;

								}



							}


							if (overlap_p == true) {
								break;
							}

						}


					}


					if (overlap_p == true) {


						break;

					}


				}
			}

			if (overlap_p == false) {
				acceptance = false;

				nanorodtry.updatepoints();

				//energy with pad:
				if ((nanorodtry.center[0] < phb + cutoff2) && (nanorodtry.center[0] > plb - cutoff2)) {

					for (int k = 0; k < nanorodtry.numSphere; k++) {
						nanorodtry.engWthPad += nanorodtry.energySpherePadVDW(nanorodtry.nanospheres[k].centerSphere[0], z0, nanorodtry.coreRad, plb, phb, nanorodtry.coatThick, padLgth_vlu);
					}
				}


				//energy with others
				double dist1;
				double dist2;
				double dist3;
				double dist4;
				double dist5;
				double dist6;
				double dist7;
				int k, l;
				double helpScaler = 0;
#pragma omp parallel for reduction(+: helpScaler) private(k,l,dist1,dist2,dist3,dist4,dist5,dist6,dist7)
				for (int i = 0; i < numrod; i++) {

					if (i != randi) {

						double rodDist = (Dist(nanorodtry.center[0], nanorodtry.center[1], nanorod[i].center[0], nanorod[i].center[1]));

						if (rodDist > cutoff1) {
							//do nothing
						}
						else {

							nanorod[i].updatepoints();// come back to this
							dist1 = Dist(nanorodtry.endCenter1X, nanorodtry.endCenter1Y, nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
							dist2 = Dist(nanorodtry.endCenter2X, nanorodtry.endCenter2Y, nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
							dist3 = Dist(nanorodtry.endCenter1X, nanorodtry.endCenter1Y, nanorod[i].endCenter2X, nanorod[i].endCenter2Y);
							dist4 = Dist(nanorodtry.endCenter2X, nanorodtry.endCenter2Y, nanorod[i].endCenter2X, nanorod[i].endCenter2Y);


							if (dist1 < cutoffInEE) {
								helpScaler += nanorodtry.ESt_ER_ER(dist1);
							}
							if (dist2 < cutoffInEE) {
								helpScaler += nanorodtry.ESt_ER_ER(dist2);
							}
							if (dist3 < cutoffInEE) {
								helpScaler += nanorodtry.ESt_ER_ER(dist3);
							}
							if (dist4 < cutoffInEE) {
								helpScaler += nanorodtry.ESt_ER_ER(dist4);
							}

							for (k = 0; k < nanorodtry.numSphere; k++) {
								dist5 = Dist(nanorodtry.nanospheres[k].centerSphere[0], nanorodtry.nanospheres[k].centerSphere[1], nanorod[i].endCenter1X, nanorod[i].endCenter1Y);
								dist6 = Dist(nanorodtry.nanospheres[k].centerSphere[0], nanorodtry.nanospheres[k].centerSphere[1], nanorod[i].endCenter2X, nanorod[i].endCenter2Y);

								if (dist5 < cutoffInSE) {
									helpScaler += nanorodtry.ESt_SR_ER(dist5);
								}
								if (dist6 < cutoffInSE) {
									helpScaler += nanorodtry.ESt_SR_ER(dist6);
								}

								for (l = 0; l < nanorod[i].numSphere; l++) {
									dist7 = Dist(nanorodtry.nanospheres[k].centerSphere[0], nanorodtry.nanospheres[k].centerSphere[1], nanorod[i].nanosphere[l].centerSphere[0], nanorod[i].nanosphere[l].centerSphere[1]);
									if (dist7 <= 385) {
										dist7 = 385;
									}

									if (dist7 < cutoffInSS) {
										helpScaler += nanorodtry.ESt_SR_SR(dist7) + nanorodtry.energySphereSpherVDW(dist7, nanorodtry.coreRad, nanorodtry.coatThick);//nanorodtry.engWthOthr
									}


								}



							}



						}
					}
				}
				nanorodtry.engWthOthr = helpScaler;
				nanorodtry.updateEnergy();

				if (nanorodtry.engTot <= nanorod[randi].engTot) {
					acceptance = true;
				}
				if (nanorodtry.engTot > nanorod[randi].engTot) {
					double p = exp((nanorod[randi].engTot - nanorodtry.engTot) / (kBoltzmann*Temp));
					double rad = (double)fRand(0., 1.);
					if (p >= rad) {
						acceptance = true;


					}
				}
				if (acceptance == true) {
					nanorod[randi].center[0] = nanorodtry.center[0];
					nanorod[randi].center[1] = nanorodtry.center[1];
					nanorod[randi].theta = nanorodtry.theta;



					for (int k = 0; k < nanorod[randi].numSphere; k++) {// finds the sphere positions for each 

						nanorod[randi].nanosphere[k].centerSphere[0] = nanorod[randi].findingCenters(k).x;
						nanorod[randi].nanosphere[k].centerSphere[1] = nanorod[randi].findingCenters(k).y;

					}
					nanorod[randi].engWthOthr = nanorodtry.engWthOthr;// accepting the values
					nanorod[randi].engWthPad = nanorodtry.engWthPad;
					nanorod[randi].engTot = nanorodtry.engTot;



					

				}



			}
			if (acceptance == false || overlap_p == true) {
				j--;
			}

			nanorodtry.engWthOthr = 0.;//setting values to zero again for future use..nanorodtry is agent carrying potential monte carlo step
			nanorodtry.engWthPad = 0.;
			nanorodtry.updateEnergy();








			/// printing the results

			//printing the spheres coordinates
			if (j % printSphereCoor == 0) {


				out.open(generatefilename((j / printSphereCoor), "Spherecoordinates"));
				int helpNum = 0;
				out << "{";

				for (int p = 0; p < numrod; p++) {




					for (int i = 0; i < nanorod[p].numSphere; i++) {
						helpNum++;// for removing one","

						if (helpNum != 1) {

							out << ",{" << nanorod[p].nanosphere[i].centerSphere[0] << ",";
							out << nanorod[p].nanosphere[i].centerSphere[1] << ",0}";
						}

						else {
							out << "{" << nanorod[p].nanosphere[i].centerSphere[0] << ",";
							out << nanorod[p].nanosphere[i].centerSphere[1] << ",0}";

						}
					}



				}
				out << "}";
				out.close();
			}



			//finding energy at specific time
			if (j % stepEng == 0) {
				energyTot[j / stepEng] = 0;
				double helpScaler = 0;

#pragma omp parallel for reduction(+: helpScaler) 
				for (int i = 0; i < numrod; i++) {
					helpScaler += nanorod[i].engWthOthr / 2 / e0 + nanorod[i].engWthPad / e0;
				}
				energyTot[j / stepEng] = helpScaler;
			}


			//printing energy
			if (j % printEng == 0) {

				out.open(generatefilename((j / printEng), "TotalEnergy"));



				out << "{";
				for (int f = 0; f <= j / stepEng; f++) {

					out << energyTot[f];
					if (f != j / stepEng) {
						out << ",";
					}
				}
				out << "}";
				out.close();
			}



			if (j%printCyl == 0) {

				out.open(generatefilename((j / printCyl), "Cylider Coordinates"));

				out << "{";
				int helpNum = 0;
				for (int p = 0; p < numrod; p++) {

					nanorod[p].updatepoints();
					if ((nanorod[p].center[0] > (nanorod[p].length / 2.)) && (nanorod[p].center[1] > (nanorod[p].length / 2.)) && (nanorod[p].center[0] < (lx - nanorod[p].length / 2.)) && (nanorod[p].center[1] < (ly - nanorod[p].length / 2.))) {
						helpNum++;
						if (helpNum != 1) {
							out << ",{Yellow, Cylinder[{ {" << nanorod[p].endCenter1X << ", " << nanorod[p].endCenter1Y << ", " << z0 << "}, { " << nanorod[p].endCenter2X << ", " << nanorod[p].endCenter2Y << ", " << z0 << "}}, 145],Opacity[0.2],Cylinder[{ {" << nanorod[p].endCenter1X << ", " << nanorod[p].endCenter1Y << ", " << z0 << "}, { " << nanorod[p].endCenter2X << ", " << nanorod[p].endCenter2Y << ", " << z0 << "}}, 170]}";
						}
						else {
							out << "{Yellow, Cylinder[{ {" << nanorod[p].endCenter1X << ", " << nanorod[p].endCenter1Y << ", " << z0 << "}, { " << nanorod[p].endCenter2X << ", " << nanorod[p].endCenter2Y << ", " << z0 << "}}, 145],Opacity[0.2],Cylinder[{ {" << nanorod[p].endCenter1X << ", " << nanorod[p].endCenter1Y << ", " << z0 << "}, { " << nanorod[p].endCenter2X << ", " << nanorod[p].endCenter2Y << ", " << z0 << "}}, 170]}";
						}


					}

				}

				out << "}";
				out.close();





			}




			if (j > enoughStepOrder) {


				if (j % calOrder == 0) {

					double helpScaler = 0;
					int jjj = 0;
					int nn = 0;

#pragma omp parallel for  private(jjj,nn,counter,helpScaler)
					for (int iii = 0; iii < xmeshn; iii++) {
						for (jjj = 0; jjj < ymeshn; jjj++) {
							oriordprm[iii][jjj][j / calOrder] = 0;
							helpScaler = 0;
							counter = 0.00000000001;
							for (nn = 0; nn < numrod; nn++) {
								if (nanorod[nn].center[0]<(iii*lx / xmeshn) + cutfx && nanorod[nn].center[0]>(iii*lx / xmeshn) - cutfx
									/*&& ReCenter[nn][1][0] > (jjj *ly / ymeshn) - cutfy && ReCenter[nn][1][0] < (jjj*ly / ymeshn) + cutfy*/) {

									helpScaler += (3 * cos(nanorod[nn].theta)*cos(nanorod[nn].theta) - 1) / 2.;//oriordprm[iii][jjj][j / calOrder]
									counter++;

								}
							}
							helpScaler = helpScaler / counter;
							oriordprm[iii][jjj][j / calOrder] = helpScaler;
						}
					}
				}

				if (j% printOrder == 0) {

					out.open(generatefilename((j / printOrder), "Order Parameter"));

					out << "{";


					for (int iii = 0; iii < xmeshn; iii++) {
						for (int jjj = 0; jjj < ymeshn; jjj++) {
							oriordprmAvg[iii][jjj] = 0;
							for (int kkk = 0; kkk < (j / calOrder); kkk++) {
								oriordprmAvg[iii][jjj] += oriordprm[iii][jjj][kkk];
							}
							oriordprmAvg[iii][jjj] = oriordprmAvg[iii][jjj] / double(j / calOrder);


							if (iii != 0) {
								out << ",{" << iii << jjj << "," << oriordprmAvg[iii][jjj] << "}";
							}
							if (iii == 0) {
								out << "{" << iii << "," << oriordprmAvg[iii][jjj] << "}";
							}

						}
					}

					out << "}";
					out.close();
				}



			}


			//calculating total order parameter

			if (j % calTotOrder == 0) {
				totOrder[j / calTotOrder] = 0.;
				double helpScaler = 0.;
				double helpnum = 0.;
#pragma omp parallel for  reduction(+:helpScaler)
				for (int p = 0; p < numrod; p++) {
					if (nanorod[p].center[0]<phb && nanorod[p].center[0]>plb) {
						helpnum++;
						helpScaler += (3. * cos(nanorod[p].theta)*cos(nanorod[p].theta) - 1.) / 2.;
					}
				}
				helpScaler = helpScaler / (double(helpnum));
				totOrder[j / calTotOrder] = helpScaler;
				//finding the avg up untill that point
				totOrderAvg[j / calTotOrder] = 0;
				double helpnum2 = 0.;
				double helpScaler2 = 0;
#pragma omp parallel for reduction(+:helpScaler2) reduction (+:helpnum2)
				for (int t = (j / calTotOrder); t >= (j / (calTotOrder)-50000); t--) {
					if (t >= 0) {
						helpnum2++;
						helpScaler2 += totOrder[t];
					}
				}

				if (helpnum2 != 0) {
					totOrderAvg[j / calTotOrder] = helpScaler2;
					totOrderAvg[j / calTotOrder] = totOrderAvg[j / calTotOrder] / helpnum2;
				}
			}
			//printing total order parameter
			if (j % printTotOrder == 0) {
				out.open(generatefilename((j / printTotOrder), "TotalOrder"));
				out << "{";
				for (int f = 0; f <= j / printTotOrder; f++) {
					out << totOrderAvg[f];
					if (f != j / printTotOrder) {
						out << ",";
					}
				}
				out << "}";
				out.close();
			}

			if (j % 1000 == 0) {
				acceptanceRate = (double(j)) / (double(total));
				//cout << acceptanceRate << "\n";
			}
		


		if (j % printCoor == 0) {
			out.open(generatefilename((j / printCoor), "Coordinates"));



			for (int p = 0; p < numrod; p++) {

				out << nanorod[p].center[0] << "       " << nanorod[p].center[1] << "       " << nanorod[p].theta << "       " << nanorod[p].length << "       " << nanorod[p].coreRadf << "       " << nanorod[p].coatThickf << "       " << nanorod[p].rodNum << "       \n";
			}

			out.close();

		}
	}

	






	int in;
	cin >> in;

	return 0;


}


