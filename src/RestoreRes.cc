#include <cmath>
#include <vector>
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../utils/ProgressBar.cc"

void RestoreMass(const string, const double, const double, TH2D *, const int);
double getMass(const double, const double, const double, const double, const double);
bool checkEnergy(double *, double *, const double, const double, const double);

int RestoreRes()
{
	const double m1 = Mass.proton;
	const double m2 = Mass.kaon;

	string output_file_name = "../output/Resonances/Lambda1520.root";

	CheckOutputFile(output_file_name);

	TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
	TH2D *mass_distr = new TH2D("mass_distr", "mass_distr", 100, 0, 10, 4000, 0, 8);

	//RestoreMass("../data/Resonances/Kstar.root", m1, m2, mass_distr, 1000);
	RestoreMass("../data/Resonances/UnindentPPi1169.root", m1, m2, mass_distr, 100);
	RestoreMass("../data/Resonances/Lambda1520.root", m1, m2, mass_distr, 10);
	//RestoreMass("../data/Resonances/Lambda.root", m1, m2, mass_distr, 1000);

	mass_distr->Draw();
	mass_distr->Write();
	
	return 0;
}

void RestoreMass(const string input_file_name, const double m1, const double m2, TH2D *mass_distr, const int weight = 1)
{
	cout << "Restoring from file " << input_file_name << endl;

	CheckInputFile(input_file_name);
	TFile *input_file = new TFile(input_file_name.c_str(), "READ");

	TTree *D1 = (TTree*) input_file->Get("D1");
	TTree *D2 = (TTree*) input_file->Get("D2");

	double p1[3], p2[3];
	double mom1, mom2;
	double momentum, mass, pT;

	double res_num = D1->GetEntries();

	float progress = 0;
	for (double counter = 0; counter < res_num; counter++)
	{
		if (counter/res_num >= progress)
		{
	
			ProgressBar(progress);
			progress+=0.01;

		}

		D1->GetEntry(static_cast<unsigned long>(counter));
		D2->GetEntry(static_cast<unsigned long>(counter));

		p1[0] = D1->GetLeaf("px")->GetValue();	
		p1[1] = D1->GetLeaf("py")->GetValue();	
		p1[2] = D1->GetLeaf("pz")->GetValue();	

		p2[0] = D2->GetLeaf("px")->GetValue();	
		p2[1] = D2->GetLeaf("py")->GetValue();	
		p2[2] = D2->GetLeaf("pz")->GetValue();	

		//momentum for the next 3 lines are squared

		momentum = pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2) + pow(p1[2]+p2[2], 2);

		mom1 = pow(p1[0], 2) + pow(p1[1], 2) + pow(p1[2], 2);
		mom2 = pow(p2[0], 2) + pow(p2[1], 2) + pow(p2[2], 2);

		mass = getMass(mom1, mom2, momentum, m1, m2);

		double pT = sqrt(pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2));
		if (pT < 0.6) continue;

		if (!checkEnergy(p1, p2, mass, m1, m2)) continue;

		mass_distr->Fill(pT, mass, weight*2);
	}

	ProgressBar(1);

	cout << endl;

	delete D1, D2;
}

double getMass(const double mom1, const double mom2, const double momentum, const double m1, const double m2)
{
	double e1 = sqrt(mom1 + m1*m1);
	double e2 = sqrt(mom2 + m2*m2);

	double mass = sqrt(pow(e1+e2, 2) - momentum);

	return mass;
}

bool checkEnergy(double *p1, double *p2, const double mass, const double m1, const double m2)
{
	double energy, e1, e2, vel, gamma, phi;
	double mom = 0, mom1 = 0, mom2 = 0;
	const double pi = 3.14159265359;
	double p1_tmp[3], p2_tmp[3];

	for (int count = 0; count < 3; count++)
	{
		double p;

		p = p1[count]+p2[count];

		p1_tmp[count] = p1[count];
		p2_tmp[count] = p2[count];

		mom1 += p1[count]*p1[count];
		mom2 += p2[count]*p2[count];

		mom += p*p;
	}
	
	e1 = sqrt(mom1 + m1*m1);
	e2 = sqrt(mom2 + m2*m2);
	energy = e1 + e2;

	vel = sqrt(mom)/energy;
	gamma = energy/mass;

	//transforming the cortesian basis for momentum to be p=px
	phi = -atan((p1[1]+p2[1])/(p1[0]+p2[0]));

	p1[0] = p1_tmp[0]*cos(phi) - p1_tmp[1]*sin(phi);
	p2[0] = p2_tmp[0]*cos(phi) - p2_tmp[1]*sin(phi);

	p1[1] = p1_tmp[0]*sin(phi) + p1_tmp[1]*cos(phi);
	p2[1] = p2_tmp[0]*sin(phi) + p2_tmp[1]*cos(phi);

	p1_tmp[0] = p1[0];
	p2_tmp[0] = p2[0];

	phi = -atan((p1[2]+p2[2])/(p1[0]+p2[0]));

	p1[0] = p1_tmp[0]*cos(phi) - p1_tmp[2]*sin(phi);
	p2[0] = p2_tmp[0]*cos(phi) - p2_tmp[2]*sin(phi);

	p1[2] = p1_tmp[0]*sin(phi) + p1_tmp[2]*cos(phi);
	p2[2] = p2_tmp[0]*sin(phi) + p2_tmp[2]*cos(phi);

	e1 = gamma*(e1 - vel*abs(p1[0]));
	e2 = gamma*(e2 - vel*abs(p2[0]));

	energy = e1 + e2;

	double e1_true = (mass*mass + m1*m1 - m2*m2)/(2*mass);
	double e2_true = (mass*mass + m2*m2 - m1*m1)/(2*mass);

	if (abs(e1 - e1_true) < 0.001 && abs (e2 - e2_true) < 0.001) {
		//cout << mass << " " << energy - mass << " " << e1 - e1_true << " " << e2 - e2_true << endl;
		return true;
	}

	return false;
}
