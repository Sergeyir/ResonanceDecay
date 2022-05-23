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
bool checkMom(double *, double *, const double);

int RestoreRes()
{
	const double m1 = Mass.proton;
	const double m2 = Mass.pion;

	string output_file_name = "../output/Resonances/Lambda1520.root";

	CheckOutputFile(output_file_name);

	TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
	TH2D *mass_distr = new TH2D("mass_distr", "mass_distr", 100, 0, 10, 4000, 0, 8);

	//RestoreMass("../data/Resonances/Kstar.root", m1, m2, mass_distr, 1);
	//RestoreMass("../data/Resonances/UnindentPiK1169.root", m1, m2, mass_distr, 100);
	RestoreMass("../data/Resonances/Lambda1520.root", m1, m2, mass_distr, 1);
	//RestoreMass("../data/Resonances/Lambda.root", m1, m2, mass_distr, 100);

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

		if (!checkMom(p1, p2, mass)) continue;

		double pT = sqrt(pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2));

		mass_distr->Fill(pT, mass, weight);
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

bool checkMom(double *p1, double *p2, const double mass)
{
	double mom = 0, energy, vel, phi, theta;
	double mom1 = 0, mom2 = 0;
	const double pi = 3.14159265359;
	double p1_tmp[3], p2_tmp[3];

	for (int count = 0; count < 3; count++)
	{
		double p;

		p = p1[count]+p2[count];

		p1_tmp[count] = p1[count];
		p2_tmp[count] = p2[count];

		mom += p*p;

		mom1 += p1[count]*p1[count];
		mom2 += p2[count]*p2[count];
	}
	
	mom1 = sqrt(mom1);
	mom2 = sqrt(mom2);

	phi = atan((p1[1]+p2[1])/(p1[0]+p2[0]));
	theta = acos((p1[2]+p2[2])/sqrt(mom));

	energy = sqrt(mom + mass*mass);

	vel = sqrt(mom)/energy;

	p1[1] = p1_tmp[1]*cos(-theta) + p1_tmp[2]*sin(-theta);
	p2[1] = p2_tmp[1]*cos(-theta) + p2_tmp[2]*sin(-theta);

	p1[2] = p1_tmp[1]*sin(-theta) + p1_tmp[2]*cos(-theta);
	p2[2] = p2_tmp[1]*sin(-theta) + p2_tmp[2]*cos(-theta);

	p1_tmp[1] = p1[1];
	p2_tmp[1] = p2[1];

	p1[0] = p1_tmp[0]*cos(-phi) - p1_tmp[1]*sin(-phi);
	p2[0] = p2_tmp[0]*cos(-phi) - p2_tmp[1]*sin(-phi);

	p1[1] = p1_tmp[0]*sin(-phi) + p1_tmp[1]*cos(-phi);
	p2[1] = p2_tmp[0]*sin(-phi) + p2_tmp[1]*cos(-phi);

	p1[2] = p1[2]-vel*energy;
	p2[2] = p2[2]-vel*energy;

	mom = pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2) + pow(p1[2]+p2[2], 2);

	cout << mom << endl;

	if (mom < 0.1) return true;

	return false;
}
