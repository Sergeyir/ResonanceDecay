#include <cmath>
#include "TRandom.h"
#include "TFile.h"
#include "TH2.h"
#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../utils/ProgressBar.cc"
#include "../func/Tsallis.cc"
#include "TTree.h"
#include "TFile.h"

int DecayRes() {
	
	//parameters
	double part_number = 1E6;
	const int seed = 1;

	double m1 = Mass.kaon;
	double m2 = Mass.pion;

	double mean = Mass.Kstar;
	double width = Width.Kstar;
	//double mean = 1.169;
	//double width = 0.024;

	string output_file_name = "../data/Resonances/Kstar.root";
	CheckOutputFile(output_file_name);
	TFile *output = new TFile(output_file_name.c_str(), "RECREATE");
	
	const double pi = 3.14159265359;

	TRandom *rand = new TRandom(seed);
	
	double mass, momentum, energy;
	double theta, phi;
	double mom1, mom2, e1, e2;
	double p1[3], p2[3];

	double progress = 0.;

	TTree *D1 = new TTree("D1", "Daughter1");
	TTree *D2 = new TTree("D2", "Daughter2");

	D1->Branch("px", &p1[0]);
	D1->Branch("py", &p1[1]);
	D1->Branch("pz", &p1[2]);

	D2->Branch("px", &p2[0]);
	D2->Branch("py", &p2[1]);
	D2->Branch("pz", &p2[2]);

	TH2D *mass_distr = new TH2D("mass_distr", "mass_distr", 100, 0, 10, 4000, 0, 8);
	mass_distr->SetDefaultSumw2();

	for (double counter = 0; counter < part_number; counter++)
	{
		if (counter/part_number >= progress) {
			ProgressBar(progress);
			progress += 0.01;
		}

		mass = rand->BreitWigner(mean, width);
		if (mass < m1 + m2) continue;
		
		momentum = Tsallis(mass, (unsigned int) counter + seed + mass);

		energy = sqrt(mass*mass + momentum*momentum);

		e1 = (mass*mass + m1*m1 - m2*m2)/(2*mass);
		e2 = (mass*mass + m2*m2 - m1*m1)/(2*mass);

		const double vel = momentum/energy;
		const double gamma = energy/mass;

		theta = rand->Uniform(pi/2);
		phi = rand->Uniform(pi);

		mom1 = sqrt(e1*e1 - m1*m1);
		mom2 = sqrt(e2*e2 - m2*m2);

		p1[0] = mom1*cos(phi)*sin(theta);
		p2[0] = -mom2*cos(phi)*sin(theta);

		p1[1] = mom1*sin(phi)*sin(theta);
		p2[1] = -mom2*sin(phi)*sin(theta);

		p1[2] = mom1*cos(theta);
		p2[2] = -mom2*cos(theta);

		p1[2] = gamma*(p1[2] + vel*e1);
		p2[2] = gamma*(p2[2] + vel*e2);

		mom1 = 0;
		mom2 = 0;
		for (auto p : p1) mom1 += p*p;
		for (auto p : p2) mom2 += p*p;

		e1 = sqrt(m1*m1 + mom1);
		e2 = sqrt(m2*m2 + mom2);

		//computing spherical symmetry angles of a direction of the decayed particle momentum
		phi = rand->Uniform(pi);

		double p1_temp = p1[0];
		double p2_temp = p2[0];

		//transforming the basis by rotating it around z and then x axis
		p1[0] = p1_temp*cos(phi) - p1[1]*sin(phi);
		p2[0] = p2_temp*cos(phi) - p2[1]*sin(phi);

		p1[1] = p1_temp*sin(phi) + p1[1]*cos(phi);
		p2[1] = p2_temp*sin(phi) + p2[1]*cos(phi);

		p1_temp = p1[1];
		p2_temp = p2[1];

		phi = rand->Uniform(pi);

		p1[1] = p1_temp*cos(phi) - p1[2]*sin(phi);
		p2[1] = p2_temp*cos(phi) - p2[2]*sin(phi);

		p1[2] = p1_temp*sin(phi) + p1[2]*cos(phi);
		p2[2] = p2_temp*sin(phi) + p2[2]*cos(phi);

		mom1 = 0;
		mom2 = 0;
		for (auto p : p1) mom1 += p*p;
		for (auto p : p2) mom2 += p*p;

		e1 = sqrt(m1*m1 + mom1);
		e2 = sqrt(m2*m2 + mom2);

		double pT = sqrt(pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2));

		mass_distr->Fill(pT, mass);

		D1->Fill();
		D2->Fill();
	}

	ProgressBar(1);
	cout << endl;

	D1->Write();
	D2->Write();
	mass_distr->Write();

	return 0;
}
