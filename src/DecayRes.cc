#include <cmath>
#include "TRandom.h"
#include "TFile.h"
#include "TH2.h"
#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../lib/ProgressBar.h"
#include "../func/Tsallis.cc"
#include "TTree.h"
#include "TFile.h"

const double pi = 3.14159265359;

struct
{
	const float pion = 0.134957; 
	const float kaon = 0.493677; 
	const float proton = 0.938272; 
	const float electron = 0.134957; 
} Mass;

void Init(std::string, const float, const float, const float, const float, const int);

int DecayRes()
{
	//light unflavored mesons
	Init("rho770", 0.775, 149E-3, Mass.pion, Mass.pion);
	Init("omega782", 0.782, 8.68E-3, Mass.pion, Mass.pion);
	Init("phi1020", 1.019, 4.29-3, Mass.kaon, Mass.kaon);
	Init("f2_1270", 1.2755, 186.7E-3, Mass.pion, Mass.pion);
	Init("a2_1320", 1.3182, 105E-3, Mass.kaon, Mass.kaon);
	Init("f2_1430", 1.43, 46E-3, Mass.kaon, Mass.kaon);
	Init("a0_1450", 1.474, 265E-3, Mass.kaon, Mass.kaon);
	Init("f0_1500_pipi", 1.506, 112E-3, Mass.pion, Mass.pion);
	Init("f0_1500_kk", 1.506, 112E-3, Mass.kaon, Mass.kaon);
	Init("f2p_1525_pipi", 1.5174, 86E-3, Mass.pion, Mass.pion);
	Init("f2p_1525_kk", 1.5174, 86E-3, Mass.kaon, Mass.kaon);
	
	return 0;
}

void Init(std::string name, const float mean, const float sigma, const float m1, const float m2, const int part_number = 1E5, const int seed = 1)
{
	
	std::string output_file_name = "../data/" + name + ".root";
	CheckOutputFile(output_file_name);
	TFile *output = new TFile(output_file_name.c_str(), "RECREATE");

	std::cout << "Writing data in file " << output_file_name << endl;

	TRandom *rand = new TRandom(seed);
	
	double mass, momentum, energy;
	double theta, phi;
	double mom1, mom2, e1, e2;
	double p1[3], p2[3];

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

	for (int counter = 0; counter < part_number; counter++)
	{
		ProgressBar.Block2((float) (counter+1.)/part_number);

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

	D1->Write();
	D2->Write();
	mass_distr->Write();


}
