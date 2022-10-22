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

void Init(std::string, const float, const float, const float, const float, const int, const int = 1, bool = true);

int DecayRes()
{
	const int part_number = 1e4;
	//light unflavored mesons
	Init("rho770", 0.775, 149E-3, Mass.pion, Mass.pion, part_number);
	Init("omega782", 0.782, 8.68E-3, Mass.pion, Mass.pion, part_number);
	Init("phi1020", 1.019, 4.29-3, Mass.kaon, Mass.kaon, part_number);
	Init("f2_1270", 1.2755, 186.7E-3, Mass.pion, Mass.pion, part_number);
	Init("a2_1320", 1.3182, 105E-3, Mass.kaon, Mass.kaon, part_number);
	Init("f2_1430", 1.43, 46E-3, Mass.kaon, Mass.kaon, part_number);
	Init("a0_1450", 1.474, 265E-3, Mass.kaon, Mass.kaon, part_number);
	Init("f0_1500_pipi", 1.506, 112E-3, Mass.pion, Mass.pion, part_number);
	Init("f0_1500_kk", 1.506, 112E-3, Mass.kaon, Mass.kaon, part_number);
	Init("f2p_1525_pipi", 1.5174, 86E-3, Mass.pion, Mass.pion, part_number);
	Init("f2p_1525_kk", 1.5174, 86E-3, Mass.kaon, Mass.kaon, part_number);
	Init("rho3_1690_pipi", 1.688, 161E-3, Mass.pion, Mass.pion, part_number);
	Init("rho3_1690_kk", 1.688, 161E-3, Mass.kaon, Mass.kaon, part_number);
	Init("a2_1700", 1.698, 265E-3, Mass.kaon, Mass.kaon, part_number);
	Init("f2_2000", 1.996, 312E-3, Mass.pion, Mass.pion, part_number);
	Init("rho_2000", 2., 244E-3, Mass.pion, Mass.pion, part_number);
	Init("f0_2060", 2.050, 80E-3, Mass.pion, Mass.pion, part_number);
	Init("f4_2050_pipi", 2.018, 237E-3, Mass.pion, Mass.pion, part_number);
	Init("f4_2050_kk", 2.018, 237E-3, Mass.kaon, Mass.kaon, part_number);
	Init("f2_2240", 2.240, 241E-3, Mass.pion, Mass.pion, part_number);
	Init("f6_2510", 2.465, 255E-3, Mass.pion, Mass.pion, part_number);

	//strange mesons
	Init("KS", 0.8954, 2E-3, Mass.pion, Mass.pion, part_number);
	Init("Kstar700", 0.845, 486E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar892", 0.892, 51.4E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar1410", 1.414, 232E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar0_1430", 1.425, 270E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar2_1430", 1.427, 100E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar_1680", 1.718, 322E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar3_1780", 1.779, 161E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar0_1950", 1.944, 100E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar4_2045", 2.048, 199E-3, Mass.kaon, Mass.pion, part_number);
	Init("Kstar5_2380", 2.382, 178E-3, Mass.kaon, Mass.pion, part_number);

	//charmed mesons
	Init("D0", 1.865, 2E-3, Mass.kaon, Mass.pion, part_number);

	//bottom mesons
	Init("B0_kpi", 5.279, 2E-3, Mass.kaon, Mass.pion, part_number);
	Init("B0_kk", 5.279, 2E-3, Mass.kaon, Mass.kaon, part_number);
	Init("B0_kp", 5.279, 2E-3, Mass.kaon, Mass.proton, part_number);

	//c bar mesons
	Init("etac1s", 2.983, 32E-3, Mass.proton, Mass.proton, part_number);
	Init("jpsi1s", 3.096, 92.6E-3, Mass.electron, Mass.electron, part_number);

	//b bar mesons
	Init("upsilon1s", 9.4603, 54.02E-3, Mass.electron, Mass.electron, part_number);
	Init("upsilon2s", 10.023, 331.5E-3, Mass.electron, Mass.electron, part_number);
	Init("upsilon3s", 10.335, 331.5E-3, Mass.electron, Mass.electron, part_number);

	//N baryons
	Init("N1440", 1.44, 325E-3, Mass.proton, Mass.pion, part_number);
	Init("N1520", 1.515, 110E-3, Mass.proton, Mass.pion, part_number);
	Init("N1535", 1.52, 150E-3, Mass.proton, Mass.pion, part_number);
	Init("N1650", 1.65, 125E-3, Mass.proton, Mass.pion, part_number);
	Init("N1675", 1.67, 150E-3, Mass.proton, Mass.pion, part_number);
	Init("N1680", 1.68, 122.5E-3, Mass.proton, Mass.pion, part_number);
	Init("N1700", 1.725, 200E-3, Mass.proton, Mass.pion, part_number);
	Init("N1710", 1.72, 140E-3, Mass.proton, Mass.pion, part_number);
	Init("N1720", 1.725, 275E-3, Mass.proton, Mass.pion, part_number);
	Init("N1860", 1.928, 376E-3, Mass.proton, Mass.pion, part_number);
	Init("N1880", 1.880, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("N1895", 1.885, 190E-3, Mass.proton, Mass.pion, part_number);
	Init("N1900", 1.920, 220E-3, Mass.proton, Mass.pion, part_number);
	Init("N2000", 2., 300E-3, Mass.proton, Mass.pion, part_number);
	Init("N1990", 2.025, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("N2060", 2.115, 375E-3, Mass.proton, Mass.pion, part_number);
	Init("N2100", 2.1, 260E-3, Mass.proton, Mass.pion, part_number);
	Init("N2120", 2.11, 310E-3, Mass.proton, Mass.pion, part_number);
	Init("N2190", 2.18, 400E-3, Mass.proton, Mass.pion, part_number);
	Init("N2200", 2.25, 375E-3, Mass.proton, Mass.pion, part_number);
	Init("N2250", 2.285, 450E-3, Mass.proton, Mass.pion, part_number);
	Init("N2600", 2.55, 650E-3, Mass.proton, Mass.pion, part_number);
	Init("N2700", 2.612, 350E-3, Mass.proton, Mass.pion, part_number);

	//delta baryon
	Init("Delta1232", 1.2313, 117E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1600", 1.57, 250E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1620", 1.59, 130E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1700", 1.720, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1900", 1.88, 250E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1905", 1.855, 335E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1910", 1.9, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1920", 1.92, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1930", 1.95, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1940", 2., 400E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta1950", 1.9275, 285E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2000", 2.015, 500E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2150", 2.15, 200E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2200", 2.2, 350E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2300", 2.3, 360E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2350", 2.35, 350E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2390", 2.375, 300E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2400", 2.45, 500E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2420", 2.5, 500E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2750", 2.8, 350E-3, Mass.proton, Mass.pion, part_number);
	Init("Delta2950", 2.99, 330E-3, Mass.proton, Mass.pion, part_number);

	
	//lambda baryons
	Init("Lambda", 1.115, 2E-3, Mass.proton, Mass.pion, part_number);
	Init("Lambda1520", 1.5195, 15.6E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1600", 1.6, 200E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1670", 1.67, 30E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1690", 1.685, 70E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1710", 1.713, 180E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1800", 1.8, 200E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1810", 1.79, 110E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1820", 1.82, 80E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1830", 1.825, 80E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda1890", 1.89, 120E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2000", 2., 275E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2050", 2.056, 493E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2070", 2.070, 370E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2080", 2.082, 181E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2100", 2.1, 175E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2110", 2.09, 250E-3, Mass.proton, Mass.kaon, part_number);
	Init("Lambda2350", 2.355, 175E-3, Mass.proton, Mass.kaon, part_number);

	//sigma baryons
	Init("Sigma1620", 1.625, 70E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1660", 1.66, 150E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1670", 1.675, 70E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1750", 1.75, 150E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1775", 1.775, 120E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1780", 1.78, 150E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1880", 1.880, 150E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1900", 1.925, 165E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1910", 1.910, 225E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1915", 1.918, 120E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma1940", 1.940, 300E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2010", 2.005, 178E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2030", 2.038, 225E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2100", 2.1, 260E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2110", 2.105, 313E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2230", 2.24, 345E-3, Mass.proton, Mass.kaon, part_number);
	Init("Sigma2250", 2.245, 105E-3, Mass.proton, Mass.kaon, part_number);
	
	return 0;
}

void Init(std::string part_name, const float mean, const float sigma, const float m1, const float m2, const int part_number, const int seed = 1, bool do_antipart = true)
{
	std::string output_file_name = "../data/" + part_name + ".root";
	CheckOutputFile(output_file_name);
	TFile *output = new TFile(output_file_name.c_str(), "RECREATE");
	std::cout OutputColor.BoldGren << "Info: " << OutputColor.reset << "Generating " << part_number << " particles with the name: " << part_name << endl;

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

		mass = rand->BreitWigner(mean, sigma);
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

	if (do_antipart == true && m1 != m2)
	{
		std::string antipart_name = "anti" + part_name;
		Init(antipart_name, mean, sigma, m2, m1, part_number, rand->Integer(10000), false);
	}

	delete rand;
	delete D1;
	delete D2;
	delete mass_distr;
	delete output;
}
