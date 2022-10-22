#include <cmath>
#include <vector>
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../lib/ProgressBar.h"

struct
{
	std::vector<std::string> part_name;
	std::vector<float> m1, m2, weight;
} Par;

void AddEntry(std::string, float, float, float, bool = true);

void Init(TH2F *);
	
void RestoreMass(const string, const float, const float, TH2F *, const float = 1);

double getMass(const float, const float, const double, const double, const double);

int RestoreRes()
{
	const double m1 = Mass.proton;
	const double m2 = Mass.kaon;

	string output_file_name = "../output/PK.root";

	CheckOutputFile(output_file_name);

	TH2F *mass_distr = new TH2F("mass_distr", "mass_distr", 100, 0, 10, 6000, 0, 12);

	AddEntry("rho770", m1, m2, 1, 0);
	AddEntry("omega782", m1, m2, 0.015, 0);
	AddEntry("phi1020", m1, m2, 1, 0);
	AddEntry("f2_1270", m1, m2, 0.84, 0);
	AddEntry("a2_1320", m1, m2, 0.049, 0);
	AddEntry("a0_1450", m1, m2, 1, 0);
	AddEntry("f0_1500_pipi", m1, m2, 0.08, 0);
	AddEntry("f0_1500_kk", m1, m2, 0.345, 0);
	AddEntry("f2p_1525_pipi", m1, m2, 0.083, 0);
	AddEntry("f2p_1525_kk", m1, m2, 0.876, 0);
	AddEntry("rho3_1690_pipi", m1, m2, 0.236, 0);
	AddEntry("rho3_1690_kk", m1, m2, 0.0158, 0);
	AddEntry("a2_1700", m1, m2, 0.019, 0);
	AddEntry("f2_2000", m1, m2, 1, 0);
	AddEntry("rho_2000", m1, m2, 1, 0);
	AddEntry("f0_2060", m1, m2, 1, 0);
	AddEntry("f4_2050_pipi", m1, m2, 0.17, 0);
	AddEntry("f4_2050_kk", m1, m2, 0.068, 0);
	AddEntry("f2_2240", m1, m2, 1, 0);
	AddEntry("f6_2510", m1, m2, 0.06, 0);
	
	AddEntry("KS", m1, m2, 0.692, 0);
	AddEntry("Kstar700", m1, m2, 1);
	AddEntry("Kstar892", m1, m2, 1);
	AddEntry("Kstar1410", m1, m2, 0.066);
	AddEntry("Kstar0_1430", m1, m2, 0.93);
	AddEntry("Kstar2_1430", m1, m2, 0.499);
	AddEntry("Kstar_1680", m1, m2, 0.387);
	AddEntry("Kstar3_1780", m1, m2, 0.31);
	AddEntry("Kstar0_1950", m1, m2, 0.31);
	AddEntry("Kstar4_2045", m1, m2, 0.099);
	AddEntry("Kstar5_2380", m1, m2, 0.061);

	//charmed mesons
	AddEntry("D0", m1, m2, 0.003947);

	//bottom mesons
	//weight values for B0 are fiducial
	AddEntry("B0_kpi", m1, m2, 0.05);
	AddEntry("B0_kk", m1, m2, 0.5, 0);
	AddEntry("B0_kp", m1, m2, 0.05);

	//c bar mesons
	AddEntry("etac1s", m1, m2, 0.0144E-3, 0);
	AddEntry("jpsi1s", m1, m2, 0.0597, 0);

	//b bar mesons
	AddEntry("upsilon1s", m1, m2, 0.0238E-3, 0);
	AddEntry("upsilon2s", m1, m2, 0.0191E-3, 0);
	AddEntry("upsilon3s", m1, m2, 0.0191E-3, 0);

	//N baryons
	AddEntry("N1440", m1, m2, 0.15);
	AddEntry("N1520", m1, m2, 0.15);
	AddEntry("N1535", m1, m2, 0.1);
	AddEntry("N1650", m1, m2, 0.15);
	AddEntry("N1650", m1, m2, 0.1);
	AddEntry("N1680", m1, m2, 0.16);
	AddEntry("N1700", m1, m2, 0.06);
	AddEntry("N1710", m1, m2, 0.063);
	AddEntry("N1720", m1, m2, 0.06);
	AddEntry("N1860", m1, m2, 0.06);
	AddEntry("N1880", m1, m2, 0.07);
	AddEntry("N1895", m1, m2, 0.05);
	AddEntry("N1900", m1, m2, 0.05);
	AddEntry("N2000", m1, m2, 0.025);
	AddEntry("N1990", m1, m2, 0.015);
	AddEntry("N2060", m1, m2, 0.05);
	AddEntry("N2100", m1, m2, 0.1);
	AddEntry("N2120", m1, m2, 0.05);
	AddEntry("N2190", m1, m2, 0.075);
	AddEntry("N2200", m1, m2, 0.0113);
	AddEntry("N2250", m1, m2, 0.05);
	AddEntry("N2600", m1, m2, 0.02);
	AddEntry("N2700", m1, m2, 0.03);

	//delta baryons
	AddEntry("Delta1232", m1, m2, 0.5);
	AddEntry("Delta1600", m1, m2, 0.08);
	AddEntry("Delta1620", m1, m2, 0.15);
	AddEntry("Delta1700", m1, m2, 0.0375);
	AddEntry("Delta1900", m1, m2, 0.04);
	AddEntry("Delta1905", m1, m2, 0.07);
	AddEntry("Delta1910", m1, m2, 0.1);
	AddEntry("Delta1920", m1, m2, 0.0625);
	AddEntry("Delta1930", m1, m2, 0.05);
	AddEntry("Delta1940", m1, m2, 0.05);
	AddEntry("Delta1950", m1, m2, 0.1);
	AddEntry("Delta2000", m1, m2, 0.035);
	AddEntry("Delta2150", m1, m2, 0.04);
	AddEntry("Delta2200", m1, m2, 0.025);
	AddEntry("Delta2300", m1, m2, 0.0225);
	AddEntry("Delta2350", m1, m2, 0.085);
	AddEntry("Delta2390", m1, m2, 0.0375);
	AddEntry("Delta2400", m1, m2, 0.04);
	AddEntry("Delta2420", m1, m2, 0.0375);
	AddEntry("Delta2750", m1, m2, 0.02);
	AddEntry("Delta2950", m1, m2, 0.02);

	//lambda baryons
	AddEntry("Lambda", m1, m2, 0.639);
	AddEntry("Lambda1520", m1, m2, 0.225);
	AddEntry("Lambda1600", m1, m2, 0.1125);
	AddEntry("Lambda1670", m1, m2, 0.125);
	AddEntry("Lambda1690", m1, m2, 0.125);
	AddEntry("Lambda1710", m1, m2, 0.115);
	AddEntry("Lambda1800", m1, m2, 0.1875);
	AddEntry("Lambda1810", m1, m2, 0.1);
	AddEntry("Lambda1820", m1, m2, 0.15);
	AddEntry("Lambda1890", m1, m2, 0.15);
	AddEntry("Lambda2000", m1, m2, 0.0675);
	AddEntry("Lambda2050", m1, m2, 0.0475);
	AddEntry("Lambda2070", m1, m2, 0.03);
	AddEntry("Lambda2080", m1, m2, 0.02);
	AddEntry("Lambda2100", m1, m2, 0.075);
	AddEntry("Lambda2110", m1, m2, 0.075);
	AddEntry("Lambda2350", m1, m2, 0.03);

	//sigma baryons
	AddEntry("Sigma1620", m1, m2, 0.175);
	AddEntry("Sigma1660", m1, m2, 0.025);
	AddEntry("Sigma1670", m1, m2, 0.02);
	AddEntry("Sigma1750", m1, m2, 0.045);
	AddEntry("Sigma1775", m1, m2, 0.1);
	AddEntry("Sigma1780", m1, m2, 0.01);
	AddEntry("Sigma1880", m1, m2, 0.1);
	AddEntry("Sigma1900", m1, m2, 0.225);
	AddEntry("Sigma1910", m1, m2, 0.0125);
	AddEntry("Sigma1915", m1, m2, 0.0125);
	AddEntry("Sigma1940", m1, m2, 0.065);
	AddEntry("Sigma2010", m1, m2, 0.035);
	AddEntry("Sigma2030", m1, m2, 0.1);
	AddEntry("Sigma2100", m1, m2, 0.04);
	AddEntry("Sigma2110", m1, m2, 0.065);
	AddEntry("Sigma2230", m1, m2, 0.03);
	AddEntry("Sigma2250", m1, m2, 0.03);

	Init(mass_distr);

	mass_distr->Draw();
	
	TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
	mass_distr->Write();
	
	return 0;
}

void AddEntry(std::string part_name, const float m1, const float m2, const float channel_fract, bool do_antipart = true)
{
	std::string input_file_name = "../data/" + part_name + ".root";
	CheckInputFile(input_file_name);
	
	Par.part_name.push_back(part_name);
	Par.m1.push_back(m1);
	Par.m2.push_back(m2);
	

	//adding antiparticles
	if (do_antipart == true && m1 != m2)
	{
		Par.weight.push_back(channel_fract);
		std::string antipart_name = "anti" + part_name;
		AddEntry(antipart_name, m1, m2, channel_fract, false);
	}
	else
	{
		Par.weight.push_back(channel_fract*2.);
	}
}

void Init(TH2F *mass_distr)
{
	for (int count = 0; count < Par.m1.size(); count++)
	{
		RestoreMass(Par.part_name[count], Par.m1[count], Par.m2[count], mass_distr);
	}
}

void RestoreMass(std::string part_name, const float m1, const float m2, TH2F *mass_distr, const float weight = 1)
{
	cout << "Restoring " << part_name << endl;

	std::string input_file_name = "../data/" + part_name + ".root";
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
		ProgressBar.Block((float) (counter + 1.)/res_num);

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

		//if (checkEnergy(p1, p2, mass, m1, m2)) continue;

		mass_distr->Fill(pT, mass, weight*2);
	}

	delete input_file, D1, D2;
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
