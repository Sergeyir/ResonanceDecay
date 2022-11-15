#include <cmath>
#include <vector>

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

#include "../lib/ErrorHandler.h"
#include "../lib/ProgressBar.h"
#include "../lib/Particles.h"

struct
{
	std::vector<std::string> part_name;
	std::vector<double> m1, m2, weight;
} Par;

void AddEntry(std::string, const double, const double, const double = 1, bool = 1);

void RestoreMass(const string, const double, const double, TH2F *, TFile *, double = 1);

double GetMass(const double, const double, const double, const double, const double);

int RestoreRes()
{
	gStyle->SetPalette(kVisibleSpectrum);
	const double m1 = Mass.kaon;
	const double m2 = Mass.kaon;

	std::string channel = "nopid";

	string output_file_name = "../output/" + channel + ".root";

	CheckOutputFile(output_file_name);

	TH2F *mass_distr = new TH2F("mass_distr", "mass_distr", 100, 0, 10, 6000, 0, 12);

	//unindentified detector anomalies --------//
	//-----------------------------------------//

	//long living particles -------------------//
	
	//Any channel can be observed in the detector for that reason
	//even if it is not the channel the decay was

	if (channel == "pipi" || channel == "nopid")
	{
		//strange mesons
		AddEntry("KS", m1, m2);
	}
	else
	{
		AddEntry("KS", m1, m2, 0.1);
	}
	
	if (channel == "kpi" || channel == "pik" || channel == "nopid")
	{
		//charmed mesons
		AddEntry("D0", m1, m2);
		//bottom mesons
		AddEntry("B0_kpi", m1, m2);
		//lambda baryons
		AddEntry("Lambda", m1, m2);
	}
	else
	{
		AddEntry("D0", m1, m2, 0.1);
		AddEntry("B0_kpi", m1, m2, 0.1);
		AddEntry("Lambda", m1, m2, 0.1);
	}
	
	if (channel == "kk" || channel == "nopid")
	{
		//bottom mesons
		AddEntry("B0_kk", m1, m2);
	}
	else
	{
		AddEntry("B0_kk", m1, m2, 0.1);
	}
	if (channel == "kp" || channel == "nopid")
	{
		//bottom mesons
		AddEntry("B0_kp", m1, m2);
	}
	else
	{
		//bottom mesons
		AddEntry("B0_kp", m1, m2, 0.1);
	}
	//-----------------------------------------//

	//pipi channel ----------------------------//
	if (channel == "pipi" || channel == "nopid")
	{
		//light unflavored mesons
		AddEntry("rho770", m1, m2);
		AddEntry("omega782", m1, m2);
		AddEntry("f2_1270", m1, m2);
		AddEntry("f0_1500_pipi", m1, m2);
		AddEntry("f2p_1525_pipi", m1, m2);
		AddEntry("f2_2000", m1, m2);
		AddEntry("rho_2000", m1, m2);
		AddEntry("f0_2060", m1, m2);
		AddEntry("f4_2050_pipi", m1, m2);
		AddEntry("rho3_1690_pipi", m1, m2);
		AddEntry("f2_2240", m1, m2);
		AddEntry("f6_2510", m1, m2);
	}
	//-----------------------------------------//

	//kpi channel -----------------------------//
	if (channel == "kpi" || channel == "pik" || channel == "nopid")
	{
		//strange mesons
		AddEntry("Kstar700", m1, m2);
		AddEntry("Kstar892", m1, m2);
		AddEntry("Kstar1410", m1, m2);
		AddEntry("Kstar0_1430", m1, m2);
		AddEntry("Kstar2_1430", m1, m2);
		AddEntry("Kstar_1680", m1, m2);
		AddEntry("Kstar3_1780", m1, m2);
		AddEntry("Kstar0_1950", m1, m2);
		AddEntry("Kstar4_2045", m1, m2);
		AddEntry("Kstar5_2380", m1, m2);
	}
	//-----------------------------------------//

	//kk channel ------------------------------//
	if (channel == "kk" || channel == "nopid")
	{
		//light unflavored mesons
		AddEntry("phi1020", m1, m2);
		AddEntry("a2_1320", m1, m2);
		AddEntry("a0_1450", m1, m2);
		AddEntry("f0_1500_kk", m1, m2);
		AddEntry("f2p_1525_kk", m1, m2);
		AddEntry("rho3_1690_kk", m1, m2);
		AddEntry("a2_1700", m1, m2);
		AddEntry("f4_2050_kk", m1, m2);
	}
	//-----------------------------------------//

	//ppi channel -----------------------------//
	if (channel == "ppi" || channel == "pip" || channel == "nopid")
	{
		//delta baryons
		AddEntry("Delta1232", m1, m2);
		AddEntry("Delta1600", m1, m2);
		AddEntry("Delta1620", m1, m2);
		AddEntry("Delta1700", m1, m2);
		AddEntry("Delta1900", m1, m2);
		AddEntry("Delta1905", m1, m2);
		AddEntry("Delta1910", m1, m2);
		AddEntry("Delta1920", m1, m2);
		AddEntry("Delta1930", m1, m2);
		AddEntry("Delta1940", m1, m2);
		AddEntry("Delta1950", m1, m2);
		AddEntry("Delta2000", m1, m2);
		AddEntry("Delta2150", m1, m2);
		AddEntry("Delta2200", m1, m2);
		AddEntry("Delta2300", m1, m2);
		AddEntry("Delta2350", m1, m2);
		AddEntry("Delta2390", m1, m2);
		AddEntry("Delta2400", m1, m2);
		AddEntry("Delta2420", m1, m2);
		AddEntry("Delta2750", m1, m2);
		AddEntry("Delta2950", m1, m2);
	}
	//-----------------------------------------//
	
	//pk channel ------------------------------//
	if (channel == "pk" || channel == "kp" || channel == "nopid")
	{
		//N baryons
		AddEntry("N1440", m1, m2);
		AddEntry("N1520", m1, m2);
		AddEntry("N1535", m1, m2);
		AddEntry("N1650", m1, m2);
		AddEntry("N1650", m1, m2);
		AddEntry("N1680", m1, m2);
		AddEntry("N1700", m1, m2);
		AddEntry("N1710", m1, m2);
		AddEntry("N1720", m1, m2);
		AddEntry("N1860", m1, m2);
		AddEntry("N1880", m1, m2);
		AddEntry("N1895", m1, m2);
		AddEntry("N1900", m1, m2);
		AddEntry("N2000", m1, m2);
		AddEntry("N1990", m1, m2);
		AddEntry("N2060", m1, m2);
		AddEntry("N2100", m1, m2);
		AddEntry("N2120", m1, m2);
		AddEntry("N2190", m1, m2);
		AddEntry("N2200", m1, m2);
		AddEntry("N2250", m1, m2);
		AddEntry("N2600", m1, m2);
		AddEntry("N2700", m1, m2);
		
		//lambda baryons
		AddEntry("Lambda1520", m1, m2);
		AddEntry("Lambda1600", m1, m2);
		AddEntry("Lambda1670", m1, m2);
		AddEntry("Lambda1690", m1, m2);
		AddEntry("Lambda1710", m1, m2);
		AddEntry("Lambda1800", m1, m2);
		AddEntry("Lambda1810", m1, m2);
		AddEntry("Lambda1820", m1, m2);
		AddEntry("Lambda1830", m1, m2);
		AddEntry("Lambda1890", m1, m2);
		AddEntry("Lambda2000", m1, m2);
		AddEntry("Lambda2050", m1, m2);
		AddEntry("Lambda2070", m1, m2);
		AddEntry("Lambda2080", m1, m2);
		AddEntry("Lambda2100", m1, m2);
		AddEntry("Lambda2110", m1, m2);
		AddEntry("Lambda2350", m1, m2);
		
		//sigma baryons
		AddEntry("Sigma1620", m1, m2);
		AddEntry("Sigma1660", m1, m2);
		AddEntry("Sigma1670", m1, m2);
		AddEntry("Sigma1750", m1, m2);
		AddEntry("Sigma1775", m1, m2);
		AddEntry("Sigma1780", m1, m2);
		AddEntry("Sigma1880", m1, m2);
		AddEntry("Sigma1900", m1, m2);
		AddEntry("Sigma1910", m1, m2);
		AddEntry("Sigma1915", m1, m2);
		AddEntry("Sigma1940", m1, m2);
		AddEntry("Sigma2010", m1, m2);
		AddEntry("Sigma2030", m1, m2);
		AddEntry("Sigma2100", m1, m2);
		AddEntry("Sigma2110", m1, m2);
		AddEntry("Sigma2230", m1, m2);
		AddEntry("Sigma2250", m1, m2);
	}
	//-----------------------------------------//

	//pp channel ------------------------------//
	if (channel == "pp" || channel == "nopid")
	{
		//c bar mesons
		AddEntry("etac1s", m1, m2);
	}
	//-----------------------------------------//

	//ee channel ------------------------------//
	if (channel == "ee" || channel == "nopid")
	{
		//c bar mesons
		AddEntry("jpsi1s", m1, m2);
		
		//b bar mesons
		AddEntry("upsilon1s", m1, m2);
		AddEntry("upsilon2s", m1, m2);
		AddEntry("upsilon3s", m1, m2);
	}
	//-----------------------------------------//

	TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
	
	for (int count = 0; count < Par.m1.size(); count++)
	{
		cout << OutputColor.bold_green << "[" << count+1 << " out of " << Par.m1.size() << "] " << OutputColor.reset << "Restoring ";
		RestoreMass(Par.part_name[count], Par.m1[count], Par.m2[count], mass_distr, output_file, Par.weight[count]);
	}

	mass_distr->Draw("COL");
	mass_distr->Write();
	
	return 0;
}

void AddEntry(std::string part_name, const double m1, const double m2, const double weight = 1, bool do_antipart = true)
{
	std::string input_file_name = "../data/" + part_name + ".root";
	CheckInputFile(input_file_name);
	
	Par.part_name.push_back(part_name);
	Par.m1.push_back(m1);
	Par.m2.push_back(m2);
	
	Par.weight.push_back(weight);

	//adding antiparticles
	if (do_antipart == true)
	{
		std::string antipart_name = "A" + part_name;
		AddEntry(antipart_name, m2, m1, weight, false);
	}
}

void RestoreMass(std::string part_name, const double m1, const double m2, TH2F *mass_distr, TFile *output_file, double weight = 1)
{
	std::string input_file_name = "../data/" + part_name + ".root";
	TFile *input_file = new TFile(input_file_name.c_str(), "READ");

	TTree *D1 = (TTree*) input_file->Get("D1");
	TTree *D2 = (TTree*) input_file->Get("D2");

	unsigned const int res_num = D1->GetEntries();

	std::cout << res_num << " particles with the name: " << part_name << std::endl;

	double p1[3], p2[3];
	double mom1, mom2;
	double momentum, mass, pT;

	double progress = 0;

	TH1F *stat = (TH1F*) input_file->Get("stat");

	weight *= stat->GetBinContent(5);

	TH2F part_distr = TH2F(part_name.c_str(), part_name.c_str(), 100, 0, 10, 6000, 0, 12);
	
	for (double counter = 0; counter < res_num; counter++)
	{
		ProgressBar.Block((double) (counter + 1.)/res_num);

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

		mass = GetMass(mom1, mom2, momentum, m1, m2);

		double pT = sqrt(pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2));
		if (pT < 0.6) continue;

		mass_distr->Fill(pT, mass, weight);
		part_distr.Fill(pT, mass, weight);
	}
	
	input_file->Close();

	output_file->cd();
	part_distr.Write();

	delete input_file;	
}

double GetMass(const double mom1, const double mom2, const double momentum, const double m1, const double m2)
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
