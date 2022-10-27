#include <cmath>

#include "TRandom.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

#include "../lib/Particles.h"
#include "../lib/ErrorHandler.h"
#include "../lib/ProgressBar.h"
#include "../func/Tsallis.h"
#include "../lib/Time.h"

const double ptmin = 0.3;
const double ptmax = 8.;
const double pi = 3.14159265359;

struct
{
	std::vector<std::string> name;
	std::vector<double> mean, sigma, m1, m2, mode;
	std::vector<int> number, seed;
	std::vector<std::string> channel;
} Part;

void PrintSeparator(std::string left_edge = "/", std::string body = "-", std::string right_edge = "/")
{
		std::cout << left_edge;
		for (int i = 0; i < 84; i++) std::cout << body;
		std::cout << right_edge << std::endl;
}

void AddEntry(std::string, const double, const double, const double, const double, const int, const double = 1, const int = 1, bool = true);

void Init(std::string, const double, const double, const double, const double, const int, const double, const int);

int DecayRes()
{
	const int part_number = 1e5;

	//light unflavored mesons
	AddEntry("rho770", 0.775, 149E-3, Mass.pion, Mass.pion, part_number, 1);
	AddEntry("omega782", 0.782, 8.68E-3, Mass.pion, Mass.pion, part_number, 0.015);
	AddEntry("phi1020", 1.019, 4.29E-3, Mass.kaon, Mass.kaon, part_number, 1);
	AddEntry("f2_1270", 1.2755, 186.7E-3, Mass.pion, Mass.pion, part_number, 0.84);
	AddEntry("a2_1320", 1.3182, 105E-3, Mass.kaon, Mass.kaon, part_number, 0.049);
	AddEntry("a0_1450", 1.474, 265E-3, Mass.kaon, Mass.kaon, part_number, 1);
	AddEntry("f0_1500_pipi", 1.506, 112E-3, Mass.pion, Mass.pion, part_number, 0.08);
	AddEntry("f0_1500_kk", 1.506, 112E-3, Mass.kaon, Mass.kaon, part_number, 0.345);
	AddEntry("f2p_1525_pipi", 1.5174, 86E-3, Mass.pion, Mass.pion, part_number, 0.083);
	AddEntry("f2p_1525_kk", 1.5174, 86E-3, Mass.kaon, Mass.kaon, part_number, 0.876);
	AddEntry("rho3_1690_pipi", 1.688, 161E-3, Mass.pion, Mass.pion, part_number, 0.236);
	AddEntry("rho3_1690_kk", 1.688, 161E-3, Mass.kaon, Mass.kaon, part_number, 0.0158);
	AddEntry("a2_1700", 1.698, 265E-3, Mass.kaon, Mass.kaon, part_number, 0.019);
	AddEntry("f2_2000", 1.996, 312E-3, Mass.pion, Mass.pion, part_number, 1);
	AddEntry("rho_2000", 2., 244E-3, Mass.pion, Mass.pion, part_number, 1);
	AddEntry("f0_2060", 2.050, 80E-3, Mass.pion, Mass.pion, part_number, 1);
	AddEntry("f4_2050_pipi", 2.018, 237E-3, Mass.pion, Mass.pion, part_number, 0.17);
	AddEntry("f4_2050_kk", 2.018, 237E-3, Mass.kaon, Mass.kaon, part_number, 0.068);
	AddEntry("f2_2240", 2.240, 241E-3, Mass.pion, Mass.pion, part_number, 1);
	AddEntry("f6_2510", 2.465, 255E-3, Mass.pion, Mass.pion, part_number, 0.06);

	//strange mesons
	AddEntry("KS", 0.8954, 2E-3, Mass.pion, Mass.pion, part_number, 0.692);
	AddEntry("Kstar700", 0.845, 486E-3, Mass.kaon, Mass.pion, part_number, 1);
	AddEntry("Kstar892", 0.892, 51.4E-3, Mass.kaon, Mass.pion, part_number, 1);
	AddEntry("Kstar1410", 1.414, 232E-3, Mass.kaon, Mass.pion, part_number, 0.066);
	AddEntry("Kstar0_1430", 1.425, 270E-3, Mass.kaon, Mass.pion, part_number, 0.93);
	AddEntry("Kstar2_1430", 1.427, 100E-3, Mass.kaon, Mass.pion, part_number, 0.499);
	AddEntry("Kstar_1680", 1.718, 322E-3, Mass.kaon, Mass.pion, part_number, 0.387);
	AddEntry("Kstar3_1780", 1.779, 161E-3, Mass.kaon, Mass.pion, part_number, 0.31);
	AddEntry("Kstar0_1950", 1.944, 100E-3, Mass.kaon, Mass.pion, part_number, 0.31);
	AddEntry("Kstar4_2045", 2.048, 199E-3, Mass.kaon, Mass.pion, part_number, 0.099);
	AddEntry("Kstar5_2380", 2.382, 178E-3, Mass.kaon, Mass.pion, part_number, 0.061);

	//charmed mesons
	AddEntry("D0", 1.865, 2E-3, Mass.kaon, Mass.pion, part_number, 0.003947);

	//bottom mesons
	AddEntry("B0_kpi", 5.279, 2E-3, Mass.kaon, Mass.pion, part_number, 0.05);
	AddEntry("B0_kk", 5.279, 2E-3, Mass.kaon, Mass.kaon, part_number, 0.5);
	AddEntry("B0_kp", 5.279, 2E-3, Mass.kaon, Mass.proton, part_number, 0.1);

	//c bar mesons
	AddEntry("etac1s", 2.983, 32E-3, Mass.proton, Mass.proton, part_number, 0.0144E-3);
	AddEntry("jpsi1s", 3.096, 92.6E-3, Mass.electron, Mass.electron, part_number, 0.0597);

	//b bar mesons
	AddEntry("upsilon1s", 9.4603, 54.02E-3, Mass.electron, Mass.electron, part_number, 0.0238E-3);
	AddEntry("upsilon2s", 10.023, 331.5E-3, Mass.electron, Mass.electron, part_number, 0.0191E-3);
	AddEntry("upsilon3s", 10.335, 331.5E-3, Mass.electron, Mass.electron, part_number, 0.0191E-3);

	//N baryons
	AddEntry("N1440", 1.44, 325E-3, Mass.proton, Mass.pion, part_number, 0.15);
	AddEntry("N1520", 1.515, 110E-3, Mass.proton, Mass.pion, part_number, 0.15);
	AddEntry("N1535", 1.52, 150E-3, Mass.proton, Mass.pion, part_number, 0.1);
	AddEntry("N1650", 1.65, 125E-3, Mass.proton, Mass.pion, part_number, 0.15);
	AddEntry("N1675", 1.67, 150E-3, Mass.proton, Mass.pion, part_number, 0.1);
	AddEntry("N1680", 1.68, 122.5E-3, Mass.proton, Mass.pion, part_number, 0.16);
	AddEntry("N1700", 1.725, 200E-3, Mass.proton, Mass.pion, part_number, 0.06);
	AddEntry("N1710", 1.72, 140E-3, Mass.proton, Mass.pion, part_number, 0.063);
	AddEntry("N1720", 1.725, 275E-3, Mass.proton, Mass.pion, part_number, 0.06);
	AddEntry("N1860", 1.928, 376E-3, Mass.proton, Mass.pion, part_number, 0.06);
	AddEntry("N1880", 1.880, 300E-3, Mass.proton, Mass.pion, part_number, 0.07);
	AddEntry("N1895", 1.885, 190E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("N1900", 1.920, 220E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("N2000", 2., 300E-3, Mass.proton, Mass.pion, part_number, 0.025);
	AddEntry("N1990", 2.025, 300E-3, Mass.proton, Mass.pion, part_number, 0.015);
	AddEntry("N2060", 2.115, 375E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("N2100", 2.1, 260E-3, Mass.proton, Mass.pion, part_number, 0.1);
	AddEntry("N2120", 2.11, 310E-3, Mass.proton, Mass.pion, part_number, 0.005);
	AddEntry("N2190", 2.18, 400E-3, Mass.proton, Mass.pion, part_number, 0.075);
	AddEntry("N2200", 2.25, 375E-3, Mass.proton, Mass.pion, part_number, 0.0113);
	AddEntry("N2250", 2.285, 450E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("N2600", 2.55, 650E-3, Mass.proton, Mass.pion, part_number, 0.02);
	AddEntry("N2700", 2.612, 350E-3, Mass.proton, Mass.pion, part_number, 0.03);

	//delta baryon
	AddEntry("Delta1232", 1.2313, 117E-3, Mass.proton, Mass.pion, part_number, 0.5);
	AddEntry("Delta1600", 1.57, 250E-3, Mass.proton, Mass.pion, part_number, 0.08);
	AddEntry("Delta1620", 1.59, 130E-3, Mass.proton, Mass.pion, part_number, 0.15);
	AddEntry("Delta1700", 1.720, 300E-3, Mass.proton, Mass.pion, part_number, 0.0375);
	AddEntry("Delta1900", 1.88, 250E-3, Mass.proton, Mass.pion, part_number, 0.04);
	AddEntry("Delta1905", 1.855, 335E-3, Mass.proton, Mass.pion, part_number, 0.07);
	AddEntry("Delta1910", 1.9, 300E-3, Mass.proton, Mass.pion, part_number, 0.1);
	AddEntry("Delta1920", 1.92, 300E-3, Mass.proton, Mass.pion, part_number, 0.0625);
	AddEntry("Delta1930", 1.95, 300E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("Delta1940", 2., 400E-3, Mass.proton, Mass.pion, part_number, 0.05);
	AddEntry("Delta1950", 1.9275, 285E-3, Mass.proton, Mass.pion, part_number, 0.1);
	AddEntry("Delta2000", 2.015, 500E-3, Mass.proton, Mass.pion, part_number, 0.035);
	AddEntry("Delta2150", 2.15, 200E-3, Mass.proton, Mass.pion, part_number, 0.04);
	AddEntry("Delta2200", 2.2, 350E-3, Mass.proton, Mass.pion, part_number, 0.025);
	AddEntry("Delta2300", 2.3, 360E-3, Mass.proton, Mass.pion, part_number, 0.0225);
	AddEntry("Delta2350", 2.35, 350E-3, Mass.proton, Mass.pion, part_number, 0.085);
	AddEntry("Delta2390", 2.375, 300E-3, Mass.proton, Mass.pion, part_number, 0.0375);
	AddEntry("Delta2400", 2.45, 500E-3, Mass.proton, Mass.pion, part_number, 0.04);
	AddEntry("Delta2420", 2.5, 500E-3, Mass.proton, Mass.pion, part_number, 0.0375);
	AddEntry("Delta2750", 2.8, 350E-3, Mass.proton, Mass.pion, part_number, 0.02);
	AddEntry("Delta2950", 2.99, 330E-3, Mass.proton, Mass.pion, part_number, 0.02);

	
	//lambda baryons
	AddEntry("Lambda", 1.115, 2E-3, Mass.proton, Mass.pion, part_number, 0.639);
	AddEntry("Lambda1520", 1.5195, 15.6E-3, Mass.proton, Mass.kaon, part_number, 0.225);
	AddEntry("Lambda1600", 1.6, 200E-3, Mass.proton, Mass.kaon, part_number, 0.1125);
	AddEntry("Lambda1670", 1.67, 30E-3, Mass.proton, Mass.kaon, part_number, 0.125);
	AddEntry("Lambda1690", 1.685, 70E-3, Mass.proton, Mass.kaon, part_number, 0.125);
	AddEntry("Lambda1710", 1.713, 180E-3, Mass.proton, Mass.kaon, part_number, 0.115);
	AddEntry("Lambda1800", 1.8, 200E-3, Mass.proton, Mass.kaon, part_number, 0.1875);
	AddEntry("Lambda1810", 1.79, 110E-3, Mass.proton, Mass.kaon, part_number, 0.1);
	AddEntry("Lambda1820", 1.82, 80E-3, Mass.proton, Mass.kaon, part_number, 0.15);
	AddEntry("Lambda1830", 1.825, 80E-3, Mass.proton, Mass.kaon, part_number, 0.015);
	AddEntry("Lambda1890", 1.89, 120E-3, Mass.proton, Mass.kaon, part_number, 0.15);
	AddEntry("Lambda2000", 2., 275E-3, Mass.proton, Mass.kaon, part_number, 0.0675);
	AddEntry("Lambda2050", 2.056, 493E-3, Mass.proton, Mass.kaon, part_number, 0.0475);
	AddEntry("Lambda2070", 2.070, 370E-3, Mass.proton, Mass.kaon, part_number, 0.03);
	AddEntry("Lambda2080", 2.082, 181E-3, Mass.proton, Mass.kaon, part_number, 0.02);
	AddEntry("Lambda2100", 2.1, 175E-3, Mass.proton, Mass.kaon, part_number, 0.075);
	AddEntry("Lambda2110", 2.09, 250E-3, Mass.proton, Mass.kaon, part_number, 0.075);
	AddEntry("Lambda2350", 2.355, 175E-3, Mass.proton, Mass.kaon, part_number, 0.03);

	//sigma baryons
	AddEntry("Sigma1620", 1.625, 70E-3, Mass.proton, Mass.kaon, part_number, 0.175);
	AddEntry("Sigma1660", 1.66, 150E-3, Mass.proton, Mass.kaon, part_number, 0.025);
	AddEntry("Sigma1670", 1.675, 70E-3, Mass.proton, Mass.kaon, part_number, 0.02);
	AddEntry("Sigma1750", 1.75, 150E-3, Mass.proton, Mass.kaon, part_number, 0.045);
	AddEntry("Sigma1775", 1.775, 120E-3, Mass.proton, Mass.kaon, part_number, 0.1);
	AddEntry("Sigma1780", 1.78, 150E-3, Mass.proton, Mass.kaon, part_number, 0.01);
	AddEntry("Sigma1880", 1.880, 150E-3, Mass.proton, Mass.kaon, part_number, 0.1);
	AddEntry("Sigma1900", 1.925, 165E-3, Mass.proton, Mass.kaon, part_number, 0.225);
	AddEntry("Sigma1910", 1.910, 225E-3, Mass.proton, Mass.kaon, part_number, 0.0125);
	AddEntry("Sigma1915", 1.918, 120E-3, Mass.proton, Mass.kaon, part_number, 0.0125);
	AddEntry("Sigma1940", 1.940, 300E-3, Mass.proton, Mass.kaon, part_number, 0.065);
	AddEntry("Sigma2010", 2.005, 178E-3, Mass.proton, Mass.kaon, part_number, 0.035);
	AddEntry("Sigma2030", 2.038, 225E-3, Mass.proton, Mass.kaon, part_number, 0.1);
	AddEntry("Sigma2100", 2.1, 260E-3, Mass.proton, Mass.kaon, part_number, 0.04);
	AddEntry("Sigma2110", 2.105, 313E-3, Mass.proton, Mass.kaon, part_number, 0.065);
	AddEntry("Sigma2230", 2.24, 345E-3, Mass.proton, Mass.kaon, part_number, 0.03);
	AddEntry("Sigma2250", 2.245, 105E-3, Mass.proton, Mass.kaon, part_number, 0.03);
	
	unsigned int sum = 0;

	unsigned int done = 0;

	std::cout << "Particles that will be generated: " << std::endl << std::endl; 
	
	std::cout << "	Name	mean(GeV)	sigma(GeV)	channel	mode" << std::endl; 
	
	for (int i = 0; i < Part.m1.size(); i++)
	{
		std::cout << "	" << Part.name[i] << "	" << Part.mean[i] << "	" << Part.sigma[i] << "	" << Part.channel[i] << "	" << Part.mode[i] << std::endl;
		sum += Part.number[i];
	}

	std::cout << std::endl;
	
	std::cout << sum << " particles will be generated" << std::endl << std::endl;

	
	for (int i = 0; i < Part.m1.size(); i++)
	{
		PrintSeparator("_", "_", "/");
		
		std::cout << OutputColor.bold_green << "[" << i + 1 << " out of " << Part.m1.size() << "]" << OutputColor.reset << " Generating " << Part.number[i] << " particles: " << Part.name[i] << "->" << Part.channel[i] << std::endl;

		chrono_t start = std::chrono::high_resolution_clock::now();
		
		Init(Part.name[i], Part.mean[i], Part.sigma[i], Part.m1[i], Part.m2[i], Part.number[i], Part.mode[i], Part.seed[i]);

		std::cout << std::endl;
	
		chrono_t stop = std::chrono::high_resolution_clock::now();

		done += Part.number[i];
		
		std::cout << "At ";
		Time.PrintCurrent();
		std::cout << "Processing time: ";
		Time.PrintDuration(start, stop);
		std::cout << "Expected remaining time: ";	
		Time.PrintRemaining(start, stop, (double) sum - done, Part.number[i]);
	}
	
	return 0;
}

void AddEntry(std::string part_name, const double mean, const double sigma, const double m1, const double m2, const int part_number, const double mode = 1, const int seed = 1, bool do_antipart = true)
{
	if (m1 == Mass.pion && m2 == Mass.pion) Part.channel.push_back("pi+pi");
	else if (m1 == Mass.kaon && m2 == Mass.pion) Part.channel.push_back("k+pi");
	else if (m2 == Mass.kaon && m1 == Mass.pion) Part.channel.push_back("pi+k");
	else if (m1 == Mass.proton && m2 == Mass.pion) Part.channel.push_back("p+pi");
	else if (m1 == Mass.kaon && m2 == Mass.kaon) Part.channel.push_back("k+k");
	else if (m2 == Mass.proton && m1 == Mass.pion) Part.channel.push_back("pi+p");
	else if (m1 == Mass.proton && m2 == Mass.kaon) Part.channel.push_back("p+k");
	else if (m2 == Mass.proton && m1 == Mass.kaon) Part.channel.push_back("k+p");
	else if (m1 == Mass.proton && m2 == Mass.proton) Part.channel.push_back("p+p");
	else if (m1 == Mass.electron && m2 == Mass.electron) Part.channel.push_back("e+e");
	else Part.channel.push_back("no");
	
	Part.name.push_back(part_name);
	Part.mean.push_back(mean);
	Part.sigma.push_back(sigma);
	Part.m1.push_back(m1);
	Part.m2.push_back(m2);
	Part.number.push_back(part_number);
	Part.seed.push_back(seed);
	Part.mode.push_back(mode);

	if (do_antipart == true)
	{
		std::string antipart_name = "anti" + part_name;
		AddEntry(antipart_name, mean, sigma, m2, m1, part_number, mode, seed, false);
	}
}

void Init(std::string part_name, const double mean, const double sigma, const double m1, const double m2, const int part_number, const double mode, const int seed)
{
	std::string output_file_name = "../data/" + part_name + ".root";
	CheckOutputFile(output_file_name);
	TFile *output = new TFile(output_file_name.c_str(), "RECREATE");

	TRandom *rand = new TRandom(seed);
	
	double p1[3], p2[3];

	TTree *D1 = new TTree("D1", "Daughter1");
	TTree *D2 = new TTree("D2", "Daughter2");

	TH1F *stat = new TH1F("stat", "stat", 10, 0, 10);
	
	stat->SetBinContent(1, mean);
	stat->SetBinContent(2, sigma);
	stat->SetBinContent(3, m1);
	stat->SetBinContent(4, m2);
	stat->SetBinContent(5, mode*exp(1./mean));

	D1->Branch("px", &p1[0]);
	D1->Branch("py", &p1[1]);
	D1->Branch("pz", &p1[2]);

	D2->Branch("px", &p2[0]);
	D2->Branch("py", &p2[1]);
	D2->Branch("pz", &p2[2]);

	TH2D *mass_distr = new TH2D("mass_distr", "mass_distr", 100, 0, 10, 6000, 0, 12);
	mass_distr->SetDefaultSumw2();

	for (int i = 0; i < part_number; i++)
	{
		ProgressBar.Block2((double) (i+1.)/part_number);

		double mass = rand->BreitWigner(mean, sigma);
		if (mass < m1 + m2) continue;
		
		double momentum = Tsallis.GetMom(mass, (unsigned int) i + seed + mass, ptmin, ptmax);
		//momentum = rand->Uniform(0.3, 8.); 
		
		double energy = sqrt(mass*mass + momentum*momentum);

		double e1 = (mass*mass + m1*m1 - m2*m2)/(2*mass);
		double e2 = (mass*mass + m2*m2 - m1*m1)/(2*mass);

		double vel = momentum/energy;
		double gamma = energy/mass;

		double mom1 = sqrt(e1*e1 - m1*m1);
		double mom2 = -1.*sqrt(e2*e2 - m2*m2);

		double theta = rand->Uniform(pi);
		double phi = rand->Uniform(2.*pi);

		p1[0] = mom1*cos(phi)*sin(theta);
		p2[0] = mom2*cos(phi)*sin(theta);

		p1[1] = mom1*sin(phi)*sin(theta);
		p2[1] = mom2*sin(phi)*sin(theta);

		p1[2] = mom1*cos(theta);
		p2[2] = mom2*cos(theta);
		
		p1[2] = gamma*(p1[2] + vel*e1);
		p2[2] = gamma*(p2[2] + vel*e2);

		mom1 = 0;
		mom2 = 0;
		
		for (double p : p1) mom1 += p*p;
		for (double p : p2) mom2 += p*p;
		
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
		
		for (double p : p1) mom1 += p*p;
		for (double p : p2) mom2 += p*p;

		e1 = sqrt(m1*m1 + mom1);
		e2 = sqrt(m2*m2 + mom2);

		double pt = sqrt(pow(p1[0]+p2[0], 2) + pow(p1[1]+p2[1], 2));
		if (pt < 0.3) continue;
		
		mass_distr->Fill(pt, mass);
		
		D1->Fill();
		D2->Fill();
	}

	D1->Write();
	D2->Write();
	mass_distr->Write();
	stat->Write();

	delete rand;
	delete D1;
	delete D2;
	delete mass_distr;
	delete stat;
	delete output;
}
