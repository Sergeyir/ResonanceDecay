#pragma once

#include <iostream>

#include "TMath.h"

#include "../lib/Tools.h"
#include "../lib/Colors.h"
#include "../lib/ErrorHandler.h"

struct
{
	Double_t GetMom(Double_t m, unsigned int seed, const float mom_min, const float mom_max)
	{
		double norm = -1;
		Tools_t<double> tools;
		//calculating norm constant
		for (float p = mom_min; p <= mom_max; p += 0.01)
		{
			norm = tools.GetMax((0.0327538/pow((1+sqrt(m*m+p*p)*2.14),6.6667)-0.0385338/pow((1+sqrt(m*m+p*p)*2.14),5.6667)), norm);
		}
		
		TRandom *rand = new TRandom(seed);
		
		//integrating the distribution for specific mass to calculate constant
	
		float mom;
		bool check = true;
		
		while (check) {
		
			//the interval for momentum does not start from 0 to increase statistics for larger momentum
			mom = rand->Uniform(mom_min, mom_max);
			
			const float prob = (0.0327538/pow((1+sqrt(m*m+mom*mom)*2.14),6.6667)-0.0385338/pow((1+sqrt(m*m+mom*mom)*2.14),5.6667))/norm;
			if (prob >= rand->Uniform(0, 1)) check = false;
		}

		delete rand;
		return mom;
	}
} Tsallis;
