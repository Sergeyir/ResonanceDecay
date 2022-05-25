#include "TMath.h"
  
  
   Double_t Tsallis(Double_t m, int seed)
{
	TRandom *r = new TRandom(seed);
	
	Double_t p, prob = 0, norm = -106254.117, prob1;
	int i = 0;
	m = m/1000;
	
	//integrating the distribution for specific mass to calculate constant
	
	for (p = 0; p <= 10; p = p + 0.01) {
	
		prob = prob + norm*(0.0327538/pow((1+sqrt(m*m+p*p)*2.14),6.6667)-0.0385338/pow((1+sqrt(m*m+p*p)*2.14),5.6667));
		
	}
	
	//calculating constant
	
	norm = norm/prob;
	
	while (i != 1) {
	
		//the interval for momentum does not start from 0 to increase statistics for larger momentum
	
		p = r->Uniform(0.6, 8);
		
		prob = norm*(0.0327538/pow((1+sqrt(m*m+p*p)*2.14),6.6667)-0.0385338/pow((1+sqrt(m*m+p*p)*2.14),5.6667));
		if (prob >= r->Uniform(0,1)) i = 1;
		
	}
	
	return p;
}
