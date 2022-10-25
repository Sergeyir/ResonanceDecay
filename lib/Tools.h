#pragma once

template <class T> 
struct Tools_t
{
	T GetMax(T val1, T val2, T val3)
	{
		if (val1 > val2 && val1 > val3) return val1;
		else if (val2 > val3) return val2;
		else return val3;
	}

	T GetMin(T val1, T val2, T val3)
	{
		if (val1 < val2 && val1 < val3) return val1;
		else if (val2 < val3) return val2;
		else return val3;
	}

	T GetMax(T val1, T val2)
	{
		if (val1 > val2) return val1;
		else return val2;
	}

	T GetMin(T val1, T val2)
	{
		if (val1 < val2) return val1;
		else return val2;
	}
};
