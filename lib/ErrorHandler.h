#pragma once

#include <iostream>
#include <string>
#include <fstream>

void PrintError()
{
	printf("\033[1m\033[31mError: ");
	printf("\033[0m");
	exit(1);
}


void PrintErrorWithoutExit()
{
	printf("\033[1m\033[31mError: ");
	printf("\033[0m");
}

void PrintError(string message)
{
	printf("\033[1m\033[31mError: ");
	printf("\033[0m");
	cout << message << endl;
	exit(1);
}

void PrintWarning()
{
	printf("\033[1m\033[33mWarning: ");
	printf("\033[0m");
}

void PrintWarning(string message)
{
	printf("\033[1m\033[33mWarning: ");
	printf("\033[0m");
	cout << message << endl;
}

void CheckInputFile(string name)
{
	ifstream file(name.c_str());

	if(!file.is_open())
	{
		PrintErrorWithoutExit();
		cout << "File " << name << " not found" << endl;
		exit(1);
	}
}

void CheckOutputFile(string name)
{
	ofstream file(name.c_str());

	if(!file.is_open())
	{
		PrintErrorWithoutExit();
		cout << "File " << name << " cannot be created" << endl;
		exit(1);
	}
}
