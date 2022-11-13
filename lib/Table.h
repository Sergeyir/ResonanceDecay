#pragma once

#include "OutputTool.h"
#include "ErrorHandler.h"

template<typename... Ts>
class Table
{
	private:
	
	int row_length;

	public:
	
	Table() 
	{
		row_length = OutputTool.length;
	};

	void Begin(std::string name)
	{
		OutputTool.PrintSimpleSeparator("|", "=", "|");
		OutputTool.PrintSeparator(name, OutputColor.green, "|", " ", "|");
		OutputTool.PrintSimpleSeparator("|", "=", "|");
	}

	void PrintHeader(Ts... args)
	{
		constexpr int size = sizeof...(args);

		std::string dummy[size] = {(std::string) args...};
		
		float cell_size = (float) (row_length - 3*size - 1)/size;

		std::cout << "| ";
		for (int i = 0; i < size; i++)
		{
			std::cout << dummy[i];
			if (i < size - 1)
			{	
				int new_cell_size = row_length - (int) cell_size*(size-1);
				for (int j = 0; j < cell_size - dummy[i].length(); j++)
				{
					std::cout << " ";
				}
				std::cout << " | ";
			}
			else
			{
				for (int j = 0; j < cell_size - dummy[i].length() - 1; j++)
				{
					std::cout << " ";
				}
			}
		}
		std::cout << " |" << std::endl;

		OutputTool.PrintSimpleSeparator();
	}

	void PrintRow(Ts... args)
	{
		constexpr int size = sizeof...(args);
		
		std::string dummy[size] = {args...};
		
		float cell_size = (float) (row_length - 3*size - 1)/size;

		std::cout << "| ";
		for (int i = 0; i < size; i++)
		{
			std::cout << dummy[i];
			if (i < size - 1)
			{	
				int new_cell_size = row_length - (int) cell_size*(size-1);
				for (int j = 0; j < cell_size - dummy[i].length(); j++)
				{
					std::cout << " ";
				}
				std::cout << " | ";
			}
			else
			{
				for (int j = 0; j < cell_size - dummy[i].length() - 1; j++)
				{
					std::cout << " ";
				}
			}
		}
		std::cout << " |" << std::endl;
	}

	void End()
	{
		OutputTool.PrintSimpleSeparator("|", "=", "|");
	}

};
