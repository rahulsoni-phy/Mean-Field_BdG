#ifndef Parameters_DL_class
#define Parameters_DL_class
#include "tensor.h"

class Parameters_DL{

	public:
		int lx, ly, ns;
		int n_orbs;
		double OnsiteE,Delta;
		double lambda_RSOC, Bmag, mu;

		double hopping;
		double Temperature,beta;

		int IterMax,RandomSeed;
		double Convergence_Error;

		void Initialize(string inputfile_);
		double matchstring(string file, string match);
		string matchstring2(string file, string match);

};

void Parameters_DL::Initialize(string inputfile_){


	cout << "____________________________________" << endl;
	cout << "Reading the inputfile for specific model: " << inputfile_ << endl;
	cout << "____________________________________" << endl;

	lx = int(matchstring(inputfile_, "Xsite"));
	ly = int(matchstring(inputfile_, "Ysite"));
	n_orbs = int(matchstring(inputfile_, "N_Orbs"));

	hopping = double(matchstring(inputfile_, "Hopping"));
	lambda_RSOC =matchstring(inputfile_, "lambda_RSOC");
	Bmag = double(matchstring(inputfile_, "B_Field"));


	Delta = double(matchstring(inputfile_, "Delta_"));
	mu = double(matchstring(inputfile_, "MU"));
	OnsiteE = double(matchstring(inputfile_, "Onsite_E"));

	IterMax = int(matchstring(inputfile_,"No_of_SelfConsistency_iters"));
	Convergence_Error=matchstring(inputfile_,"Convergence_Error");
	RandomSeed = matchstring(inputfile_,"RandomSeed");
	Temperature = matchstring(inputfile_,"Temperature");
	beta=(1.0/Temperature);
	
	ns = lx * ly;
	cout << "TotalNumberOf Unit cells = " << ns << endl;
	cout << "____________________________________" << endl;
}

double Parameters_DL::matchstring(string file, string match){

	string test;
	string line;
	ifstream readFile(file);
	double amount;
	bool pass = false;
	while (std::getline(readFile, line))
	{
		std::istringstream iss(line);
		if (std::getline(iss, test, '=') && pass == false)
		{
			// ---------------------------------
			if (iss >> amount && test == match)
			{
				// cout << amount << endl;
				pass = true;
			}
			else
			{
				pass = false;
			}
			// ---------------------------------
			if (pass)
				break;
		}
	}
	if (pass == false)
	{
		string errorout = match;
		errorout += "= argument is missing in the input file!";
		throw std::invalid_argument(errorout);
	}
	cout << match << " = " << amount << endl;
	return amount;
}

string Parameters_DL::matchstring2(string file, string match){


	string line;
	ifstream readFile(file);
	string amount;
	int offset;

	if (readFile.is_open())
	{
		while (!readFile.eof())
		{
			getline(readFile, line);

			if ((offset = line.find(match, 0)) != string::npos)
			{
				amount = line.substr(offset + match.length() + 1);
			}
		}
		readFile.close();
	}
	else
	{
		cout << "Unable to open input file while in the Parameters_DL class." << endl;
	}

	cout << match << " = " << amount << endl;
	return amount;
}

#endif
