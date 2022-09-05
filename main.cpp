#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <random>
#include <complex>
#include <cmath>
#include <cassert>
using namespace std;

#include "src/Matrix.h"
#include "src/tensor.h"
#include "random"

//Model files
#include "src/Parameters_DL.h"
#include "src/Coordinates_DL.h"
#include "src/Kspace_calculation_DL_SC.h"

int main(int argc, char *argv[]) {

	string ex_string_original =argv[0];
	string ex_string;
	//ex_string.substr(ex_string_original.length()-5);
	ex_string=ex_string_original.substr (2);
	cout<<"'"<<ex_string<<"'"<<endl;

	if(ex_string=="k_space_SelfConsistency"){
		string model_inputfile = argv[1];

		if (argc<2) { throw std::invalid_argument("USE:: executable inputfile"); }

		Parameters_DL Parameters_DL_;
		Parameters_DL_.Initialize(model_inputfile);

		Coordinates_DL Coordinates_DL_(Parameters_DL_.lx, Parameters_DL_.ly, Parameters_DL_.n_orbs);

		mt19937_64 Generator_(Parameters_DL_.RandomSeed);
		Kspace_calculation_DL Kspace_calculation_DL_(Parameters_DL_, Coordinates_DL_,Generator_);

		Kspace_calculation_DL_.SelfConsistency();

	}

}
