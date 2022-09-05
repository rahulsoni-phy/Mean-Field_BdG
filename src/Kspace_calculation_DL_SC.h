#include <algorithm>
#include <functional>
#include <math.h>
#include "tensor.h"
#include "Parameters_DL.h"
#include "Coordinates_DL.h"
#include "random"
#include "Matrix.h"
#define PI acos(-1.0)

#ifndef Kspace_calculation_DL_class
#define Kspace_calculation_DL_class

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *, 
		std::complex<double> *,int *, double *, int *);

class Kspace_calculation_DL{
	public:
		Kspace_calculation_DL(Parameters_DL &Parameters__, Coordinates_DL &Coordinates__, mt19937_64& Generator1__ )
			: Parameters_(Parameters__), Coordinates_(Coordinates__), Generator1_(Generator1__)
		{
			Initialize();
		}

		void Initialize();
		void Diagonalize(char option);
		double random1();
		void SelfConsistency();
		void Create_Kspace_Spectrum();
		void Arranging_spectrum();
		//		double chemicalpotential(double Particles);
		void Get_new_OPs_and_error();
		//		void Get_Energies();
		void Get_Bands();
		void Calculate_ChernNumbers();
		void Create_Kspace_Spectrum_in_double_brillouin_zone();
		void Calculate_ChernNumbers_in_double_brillouin_zone();

		mt19937_64 &Generator1_;
		uniform_real_distribution<double> dis1_;
		Parameters_DL &Parameters_;
		Coordinates_DL &Coordinates_;
		int lx_, ly_, ncells_, n_orbs_;
		Matrix<complex<double>> Ham_;
		vector<double> eigs_;
		Mat_2_Complex_doub Eigvectors_;
		Mat_1_doub Eigenvalues_;
		Mat_2_doub Eigenvalues_sector;
		Mat_2_Complex_doub Eigvectors_saved;
		Mat_1_doub Eigenvalues_saved;
		Mat_1_doub Kx_values, Ky_values;

		Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
		Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

		Mat_3_Complex_doub OPs_,OPs_new_;

		double mu_;
		double OP_error_;
		double E_quant, E_class;

};


void Kspace_calculation_DL::Calculate_ChernNumbers_in_double_brillouin_zone(){

}


void Kspace_calculation_DL::Calculate_ChernNumbers(){

}


void Kspace_calculation_DL::Get_Bands(){

	string File_Out_Bands;
	File_Out_Bands = "Bands.txt";
	ofstream file_out_bands(File_Out_Bands.c_str());
	int k_index;
	int n1, n2;

	Mat_1_intpair k_path;
	k_path.clear();
	Mat_1_intpair k_path2;
	k_path2.clear();
	pair_int temp_pair;


	// ---k_path---------
	//--------\Gamma to K----------------
	n1=0;
	n2=0;
	while (n2<=int(Parameters_.ly/3))
	{
		temp_pair.first = n1;
		temp_pair.second = n2;
		k_path.push_back(temp_pair);
		n2++;
		n1=int((2*Parameters_.lx/Parameters_.ly)*n2);
	}
	//----------------------------------

	//--------K to M-----------------
	n2=int(Parameters_.ly/3);
	n1=int((2*Parameters_.lx/Parameters_.ly)*n2);
	n2++;
	n1--;
	while (n1>=int(Parameters_.lx/2))
	{
		temp_pair.first = n1;
		temp_pair.second = n2;
		k_path.push_back(temp_pair);
		n2++;
		n1--;
	}
	//----------------------------------
	//--------M to \Gamma[with one extra point,
	//                  because in gnuplot use "set pm3d corners2color c1"
	//                  ]-----------------
	n1=int(Parameters_.lx/2);
	n2=int(Parameters_.ly/2);
	n2--;
	n1--;
	while (n1>=0)
	{
		temp_pair.first = n1;
		temp_pair.second = n2;
		k_path.push_back(temp_pair);
		n2--;
		n1--;
	}

	temp_pair.first = 0;
	temp_pair.second = 0;
	k_path.push_back(temp_pair);

	//----------------------------------
	cout<<"PRINTING PATH"<<endl;
	for (int k_point = 0; k_point < k_path.size(); k_point++)
	{
		cout<<k_path[k_point].first<< "   "<<k_path[k_point].second<<endl;
	}
	//----k_path done-------

	for (int k_point = 0; k_point < k_path.size(); k_point++)
	{
		k_index=Coordinates_.Ncell(k_path[k_point].first, k_path[k_point].second);

		file_out_bands<<k_point<<"   ";
		for(int band=0;band<12;band++){
			file_out_bands<<Eigenvalues_[12*k_index + band]<<"   ";
		}
		file_out_bands<<endl;
	}

}

void Kspace_calculation_DL::Arranging_spectrum(){

	Eigvectors_saved = Eigvectors_;
	Eigenvalues_saved = Eigenvalues_;

	double value_;
	Mat_1_Complex_doub Vec_temp;
	Vec_temp.resize(12);
	for(int i=0;i<Eigenvalues_.size();i++){
		for(int j=i+1;j<Eigenvalues_.size();j++){

			if(Eigenvalues_[j]<Eigenvalues_[i]){
				value_=Eigenvalues_[i];

				for(int comp=0;comp<12;comp++){
					Vec_temp[comp]=Eigvectors_[i][comp];
				}
				Eigenvalues_[i]=Eigenvalues_[j];

				for(int comp=0;comp<12;comp++){
					Eigvectors_[i][comp]=Eigvectors_[j][comp];
				}
				Eigenvalues_[j]=value_;

				for(int comp=0;comp<12;comp++){
					Eigvectors_[j][comp]=Vec_temp[comp];
				}
			}
		}
	}
}

double Kspace_calculation_DL::random1(){

	return dis1_(Generator1_);

}

void Kspace_calculation_DL::Initialize(){

	ly_ = Parameters_.ly;
	lx_ = Parameters_.lx;
	ncells_ = lx_ * ly_;
	n_orbs_ = Parameters_.n_orbs;
	int space = 2 * ncells_ * n_orbs_;

	int ind1,ind2;


	Ham_.resize(12, 12);
	Eigvectors_.resize(ncells_*12);
	Eigenvalues_.resize(ncells_*12);

	Eigenvalues_sector.resize(ncells_);
	for(int i=0;i<ncells_;i++){
		Eigenvalues_sector[i].resize(12);
	}

	for(int i=0;i<ncells_*12;i++){
		Eigvectors_[i].resize(12); //Eigenvector number
	}

	Kx_values.resize(ncells_);
	Ky_values.resize(ncells_);

	OPs_.resize(ncells_);
	OPs_new_.resize(ncells_);

	for(int i=0;i<ncells_;i++){
		OPs_[i].resize(6);
		OPs_new_[i].resize(6);
		for(int m=0;m<6;m++){
			OPs_[i][m].resize(6);
			OPs_new_[i][m].resize(6);
		}
	}

	for(int i=0;i<ncells_;i++){
		for(int alpha=0;alpha<3;alpha++){
			for(int s=0;s<2;s++){
				ind1 = alpha + 3*s;

				for(int beta=0;beta<3;beta++){
					for(int sp=0;sp<2;sp++){
						ind2 = beta + 3*sp;

						if(abs(alpha-beta)==1){
							OPs_[i][ind1][ind2]= complex<double>(random1(), random1());
						}
						else{
							OPs_[i][ind1][ind2]= zero_complex;
						}
					}
				}
			}
		}
	}


}




void Kspace_calculation_DL::Create_Kspace_Spectrum(){

	int k_index,ind1,ind2;
	double k1x_val, k2x_val, k1y_val, k2y_val;
	int HS_=6;

	Mat_2_Complex_doub Gamma_k_,Gamma_k_plus_,Gamma_k_minus_; //A good change
	Gamma_k_.resize(lx_);	Gamma_k_plus_.resize(lx_); Gamma_k_minus_.resize(lx_);
	for(int k1=0;k1<lx_;k1++){
		Gamma_k_[k1].resize(ly_);
		Gamma_k_plus_[k1].resize(ly_);
		Gamma_k_minus_[k1].resize(ly_);
	}

	for(int k1=0;k1<lx_;k1++){
		for(int k2=0;k2<ly_;k2++){

			k_index = Coordinates_.Ncell(k1,k2);

			Ham_.fill(0.0);

			Gamma_k_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*((2.0*PI*k1)/(1.0*lx_))) 
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)))) );

			Gamma_k_plus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) + (2.0*PI/3.0))) 
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) + (4.0*PI/3.0))) );

			Gamma_k_minus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) - (2.0*PI/3.0)))
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) - (4.0*PI/3.0))) );


			//Creating Lower Triangular Hamiltonian
			//Adding Chemical Potential, Zeeman, and Onsite Energy
			for(int s=0;s<2;s++){
				for(int alpha=0;alpha<3;alpha++){
					ind1=alpha+3*s;
					if(s==0){
						Ham_(ind1,ind1) += one_complex*( -1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ + ind1,HS_ + ind1) += one_complex*( 1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(s==1){
						Ham_(ind1,ind1) += one_complex*( 1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ + ind1,HS_ + ind1) += one_complex*( -1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(alpha==1){
						Ham_(ind1,ind1) += -1.0*Parameters_.OnsiteE*one_complex;
						Ham_(HS_ + ind1,HS_ + ind1) += 1.0*Parameters_.OnsiteE*one_complex;
					}
				}
			}

			//Adding Hopping Term
			for(int ind=0;ind<6;ind++){
				if(ind%3!=0){
					Ham_(ind,ind-1) += -1.0*Parameters_.hopping*Gamma_k_[k1][k2];
					Ham_(HS_ + ind,HS_ + ind-1) += 1.0*Parameters_.hopping*Gamma_k_[k1][k2];
				}
			}

			//Adding Rashba-SOC Term
			Ham_(3,1) += -iota_complex*Parameters_.lambda_RSOC*conj(Gamma_k_minus_[k1][k2]);
			Ham_(4,0) += iota_complex*Parameters_.lambda_RSOC*Gamma_k_plus_[k1][k2];
			Ham_(4,2) += iota_complex*Parameters_.lambda_RSOC*conj(Gamma_k_minus_[k1][k2]);
			Ham_(5,1) += -iota_complex*Parameters_.lambda_RSOC*Gamma_k_plus_[k1][k2];

			Ham_(HS_ +3,HS_ +1) += -iota_complex*Parameters_.lambda_RSOC*conj(Gamma_k_plus_[k1][k2]);
			Ham_(HS_ +4,HS_ +0) += iota_complex*Parameters_.lambda_RSOC*Gamma_k_minus_[k1][k2];
			Ham_(HS_ +4,HS_ +2) += iota_complex*Parameters_.lambda_RSOC*conj(Gamma_k_plus_[k1][k2]);
			Ham_(HS_ +5,HS_ +1) += -iota_complex*Parameters_.lambda_RSOC*Gamma_k_minus_[k1][k2];

			//Pairing: Onsite Term (M_Delta Term)
			for(int alpha=0;alpha<3;alpha++){
				for(int s=0;s<2;s++){
					ind1 = alpha + 3*s;
					for(int beta=0;beta<3;beta++){
						for(int sp=0;sp<2;sp++){
							ind2 = beta + 3*sp;

							Ham_(HS_+ind1,ind2) += 1.0*conj( OPs_[k_index][ind2][ind1] );

						}
					}
				}
			}

			//Creating Upper Triangular Hamiltonian (without diagonal)
			for(int row=0;row<2*HS_;row++){
				for(int col=row+1;col<2*HS_;col++){
					Ham_(row,col) = conj(Ham_(col,row));
				}
			}

			char Dflag='V'; //Since we are getting Eigenvectors
			Diagonalize(Dflag);

			for(int row=0;row<2*HS_;row++){
				Eigenvalues_[12*k_index + row]=eigs_[row]; //Each k=(k1,k2) have 12 Eigens
				Eigenvalues_sector[k_index][row]=eigs_[row];
				for(int col=0;col<2*HS_;col++){
					Eigvectors_[12*k_index + col][row]=Ham_(row,col);
				}
			}
		}
	}
}

void Kspace_calculation_DL::Get_new_OPs_and_error(){

	int k_index, ind1, ind2;
	for(int k1=0;k1<lx_;k1++){
		for(int k2=0;k2<ly_;k2++){
			k_index = Coordinates_.Ncell(k1,k2);

			for(int alpha=0;alpha<3;alpha++){
				for(int s=0;s<2;s++){
					ind1 = alpha + 3*s;

					for(int beta=0;beta<3;beta++){
						for(int sp=0;sp<2;sp++){
							ind2 = beta + 3*sp;

							if( abs(alpha - beta)==1 ){
								for(int n=0;n<6;n++){
									OPs_new_[k_index][ind1][ind2] += ()*(1.0/( exp((Eigenvalues_sector[k_index][n])*Parameters_.beta ) + 1.0)) + 
										()*(1.0 - 1.0/( exp((Eigenvalues_sector[k_index][n])*Parameters_.beta ) + 1.0));
								}
							}

							else{
								OPs_new_[k_index][ind1][ind2] = zero_complex;
							}
						}
					}
				}
			}
		}
	}

	//Error in Order Params:
	OP_error_=0.0;
	for(int k_ind=0;k_ind<ncells_;k_ind++){
		for(int ind1=0;ind1<6;ind1++){
			for(int ind2=0;ind2<6;ind2++){
				OP_error_ += abs((OPs_[k_ind][ind1][ind2] - OPs_new_[k_ind][ind1][ind2])*conj( OPs_[k_ind][ind1][ind2] - OPs_new_[k_ind][ind1][ind2] ));
			}
		}
	}
	OP_error_ = sqrt(OP_error_);
}


void Kspace_calculation_DL::SelfConsistency(){

	string File_out_progress;
	File_out_progress = "output_Kspace_SelfConsistency.txt";
	ofstream file_out_progress(File_out_progress.c_str());

	cout<<"error targetted = "<<Parameters_.Convergence_Error<<endl;
	cout<<"Max iterations = "<<Parameters_.IterMax<<endl;

	OP_error_=10.0;
	int iter=0;
	while( (OP_error_>=Parameters_.Convergence_Error) && (iter<=Parameters_.IterMax) ){
		Create_Kspace_Spectrum();
		Get_new_OPs_and_error();

		file_out_progress<<iter<<"	"<<OP_error_<<endl;

		for(int i=0;i<ncells_;i++){
			for(int ind1=0;ind1<6;ind1++){
				for(int ind2=0;ind2<6;ind2++){
					OPs_[i][ind1][ind2] = OPs_new_[i][ind1][ind2];
				}
			}
		}

		iter++;
	}

	Create_Kspace_Spectrum();
	Get_Bands();
	Calculate_ChernNumbers();
	cout<<"Chern numbers after doubling brillouin zone"<<endl;
	Create_Kspace_Spectrum_in_double_brillouin_zone();
	cout<<"Spectrum created"<<endl;
	Calculate_ChernNumbers_in_double_brillouin_zone();
}


void Kspace_calculation_DL::Diagonalize(char option){

	char jobz=option;
	char uplo='L';
	int n=Ham_.n_row();
	int lda=Ham_.n_col();
	vector<complex<double>> work(3);
	vector<double> rwork(3*n -2);
	int info;
	int lwork= -1;

	eigs_.resize(Ham_.n_row());
	fill(eigs_.begin(),eigs_.end(),0);

	// query:
	zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

	//lwork = int(real(work[0]))+1;
	lwork = int((work[0].real()));
	work.resize(lwork);

	// real work:
	zheev_(&jobz,&uplo,&n,&(Ham_(0,0)),&lda,&(eigs_[0]),&(work[0]),&lwork,&(rwork[0]),&info);

	if (info!=0) {
		std::cerr<<"info="<<info<<"\n";
		perror("diag: zheev: failed with info!=0.\n");
	}

	// Ham_.print();
	//
	// for(int i=0;i<eigs_.size();i++){
	//	cout<<eigs_[i]<<endl;
	//}

}


#endif
