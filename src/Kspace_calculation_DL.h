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
		//		double random1();
		void SelfConsistency();
		void Create_Kspace_Spectrum();
		void Arranging_spectrum();
		//		double chemicalpotential(double Particles);
		//		void Get_new_OPs_and_error();
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
		Mat_2_Complex_doub Eigvectors_saved;
		Mat_1_doub Eigenvalues_saved;
		Mat_1_doub Kx_values, Ky_values;

		Matrix<double> Tx, Ty, Tpxpy, Tpxmy;
		Matrix<double> IntraCell_Hopp, InterCell_px_Hopp, InterCell_py_Hopp, InterCell_pxmy_Hopp;

		Mat_1_Complex_doub OPs_, OPs_new_;
		double mu_;
		double OP_error_;
		double E_quant, E_class;

};

void Kspace_calculation_DL::Calculate_ChernNumbers_in_double_brillouin_zone(){


	//cout<<"here"<<endl;
	Matrix<complex<double>> F_mat; //F1, F2, F3, F4, F5;
	F_mat.resize(12, 4*lx_*ly_);
	//cout<<"here"<<endl;

	complex<double> Ux_k, Uy_k, Ux_kpy, Uy_kpx;
	vector<complex<double>> F_bands;
	F_bands.resize(12);
	vector<complex<double>> Chern_num;
	Chern_num.resize(12);


	for (int band = 0; band < 12; band++)
	{
		//cout<<"here"<<endl;
		string file_Fk="Fk_band"+to_string(band)+"_double_brillouin_zone.txt";
		ofstream fl_Fk_out(file_Fk.c_str());
		fl_Fk_out<<"#nx  ny  tilde_F(nx,ny).real()  tilde_F(nx,ny).imag()  ArgofLog.real()  ArgofLog.imag()"<<endl;
		fl_Fk_out<<"#Extra momentum point for pm3d corners2color c1"<<endl;

		F_bands[band] = 0.0;
		for (int nx = 0; nx < 2*lx_; nx++)
		{
			for (int ny = 0; ny < 2*ly_; ny++)
			{

				int n = nx + 2*lx_*ny;
				int n_left, n_right, nx_left, ny_left, nx_right, ny_right;

				//U1_k
				Ux_k = 0;
				n_left = n;
				nx_right = (nx + 1) % (2*lx_);
				ny_right = ny;
				n_right = nx_right + 2*lx_*ny_right;
				for (int comp = 0; comp < 12; comp++)
				{
					Ux_k +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Ux_k = Ux_k * (1.0 / abs(Ux_k));

				//U2_kpx
				Uy_kpx = 0;
				nx_left = (nx + 1) % (2*lx_);
				ny_left = ny;
				n_left = nx_left + 2*lx_*ny_left;
				nx_right = nx_left;
				ny_right = (ny_left + 1) % (2*ly_);
				n_right = nx_right + 2*lx_*ny_right;
				for (int comp = 0; comp < 12; comp++)
				{
					Uy_kpx +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Uy_kpx = Uy_kpx * (1.0 / abs(Uy_kpx));

				//U1_kpy
				Ux_kpy = 0;
				nx_left = nx;
				ny_left = (ny + 1) % (2*ly_);
				n_left = nx_left + 2*lx_*ny_left;
				nx_right = (nx_left + 1) % (2*lx_);
				ny_right = ny_left;
				n_right = nx_right + 2*lx_*ny_right;
				for (int comp = 0; comp < 12; comp++)
				{
					Ux_kpy +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Ux_kpy = Ux_kpy * (1.0 / abs(Ux_kpy));

				//U2_k
				Uy_k = 0;
				nx_left = nx;
				ny_left = ny;
				n_left = nx_left + 2*lx_*ny_left;
				nx_right = nx_left;
				ny_right = (ny_left + 1) % (2*ly_);
				n_right = nx_right + 2*lx_*ny_right;
				for (int comp = 0; comp < 12; comp++)
				{
					Uy_k +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Uy_k = Uy_k * (1.0 / abs(Uy_k));

				// Calculating tilde F12
				F_mat(band, n) = log(Ux_k *
						Uy_kpx *
						conj(Ux_kpy) * conj(Uy_k));

				F_bands[band] += F_mat(band, n);


				fl_Fk_out.precision(10);

				fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
					"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
					"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;




				if(abs((abs(F_mat(band, n).imag()) - PI))<0.0000001){
					cout<<ny<<"  "<<nx<<"  gives Pi for band"<< band <<endl;
					// assert (abs((abs(F_mat(band, n).imag()) - M_PI))>0.0000001);
				}


				if(ny==2*ly_-1){//For pm3d corners2color c1
					fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
				}
			}
			if(nx==2*lx_-1){//For pm3d corners2color c1
				fl_Fk_out<<endl;
				for(int ny_=0;ny_<2*ly_;ny_++){
					int n_ = nx + 2*lx_* ny_;
					fl_Fk_out<<nx<<"  "<<ny_<<"  "<<F_mat(band, n_).real()<<"  "<<F_mat(band, n_).imag()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
				}
			}
			fl_Fk_out<<endl;
		}


		Chern_num[band] = (-1.0 * iota_complex / (2 * PI*4)) * F_bands[band];
		// Chern_num_orgnl[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands_orgnl[band];
		fl_Fk_out<<"#Chern no*2pi*Iota= "<<F_bands[band].real()<<"  "<<F_bands[band].imag()<<endl;
		cout << "tilde Chern number [" << band << "] = " << Chern_num[band].real() << "        " << Chern_num[band].imag() << endl;
		//  cout << "Chern number [" << band << "] = " << Chern_num_orgnl[band].real() << " " << Chern_num_orgnl[band].imag() << endl;

	}

}


void Kspace_calculation_DL::Calculate_ChernNumbers(){



	Matrix<complex<double>> F_mat; //F1, F2, F3, F4, F5;
	F_mat.resize(12, lx_*ly_);

	complex<double> Ux_k, Uy_k, Ux_kpy, Uy_kpx;
	vector<complex<double>> F_bands;
	F_bands.resize(12);
	vector<complex<double>> Chern_num;
	Chern_num.resize(12);

	vector<complex<double>> F_bands_orgnl;
	F_bands_orgnl.resize(12);
	vector<complex<double>> Chern_num_orgnl;
	Chern_num_orgnl.resize(12);
	for (int band = 0; band < 12; band++)
	{
		string file_Fk="Fk_band"+to_string(band)+".txt";
		ofstream fl_Fk_out(file_Fk.c_str());
		fl_Fk_out<<"#nx  ny  tilde_F(nx,ny).real()  tilde_F(nx,ny).imag()  ArgofLog.real()  ArgofLog.imag()"<<endl;
		fl_Fk_out<<"#Extra momentum point for pm3d corners2color c1"<<endl;

		//        string file_Fk_orgnl="Fk_original_band"+to_string(band)+".txt";
		//        ofstream fl_Fk_orgnl_out(file_Fk_orgnl.c_str());
		//        fl_Fk_orgnl_out<<"#nx  ny  F(nx,ny).real()*(2pi/Lx)*(2pi/Ly)  F(nx,ny).imag()*(2pi/Lx)*(2pi/Ly)"<<endl;

		F_bands[band] = 0.0;
		F_bands_orgnl[band] = 0.0;
		for (int nx = 0; nx < lx_; nx++)
		{
			for (int ny = 0; ny < ly_; ny++)
			{
				int n = Coordinates_.Ncell(nx,ny);
				int n_left, n_right, nx_left, ny_left, nx_right, ny_right;

				//U1_k
				Ux_k = 0;
				n_left = n;
				nx_right = (nx + 1) % lx_;
				ny_right = ny;
				n_right = Coordinates_.Ncell(nx_right, ny_right);
				for (int comp = 0; comp < 12; comp++)
				{
					Ux_k +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Ux_k = Ux_k * (1.0 / abs(Ux_k));

				//U2_kpx
				Uy_kpx = 0;
				nx_left = (nx + 1) % lx_;
				ny_left = ny;
				n_left = Coordinates_.Ncell(nx_left, ny_left);
				nx_right = nx_left;
				ny_right = (ny_left + 1) % ly_;
				n_right = Coordinates_.Ncell(nx_right, ny_right);
				for (int comp = 0; comp < 12; comp++)
				{
					Uy_kpx +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Uy_kpx = Uy_kpx * (1.0 / abs(Uy_kpx));

				//U1_kpy
				Ux_kpy = 0;
				nx_left = nx;
				ny_left = (ny + 1) % ly_;
				n_left = Coordinates_.Ncell(nx_left, ny_left);
				nx_right = (nx_left + 1) % lx_;
				ny_right = ny_left;
				n_right = Coordinates_.Ncell(nx_right, ny_right);
				for (int comp = 0; comp < 12; comp++)
				{
					Ux_kpy +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Ux_kpy = Ux_kpy * (1.0 / abs(Ux_kpy));

				//U2_k
				Uy_k = 0;
				nx_left = nx;
				ny_left = ny;
				n_left = Coordinates_.Ncell(nx_left, ny_left);
				nx_right = nx_left;
				ny_right = (ny_left + 1) % ly_;
				n_right = Coordinates_.Ncell(nx_right, ny_right);
				for (int comp = 0; comp < 12; comp++)
				{
					Uy_k +=
						conj(Eigvectors_[12*n_left +  band][comp]) *
						Eigvectors_[12*n_right + band][comp];
				}
				Uy_k = Uy_k * (1.0 / abs(Uy_k));

				// Calculating tilde F12
				F_mat(band, n) = log(Ux_k *
						Uy_kpx *
						conj(Ux_kpy) * conj(Uy_k));

				F_bands[band] += F_mat(band, n);


				fl_Fk_out.precision(10);

				fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
					"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
					"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;




				if(abs((abs(F_mat(band, n).imag()) - PI))<0.0000001){
					cout<<ny<<"  "<<nx<<"  gives Pi for band"<< band <<endl;
					// assert (abs((abs(F_mat(band, n).imag()) - M_PI))>0.0000001);
				}


				if(ny==ly_-1){//For pm3d corners2color c1
					fl_Fk_out<<nx<<"  "<<ny<<"  "<<F_mat(band, n).real()<<"  "<<F_mat(band, n).imag()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
				}
			}
			if(nx==lx_-1){//For pm3d corners2color c1
				fl_Fk_out<<endl;
				for(int ny_=0;ny_<ly_;ny_++){
					int n_ = nx + lx_*ny_;
					fl_Fk_out<<nx<<"  "<<ny_<<"  "<<F_mat(band, n_).real()<<"  "<<F_mat(band, n_).imag()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).real()<<
						"  "<<(Ux_k*Uy_kpx*conj(Ux_kpy)*conj(Uy_k)).imag()<<endl;
				}
			}
			fl_Fk_out<<endl;
		}


		Chern_num[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands[band];
		// Chern_num_orgnl[band] = (-1.0 * iota_complex / (2 * PI)) * F_bands_orgnl[band];
		fl_Fk_out<<"#Chern no*2pi*Iota= "<<F_bands[band].real()<<"  "<<F_bands[band].imag()<<endl;
		cout << "tilde Chern number [" << band << "] = " << Chern_num[band].real() << "        " << Chern_num[band].imag() << endl;
		//  cout << "Chern number [" << band << "] = " << Chern_num_orgnl[band].real() << " " << Chern_num_orgnl[band].imag() << endl;

	}

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

void Kspace_calculation_DL::Initialize(){

	ly_ = Parameters_.ly;
	lx_ = Parameters_.lx;
	ncells_ = lx_ * ly_;
	n_orbs_ = Parameters_.n_orbs;
	int space = 2 * ncells_ * n_orbs_;

	Ham_.resize(12, 12);
	Eigvectors_.resize(ncells_*12);
	Eigenvalues_.resize(ncells_*12);

	for(int i=0;i<ncells_*12;i++){
		Eigvectors_[i].resize(12); //Eigenvector number
	}

	Kx_values.resize(ncells_);
	Ky_values.resize(ncells_);

	OPs_.resize(3);
	OPs_new_.resize(3);

	for(int i=0;i<OPs_.size();i++){
		if(i==1){
			OPs_[i] = one_complex*Parameters_.Delta2;
		}else{
			OPs_[i] = one_complex*Parameters_.Delta1;
		}
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

void Kspace_calculation_DL::Create_Kspace_Spectrum_in_double_brillouin_zone(){

	int k_index;
        double k1x_val, k2x_val, k1y_val, k2y_val;
        int HS_=6;

        Mat_2_Complex_doub Gamma_k_,Gamma_k_plus_,Gamma_k_minus_; //A good change
        Gamma_k_.resize(2*lx_);   Gamma_k_plus_.resize(2*lx_); Gamma_k_minus_.resize(2*lx_);
        for(int k1=0;k1<2*lx_;k1++){
                Gamma_k_[k1].resize(2*ly_);
                Gamma_k_plus_[k1].resize(2*ly_);
                Gamma_k_minus_[k1].resize(2*ly_);
        }

	Eigvectors_.resize(ncells_*12*4);
	Eigenvalues_.resize(ncells_*12*4);
	for(int i=0;i<ncells_*12*4;i++){
		Eigvectors_[i].resize(12); //Eigenvector number
	}

	Kx_values.resize(ncells_*4);
	Ky_values.resize(ncells_*4);


	for(int k1=0;k1<2*lx_;k1++){
		k1x_val = (2.0*PI*k1)/(1.0*lx_);
		k1y_val = ((2.0*PI*k1)/(1.0*lx_))*(-1.0/sqrt(3));

		for(int k2=0;k2<2*ly_;k2++){
			k2x_val = 0.0;
			k2y_val = ((2.0*PI*k2)/(1.0*ly_))*(2.0/sqrt(3));

			k_index = k1 + 2*lx_*k2;
			Kx_values[k_index]=k1x_val + k2x_val;
			Ky_values[k_index]=k1y_val + k2y_val;


			Ham_.fill(0.0);

			Gamma_k_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*((2.0*PI*k1)/(1.0*lx_)))
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)))) );

			Gamma_k_plus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) + (2.0*PI/3.0)))
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) + (4.0*PI/3.0))) );

			Gamma_k_minus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) - (2.0*PI/3.0)))
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) - (4.0*PI/3.0))) );


			//Creating Lower Triangular Hamiltonian
			//Adding Chemical Potential, Zeeman, and Onsite Energy
			for(int spin=0;spin<2;spin++){
				for(int orb=0;orb<3;orb++){
					if(spin==0){
						Ham_(3*spin+orb,3*spin+orb) += one_complex*( -1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += one_complex*( 1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(spin==1){
						Ham_(3*spin+orb,3*spin+orb) += one_complex*( 1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += one_complex*( -1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(orb==1){
						Ham_(3*spin+orb,3*spin+orb) += -1.0*Parameters_.OnsiteE*one_complex;
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += 1.0*Parameters_.OnsiteE*one_complex;
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
			Ham_(6,3) += -1.0*OPs_[0];
			Ham_(7,4) += -1.0*OPs_[1];
			Ham_(8,5) += -1.0*OPs_[2];
			Ham_(9,0)  += 1.0*OPs_[0];
			Ham_(10,1) += 1.0*OPs_[1];
			Ham_(11,2) += 1.0*OPs_[2];

			//Creating Upper Traingular Hamiltonian (without diagonal)
			for(int row=0;row<2*HS_;row++){
				for(int col=row+1;col<2*HS_;col++){
					Ham_(row,col) = conj(Ham_(col,row));
				}
			}

			char Dflag='V'; //Since we are getting Eigenvectors
			Diagonalize(Dflag);

			for(int row=0;row<2*HS_;row++){
				Eigenvalues_[12*k_index + row]=eigs_[row]; //Each k=(k1,k2) have 12 Eigens
				for(int col=0;col<2*HS_;col++){
					Eigvectors_[12*k_index + col][row]=Ham_(row,col);
				}
			}
		}
	}
}




void Kspace_calculation_DL::Create_Kspace_Spectrum(){

	int k_index;
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
		k1x_val = (2.0*PI*k1)/(1.0*lx_);
		k1y_val = ((2.0*PI*k1)/(1.0*lx_))*(-1.0/sqrt(3));

		for(int k2=0;k2<ly_;k2++){
			k2x_val = 0.0;
			k2y_val = ((2.0*PI*k2)/(1.0*ly_))*(2.0/sqrt(3));

			k_index = Coordinates_.Ncell(k1,k2);
			Kx_values[k_index]=k1x_val + k2x_val;
			Ky_values[k_index]=k1y_val + k2y_val;

			Ham_.fill(0.0);

			Gamma_k_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*((2.0*PI*k1)/(1.0*lx_))) 
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)))) );

			Gamma_k_plus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) + (2.0*PI/3.0))) 
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) + (4.0*PI/3.0))) );

			Gamma_k_minus_[k1][k2] = one_complex*( 1.0 + exp(iota_complex*(((2.0*PI*k1)/(1.0*lx_)) - (2.0*PI/3.0)))
					+ exp(iota_complex*(((2.0*PI*k2)/(1.0*ly_)) - (4.0*PI/3.0))) );


			//Creating Lower Triangular Hamiltonian
			//Adding Chemical Potential, Zeeman, and Onsite Energy
			for(int spin=0;spin<2;spin++){
				for(int orb=0;orb<3;orb++){
					if(spin==0){
						Ham_(3*spin+orb,3*spin+orb) += one_complex*( -1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += one_complex*( 1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(spin==1){
						Ham_(3*spin+orb,3*spin+orb) += one_complex*( 1.0*Parameters_.Bmag - 1.0*Parameters_.mu );
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += one_complex*( -1.0*Parameters_.Bmag + 1.0*Parameters_.mu );
					}
					if(orb==1){
						Ham_(3*spin+orb,3*spin+orb) += -1.0*Parameters_.OnsiteE*one_complex;
						Ham_(HS_ +3*spin +orb,HS_ +3*spin +orb) += 1.0*Parameters_.OnsiteE*one_complex;
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
			Ham_(6,3) += -1.0*OPs_[0];
			Ham_(7,4) += -1.0*OPs_[1];
			Ham_(8,5) += -1.0*OPs_[2];
			Ham_(9,0)  += 1.0*OPs_[0];
			Ham_(10,1) += 1.0*OPs_[1];
			Ham_(11,2) += 1.0*OPs_[2];

			//Creating Upper Traingular Hamiltonian (without diagonal)
			for(int row=0;row<2*HS_;row++){
				for(int col=row+1;col<2*HS_;col++){
					Ham_(row,col) = conj(Ham_(col,row));
				}
			}

			char Dflag='V'; //Since we are getting Eigenvectors
			Diagonalize(Dflag);

			for(int row=0;row<2*HS_;row++){
				Eigenvalues_[12*k_index + row]=eigs_[row]; //Each k=(k1,k2) have 12 Eigens
				for(int col=0;col<2*HS_;col++){
					Eigvectors_[12*k_index + col][row]=Ham_(row,col);
				}
			}
		}
	}
}



void Kspace_calculation_DL::SelfConsistency(){

	/*	string File_Out_progress;
		File_Out_progress = "output_Kspace_SelfConsistency.txt";
		ofstream file_out_progress(File_Out_progress.c_str());

		cout<<"error targetted = "<<Parameters_.Convergence_Error<<endl;
		cout<<"Max iterations = "<<Parameters_.IterMax<<endl;

		OP_error_=10.0;
		int iter=0;
		while( (OP_error_>=Parameters_.Convergence_Error) && (iter<=Parameters_.IterMax)){
		Create_Kspace_Spectrum();
		Arranging_spectrum();
		mu_=Parameters_.mu; //Need to be added
		Get_new_OPs_and_error();
		Get_Energies();

		file_out_progress<<iter<<"   "<<OP_error_<<"   "<<E_class<<"   "<<E_quant<<"    ";
		for(int OP_no=0;OP_no<3;OP_no++){
		file_out_progress<<OPs_[OP_no].real()<<"    "<<OPs_[OP_no].imag()<<"    ";
		}

		file_out_progress<<endl;

	//REQUIRED REVIEW
	for(int i=0;i<OPs_.size();i++){
	OPs_[i]=Parameters_.alpha_OP*OPs_[i] + (1.0-Parameters_.alpha_OP)*OPs_new_[i];
	}

	iter++;
	}
	*/

	cout<<"mu = "<<mu_<<endl;
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
