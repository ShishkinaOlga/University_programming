#include "Electro.h"


typedef CVector<3> CV;
CV vars_el[PhMed::nMedia]; //матрица значений из ForwardProblem.h
double m_El;//эффективная пористость
double *Rr; // dr

//вектора для сжатого хранения матрицы коэффициентов 
Vec_INT  *ija_p;
Vec_DP *sa_p;

Electro::Electro()
{
	Read_Write();

	AO = new double[num_zonde];
	AN = new double[num_zonde];
	K_norm = new double[num_zonde];
	K_grad = new double[num_zonde];
	nz1 = new int[num_A];
	z_A = new double[num_A];
	z_O = new double[(num_zonde)*(num_A)];
	z_N = new double[(num_zonde)*(num_A)];
	z_M = new double[(num_zonde)*(num_A)];
	nz_O = new int[(num_zonde)*(num_A)];
	nz_N = new int[(num_zonde)*(num_A)];
	nz_M = new int[(num_zonde)*(num_A)];


	ro_calc = new double[(num_zonde)*(num_A)];
	depth_k = new double[num_A];
	spatial_dr = new double[NumCells - 1];
	distance = new double[NumCells];
	ksi = new double[NumCells];
	d_ksi = new double[NumCells - 1];
	dz = new double[nz - 1];
	z = new double[nz];
	elect_cond = new double[nz*NumCells];
	elect_cond_nz_cent = new double[(nz - 1)*NumCells]; 
	elect_cond_nx_cent = new double[nz*NumCells]; 
	betta_cent = new double[(nz - 1)*NumCells];
	alfa_cent = new double[nz*NumCells];
	zu = new double[NumCells];
	du = new double[nz];
	Ci = new double[nz*NumCells];
	Ei = new double[nz*NumCells];
	Wi = new double[nz*NumCells];
	Ni = new double[nz*NumCells];
	Si = new double[nz*NumCells];
	U = new double[(nz + 1)*(NumCells + 1)];//
	x_m = new double[(nz + 1)*(NumCells + 1)];//

	ija = new Vec_INT(npsize + 1);//
	sa = new Vec_DP(npsize + 1);//

	int t1=ija->size();
	int t2=sa->size();

}


Electro::~Electro()
{
	delete[] spatial_dr;
	delete[] distance;
	delete[] ksi;
	delete[] d_ksi;
	delete[] dz;
	delete[] z;
	delete[] elect_cond;
	delete[] elect_cond_nz_cent;
	delete[] elect_cond_nx_cent;
	delete[] AM;
	delete[] MN;
	delete[] AO;
	delete[] AN;
	delete[] depth_A;
	delete[] K_norm ;
	delete[] K_grad ;
	delete[] nz1 ;
	delete[] z_A ;
	delete[] z_O ;
	delete[] z_N ;
	delete[] z_M ;
	delete[] nz_O;
	delete[] nz_N;
	delete[] nz_M;
	delete[] betta_cent;
	delete[] alfa_cent;
	delete[] zu;
	delete[] du;
	delete[] Ci;
	delete[] Ei;
	delete[] Wi;
	delete[] Ni;
	delete[] Si;
	delete[] U;
	delete[] x_m;
	if (ija != NULL)
	{
		delete ija;
		ija = NULL;
	}
	if (sa != NULL) 
	{
		delete sa;
		sa = NULL;
	}

	delete[] ro_calc;
	delete[] depth_k;

}

