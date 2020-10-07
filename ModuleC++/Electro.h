#pragma once
#include "PhMed.h"
#include"Vector.h"
#include"GridBlock.h"
#include <cmath>
#include <stdio.h>
#include "Stdafx.h"
#include "NR/nr.h"
#include <iomanip>
#pragma comment(linker, "/STACK:16777216")
#include <deque>
#include <iostream>
#include <fstream>
	



using namespace std;


class Electro
{
	
public:
	Electro();
	~Electro();
	
	const int nz = 320;// количество ячеек по оси Z (изменив, нужно изменить массив в функции Numerical_realization())
	const int nz_clay = 100; // количество ячеек, отведенных на глину (добавляется сверху и снизу пласта)
	const int NR = 125;  // количество ячеек по оси X --> зависит от разбиения в гидродинамике   
	const int nw = 9; //количество ячеек отведенных на скважину (от зонда до стенки) по оси X
	const int NumCells = NR + nw; 
	const int NP = (nz - 2) * (NumCells - 1); //общее количество ячеек без границ
	const int npsize = 5 * NP - 2 * NumCells;//значение npsize - количество ненулевых элементов в пятидиагональной матрице коэффициентов



	//для теста
	bool Test = true; //тестовый режим - true, расчёт по ГДИС - false


	const int NumLayersX = 1;//однослойные палетки - 1, двухслойные - 2, трехслойные - 3
	const int NumLayersZ = 2;//количество разнопроводимых слоев в пласте
	//

	double rw = 0.1; // in m//радиус скважины
	double h = 40.0; // in m//толщина пласта (с глиной)
	double top = 2000.0; // the first layer top


	double *spatial_dr;
	double *distance;

	double *ksi; //логарифмическая координата (log(r))
	double *d_ksi;

	double *dz; // sixe of block (z-axis)
	double *z;
	
	double ro_clay = 2.0; // Om*m//УЭС глины
	double ro_res = 100.0; // Om*m//УЭС резервуара
	double ro_ext = 100.0; // Om*m
	double ro_mud; // Om*m //УЭС бур.раствора (закачиваемой воды)

	double *elect_cond; //reservoir conductivity [Om-1]
	//double *elect_cond_nz;
	double *elect_cond_nz_cent;
	double *elect_cond_nx_cent;


	// Initial data
	bool Existing_clay = false;
	bool Pointed_source = true; // true- poited source; false - distributed source;
	bool stand_tool = false;  //   true - normal zonde; false - lateral zonde
	bool Source_b = true; // true - source is controlled by electrical force; false - source is controlled by potential;
	
	// размеры зондов в м
	double *AM;
	double *MN;
	double *AO;
	double *AN;
	//

	int num_zonde;//количество зондов
	int num_A;//количество точек измерений (положений точки А)
	double *depth_A; //глубина положения точки А = top + depth_A
	double *ro_calc;
	double *depth_k;

	double R_zonde = 0.038; // in m //радиус зонда
	double I = 0.5; // сила тока в амперах

	double k_stan = 2; // - multiple for investigation radius of standard log tool(from 2 to 3)
	double k_grad = 1; // - multiple for investigation radius of standard log tool(1)

	double *K_norm;
	double *K_grad;

	int *nz1;
	double *z_A;
	double *z_O;
	double *z_M;
	double *z_N;
	int *nz_N;
	int *nz_O;
	int *nz_M;

	
	double *alfa_cent;
	double *betta_cent;
	double *zu;
	double *du;

	double *Ci;
	double *Wi;
	double *Ei;
	double *Ni;
	double *Si;

	int iter;
	double err;
	int itmax = 5000; // max количество итераций - сильно влияет на результат!!!!!!! (меньше 100 лучше не допускать)
	int itol = 1; // варьируется в пределах 1-3; уменьшение itol => сглаживание к концу значений
	double tol = 0.00001; // толеранс

	double U_b = 0.0;
	double U_ext = 0;
	double U_int = 30.0;
	double U_bound = 0.0;/////////

	Vec_INT *ija;
	Vec_DP *sa;

	//double *b;			// compresource
	double *U; // Potential of all zone (with conditions)
	double *x_m;

	void Electro::Resistivity()
	{    
		extern double *Rr;

		for (int lineIdx = 0; lineIdx < nw; lineIdx++)
		{
			spatial_dr[lineIdx] = (rw - R_zonde) / nw;
		}
		for (int lineIdx = nw; lineIdx < (NumCells - 1); lineIdx++)
		{
			int ind = lineIdx - nw;
			spatial_dr[lineIdx] = Rr[ind]/100;
		}


		if (Test==true)
		{ 
			Test_Resist();
		}
		else
		{
			Darley_Archie();
		}

		NewGrid();
		Intialization();
		Source();
		Numerical_realization();
		
	}

	void Electro::Read_Write()
	{ 
		num_zonde = 0;
		FILE *f_grid_zonde;
		char buf_zonde[1024];
		fopen_s(&f_grid_zonde, "ZondeAM.dat", "r");
		fgets(buf_zonde, 1024, f_grid_zonde);
		for (;;)
		{
			if (fgets(buf_zonde, 1024, f_grid_zonde) == NULL)
				break;
			num_zonde++;
		}
		fclose(f_grid_zonde);

		num_A = 0;
		FILE *f_grid_A;
		char buf_A[1024];
		fopen_s(&f_grid_A, "location_A.dat", "r");
		fgets(buf_A, 1024, f_grid_A);
		for (;;)
		{
			if (fgets(buf_A, 1024, f_grid_A) == NULL)
				break;
			num_A++;
		}
		fclose(f_grid_A);


		AM = new double[(num_zonde)];
		MN = new double[(num_zonde)];
		depth_A = new double[(num_A)];

		// запись в массив сетки
		FILE *f_AM;
		char buf_AM[1024];
		

		fopen_s(&f_AM, "ZondeAM.dat", "r");
		fgets(buf_AM, 1024, f_AM);
		for (int k = 0; k < num_zonde; k++)
		{
			fgets(buf_AM, 1024, f_AM);
			sscanf(buf_AM, "%lg", &AM[k]);
		}
		fclose(f_AM);

		FILE *f_MN;
		char buf_MN[1024];


		fopen_s(&f_MN, "ZondeMN.dat", "r");
		fgets(buf_MN, 1024, f_MN);
		for (int k = 0; k < num_zonde; k++)
		{
			fgets(buf_MN, 1024, f_MN);
			sscanf(buf_MN, "%lg", &MN[k]);
		}
		fclose(f_MN);


		fopen_s(&f_grid_A, "location_A.dat", "r");
		fgets(buf_A, 1024, f_grid_A);
		for (int q = 0; q < num_A; q++)
		{
			fgets(buf_A, 1024, f_grid_A);
			sscanf(buf_A, "%lg", &depth_A[q]);
		}
		fclose(f_grid_A);

	}

	void Electro::Darley_Archie()
	{
		int a = 1;
		int m_const = 2;
		int n_const = 2;
		int temperature = 50; // reservoir temperature
		double ro_w = 1.228; //density of water (g/l)
		double Swc = 0.2;
		double Sw;
		double poro;
		double Cw;
		double wat_Resist, res_Resist;
		extern double m_El;
		typedef CVector<3> CV1;
		extern CV1 vars_el[PhMed::nMedia];

		for (int j = 0; j < nz; j++)
		{
			for (int lineIdx = nw; lineIdx < NumCells; lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				int indVar = lineIdx - nw;

				Sw = vars_el[PhMed::fracture][indVar][PhMed::w_sat] * (1 - Swc) + Swc;
				poro = m_El / (1 - Swc);
				Cw = vars_el[PhMed::fracture][indVar][PhMed::s_con] * ro_w * pow(10, 6); //convert from (g/g) to (ppm)

				if (Cw != 0)
				{
					wat_Resist = (0.0123 + 3647.5 / pow(Cw, 0.955))*(82.0 / (1.8*temperature + 39.0));
					res_Resist = a * pow(poro, (-m_const)) * wat_Resist * pow(Sw, (-n_const));
					elect_cond[indX] = 1.0 / res_Resist; // проводимость


					if (indX == nw)
					{
						ro_mud = res_Resist; 

					}
				}
				else
				{
					res_Resist = ro_res;
					elect_cond[indX] = 1.0 / res_Resist;
				}
			}
		}

		for (int j = 0; j < nz; j++)
		{
			for (int lineIdx = 0; lineIdx < nw; lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				elect_cond[indX] = 1 / ro_mud;
			}
		}

		for (int j = 0; j < nz_clay; j++)
		{
			for (int lineIdx = nw; lineIdx < NumCells; lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				elect_cond[indX] = 1 / ro_clay;
			}
		}

		for (int j = nz - nz_clay; j < nz; j++)
		{
			for (int lineIdx = nw; lineIdx < NumCells; lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				elect_cond[indX] = 1 / ro_clay;
			}
		}
	}
	void Electro::Test_Resist()
	{
		double ro_Res, ro_Well;
		double *ro_f;
		
		ro_f = new double[NumLayersZ];

		/*FILE *f_Resist;
		char buf_Resist[1024];
		fopen_s(&f_Resist, "C:\\Users\\Olga\\Desktop\\Работа ИПНГ\\новая схема\\Реализация (новая схема)\\Input\\Resist_Test.in", "r");
		fgets(buf_Resist, 1024, f_Resist);
		for (int k = 0; k < NumLayersZ; k++)
		{
			fgets(buf_Resist, 1024, f_Resist);
			sscanf(buf_Resist, "%lg", &ro_f[k]);
		}
		fclose(f_Resist);

		cout << ro_f[0] << "    " << ro_f[1] << endl;*/
		ro_f[0] = 1;
		ro_f[1] = 5;

/*
		ro_Res = ro_f[0];*/

		if (NumLayersX == 1)
		{
			if (NumLayersZ == 2)
			{
				int nz_1 = 150; //ячеек на первый слой (толщина в метрах * 4)
				for (int j = 0; j < nz_1; j++)
				{
					for (int lineIdx = 0; lineIdx < NumCells; lineIdx++)
					{
						int indX = j * NumCells + lineIdx;
						elect_cond[indX] = 1 / ro_f[0];
					}
				}
				for (int j = nz_1; j < nz; j++)
				{
					for (int lineIdx = 0; lineIdx < NumCells; lineIdx++)
					{
						int indX = j * NumCells + lineIdx;
						elect_cond[indX] = 1 / ro_f[1];
					}
				}
			}
		}
		
		delete[] ro_f;
		/*if (NumLayersX == 2)
		{
			ro_Well = ro_f[1];

			for (int j = 0; j <= nzi; j++)
			{
				for (int lineIdx = 0; lineIdx <= nw; lineIdx++)
				{
					int indX = j * (NumCells + 1) + lineIdx;
					elect_cond[indX] = 1 / ro_Well;
				}
			}
			for (int j = 0; j <= nzi; j++)
			{
				for (int lineIdx = (nw + 1); lineIdx <= NumCells; lineIdx++)
				{
					int indX = j * (NumCells + 1) + lineIdx;
					elect_cond[indX] = 1 / ro_Res;
				}
			}
		}
		if (NumLayersX == 3)
		{}*/
	}

	void Electro::NewGrid()
	{
		distance[0] = R_zonde;
		for (int i = 0; i < (NumCells - 1); i++)
		{
			distance[i + 1] = distance[i] + spatial_dr[i];
		}

		ksi[0] = log(distance[0] / distance[NumCells - 1]);
		for (int i = 0; i < (NumCells - 1); i++)
		{
			ksi[i + 1] = log(distance[i + 1] / distance[NumCells - 1]);
			d_ksi[i] = ksi[i + 1] - ksi[i];
		}


		for (int j = 0; j < (nz - 1); j++)
		{
			for (int lineIdx = 0; lineIdx < NumCells; lineIdx++)
			{
				int indZ = (j + 1) * NumCells + lineIdx;
				int indX = j * NumCells + lineIdx;
				elect_cond_nz_cent[indX] = (2 * elect_cond[indZ] * elect_cond[indX]) / (elect_cond[indZ] + elect_cond[indX]);
			}
		}

		for (int j = 0; j < nz; j++)
		{
			for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
			{
				int indZ = j * NumCells + lineIdx + 1;
				int indX = j * NumCells + lineIdx;
				elect_cond_nx_cent[indX] = (2 * elect_cond[indZ] * elect_cond[indX]) / (elect_cond[indZ] + elect_cond[indX]);
			}
		}


	}	
	void Electro::Intialization()
	{
		for (int k = 0; k < num_zonde; k++)
		{
			if (stand_tool == true)
			{
				//AO=AM/2;
				AO[k] = AM[k];
				AN[k] = AM[k] + MN[k];
				K_norm[k] = 4 * 3.14 * AO[k]; //AM[k];??

			}
			else
			{
				AO[k] = AM[k] + MN[k] / 2;
				AN[k] = AM[k] + MN[k];
				K_grad[k] = 4 * 3.14 * AM[k] * AN[k] / MN[k];
			}
		}	
	}
	void Electro::Source()
	{
		z[0] = top;
		for (int j = 0; j < (nz - 1); j++)
		{
			dz[j] = h / nz;
			z[j + 1] = z[j] + dz[j];
		}
		
		if (Pointed_source == true)
		{
			for (int q = 0; q < num_A; q++)
			{
				for (int k = 0; k < num_zonde; k++)
				{
					int ind_k;
					ind_k = k * (num_A) + q;
					z_A[q] = top + depth_A[q];
					z_O[ind_k] = z_A[q] + AO[k];
					z_M[ind_k] = z_A[q] + AM[k];
					z_N[ind_k] = z_A[q] + AN[k];
				}
			}

			for (int q = 0; q < num_A; q++)
			{
				nz1[q] = 0;
				for (int k = 0; k < num_zonde; k++)
				{
					int ind_k;
					ind_k = k * (num_A) + q;
					for (int j = 0; j < nz; j++)//////////////////////////////////????????????????????????????
					{
						if ((z[j] <= z_A[q]) && (z_A[q] <= z[j + 1]))
						{
							nz1[q] = j;
							
						}
					}
					for (int j = 0; j < nz; j++)
					{
						if ((z[j] <= z_O[ind_k]) && (z_O[ind_k] <= z[j + 1]))
						{
							nz_O[ind_k] = j;							
						}
					}

					for (int j = 0; j < nz; j++)
					{
						
						if ((z[j] <= z_N[ind_k]) && (z_N[ind_k] <= z[j + 1]))
						{
							nz_N[ind_k] = j;
						}
					}


					for (int j = 0; j < nz; j++)
					{
						//int ind_k=k*(num_A+1)+q;
						if ((z[j] <= z_M[ind_k]) && (z_M[ind_k] <= z[j + 1]))
						{
							nz_M[ind_k] = j;
							if (nz_M[ind_k] == nz1[q])
							{
								nz_M[ind_k] = nz1[q] + 1;
							}
							if (nz_M[ind_k] == nz_N[ind_k])
							{
								nz_N[ind_k] = nz_M[ind_k] + 1;
							}
						}
					}
				}
				cout << nz1[q] << "       ";
				
			}
			cout << endl;
		}
	}

	void Electro::Numerical_realization()
	{
		Vec_DP bi(NP), x(NP);

		// Calculation of alfa and betta coefficients
		for (int j = 0; j < (nz - 1); j++)
		{
			for (int lineIdx = 0; lineIdx < NumCells; lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				betta_cent[indX] = elect_cond_nz_cent[indX] / dz[j];
				//cout << elect_cond_nz_cent[indX] << endl;
			}
		}

		for (int j = 0; j < nz; j++)
		{
			for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
			{
				int indX = j * NumCells + lineIdx;
				alfa_cent[indX] = elect_cond_nx_cent[indX] / d_ksi[lineIdx];
			}
		}

		for (int i = 0; i < (NumCells - 1); i++)
		{
			zu[i] = 1 / (pow(distance[i + 1], 2) - pow(distance[i], 2));
		}

		for (int j = 0; j < (nz - 1); j++)
		{
			du[j] = 1 / (2 * dz[j]);

		}

		/*for (int i = 0; i < NumCells; i++)
		{
			int j = 1;
			zu[0] = 2*exp(-2 * ksi[0])  / (d_ksi[1] + d_ksi[0]);
			zu[i + 1] = 2*exp(-2 * ksi[i]) / (d_ksi[i + 1] + d_ksi[i]);
			du[0] = 2 / (dz[0] + dz[1]);
			du[i + 1] = 2 / (dz[j] + dz[j + 1]);
		}*/

				// FULLING OF THE INTERNAL POINTS
		

		for (int q = 0; q < num_A; q++)
		{
			for (int j = 0; j < nz; j++)
			{
				for (int i = 0; i < NumCells; i++)
				{
					int IndU_int = j * NumCells + i;
					x_m[IndU_int] = U_bound;
				}
			}





			for (int j = 0; j < nz; j++)
			{
				for (int i = 0; i < NumCells; i++)
				{
					int IndU = j * NumCells + i;
					U[IndU] = U_b;
				}
			}

			//// internal potential matrix
			//for (int j = 0; j < nz; j++)
			//{
			//	for (int i = 0; i < NumCells; i++)
			//	{
			//		int IndU_int = j * NumCells + i;
			//		x_m[IndU_int] = U_int;
			//	}
			//}

			for (int j = 0; j < nz; j++)
			{
				int IndU = j * NumCells + NumCells;
				U[IndU] = U_ext;
			}

			
			for (int j = 0; j < (nz - 2); j++)
			{
				for (int i = 0; i < (NumCells - 1); i++)
				{
					int IndU_int = j * (NumCells - 1) + i;
					x[IndU_int] = U_int;
				}
			}
		

			for (int j = 1; j < (nz - 1); j++)
			{
				for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
				{
					int indCentZ = (j - 1) * NumCells + lineIdx;
					int indX = j * NumCells + lineIdx;
					Ni[indX] = 0;
					Si[indX] = 0;
					Wi[indX] = 0;
					Ei[indX] = 0;
					Ci[indX] = 0;

					bi[indX] = 0;
					if ((lineIdx == 0)&&(j == nz1[q]))
					{
						bi[indX] = -I / (2 * 3.14 * dz[j] * (pow(distance[1], 2) - pow(distance[0], 2)));
					}
					

					if (j > 1)
					{
						Si[indX] = du[j] * (betta_cent[indCentZ]);
						Ci[indX] += -du[j] * (betta_cent[indCentZ]);
					}
					else
					{
						Ci[indX] += -du[j] * (betta_cent[indCentZ]);
						bi[indX] += -du[j] * (betta_cent[indCentZ]) * U[0 * NumCells + lineIdx];
					}

					if (j < (nz - 2))
					{
						Ni[indX] = du[j] * (betta_cent[indX]);
						Ci[indX] += -du[j] * (betta_cent[indX]);
					}
					else
					{
						Ci[indX] += -du[j] * (betta_cent[indX]);
						bi[indX] += -du[j] * (betta_cent[indX]) * U[(nz - 1) * NumCells + lineIdx];
					}

					if (lineIdx < (NumCells - 2))
					{
						Ei[indX] = zu[lineIdx] * (alfa_cent[indX]);
						Ci[indX] += -zu[lineIdx] * (alfa_cent[indX]);
					}
					else
					{
						Ci[indX] += -zu[lineIdx] * (alfa_cent[indX]);
						bi[indX] = bi[indX] - zu[lineIdx] * (alfa_cent[indX]) * U[j*NumCells + (NumCells - 1)];
					}

					if (lineIdx > 0)
					{
						Wi[indX] = zu[lineIdx] * (alfa_cent[indX - 1]);
						Ci[indX] += -zu[lineIdx] * (alfa_cent[indX - 1]);
					}
				}
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			// COMPRESSION OF COEFFIECENTS:
			for (int j = 1; j < (nz - 1); j++)
			{
				for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
				{
					int indXi = j * NumCells + lineIdx;
					int indX = (j - 1)*(NumCells - 1) + lineIdx;
					Ni[indX] = Ni[indXi];
					Si[indX] = Si[indXi];
					Wi[indX] = Wi[indXi];
					Ci[indX] = Ci[indXi];
					Ei[indX] = Ei[indXi];
					bi[indX] = bi[indXi];
				}
			}
			//Вывод диагоналей в консоль (для проверки)
		   /*for (int j = 1; j < (nz - 1); j++)
		   	for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
		   	{
		   		int indX = (j - 1)*(NumCells - 1) + lineIdx;
		   		cout << indX << " ";
		   		cout << Si[indX] << " ";
		   		cout << Wi[indX] << " ";
		   		cout << Ci[indX] << " ";
		   		cout << Ei[indX] << " ";
		   		cout << Ni[indX] << " ";
		   		cout << endl;
		   	}*/

			
			/*if (n > 10 || n < 1)
				throw "выход за границы массива";*/

			
			//Запись пятидиагональной матрицы в сжатом виде ( формат MSR ) ; данный алгоритм - требуемый эквивалент алгоритма функции NR::sprsin(..)
				int k;
				for (int j = 0; j < NP; j++) (*sa)[j] = Ci[j];
				(*ija)[0] = NP + 1;
				k = NP;
				const int newNC = NumCells - 1;//new NumCells
				for (int i = 0; i < NP; i++) {
					if (i == 0)
					{
						++k;
						(*sa)[k] = Ei[i];
						(*ija)[k] = i + 1;
						++k;
						(*sa)[k] = Ni[i];
						(*ija)[k] = i + newNC;

					}
					if ((i > 0) && (i < newNC))
					{
						++k;
						(*sa)[k] = Wi[i];
						(*ija)[k] = i - 1;

						++k;
						(*sa)[k] = Ei[i];
						(*ija)[k] = i + 1;

						++k;
						(*sa)[k] = Ni[i];
						(*ija)[k] = i + newNC;

					}
					if ((i >= newNC) && (i <= NP - newNC - 1))
					{
						++k;
						(*sa)[k] = Si[i];
						(*ija)[k] = i - newNC;

						++k;
						(*sa)[k] = Wi[i];
						(*ija)[k] = i - 1;

						++k;
						(*sa)[k] = Ei[i];
						(*ija)[k] = i + 1;

						++k;
						(*sa)[k] = Ni[i];
						(*ija)[k] = i + newNC;
					}
					if ((i > NP - newNC - 1) && (i < NP - 1))
					{
						++k;
						(*sa)[k] = Si[i];
						(*ija)[k] = i - newNC;

						++k;
						(*sa)[k] = Wi[i];
						(*ija)[k] = i - 1;

						++k;
						(*sa)[k] = Ei[i];
						(*ija)[k] = i + 1;
					}
					if (i == NP - 1)
					{
						++k;
						(*sa)[k] = Si[i];
						(*ija)[k] = i - newNC;

						++k;
						(*sa)[k] = Wi[i];
						(*ija)[k] = i - 1;
					}
					(*ija)[i + 1] = k + 1;
				}
		


			//Solves A · x = b
			NR::linbcg(bi, x, itol, tol, itmax, iter, err);



			// MOVED OF COEFFIECENTS:
			for (int j = 1; j < (nz - 1); j++)
			{
				for (int lineIdx = 0; lineIdx < (NumCells - 1); lineIdx++)
				{
					int indXi = j * NumCells + lineIdx;
					int indX = (j - 1)*(NumCells - 1) + lineIdx;
					x_m[indXi] = x[indX];
					U[indXi] = x_m[indXi];
				}
			}

			//Расчёт КС
			int ind_k; int l;
			for (k = 0; k < num_zonde; k++)
			{
				ind_k = k * (num_A)+q;
				if (stand_tool == false)				// if there is lateral device
				{
					int indX_M; int indX_N;
					indX_M = nz_M[ind_k] * NumCells;
					indX_N = nz_N[ind_k] * NumCells;
					ro_calc[ind_k] = K_grad[k] * (U[indX_M] - U[indX_N]) / I; // [K_grad]=m => [ro_calc] = Om*m; К = 4 * 3.14*AM[k] * AN[k] / MN[k];
					depth_k[q] = z_O[ind_k];
					cout << ro_calc[ind_k] << endl;
				}
				else
				{
					int indX_O;
					indX_O = nz_O[ind_k] * NumCells;
					ro_calc[ind_k] = K_norm[k] * (U[indX_O] - U_b) / I; // [K_norm]=m => [ro_calc] = Om*m;
					depth_k[q] = z_O[ind_k];
				}
			}
			
		}

		/*ofstream FILE;
		FILE.open("Output\\Electro.out");
		for (int ind_k = 0; ind_k < num_A; ind_k++)
		{
			FILE << ro_calc_k[;
			FILE << endl;
		}
		FILE.close();*/

	}



};


