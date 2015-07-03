#ifndef _Header_H
#define _Header_H


#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<cmath>
#include<typeinfo>
#include<iomanip>
#include <cstdio>
#include<cstdlib>
#include <algorithm>
#include"gnuplot_i.hpp"
// #include "gnuplot-iostream.h"


#define radian M_PI/180.0
#define gamma 1.4
#define inv_gamma 0.714285714286
#define inv_gamma_1 1.0/(gamma-1)
#define R 287
#define cp gamma*R/(gamma-1)
#define cv R/(gamma-1)

#define gamma_R gamma*R			// This is gamm*R used to calculate the speed of sound
// Constants defined for Sutherland's law for evaluating viscosity
#define	V_S_T  110.4								//	mu_ref=1.716e-5; T_ref=273.15;
#define	V_C1   1.45793265452e-06						// 	C1=(mu_ref/(pow(T_ref,1.5)))*(T_ref+S_T);
#define	T_S_T  194.4								//	K_ref=0.02414;	T_ref=273.15;
#define	T_C1   0.0025001353447							// 	C1=(K_ref/(pow(T_ref,1.5)))*(T_ref + S_T);
#define sleeplength 1

using namespace std;

extern vector<vector<double > >U,U_New,Eigen_Values;
extern vector<double> U_Cell,U_Cell_New,Left_Face_Average_Flux,Right_Face_Average_Flux,Pressure,Temperature,Density,Velocity,C,Internal_Energy;
extern vector<double> Left_Face_Dissipative_Flux,Right_Face_Dissipative_Flux,Cell_del_X,Cell_del_t,Net_Cell_Flux;
extern vector<double> G_P_L_Plus,G_P_L_Minus,G_P_R_Plus,G_P_R_Minus,G_C_L_Plus,G_C_L_Minus,G_C_R_Plus,G_C_R_Minus,Phi;
extern int Physical_Cells,Ghost_Cells,Total_Cells,Boundary_Cell_No,Interior_Cell_No,Test_Case;
extern int Boundary_Type,Implement_Type,Dissipation_Type,Flux_Scheme_Type;
extern double CFL,del_X,del_t,Eigen_Value_1,Eigen_Value_2,Eigen_Value_3,u_L,P_L,Rho_L,c_L,T_L,V_Mag_L,Inv_Rho,Terminating_Time;
extern double rho_error,rhou_error,rhoE_error,E_R,u,P,Rho,c,T,V_Mag,E,E_L,u_R,P_R,Rho_R,c_R,T_R,V_Mag_R;
extern double R_I_1,R_I_2,R_I_3,Min_delta_t,Min_delta_X,Sigma,Length,Lambda_Max,Lambda_Min;
extern string Solution_Op_File,Initial_Solution_Op_File,Plot_Temperature,Plot_Density,Plot_Pressure,Plot_Velocity,Plot_e,Exact_Solution_File;
extern Gnuplot g1,g2,g3,g4,g5,g;
extern vector<double> Exact_Density,Exact_Velocity,Exact_Pressure,Exact_e;
extern bool Is_Kinetic_Scheme;

// Functions being used in the code
void Allocate_Memory();
void Initalize();
// void Initialize(string &);
void Calculate_Convective_Flux();
void Calculate_Primitive_Variables();
void Calculate_Primitive_Variables(vector<double> &, int & );

//Functions of Numerical schemes in the code
void LLF_Dissipation(int &);
void PVU_Dissipation(int &);
void PVU_P_Dissipation(int &);
void MOVERS_Dissipation(int &);
void KFVS(int & );
void Steger_Warming(int & );
void Vanleer(int & );
void AUSM(int & );
void ROE(int &);
void HLL(int & );
void HLLC(int & );
void Rusanov_Dissipation(int &);

void Find_Error_and_Update();

void Super_Sonic_Inlet_Boundary_Condition(int &);
void Super_Sonic_Exit_Boundary_Condition(int & );
void Sub_Sonic_Inlet_Boundary_Condition(int &);
void Sub_Sonic_Exit_Boundary_Condition(int & );
void Exit_Boundary_Condition();

void Evaluate_Riemann_Invariants();
void Evaluate_Sigma();
void Reset_Variables();
void Solver();
void Evaluate_Average_Flux(int & );
void Evaluate_Dissipation_Flux(int & , int & );
void Evaluate_Net_Cell_Flux(int & );
void Evaluate_Sigma();
void Evaluate_Time_Evolution();
// Test cases being tested in 1D solver
void Sod_Test_Case_1();
void Sod_Test_Case_2();
void Toro_Test_Case_1();
void Toro_Test_Case_2();
void Toro_Test_Case_3();
void Toro_Test_Case_4();
void Toro_Test_Case_5();
void Steady_Shock();
void Steady_Contact();
void Calculate_Conserved_Variables(double & , double & , double & );
void Boundary_Conditions(int & ,int &);
void Write_Solution(string &);
void Plot_Solution();
void Plot_Solution(vector<vector<double> >);
void Slope_Limiter(vector<double> & );
void Min_Mod(int & ,vector<double> & ,int &);
void File_Name();
void Directory_Name();
void Read_Exact_Solution_Input(string & );
void Error_Function(double & , double & );
void Mod(double &, double &);
double & Maximum(double &,double &,double &);
double & Minimum(double &,double &,double &);
int  Sign(double & );


#endif  //#ifndef _Header_H