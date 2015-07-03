#include "1D_Euler_FVM.h"


string Solution_Op_File,Initial_Solution_Op_File,Plot_Density,Plot_Pressure,Plot_Velocity,Plot_e,Plot_Temperature,Exact_Solution_File;

bool Is_Kinetic_Scheme;

int main()
{

// Choses the boundary Condition
Boundary_Type = 0; // 0 - Super sonic boundaries 1- subsonic boundary
Implement_Type = 2; // 2 - for extrapolation
Is_Kinetic_Scheme = false;
//Represents Types of Numerical dissipation and Flux scheme being used
/*0 - LLF,
 * 1 - PVU,
 * 2 - PVU_P,
 * 3 - MOVERS,
 * 4- ROE, 
 */ 
 /* 5- KFVS,
 * 6 - Steger_Warming,
 * 7 - Vanleer
 * 8 - AUSM,
 * 9 - Roe, 
 * 10 - HLL,
 * 11 - HLLC,
 * 12 - Rusanov Flux*/
Dissipation_Type =3; 
/* 0 - Central scheme with dissipation,
 * 1-  KFVS,
 * 2 - Steger_Warming,
 * 3 - Vanleer 
 * 4-  AUSM,
 * 5 - HLL*/
Flux_Scheme_Type = 0; 

//CFL Condition
CFL = 0.1;
// Initial Conditions
Test_Case =7;		// Toro test cases taken from TORO text book page number 249

Ghost_Cells = 2;	// Number of ghost cells
Physical_Cells = 98;	// Total number of Physical cells 
Length =1.0;		// Domain Length being used

// This function initializes the domain according to the test cases
Initalize();
cout<<"Initalization done\n";
// Invokes the solver.
Read_Exact_Solution_Input(Exact_Solution_File);
Solver();
// This function writes output file in the corresponding Solution folder
Write_Solution(Solution_Op_File);
// Plots the solution in .png files in the respective solution folders
Plot_Solution();
return 0;
}