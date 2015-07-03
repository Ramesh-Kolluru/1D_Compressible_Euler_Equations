#include "1D_Euler_FVM.h"


// Boundary conditions are being implemented by Riemann Invariants as per Rohan and Venkat's Report also by simple extrapolation methods


void Boundary_Conditions(int & Boundary_Type,int & I_Type)
{
//     cout<<Boundary_Type<<endl;
    switch(Boundary_Type)
    {
        case 0:
            Super_Sonic_Exit_Boundary_Condition(I_Type);
            Super_Sonic_Inlet_Boundary_Condition(I_Type);
            break;
        case 1:
            Sub_Sonic_Exit_Boundary_Condition(I_Type);
            Sub_Sonic_Inlet_Boundary_Condition(I_Type);
            break;
        default :
            Super_Sonic_Exit_Boundary_Condition(I_Type);
            Super_Sonic_Inlet_Boundary_Condition(I_Type);
            break;
    }
}

void Super_Sonic_Inlet_Boundary_Condition(int & Implementaion_Type)
{
    Boundary_Cell_No = 0,Interior_Cell_No =1;
//     cout<<"In super sonic Inlet Boundary\t"<<Implementaion_Type<<endl;
    switch(Implementaion_Type)
    {
        case 1: 
            //Using Riemann Invariants
            break;
        case 2: // Using Simple extrapolation
            U[Interior_Cell_No][0] = U[Boundary_Cell_No][0];
            U[Interior_Cell_No][1] = U[Boundary_Cell_No][1];
            U[Interior_Cell_No][2] = U[Boundary_Cell_No][2];
	    Calculate_Primitive_Variables(U[Boundary_Cell_No], Interior_Cell_No);
            break;
        default :
            break;
    }
//     cout<<"Boundary Conditions Extrapolated\n";
    
}


void Super_Sonic_Exit_Boundary_Condition(int & Implementaion_Type)
{
    Boundary_Cell_No = Total_Cells-1,Interior_Cell_No =Total_Cells-2;
//     cout<<"In super sonic Exit Boundary\t"<<Implementaion_Type<<endl;
    switch(Implementaion_Type)
    {
        case 1: 
            //Using Riemann Invariants
            
            break;
        case 2: // Using Simple extrapolation u[i-2] = U[i-1]
            U[Boundary_Cell_No][0] = U[Interior_Cell_No][0];
            U[Boundary_Cell_No][1] = U[Interior_Cell_No][1];
            U[Boundary_Cell_No][2] = U[Interior_Cell_No][2];
	    Calculate_Primitive_Variables(U[Boundary_Cell_No], Boundary_Cell_No);
            break;
        default :
            break;
    }
    

}

void Sub_Sonic_Inlet_Boundary_Condition(int & Implementaion_Type)
{
    
}


void Sub_Sonic_Exit_Boundary_Condition(int & Implementaion_Type)
{
}