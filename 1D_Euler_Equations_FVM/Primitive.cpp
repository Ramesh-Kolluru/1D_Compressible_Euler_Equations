#include "1D_Euler_FVM.h"


void Calculate_Primitive_Variables()
{
    double Inv_Rho =0.0;
//     Reset_Variables();
    for(int Cell_No=0;Cell_No<Total_Cells;Cell_No++)
    {
        
        Rho = U[Cell_No][0];                        // Density
        Inv_Rho = 1.0/Rho;                          //1/rho;
        u = U[Cell_No][1]*Inv_Rho;                  // u = (rho*U/rho)
        V_Mag = 0.5*u*u;                            // 
        P = (gamma -1)*(U[Cell_No][2] - Rho*V_Mag); // Rho*E = (P/(gamma-1)) + (1/2)*u*u --------- Pressure is evaluated from this expression
        T = P/(R*Rho);                              // Temperature from EOS P = rho*R*T------------- Temperature
        c = sqrt(gamma*P*Inv_Rho);                 // C = sqrt(gamma*P/rho) ------------- speed of sound
       
        
        Eigen_Value_1 = fabs(u)+c;                        // u + a 
        Eigen_Value_2 = fabs(u)-c;                        // u - a
        Eigen_Value_3 = fabs(u);                          // u
        
        Density[Cell_No] = Rho;                     
        Velocity[Cell_No] = u;
        Pressure[Cell_No] = P;
        Temperature[Cell_No] = T;
        C[Cell_No] = c;
	Internal_Energy[Cell_No] = P/(Rho*(gamma - 1.0));		// Rho*Cv*T
        
        Eigen_Values[Cell_No][0] = (Eigen_Value_1);
        Eigen_Values[Cell_No][1] = (Eigen_Value_2);
        Eigen_Values[Cell_No][2] = (Eigen_Value_3);
    }
}
    

void Calculate_Primitive_Variables(vector<double> & U_1, int & Cell_Index)
{
    double Inv_Rho =0.0;
    
        Rho = U_1[0];                      // Density
        Inv_Rho = 1.0/Rho;                 // 1/rho;
        u = U_1[1]*Inv_Rho;                // u = (rho*U/rho)
        V_Mag = 0.5*u*u;                   // 
        P = (gamma-1)*(U_1[2] - Rho*V_Mag);// Rho*E = (P/(gamma-1)) - (1/2)*Rho*u*u - Pressure is evaluated from this expression
        T = P/(R*Rho);                     // Temperature from EOS P = rho*R*T------------- Temperature
        c = sqrt(gamma*P*Inv_Rho);         // C = sqrt(gamma*P/rho) ------------- speed of sound
       

//         cout<<Cell_Index<<"\t"<<Rho<<"\t"<<u<<"\t"<<P<<"\t"<<T<<"\t"<<c<<endl;
        
        Eigen_Value_1 = fabs(u)+c;                        // u + a 
        Eigen_Value_2 = fabs(u)-c;                        // u - a
        Eigen_Value_3 = fabs(u);                          // u
        
        Density[Cell_Index] = Rho;
        Velocity[Cell_Index] = u;
        Pressure[Cell_Index] = P;
        Temperature[Cell_Index] = T;
        C[Cell_Index] = c;
	Internal_Energy[Cell_Index] = P/(Rho*(gamma - 1.0));		// Rho*Cv*T
        
        U[Cell_Index][0] = U_1[0];
        U[Cell_Index][1] = U_1[1];
        U[Cell_Index][2] = U_1[2];
        
//         U_New[Cell_Index][0] = U_1[0];
//         U_New[Cell_Index][1] = U_1[1];
//         U_New[Cell_Index][2] = U_1[2];
        
        Eigen_Values[Cell_Index][0] = (Eigen_Value_1);
        Eigen_Values[Cell_Index][1] = (Eigen_Value_2);
        Eigen_Values[Cell_Index][2] = (Eigen_Value_3);
}
    

// This function evaluates the conservative vector U based on Pressure, Velocity and Density
    
void Calculate_Conserved_Variables(double & P, double & u, double & rho)
{
    U_Cell[0] = rho;				// Rho
    U_Cell[1] = rho*u;				// Rho* u
    U_Cell[2] = (1.0/(gamma-1))*P + 0.5*rho*u*u;	// Rho*E = (P/(gamma-1)) + (1/2)*Rho*u*u - Pressure is evaluated from this expression
}



