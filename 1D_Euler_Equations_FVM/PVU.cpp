#include "1D_Euler_FVM.h"


void PVU_Dissipation(int & Cell_No)
{
  
 double Beta =0.0, Term1 = 0.0,u_R_Minus,u_R_Plus,u_L_Minus,u_L_Plus,a_L,a_R;
 u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;T_L =0.0;u_L_Minus =0.0;u_L_Plus=0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;T_R =0.0;u_R_Minus =0.0;u_R_Plus=0.0;
 
//  cout<<"In PVU_Dissipation Function\n";  
  
// Dissipation at i-1/2 interface between i and i-1   
    
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    u_L_Plus = 0.5*(u_L + fabs(u_L));
    u_L_Minus = 0.5*(u_L - fabs(u_L));
    a_L = C[Cell_No - 1];

    
    G_C_L_Plus[0] = u_L_Plus*U[Cell_No - 1][0];
    G_C_L_Plus[1] = u_L_Plus*U[Cell_No - 1][1];
    G_C_L_Plus[2] = u_L_Plus*U[Cell_No - 1][2];
    
    G_C_L_Minus[0] = u_L_Minus*U[Cell_No - 1][0];
    G_C_L_Minus[1] = u_L_Minus*U[Cell_No - 1][1];
    G_C_L_Minus[2] = u_L_Minus*U[Cell_No - 1][2];

    Beta = 0.5*gamma/(a_L*a_L);     	// beta = 1.0/(2*R*T) = 0.5*gamma / a*a; a*a = gamma*R*T
    Term1 = 0.5*Rho_L/sqrt(M_PI*Beta); // Term1 = rho/(2*sqrt(pi*beta))
    
    G_P_L_Plus[0] = Term1;
    G_P_L_Plus[1] = 0.5*P_L  + Term1*u_L;
    G_P_L_Plus[2] = 0.5*P_L*u_L + (Term1/Rho_L)*(0.5*P_L + E_L);
    
    G_P_L_Minus[0] = -Term1;
    G_P_L_Minus[1] = 0.5*P_L  - Term1*u_L;
    G_P_L_Minus[2] = 0.5*P_L*u_L - (Term1/Rho_L)*(0.5*P_L + E_L);
    
    
// values at i Left state values
    
    u_R = Velocity[Cell_No];
    P_R = Pressure[Cell_No];
    Rho_R = Density[Cell_No];
    E_R =  U[Cell_No][2];
    u_R_Plus = 0.5*(u_R + fabs(u_R));
    u_R_Minus = 0.5*(u_R - fabs(u_R));
    a_R = C[Cell_No];
    
    G_C_R_Plus[0] = u_R_Plus*U[Cell_No][0];
    G_C_R_Plus[1] = u_R_Plus*U[Cell_No][1];
    G_C_R_Plus[2] = u_R_Plus*U[Cell_No][2];
    
    G_C_R_Minus[0] = u_R_Minus*U[Cell_No][0];
    G_C_R_Minus[1] = u_R_Minus*U[Cell_No][1];
    G_C_R_Minus[2] = u_R_Minus*U[Cell_No][2];

    Beta = 0.5*gamma/(a_R*a_R);
    Term1 = 0.5*Rho_R/sqrt(M_PI*Beta); 
    
    G_P_R_Plus[0] = Term1;
    G_P_R_Plus[1] = 0.5*P_R  + Term1*u_R;
    G_P_R_Plus[2] = 0.5*P_R*u_R + (Term1/Rho_R)*(0.5*P_R + E_R);
    
    G_P_R_Minus[0] = -Term1;
    G_P_R_Minus[1] = 0.5*P_R  - Term1*u_R;
    G_P_R_Minus[2] = 0.5*P_R*u_R - (Term1/Rho_R)*(0.5*P_R + E_R);

    Left_Face_Dissipative_Flux[0] = 0.5*((G_C_R_Plus[0] - G_C_L_Plus[0])-(G_C_R_Minus[0] - G_C_L_Minus[0])
                                        + (G_P_R_Plus[0] - G_P_L_Plus[0])-(G_P_R_Minus[0] - G_P_L_Minus[0]));
    Left_Face_Dissipative_Flux[1] = 0.5*((G_C_R_Plus[1] - G_C_L_Plus[1])-(G_C_R_Minus[1] - G_C_L_Minus[1])
                                        + (G_P_R_Plus[1] - G_P_L_Plus[1])-(G_P_R_Minus[1] - G_P_L_Minus[1]));
    Left_Face_Dissipative_Flux[2] = 0.5*((G_C_R_Plus[2] - G_C_L_Plus[2])-(G_C_R_Minus[2] - G_C_L_Minus[2]) 
                                        + (G_P_R_Plus[2] - G_P_L_Plus[2])-(G_P_R_Minus[2] - G_P_L_Minus[2]));
 
// Dissipation at i+1/2 interface between i and i+1   
    
// values at i Left state values
    u_L = Velocity[Cell_No];
    P_L = Pressure[Cell_No];
    Rho_L = Density[Cell_No];
    E_L =  U[Cell_No][2];
    a_L = C[Cell_No];
    
    u_L_Plus = 0.5*(u_L + fabs(u_L));
    u_L_Minus = 0.5*(u_L - fabs(u_L));
    

    G_C_L_Plus[0] = u_L_Plus*U[Cell_No][0];
    G_C_L_Plus[1] = u_L_Plus*U[Cell_No][1];
    G_C_L_Plus[2] = u_L_Plus*U[Cell_No][2];
    
    G_C_L_Minus[0] = u_L_Minus*U[Cell_No][0];
    G_C_L_Minus[1] = u_L_Minus*U[Cell_No][1];
    G_C_L_Minus[2] = u_L_Minus*U[Cell_No][2];

    Beta = 0.5*gamma/(a_L*a_L);
    Term1 = 0.5*Rho_L/sqrt(M_PI*Beta); 
    
    G_P_L_Plus[0] = Term1;
    G_P_L_Plus[1] = 0.5*P_L  + Term1*u_L;
    G_P_L_Plus[2] = 0.5*P_L*u_L + (Term1/Rho_L)*(0.5*P_L + E_L);
    
    G_P_L_Minus[0] = -Term1;
    G_P_L_Minus[1] = 0.5*P_L  - Term1*u_L;
    G_P_L_Minus[2] = 0.5*P_L*u_L - (Term1/Rho_L)*(0.5*P_L + E_L);
    
    
// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];
    
    u_R_Plus = 0.5*(u_R + fabs(u_R));
    u_R_Minus = 0.5*(u_R - fabs(u_R));
    

    G_C_R_Plus[0] = u_R_Plus*U[Cell_No + 1][0];
    G_C_R_Plus[1] = u_R_Plus*U[Cell_No + 1][1];
    G_C_R_Plus[2] = u_R_Plus*U[Cell_No + 1][2];
    
    G_C_R_Minus[0] = u_R_Minus*U[Cell_No + 1][0];
    G_C_R_Minus[1] = u_R_Minus*U[Cell_No + 1][1];
    G_C_R_Minus[2] = u_R_Minus*U[Cell_No + 1][2];

    Beta = 0.5*gamma/(a_R*a_R);
    Term1 = 0.5*Rho_R/sqrt(M_PI*Beta); 
    
    G_P_R_Plus[0] = Term1;
    G_P_R_Plus[1] = 0.5*P_R  + Term1*u_R;
    G_P_R_Plus[2] = 0.5*P_R*u_R + (Term1/Rho_R)*(0.5*P_R + E_R);
    
    G_P_R_Minus[0] = -Term1;
    G_P_R_Minus[1] = 0.5*P_R  - Term1*u_R;
    G_P_R_Minus[2] = 0.5*P_R*u_R - (Term1/Rho_R)*(0.5*P_R + E_R);

        
    Right_Face_Dissipative_Flux[0] = 0.5*((G_C_R_Plus[0] - G_C_L_Plus[0])-(G_C_R_Minus[0] - G_C_L_Minus[0]) 
                                        + (G_P_R_Plus[0] - G_P_L_Plus[0])-(G_P_R_Minus[0] - G_P_L_Minus[0]));
    Right_Face_Dissipative_Flux[1] = 0.5*((G_C_R_Plus[1] - G_C_L_Plus[1])-(G_C_R_Minus[1] - G_C_L_Minus[1])
                                        + (G_P_R_Plus[1] - G_P_L_Plus[1])-(G_P_R_Minus[1] - G_P_L_Minus[1]));
    Right_Face_Dissipative_Flux[2] = 0.5*((G_C_R_Plus[2] - G_C_L_Plus[2])-(G_C_R_Minus[2] - G_C_L_Minus[2])
                                        + (G_P_R_Plus[2] - G_P_L_Plus[2])-(G_P_R_Minus[2] - G_P_L_Minus[2]));
}