#include "1D_Euler_FVM.h"

void Mod(double & lambda, double & mod_lambda)
{
  (lambda>=0.0?mod_lambda = lambda: mod_lambda = -lambda);
}

void Steger_Warming(int & Cell_No)
{
  double lambda_1=0.0, lambda_1_Minus = 0.0,a_L=0.0,a=0.0,a_R=0.0,lambda_1_Plus=0.0,lambda_2=0.0,lambda_2_Minus=0.0,lambda_2_Plus=0.0,lambda_3=0.0,lambda_3_Minus=0.0,lambda_3_Plus=0.0,H_L=0.0,H_R=0.0,H=0.0;
  double mod_lambda_1 = 0.0, mod_lambda_2 = 0.0, mod_lambda_3 = 0.0;
 vector<double> G_1_Plus(3,0.0),G_1_Minus(3,0.0),G_2_Plus(3,0.0),G_2_Minus(3,0.0),G_3_Plus(3,0.0),G_3_Minus(3,0.0);
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;T_L =0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;T_R =0.0;
 
//  cout<<"In PVU_Dissipation Function\n";  
 
 /* Steger and Warming Scheme Flux as given in Toro book Page No 276 equation 8.67
  * 
  * lambda_1 = u-a;
  * lambda_2 = u;
  * lambda_3 = u+a
  * 
  * lambda_i_Plus = 0.5*(lambda_i + mod(lambda_i)); i = 1,2,3
  * lambda_i_Minus = 0.5*(lambda_i - mod(lambda_i));
  * 
  * mod(lambda_i) 	= lambda_i if lambda_i>=0
  * 			= -lambda_i if lambda_i<=0
  * 
  * H = (E + P )/Rho = (u*u)/2 + (a*a)/(gamma-1)
  * 
  * F_Plus[0] = (rho/(2*gamma))*(lambda_1_Plus + 2*(gamma-1)*lambda_2_Plus +lambda_3_Plus)
  * F_Plus[1] = (rho/(2*gamma))*((u-a)*lambda_1_Plus + u*2*(gamma-1)*lambda_2_Plus + (u+a)*lambda_3_Plus)
  * F_Plus[2] = (rho/(2*gamma))*((H-u*a)*lambda_1_Plus + u*u*(gamma-1)*lambda_2_Plus +(H+u*a)*lambda_3_Plus)
  * 
  * F_Minus[0] = (rho/(2*gamma))*(lambda_1_Minus + 2*(gamma-1)*lambda_2_Minus +lambda_3_Minus)
  * F_Minus[1] = (rho/(2*gamma))*((u-a)*lambda_1_Minus + u*2*(gamma-1)*lambda_2_Minus + (u+a)*lambda_3_Minus)
  * F_Minus[2] = (rho/(2*gamma))*((H-u*a)*lambda_1_Minus + u*u*(gamma-1)*lambda_2_Minus +(H+u*a)*lambda_3_Minus)
  * */
  
// Upwind flux at i-1/2 interface between i and i-1   
    
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    a_L = C[Cell_No - 1];
    H_L = (E_L + P_L)/Rho_L;
    
    lambda_1 = u_L - a_L;
    lambda_2 = u_L ;
    lambda_3 = u_L + a_L;
    
    Mod(lambda_1,mod_lambda_1);
    lambda_1_Minus = 0.5*(lambda_1 - mod_lambda_1);
    lambda_1_Plus = 0.5*(lambda_1 + mod_lambda_1);
    
    Mod(lambda_2,mod_lambda_2);
    lambda_2_Minus = 0.5*(lambda_2 - mod_lambda_2);
    lambda_2_Plus = 0.5*(lambda_2 + mod_lambda_2);
    
    Mod(lambda_3,mod_lambda_3);
    lambda_3_Minus = 0.5*(lambda_3 - mod_lambda_3);
    lambda_3_Plus = 0.5*(lambda_3 + mod_lambda_3);
    
    G_1_Plus[0] = (0.5*Rho_L/gamma)*(lambda_1_Plus + 2.0*(gamma-1)*lambda_2_Plus +lambda_3_Plus);
    G_1_Plus[1] = (0.5*Rho_L/gamma)*((u_L-a_L)*lambda_1_Plus + u_L*2.0*(gamma-1)*lambda_2_Plus + (u_L+a_L)*lambda_3_Plus);
    G_1_Plus[2] = (0.5*Rho_L/gamma)*((H_L-u_L*a_L)*lambda_1_Plus + u_L*u_L*(gamma-1)*lambda_2_Plus +(H_L+u_L*a_L)*lambda_3_Plus);
    
    G_1_Minus[0] = (0.5*Rho_L/gamma)*(lambda_1_Minus + 2.0*(gamma-1)*lambda_2_Minus +lambda_3_Minus);
    G_1_Minus[1] = (0.5*Rho_L/gamma)*((u_L-a_L)*lambda_1_Minus + u_L*2.0*(gamma-1)*lambda_2_Minus + (u_L+a_L)*lambda_3_Minus);
    G_1_Minus[2] = (0.5*Rho_L/gamma)*((H_L-u_L*a_L)*lambda_1_Minus + u_L*u_L*(gamma-1)*lambda_2_Minus +(H_L+u_L*a_L)*lambda_3_Minus);
    
    
// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    a = C[Cell_No];
    H = (E + P)/Rho;
    
    lambda_1 = u - a;
    lambda_2 = u ;
    lambda_3 = u + a;
    
    Mod(lambda_1,mod_lambda_1);
    lambda_1_Minus = 0.5*(lambda_1 - mod_lambda_1);
    lambda_1_Plus = 0.5*(lambda_1 + mod_lambda_1);
    
    Mod(lambda_2,mod_lambda_2);
    lambda_2_Minus = 0.5*(lambda_2 - mod_lambda_2);
    lambda_2_Plus = 0.5*(lambda_2 + mod_lambda_2);
    
    Mod(lambda_3,mod_lambda_3);
    lambda_3_Minus = 0.5*(lambda_3 - mod_lambda_3);
    lambda_3_Plus = 0.5*(lambda_3 + mod_lambda_3);
    
    G_2_Plus[0] = (0.5*Rho/gamma)*(lambda_1_Plus + 2.0*(gamma-1)*lambda_2_Plus +lambda_3_Plus);
    G_2_Plus[1] = (0.5*Rho/gamma)*((u-a)*lambda_1_Plus + u*2.0*(gamma-1)*lambda_2_Plus + (u+a)*lambda_3_Plus);
    G_2_Plus[2] = (0.5*Rho/gamma)*((H-u*a)*lambda_1_Plus + u*u*(gamma-1)*lambda_2_Plus +(H+u*a)*lambda_3_Plus);
    
    G_2_Minus[0] = (0.5*Rho/gamma)*(lambda_1_Minus + 2.0*(gamma-1)*lambda_2_Minus +lambda_3_Minus);
    G_2_Minus[1] = (0.5*Rho/gamma)*((u-a)*lambda_1_Minus + u*2.0*(gamma-1)*lambda_2_Minus + (u+a)*lambda_3_Minus);
    G_2_Minus[2] = (0.5*Rho/gamma)*((H-u*a)*lambda_1_Minus + u*u*(gamma-1)*lambda_2_Minus +(H+u*a)*lambda_3_Minus);

// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];
    H_R = (E_R + P_R)/Rho_R;
    
    lambda_1 = u_R - a_R;
    lambda_2 = u_R ;
    lambda_3 = u_R + a_R;
    
    Mod(lambda_1,mod_lambda_1);
    lambda_1_Minus = 0.5*(lambda_1 - mod_lambda_1);
    lambda_1_Plus = 0.5*(lambda_1 + mod_lambda_1);
    
    Mod(lambda_2,mod_lambda_2);
    lambda_2_Minus = 0.5*(lambda_2 - mod_lambda_2);
    lambda_2_Plus = 0.5*(lambda_2 + mod_lambda_2);
    
    Mod(lambda_3,mod_lambda_3);
    lambda_3_Minus = 0.5*(lambda_3 - mod_lambda_3);
    lambda_3_Plus = 0.5*(lambda_3 + mod_lambda_3);
    
    G_3_Plus[0] = (0.5*Rho_R/gamma)*(lambda_1_Plus + 2.0*(gamma-1)*lambda_2_Plus +lambda_3_Plus);
    G_3_Plus[1] = (0.5*Rho_R/gamma)*((u_R-a_R)*lambda_1_Plus + u_R*2.0*(gamma-1)*lambda_2_Plus + (u_R+a_R)*lambda_3_Plus);
    G_3_Plus[2] = (0.5*Rho_R/gamma)*((H_R-u_R*a_R)*lambda_1_Plus + u_R*u_R*(gamma-1)*lambda_2_Plus +(H_R+u_R*a_R)*lambda_3_Plus);
    
    G_3_Minus[0] = (0.5*Rho_R/gamma)*(lambda_1_Minus + 2.0*(gamma-1)*lambda_2_Minus +lambda_3_Minus);
    G_3_Minus[1] = (0.5*Rho_R/gamma)*((u_R-a_R)*lambda_1_Minus + u_R*2.0*(gamma-1)*lambda_2_Minus + (u_R+a_R)*lambda_3_Minus);
    G_3_Minus[2] = (0.5*Rho_R/gamma)*((H_R-u_R*a_R)*lambda_1_Minus + u_R*u_R*(gamma-1)*lambda_2_Minus +(H_R+u_R*a_R)*lambda_3_Minus);
    
    
//    Flux at i+1/2 interface = G+|i + G-|i+1
    Right_Face_Average_Flux[0] = (G_2_Plus[0] + G_3_Minus[0]);
    Right_Face_Average_Flux[1] = (G_2_Plus[1] + G_3_Minus[1]);
    Right_Face_Average_Flux[2] = (G_2_Plus[2] + G_3_Minus[2]);
    
    Right_Face_Dissipative_Flux[0] = 0.0;
    Right_Face_Dissipative_Flux[1] = 0.0;
    Right_Face_Dissipative_Flux[2] = 0.0;

//    Flux at i -1/2 interface = G+|i-1 + G-|i
    Left_Face_Average_Flux[0] = (G_1_Plus[0] + G_2_Minus[0]);
    Left_Face_Average_Flux[1] = (G_1_Plus[1] + G_2_Minus[1]);
    Left_Face_Average_Flux[2] = (G_1_Plus[2] + G_2_Minus[2]);
    
    Left_Face_Dissipative_Flux[0] = 0.0;
    Left_Face_Dissipative_Flux[1] = 0.0;
    Left_Face_Dissipative_Flux[2] = 0.0;
}