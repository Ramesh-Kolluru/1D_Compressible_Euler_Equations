#include "1D_Euler_FVM.h"


void Vanleer(int & Cell_No)
{
  double Mach_No_L=0.0,Mach_No_R=0.0,Mach_No=0.0,a_L=0.0,a=0.0,a_R=0.0,Term1 =0.0,Term2 = 0.0, Term3 =0.0;
 vector<double> G_1_Plus(3,0.0),G_1_Minus(3,0.0),G_2_Plus(3,0.0),G_2_Minus(3,0.0),G_3_Plus(3,0.0),G_3_Minus(3,0.0);
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;
 
//  cout<<"In Vanleer Function\n";  
 
 /* Vanleer Scheme Flux as given in Toro book Page No 277 equation 8.71,8.72,8.73,8.74,8.75
  * 
  * F_Plus[0] = ((1/4)(Rho*a)(1+M)^2)*(1)
  * F_Plus[1] = ((1/4)(Rho*a)(1+M)^2)*((2*a/gamma)*(0.5*(gamma-1))*M + 1)
  * F_Plus[2] = ((1/4)(Rho*a)(1+M)^2)*((2*a^2/(gamma^2-1))*(0.5*(gamma-1))*M + 1)^2
  * 
  * F_Minus[0] = -(1/4)(Rho*a)(1-M)^2)*1
  * F_Minus[1] = -(1/4)(Rho*a)(1-M)^2)*((2*a/(gamma-1))*(0.5*(gamma-1))*M - 1)
  * F_Minus[2] = -(1/4)(Rho*a)(1-M)^2)*((2*a^2/(gamma^2-1))*(0.5*(gamma-1))*M - 1)^2
  * */
  
// Upwind flux at i-1/2 interface between i and i-1   
    
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    a_L = C[Cell_No - 1];
    Mach_No_L = u_L/a_L;
    
//     cout<<"Left state Mach Number\t"<<Mach_No_L<<endl;
      Term1 = 0.25*(Rho_L*a_L)*(1.0 + Mach_No_L)*(1.0 + Mach_No_L);
      Term2 = (0.5*(gamma-1.0)*Mach_No_L + 1.0);
      Term3 = ((2.0*a_L*a_L)/(gamma*gamma -1.0));
      
      G_1_Plus[0] = Term1;
      G_1_Plus[1] = Term1*((2.0*a_L/gamma)*Term2);
      G_1_Plus[2] = Term1*(Term3*Term2*Term2);
      
      
      Term1 = -0.25*(Rho_L*a_L)*(1.0 - Mach_No_L)*(1.0 - Mach_No_L);
      Term2 = (0.5*(gamma-1.0)*Mach_No_L - 1);
      Term3 = ((2.0*a_L*a_L)/(gamma*gamma -1.0));
    
      G_1_Minus[0] = Term1;
      G_1_Minus[1] = Term1*((2.0*a_L/gamma)*Term2);
      G_1_Minus[2] = Term1*(Term3*Term2*Term2); 
    
    
    
    if(Mach_No_L > 1.0) 
    {

      G_1_Plus[0] = Rho_L*a_L*Mach_No_L;
      G_1_Plus[1] = Rho_L*a_L*a_L*(Mach_No_L*Mach_No_L + 1.0/gamma);
      G_1_Plus[2] = Rho_L*a_L*a_L*a_L*Mach_No_L*(0.5*Mach_No_L*Mach_No_L + 1.0/(gamma-1.0));
      
    
      G_1_Minus[0] = 0.0;
      G_1_Minus[1] = 0.0;
      G_1_Minus[2] = 0.0; 

    }
    if(Mach_No_L < -1.0)
    {
    G_1_Minus[0] = Rho_L*a_L*Mach_No_L;
    G_1_Minus[1] = Rho_L*a_L*a_L*(Mach_No_L*Mach_No_L + 1.0/gamma);
    G_1_Minus[2] = Rho_L*a_L*a_L*a_L*Mach_No_L*(0.5*Mach_No_L*Mach_No_L + 1.0/(gamma-1.0));
    
   
    G_1_Plus[0] = 0.0;
    G_1_Plus[1] = 0.0;
    G_1_Plus[2] = 0.0; 
      
    }
    
    
// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    a = C[Cell_No];
    Mach_No = u/a;
    
//     cout<<"ith Cell Mach Number\t"<<Mach_No<<endl;
    
    Term1 = 0.25*(Rho*a)*(1.0 + Mach_No)*(1.0 + Mach_No);
    Term2 = (0.5*(gamma-1.0)*Mach_No + 1.0);
    Term3 = ((2.0*a*a)/(gamma*gamma -1.0));
    
    G_2_Plus[0] = Term1;
    G_2_Plus[1] = Term1*((2.0*a/gamma)*Term2);
    G_2_Plus[2] = Term1*(Term3*Term2*Term2);
    
    
    Term1 = -0.25*(Rho*a)*(1.0 - Mach_No)*(1.0 - Mach_No);
    Term2 = (0.5*(gamma-1.0)*Mach_No - 1.0);
    Term3 = ((2.0*a*a)/(gamma*gamma -1.0));
    
    G_2_Minus[0] = Term1;
    G_2_Minus[1] = Term1*((2.0*a/gamma)*Term2);
    G_2_Minus[2] = Term1*(Term3*Term2*Term2); 
    
    
    if(Mach_No > 1.0) 
    {

    G_2_Plus[0] = Rho*a*Mach_No;
    G_2_Plus[1] = Rho*a*a*(Mach_No*Mach_No + 1.0/gamma);
    G_2_Plus[2] = Rho*a*a*a*Mach_No*(0.5*Mach_No*Mach_No + 1.0/(gamma-1.0));
    
   
    G_2_Minus[0] = 0.0;
    G_2_Minus[1] = 0.0;
    G_2_Minus[2] = 0.0; 

    }
    if(Mach_No < -1.0)
    {
    G_2_Minus[0] = Rho*a*Mach_No;
    G_2_Minus[1] = Rho*a*a*(Mach_No*Mach_No + 1.0/gamma);
    G_2_Minus[2] = Rho*a*a*a*Mach_No*(0.5*Mach_No*Mach_No + 1.0/(gamma-1.0));
    
   
    G_2_Plus[0] = 0.0;
    G_2_Plus[1] = 0.0;
    G_2_Plus[2] = 0.0; 
      
    }
    

// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];
    Mach_No_R = u_R/a_R;
    
//     cout<<"i+1 th Cell Mach Number\t"<<Mach_No_R<<endl;
    
    
    Term1 = 0.25*(Rho_R*a_R)*(1.0 + Mach_No_R)*(1.0 + Mach_No_R);
    Term2 = (0.5*(gamma-1.0)*Mach_No_R + 1.0);
    Term3 = ((2.0*a_R*a_R)/(gamma*gamma -1.0));
    
    G_3_Plus[0] = Term1;
    G_3_Plus[1] = Term1*((2.0*a_R/gamma)*Term2);
    G_3_Plus[2] = Term1*(Term3*Term2*Term2);
    
    
    Term1 = -0.25*(Rho_R*a_R)*(1.0 - Mach_No_R)*(1.0 - Mach_No_R);
    Term2 = (0.5*(gamma-1.0)*Mach_No_R - 1.0);
    Term3 = ((2.0*a_R*a_R)/(gamma*gamma -1.0));
    
    G_3_Minus[0] = Term1;
    G_3_Minus[1] = Term1*((2.0*a_R/gamma)*Term2);
    G_3_Minus[2] = Term1*(Term3*Term2*Term2); 
     
    
    
    if(Mach_No_R > 1.0) 
    {

    G_3_Plus[0] = Rho_R*a_R*Mach_No_R;
    G_3_Plus[1] = Rho_R*a_R*a_R*(Mach_No_R*Mach_No_R + 1.0/gamma);
    G_3_Plus[2] = Rho_R*a_R*a_R*a_R*Mach_No_R*(0.5*Mach_No_R*Mach_No_R + 1.0/(gamma-1.0));
    
   
    G_3_Minus[0] = 0.0;
    G_3_Minus[1] = 0.0;
    G_3_Minus[2] = 0.0; 

    }
    if(Mach_No_R < -1.0)
    {
    G_3_Minus[0] = Rho_R*a_R*Mach_No_R;
    G_3_Minus[1] = Rho_R*a_R*a_R*(Mach_No_R*Mach_No_R + 1.0/gamma);
    G_3_Minus[2] = Rho_R*a_R*a_R*a_R*Mach_No_R*(0.5*Mach_No_R*Mach_No_R + 1.0/(gamma-1.0));
    
   
    G_3_Plus[0] = 0.0;
    G_3_Plus[1] = 0.0;
    G_3_Plus[2] = 0.0; 
      
    }
    
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