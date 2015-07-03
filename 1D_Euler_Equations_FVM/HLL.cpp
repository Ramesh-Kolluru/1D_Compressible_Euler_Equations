#include "1D_Euler_FVM.h"


void HLL(int & Cell_No)
{
  double a_L=0.0,a=0.0,a_R=0.0,S_R=0.0,S_L=0.0;
 vector<double> Flux1(3,0.0),Flux2(3,0.0),Flux3(3,0.0);
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;
 
//  cout<<"In HLL Function\n";  
//  cout<<"Cell No \t"<<Cell_No<<endl;
 
 /* HLL Scheme Flux as given by 
  * 
  * F_C_Plus[0] = 
  * F_C_Plus[1] = 
  * F_C_Plus[2] = 
  * 
  * F_C_Minus[0] = 
  * F_C_Minus[1] = 
  * F_C_Minus[2] = 
  * */
  
// Upwind flux at i-1/2 interface between i and i-1   
    
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    a_L = C[Cell_No - 1];

    Flux1[0] = Rho_L*u_L;
    Flux1[1] = Rho_L*u_L*u_L + P_L;
    Flux1[2] = (E_L + P_L)*u_L;


// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    a = C[Cell_No];

    Flux2[0] = Rho*u;
    Flux2[1] = Rho*u*u + P;
    Flux2[2] = (E + P)*u;
    

// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];

    Flux3[0] = Rho_R*u_R;
    Flux3[1] = Rho_R*u_R*u_R + P_R;
    Flux3[2] = (E_R + P_R)*u_R;
    


// cout<<"M_I \t"<<(M_Plus + M_R_Minus)<<endl;
//    Flux at i+1/2 interface = G+|i + G-|i+1
    S_L = min(u-a,u_R-a_R);
    S_R = max(u+a,u_R + a_R);
    if((S_L <=0.0) and (S_R >= 0.0))
    {
      Right_Face_Average_Flux[0] = (S_R* Flux2[0] - S_L*Flux3[0] + S_L*S_R*(U[Cell_No + 1][0] - U[Cell_No][0]))/(S_R - S_L);
      Right_Face_Average_Flux[1] = (S_R* Flux2[1] - S_L*Flux3[1] + S_L*S_R*(U[Cell_No + 1][1] - U[Cell_No][1]))/(S_R - S_L);
      Right_Face_Average_Flux[2] = (S_R* Flux2[2] - S_L*Flux3[2] + S_L*S_R*(U[Cell_No + 1][2] - U[Cell_No][2]))/(S_R - S_L);
    }
    if (S_L >= 0.0)
    {
      Right_Face_Average_Flux[0] = Flux2[0];
      Right_Face_Average_Flux[1] = Flux2[1];
      Right_Face_Average_Flux[2] = Flux2[2];
    }
    if(S_R <= 0.0)
    {
      Right_Face_Average_Flux[0] = Flux3[0];
      Right_Face_Average_Flux[1] = Flux3[1];
      Right_Face_Average_Flux[2] = Flux3[2];
    }
    
    Right_Face_Dissipative_Flux[0] = 0.0;
    Right_Face_Dissipative_Flux[1] = 0.0;
    Right_Face_Dissipative_Flux[2] = 0.0;
// cout<<"M_I \t"<<(M_Minus + M_L_Plus)<<endl;
//    Flux at i -1/2 interface = G+|i-1 + G-|i
    S_L = min(u_L-a_L,u-a);
    S_R = max(u_L+a_L,u + a);

    if((S_L <=0.0) and (S_R >= 0.0))
    {
      Left_Face_Average_Flux[0] = (S_R*Flux1[0] - S_L*Flux2[0] + S_L*S_R*(U[Cell_No][0] - U[Cell_No - 1][0]))/(S_R - S_L);
      Left_Face_Average_Flux[1] = (S_R*Flux1[1] - S_L*Flux2[1] + S_L*S_R*(U[Cell_No][1] - U[Cell_No - 1][1]))/(S_R - S_L);
      Left_Face_Average_Flux[2] = (S_R*Flux1[2] - S_L*Flux2[2] + S_L*S_R*(U[Cell_No][2] - U[Cell_No - 1][2]))/(S_R - S_L);
    }
    if (S_L >= 0.0)
    {
      Left_Face_Average_Flux[0] = Flux1[0];
      Left_Face_Average_Flux[1] = Flux1[1];
      Left_Face_Average_Flux[2] = Flux1[2];
    }
    if(S_R <= 0.0)
    {
      Left_Face_Average_Flux[0] = Flux2[0];
      Left_Face_Average_Flux[1] = Flux2[1];
      Left_Face_Average_Flux[2] = Flux2[2];
    }
    
    Left_Face_Dissipative_Flux[0] = 0.0;
    Left_Face_Dissipative_Flux[1] = 0.0;
    Left_Face_Dissipative_Flux[2] = 0.0;
}