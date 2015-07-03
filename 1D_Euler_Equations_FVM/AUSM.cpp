#include "1D_Euler_FVM.h"


void AUSM(int & Cell_No)
{
  double Mach_No_L=0.0,Mach_No_R=0.0,Mach_No=0.0,a_L=0.0,a=0.0,a_R=0.0,H_L=0.0,H=0.0,H_R=0.0;
 vector<double> G_C_1(3,0.0),G_C_2(3,0.0),G_C_3(3,0.0);
 double P_Plus=0.0,P_Minus =0.0,P_L_Minus=0.0,P_L_Plus=0.0,P_R_Minus=0.0,P_R_Plus=0.0,M_L_Minus=0.0,M_R_Plus=0.0,M_Minus=0.0,M_Plus=0.0,M_L_Plus=0.0,M_R_Minus=0.0;
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;
 
//  cout<<"In AUSM Function\n";  
//  cout<<"Cell No \t"<<Cell_No<<endl;
 
 /* AUSM Scheme Flux as given by 
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
    Mach_No_L = u_L/a_L;
    H_L = (E_L + P_L)/Rho_L;

    G_C_1[0] = Rho_L*a_L;
    G_C_1[1] = Rho_L*a_L*u_L;
    G_C_1[2] = Rho_L*a_L*H_L;

//     cout<<"i-1 Mach number \t"<<Mach_No_L<<endl;

// Condition of  |M| <= 1
    if((Mach_No_L <= 1.0) and (Mach_No_L >= -1.0))
    {
//		M+ =  (1/4)*(M+1)^2
      M_L_Plus =  0.25*(Mach_No_L + 1.0)*(Mach_No_L + 1.0);
//		M- = (1/4)*(M-1)^2
      M_L_Minus = -0.25*(Mach_No_L - 1.0)*(Mach_No_L - 1.0);
//		(P/4)*(M + 1)^2(2 - M)
      P_L_Plus = 0.25*P_L*(Mach_No_L + 1.0)*(Mach_No_L + 1.0)*(2.0 - Mach_No_L);
//		(P/4)*(M - 1)^2(2 + M)
      P_L_Minus = 0.25*P_L*(Mach_No_L - 1.0)*(Mach_No_L - 1.0)*(2.0 + Mach_No_L);
    }
    else
    {
//	 M+ = 0.5*(M + |M|)
    M_L_Plus = 0.5*(Mach_No_L + fabs(Mach_No_L));
    
//	 M+ = 0.5*(M - |M|)
    M_L_Minus = 0.5*(Mach_No_L - fabs(Mach_No_L));

// 	P+ = P/2*(M + |M|)/M
    P_L_Plus = 0.5*P_L*(Mach_No_L + fabs(Mach_No_L))/Mach_No_L;
    
// 	P+ = P/2*(M - |M|)/M
    P_L_Minus = 0.5*P_L*(Mach_No_L - fabs(Mach_No_L))/Mach_No_L; 

    }
        
// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    a = C[Cell_No];
    Mach_No = u/a;
    H = (E + P)/Rho;

    G_C_2[0] = Rho*a;
    G_C_2[1] = Rho*a*u;
    G_C_2[2] = Rho*a*H;
    
// cout<<"i Mach number \t"<<Mach_No<<endl;
// Condition of  |M| <= 1
    if((Mach_No <= 1.0) and (Mach_No >= -1.0))
    {
//		M+ =  (1/4)*(M+1)^2
      M_Plus =  0.25*(Mach_No + 1.0)*(Mach_No + 1.0);
//		M- = (1/4)*(M-1)^2
      M_Minus = -0.25*(Mach_No - 1.0)*(Mach_No - 1.0);
//		(P/4)*(M + 1)^2(2 - M)
      P_Plus = 0.25*P*(Mach_No + 1.0)*(Mach_No + 1.0)*(2.0 - Mach_No);
//		(P/4)*(M - 1)^2(2 + M)
      P_Minus = 0.25*P*(Mach_No - 1.0)*(Mach_No - 1.0)*(2.0 + Mach_No);
    }
    else
    {
//	 M+ = 0.5*(M + |M|)
    M_Plus = 0.5*(Mach_No + fabs(Mach_No));
    
//	 M+ = 0.5*(M - |M|)
    M_Minus = 0.5*(Mach_No - fabs(Mach_No));

// 	P+ = P/2*(M + |M|)/M
    P_Plus = 0.5*P*(Mach_No + fabs(Mach_No))/Mach_No;
    
// 	P+ = P/2*(M - |M|)/M
    P_Minus = 0.5*P*(Mach_No - fabs(Mach_No))/Mach_No; 

    }


// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];
    Mach_No_R = u_R/a_R;
    H_R = (E_R + P_R)/Rho_R;

    G_C_3[0] = Rho_R*a_R;
    G_C_3[1] = Rho_R*a_R*u_R;
    G_C_3[2] = Rho_R*a_R*H_R;
    
//   cout<<"i+1 Mach number \t"<<Mach_No_R<<endl;

// Condition of  |M| <= 1
    if((Mach_No_R <= 1.0) and (Mach_No_R >= -1.0))
    {
//		M+ =  (1/4)*(M+1)^2
      M_R_Plus =  0.25*(Mach_No_R + 1.0)*(Mach_No_R + 1.0);
//		M- = (1/4)*(M-1)^2
      M_R_Minus = -0.25*(Mach_No_R - 1.0)*(Mach_No_R - 1.0);
//		(P/4)*(M + 1)^2(2 - M)
      P_R_Plus = 0.25*P_R*(Mach_No_R + 1.0)*(Mach_No_R + 1.0)*(2.0 - Mach_No_R);
//		(P/4)*(M - 1)^2(2 + M)
      P_R_Minus = 0.25*P_R*(Mach_No_R - 1.0)*(Mach_No_R - 1.0)*(2.0 + Mach_No_R);
    }
    else
    {
//	 M+ = 0.5*(M + |M|)
    M_R_Plus = 0.5*(Mach_No_R + fabs(Mach_No_R));
    
//	 M+ = 0.5*(M - |M|)
    M_R_Minus = 0.5*(Mach_No_R - fabs(Mach_No_R));

// 	P+ = P/2*(M + |M|)/M
    P_R_Plus = 0.5*P_R*(Mach_No_R + fabs(Mach_No_R))/Mach_No_R;
    
// 	P+ = P/2*(M - |M|)/M
    P_R_Minus = 0.5*P_R*(Mach_No_R - fabs(Mach_No_R))/Mach_No_R; 
    }

// cout<<"M_I \t"<<(M_Plus + M_R_Minus)<<endl;
//    Flux at i+1/2 interface = G+|i + G-|i+1
    Right_Face_Average_Flux[0] = (M_Plus + M_R_Minus)*0.5*(G_C_2[0] + G_C_3[0])-0.5*fabs(M_Plus + M_R_Minus)*(G_C_3[0] - G_C_2[0]);
    Right_Face_Average_Flux[1] = (M_Plus + M_R_Minus)*0.5*(G_C_2[1] + G_C_3[1])-0.5*fabs(M_Plus + M_R_Minus)*(G_C_3[1] - G_C_2[1]) + (P_Plus + P_R_Minus) ;
    Right_Face_Average_Flux[2] = (M_Plus + M_R_Minus)*0.5*(G_C_2[2] + G_C_3[2])-0.5*fabs(M_Plus + M_R_Minus)*(G_C_3[2] - G_C_2[2]);
    
    Right_Face_Dissipative_Flux[0] = 0.0;
    Right_Face_Dissipative_Flux[1] = 0.0;
    Right_Face_Dissipative_Flux[2] = 0.0;
// cout<<"M_I \t"<<(M_Minus + M_L_Plus)<<endl;
//    Flux at i -1/2 interface = G+|i-1 + G-|i
    Left_Face_Average_Flux[0] = (M_Minus + M_L_Plus)*0.5*(G_C_1[0] + G_C_2[0])  -0.5*fabs(M_Minus + M_L_Plus)*(G_C_2[0]-G_C_1[0]);
    Left_Face_Average_Flux[1] = (M_Minus + M_L_Plus)*0.5*(G_C_1[1] + G_C_2[1])  -0.5*fabs(M_Minus + M_L_Plus)*(G_C_2[1]-G_C_1[1]) + (P_Minus + P_L_Plus);
    Left_Face_Average_Flux[2] = (M_Minus + M_L_Plus)*0.5*(G_C_1[2] + G_C_2[2])  -0.5*fabs(M_Minus + M_L_Plus)*(G_C_2[2]-G_C_1[2]);
    
    Left_Face_Dissipative_Flux[0] = 0.0;
    Left_Face_Dissipative_Flux[1] = 0.0;
    Left_Face_Dissipative_Flux[2] = 0.0;
}