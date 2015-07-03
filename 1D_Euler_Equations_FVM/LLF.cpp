#include "1D_Euler_FVM.h"


void LLF_Dissipation(int & Cell_No)
{
 double Mod_Alpha = 0.0;
 
//	Cell_No = i , Cell_No +1 = i+1, Cell_No -1 = i-1
 
//  Mod_Alpha = max(Eigen_Values[Cell_No][0],Eigen_Values[Cell_No - 1][0]);
    Mod_Alpha = 0.5*(Eigen_Values[Cell_No][0]+Eigen_Values[Cell_No - 1][0]);
//  Mod_Alpha = 0.5 * del_X/ Min_delta_t;
    
    Left_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha*(U[Cell_No][0] - U[Cell_No-1][0]);
    Left_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha*(U[Cell_No][1] - U[Cell_No-1][1]);
    Left_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha*(U[Cell_No][2] - U[Cell_No-1][2]);
 
//  Mod_Alpha = max(Eigen_Values[Cell_No][0], Eigen_Values[Cell_No + 1][0]);
    Mod_Alpha = 0.5*(Eigen_Values[Cell_No][0]+Eigen_Values[Cell_No + 1][0]);
//  Mod_Alpha = 0.5 * del_X/ Min_delta_t;
 
    Right_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha*(U[Cell_No + 1][0] - U[Cell_No][0]);
    Right_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha*(U[Cell_No + 1][1] - U[Cell_No][1]);
    Right_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha*(U[Cell_No + 1][2] - U[Cell_No][2]);
}



void Rusanov_Dissipation(int & Cell_No)
{
 double Mod_Alpha = 0.0;
 vector<double> S(4,0.0);
 
 u_L=0.0; u_R=0.0; c_L=0.0;c_R=0.0;
//	Cell_No = i , Cell_No +1 = i+1, Cell_No -1 = i-1
 
    u_L = Velocity[Cell_No - 1]; c_L = C[Cell_No - 1];c_R = C[Cell_No]; u_R = Velocity[Cell_No];
    S[0] = fabs(u_L - c_L);
    S[1] = fabs(u_L + c_L);
    S[2] = fabs(u_R - c_R);
    S[3] = fabs(u_R + c_R);
  
    Mod_Alpha = *max_element(S.begin(),S.end());
    
    Left_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha*(U[Cell_No][0] - U[Cell_No-1][0]);
    Left_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha*(U[Cell_No][1] - U[Cell_No-1][1]);
    Left_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha*(U[Cell_No][2] - U[Cell_No-1][2]);
 
    u_L = Velocity[Cell_No]; c_L = C[Cell_No];c_R = C[Cell_No + 1]; u_R = Velocity[Cell_No + 1];
    S[0] = fabs(u_L - c_L);
    S[1] = fabs(u_L + c_L);
    S[2] = fabs(u_R - c_R);
    S[3] = fabs(u_R + c_R);
    Mod_Alpha = *max_element(S.begin(),S.end());

 
    Right_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha*(U[Cell_No + 1][0] - U[Cell_No][0]);
    Right_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha*(U[Cell_No + 1][1] - U[Cell_No][1]);
    Right_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha*(U[Cell_No + 1][2] - U[Cell_No][2]);
}



void MOVERS_Dissipation(int & Cell_No)
{
 double Mod_Alpha1 = 0.0,epsilon = 1e-6,Min1=0.0,Min2=0.0,Max1=0.0,Max2=0.0;//,delta=0.0,k=0.5;
 double Mod_Alpha2 = 0.0,Mod_Alpha3=0.0,F_L_1=0.0,F_R_1=0.0,F_L_2=0.0,F_R_2=0.0,F_L_3=0.0,F_R_3=0.0;
  vector<double> S(6,0.0);
 
//	Cell_No = i , Cell_No +1 = i+1, Cell_No -1 = i-1
 
 // at interface i Left state L and i+1 right state R
    Rho_L = Density[Cell_No];
    u_L = Velocity[Cell_No];
    c_L = C[Cell_No];
    P_L = Pressure[Cell_No];
    E_L = U[Cell_No][2];

    F_L_1 = Rho_L*u_L;
    F_L_2 = Rho_L*u_L*u_L + P_L;
    F_L_3 = (E_L + P_L)*u_L;
    
    Rho_R = Density[Cell_No + 1];
    u_R = Velocity[Cell_No + 1];
    c_R = C[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    E_R = U[Cell_No + 1][2];

    F_R_1 = Rho_R*u_R;
    F_R_2 = Rho_R*u_R*u_R + P_R;
    F_R_3 = (E_R + P_R)*u_R;

    S[0] = fabs(u_L - c_L);
    S[1] = fabs(u_L + c_L);
    S[2] = fabs(u_L);
    S[3] = fabs(u_R - c_R);
    S[4] = fabs(u_R + c_R);
    S[5] = fabs(u_R);

//	|alpha| = Delta F/ Delta U --- of Total Energy
    Mod_Alpha1 = fabs((F_R_1-F_L_1)/(Rho_R-Rho_L));
    Mod_Alpha2 = fabs((F_R_2-F_L_2)/((Rho_R*u_R)-(Rho_L*u_L)));
    Mod_Alpha3 = fabs((F_R_3-F_L_3)/(E_R-E_L));
//	Max = a>b ? (a>c ? a:c) :(b>c ? b : c);
    Max1 = (S[0]>S[1]) ? (S[0]>S[2]?S[0]:S[2]) : (S[1]>S[2]?S[1]:S[2]);
    Max2 = (S[3]>S[4]) ? (S[3]>S[5]?S[3]:S[5]) : (S[4]>S[5]?S[4]:S[5]);
    Lambda_Max = max(Max1, Max2);

//	Min = a<b ? (a<c ? a:c) :(b<c ? b : c);
    Min1 = (S[0]<S[1]) ? (S[0]<S[2]?S[0]:S[2]) : (S[1]<S[2]?S[1]:S[2]);
    Min2 = (S[3]<S[4]) ? (S[3]<S[5]?S[3]:S[5]) : (S[4]<S[5]?S[4]:S[5]);
    Lambda_Min = max(Min1, Min2);
    
    Mod_Alpha1 = (Mod_Alpha1 > Lambda_Max) ? Lambda_Max : (Mod_Alpha1 < Lambda_Min) ? Lambda_Min : Mod_Alpha1;
    Mod_Alpha1 = (fabs(Rho_R - Rho_L)<epsilon) ? Lambda_Min : (fabs(F_R_1-F_L_1)<epsilon) ? 0.0 :Mod_Alpha1;
    
    Mod_Alpha2 = (Mod_Alpha2 > Lambda_Max) ? Lambda_Max : (Mod_Alpha2 < Lambda_Min) ? Lambda_Min : Mod_Alpha2;
    Mod_Alpha2 = (fabs((Rho_R*u_R) - (Rho_L*u_L))<epsilon) ? Lambda_Min : (fabs(F_R_2-F_L_2)<epsilon) ? 0.0 :Mod_Alpha2;

    Mod_Alpha3 = (Mod_Alpha3 > Lambda_Max) ? Lambda_Max : (Mod_Alpha3 < Lambda_Min) ? Lambda_Min : Mod_Alpha3;
    Mod_Alpha3 = (fabs(E_R - E_L)<epsilon) ? Lambda_Min : (fabs(F_R_3-F_L_3)<epsilon) ? 0.0 :Mod_Alpha3;


    Right_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha1*(U[Cell_No + 1][0] - U[Cell_No][0]);
    Right_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha2*(U[Cell_No + 1][1] - U[Cell_No][1]);
    Right_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha3*(U[Cell_No + 1][2] - U[Cell_No][2]);
    
 // at interface i - 1  Left state L and i right state R
    Rho_L = Density[Cell_No - 1];
    u_L = Velocity[Cell_No - 1];
    c_L = C[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    E_L = U[Cell_No - 1][2];

    F_L_1 = Rho_L*u_L;
    F_L_2 = Rho_L*u_L*u_L + P_L;
    F_L_3 = (E_L + P_L)*u_L;
    
    Rho_R = Density[Cell_No];
    u_R = Velocity[Cell_No];
    c_R = C[Cell_No];
    P_R = Pressure[Cell_No];
    E_R = U[Cell_No][2];

    F_R_1 = Rho_R*u_R;
    F_R_2 = Rho_R*u_R*u_R + P_R;
    F_R_3 = (E_R + P_R)*u_R;

    S[0] = fabs(u_L - c_L);
    S[1] = fabs(u_L + c_L);
    S[2] = fabs(u_L);
    S[3] = fabs(u_R - c_R);
    S[4] = fabs(u_R + c_R);
    S[5] = fabs(u_R);

//	|alpha| = Delta F/ Delta U --- of Total Energy
    Mod_Alpha1 = fabs((F_R_1-F_L_1)/(Rho_R-Rho_L));
    Mod_Alpha2 = fabs((F_R_2-F_L_2)/((Rho_R*u_R)-(Rho_L*u_L)));
    Mod_Alpha3 = fabs((F_R_3-F_L_3)/(E_R-E_L));
//	Max = a>b ? (a>c ? a:c) :(b>c ? b : c);
    Max1 = (S[0]>S[1]) ? (S[0]>S[2]?S[0]:S[2]) : (S[1]>S[2]?S[1]:S[2]);
    Max2 = (S[3]>S[4]) ? (S[3]>S[5]?S[3]:S[5]) : (S[4]>S[5]?S[4]:S[5]);
    Lambda_Max = max(Max1, Max2);

//	Min = a<b ? (a<c ? a:c) :(b<c ? b : c);
    Min1 = (S[0]<S[1]) ? (S[0]<S[2]?S[0]:S[2]) : (S[1]<S[2]?S[1]:S[2]);
    Min2 = (S[3]<S[4]) ? (S[3]<S[5]?S[3]:S[5]) : (S[4]<S[5]?S[4]:S[5]);
    Lambda_Min = max(Min1, Min2);
    
    Mod_Alpha1 = (Mod_Alpha1 > Lambda_Max) ? Lambda_Max : (Mod_Alpha1 < Lambda_Min) ? Lambda_Min : Mod_Alpha1;
    Mod_Alpha1 = (fabs(Rho_R - Rho_L)<epsilon) ? Lambda_Min : (fabs(F_R_1-F_L_1)<epsilon) ? 0.0 :Mod_Alpha1;
    
    Mod_Alpha2 = (Mod_Alpha2 > Lambda_Max) ? Lambda_Max : (Mod_Alpha2 < Lambda_Min) ? Lambda_Min : Mod_Alpha2;
    Mod_Alpha2 = (fabs((Rho_R*u_R) - (Rho_L*u_L))<epsilon) ? Lambda_Min : (fabs(F_R_2-F_L_2)<epsilon) ? 0.0 :Mod_Alpha2;

    Mod_Alpha3 = (Mod_Alpha3 > Lambda_Max) ? Lambda_Max : (Mod_Alpha3 < Lambda_Min) ? Lambda_Min : Mod_Alpha3;
    Mod_Alpha3 = (fabs(E_R - E_L)<epsilon) ? Lambda_Min : (fabs(F_R_3-F_L_3)<epsilon) ? 0.0 :Mod_Alpha3;

    
    Left_Face_Dissipative_Flux[0] = 0.5*Mod_Alpha1*(U[Cell_No][0] - U[Cell_No-1][0]);
    Left_Face_Dissipative_Flux[1] = 0.5*Mod_Alpha2*(U[Cell_No][1] - U[Cell_No-1][1]);
    Left_Face_Dissipative_Flux[2] = 0.5*Mod_Alpha3*(U[Cell_No][2] - U[Cell_No-1][2]);

}
