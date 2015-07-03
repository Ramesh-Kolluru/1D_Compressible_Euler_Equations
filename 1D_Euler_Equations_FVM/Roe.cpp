#include "1D_Euler_FVM.h"


void ROE(int & Cell_No)
{
  double Mach_No_L=0.0,Mach_No_R=0.0,Mach_No=0.0,a_L=0.0,a=0.0,a_R=0.0,H_L=0.0,H=0.0,H_R=0.0,v_mag=0.0,Term1=0.0,delta_u_1=0.0,delta_u_2=0.0,delta_u_3=0.0;
 vector<double> K_1(3,0.0),K_2(3,0.0),K_3(3,0.0);
 double u_tilde=0.0,rho_tilde = 0.0,a_tilde =0.0,H_tilde=0.0,lambda_1=0.0,lambda_2=0.0,lambda_3=0.0,alpha_1 =0.0,alpha_2=0.0,alpha_3=0.0;
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;
 
//  cout<<"In Roe Function\n";  
//  cout<<"Cell No \t"<<Cell_No<<endl;
 
 /*
  * The dissipative flux to be determined is (1/2)*alpha_i*|lambda_i|K_i ,i = 1,2,3 in 1D
  * alpha_i are wave strengths, 
  * lambda_i = Eigen values
  * K_i - Right Eigen vectors
  * 
  * 
  */
 
 
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    a_L = C[Cell_No - 1];
    Mach_No_L = u_L/a_L;
    H_L = (E_L + P_L)/Rho_L;

        
// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    a = C[Cell_No];
    Mach_No = u/a;
    H = (E + P)/Rho;


// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    a_R = C[Cell_No + 1];
    Mach_No_R = u_R/a_R;
    H_R = (E_R + P_R)/Rho_R;

//				L	R		
//    Flux at i+1/2 interface = G+|i + G-|i+1

// Roe averaged values
    Term1 = 1.0/(sqrt(Rho) + sqrt(Rho_R));
    u_tilde = Term1*(sqrt(Rho)*u + sqrt(Rho_R)*u_R);
    v_mag = 0.5*u_tilde*u_tilde;
    H_tilde = Term1*(sqrt(Rho)*H + sqrt(Rho_R)*H_R);
    a_tilde = sqrt((gamma-1.0)*(H_tilde - v_mag));
    delta_u_1 = Rho_R - Rho;
    delta_u_2 = Rho_R*u_R - Rho*u;
    delta_u_3 = E_R - E;
    

// Roe averaged Eigen values
    lambda_1 = u_tilde - a_tilde;
    lambda_2 = u_tilde;
    lambda_3 = u_tilde + a_tilde;
// Roe averaged Right Eigen vectors
    K_1[0] =1;					K_2[0] = 1;		 	K_3[0] = 1;
    K_1[1] =u_tilde - a_tilde;			K_2[1] = u_tilde;	 	K_3[1] = u_tilde + a_tilde;
    K_1[2] =H_tilde - u_tilde*a_tilde;		K_2[2] = v_mag;			K_3[2] = H_tilde + u_tilde*a_tilde;
// Roe averaged Wave strengths
    alpha_2 = ((gamma - 1.0)/(a_tilde*a_tilde))*(delta_u_1*(H_tilde - u_tilde*u_tilde) + u_tilde*delta_u_2 - delta_u_3);
    alpha_1 = (0.5/a_tilde)*(delta_u_1*(lambda_3) - delta_u_2 - a_tilde*alpha_2);
    alpha_3 = delta_u_1 -(alpha_1 + alpha_2);
    
    Right_Face_Dissipative_Flux[0] = 0.5*(alpha_1*fabs(lambda_1)*K_1[0] + alpha_2*fabs(lambda_2)*K_2[0] + alpha_3*fabs(lambda_3)*K_3[0]);
    Right_Face_Dissipative_Flux[1] = 0.5*(alpha_1*fabs(lambda_1)*K_1[1] + alpha_2*fabs(lambda_2)*K_2[1] + alpha_3*fabs(lambda_3)*K_3[1]);
    Right_Face_Dissipative_Flux[2] = 0.5*(alpha_1*fabs(lambda_1)*K_1[2] + alpha_2*fabs(lambda_2)*K_2[2] + alpha_3*fabs(lambda_3)*K_3[2]);

//				L	R
//    Flux at i -1/2 interface = G+|i-1 + G-|i
   Term1 = 1.0/(sqrt(Rho) + sqrt(Rho_L));
    u_tilde = Term1*(sqrt(Rho)*u + sqrt(Rho_L)*u_L);
    v_mag = 0.5*u_tilde*u_tilde;
    H_tilde = Term1*(sqrt(Rho)*H + sqrt(Rho_L)*H_L);
    a_tilde = sqrt((gamma-1.0)*(H_tilde - v_mag));
    delta_u_1 = Rho - Rho_L;
    delta_u_2 = Rho*u - Rho_L*u_L;
    delta_u_3 = E - E_L ;
    

// Roe averaged Eigen values
    lambda_1 = u_tilde - a_tilde;
    lambda_2 = u_tilde;
    lambda_3 = u_tilde + a_tilde;
// Roe averaged Right Eigen vectors
    K_1[0] =1;					K_2[0] = 1;		 	K_3[0] = 1;
    K_1[1] =u_tilde - a_tilde;			K_2[1] = u_tilde;	 	K_3[1] = u_tilde + a_tilde;
    K_1[2] =H_tilde - u_tilde*a_tilde;		K_2[2] = v_mag;			K_3[2] = H_tilde + u_tilde*a_tilde;
// Roe averaged Wave strengths
    alpha_2 = ((gamma - 1.0)/(a_tilde*a_tilde))*(delta_u_1*(H_tilde - u_tilde*u_tilde) + u_tilde*delta_u_2 - delta_u_3);
    alpha_1 = (0.5/a_tilde)*(delta_u_1*(lambda_3) - delta_u_2 - a_tilde*alpha_2);
    alpha_3 = delta_u_1 -(alpha_1 + alpha_2);
    
    Left_Face_Dissipative_Flux[0] = 0.5*(alpha_1*fabs(lambda_1)*K_1[0] + alpha_2*fabs(lambda_2)*K_2[0] + alpha_3*fabs(lambda_3)*K_3[0]);
    Left_Face_Dissipative_Flux[1] = 0.5*(alpha_1*fabs(lambda_1)*K_1[1] + alpha_2*fabs(lambda_2)*K_2[1] + alpha_3*fabs(lambda_3)*K_3[1]);
    Left_Face_Dissipative_Flux[2] = 0.5*(alpha_1*fabs(lambda_1)*K_1[2] + alpha_2*fabs(lambda_2)*K_2[2] + alpha_3*fabs(lambda_3)*K_3[2]);
}